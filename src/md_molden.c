#include <md_molden.h>

#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_parse.h>
#include <core/md_arena_allocator.h>
#include <core/md_str_builder.h>

#include <md_util.h>
#include <md_system.h>
#include <md_gto.h>

#include <float.h>
#include <math.h>

/*
 * MOLDEN FILE FORMAT PARSER
 * 
 * This file implements parsing and data extraction for Molden format files.
 * The implementation follows the patterns established in md_vlx.c for consistency.
 * 
 * Molden Format Reference:
 * http://www.cmbi.ru.nl/molden/molden_format.html
 * 
 * File Structure:
 * Molden files are text-based and organized into blocks marked by [BLOCK_NAME]
 * 
 * Main Blocks:
 * - [Molden Format]: Optional version tag
 * - [Atoms]: Atomic positions (Angs or AU)
 * - [GTO]: Gaussian basis set definition
 * - [MO]: Molecular orbital data (coefficients, energies, occupancies, symmetries)
 * - [5D], [7F], [9G]: Specify spherical harmonic basis functions
 * 
 * Design Notes:
 * - Follows VLX structure with basis_set, atomic data, and orbital structures
 * - Uses arena allocator for memory management
 * - Coordinate conversion handled internally (Bohr <-> Ångström)
 * - Basis normalization follows Gaussian conventions
 */

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.5291772109029999

// Internal basis set representation (mirrors VLX structure)
typedef struct basis_set_func_t {
	uint8_t  type;          // Angular momentum quantum number (s=0, p=1, d=2, ...)
	uint8_t  param_count;   // Number of primitives in this contracted function
	uint16_t param_offset;  // Offset into parameter arrays
} basis_set_func_t;

typedef struct basis_set_basis_t {
	uint8_t  max_type;          // Maximum angular momentum for this atom type
	uint8_t  basis_func_count;  // Number of basis functions
	uint16_t basis_func_offset; // Offset into basis function array
} basis_set_basis_t;

typedef struct basis_set_t {
	str_t identifier;

	struct {
		size_t count;
		double* exponents;
		double* coefficients;
	} param;

	struct {
		size_t count;
		basis_set_func_t* data;
	} basis_func;

	// Atom basis entries indexed by atomic number
	// 0 is NULL entry, 1 = Hydrogen, 2 = Helium, etc.
	struct {
		size_t count;
		basis_set_basis_t* data;
	} atom_basis;
} basis_set_t;

// 1D data array (for energies, occupancies, etc.)
typedef struct md_molden_1d_data_t {
	size_t  size;
	double* data;
} md_molden_1d_data_t;

// 2D data array (for coefficient matrices)
typedef struct md_molden_2d_data_t {
	size_t  size[2];  // [num_rows, num_cols]
	double* data;     // Row-major storage
} md_molden_2d_data_t;

// Molecular orbital data structure
typedef struct md_molden_orbital_t {
	md_molden_2d_data_t coefficients;  // [num_ao x num_mo]
	md_molden_1d_data_t energy;        // [num_mo]
	md_molden_1d_data_t occupancy;     // [num_mo]
	str_t* symmetry;                    // [num_mo] symmetry labels (optional)
	str_t* spin;                        // [num_mo] spin labels (optional)
} md_molden_orbital_t;

// Self Consistent Field data
typedef struct md_molden_scf_t {
	md_molden_scf_type_t type;
	size_t homo_idx[2];  // HOMO indices for [alpha, beta]
	size_t lumo_idx[2];  // LUMO indices for [alpha, beta]

	md_molden_orbital_t alpha;
	md_molden_orbital_t beta;
} md_molden_scf_t;

// Main Molden data structure
typedef struct md_molden_t {
	basis_set_t basis_set;

	str_t basis_set_ident;
	md_molden_coord_unit_t coord_unit;

	size_t number_of_atoms;
	size_t number_of_alpha_electrons;
	size_t number_of_beta_electrons;

	double molecular_charge;
	size_t spin_multiplicity;

	// Arrays (length = number_of_atoms)
	dvec3_t* atom_coordinates;
	md_element_t* atomic_numbers;

	int* ao_to_atom_idx;  // Maps atomic orbitals to atom indices

	// SCF data
	md_molden_scf_t scf;

	// Atomic orbital data represented as GTOs
	md_gto_data_t gto_data;

	// Basis function format flags
	bool use_5d;   // Use 5d instead of 6d (default true for d orbitals)
	bool use_7f;   // Use 7f instead of 10f (default true for f orbitals)
	bool use_9g;   // Use 9g instead of 15g (default true for g orbitals)

	struct md_allocator_i* arena;
} md_molden_t;

// =============================
// Forward declarations
// =============================

static bool parse_molden_file(md_molden_t* molden, str_t filename);
static bool parse_atoms_block(md_molden_t* molden, md_buffered_reader_t* reader, str_t* next_block_line);
static bool parse_gto_block(md_molden_t* molden, md_buffered_reader_t* reader, str_t* next_block_line);
static bool parse_mo_block(md_molden_t* molden, md_buffered_reader_t* reader, str_t* next_block_line);
static void normalize_basis_set(basis_set_t* basis_set);
static void extract_gto_data(md_gto_data_t* out_data, const dvec3_t* atom_coordinates, const md_element_t* atomic_numbers, size_t number_of_atoms, const basis_set_t* basis_set, md_allocator_i* alloc);
static size_t extract_ao_to_atom_idx(int* out_ao_to_atom, const md_element_t* atomic_numbers, size_t number_of_atoms, const basis_set_t* basis_set);

// =============================
// Utility functions
// =============================

static int char_to_angular_momentum_type(char c) {
	switch (c) {
	case 's': case 'S': return 0;
	case 'p': case 'P': return 1;
	case 'd': case 'D': return 2;
	case 'f': case 'F': return 3;
	case 'g': case 'G': return 4;
	case 'h': case 'H': return 5;
	case 'i': case 'I': return 6;
	default: return -1;
	}
}

static inline basis_set_basis_t* basis_set_get_atom_basis(const basis_set_t* basis_set, int atomic_number) {
	if (atomic_number < (int)basis_set->atom_basis.count) {
		return basis_set->atom_basis.data + atomic_number;
	}
	return NULL;
}

static inline int compute_max_angular_momentum(const basis_set_t* basis_set, const md_element_t* atomic_numbers, size_t count) {
	ASSERT(basis_set);
	ASSERT(atomic_numbers);
	int max_angl = 0;
	for (size_t i = 0; i < count; ++i) {
		const basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_numbers[i]);
		if (atom_basis) {
			max_angl = MAX(max_angl, (int)atom_basis->max_type);
		}
	}
	return max_angl;
}

// These constants and tables are ported from VLX implementation for spherical harmonics
#define d3  3.464101615137754587
#define f5  1.581138830084189666
#define f15 7.745966692414833770
#define f3  1.224744871391589049
#define g35 4.0 * 5.916079783099616042
#define g17 4.0 * 4.183300132670377739
#define g5  4.0 * 2.236067977499789696
#define g2  4.0 * 1.581138830084189666

static const double		S_factors[] = {1.0};
static const uint8_t    S_indices[] = {0};
static const uint8_t	S_num_fac[] = {1};

static const double		P_factors[] = {1.0, 1.0, 1.0};
static const uint8_t	P_indices[] = {1, 2, 0};
static const uint8_t	P_offsets[] = {0, 1, 2};
static const uint8_t	P_num_fac[] = {1, 1, 1};

static const double		D_factors[] = {d3, d3, -1.0, -1.0, 2.0, d3, 0.5 * d3, -0.5 * d3};
static const uint8_t    D_indices[] = {1, 4, 0, 3, 5, 2, 0, 3};
static const uint8_t    D_offsets[] = {0, 1, 2, 5, 6};
static const uint8_t	D_num_fac[] = {1, 1, 3, 1, 2};

static const double		F_factors[] = {3.0 * f5, -f5, f15, 4.0 * f3, -f3, -f3, 2.0, -3.0, -3.0, 4.0 * f3, -f3, -f3, 0.5 * f15, -0.5 * f15, f5, -3.0 * f5};
static const uint8_t    F_indices[] = {1, 6, 4, 8, 1, 6, 9, 2, 7, 5, 0, 3, 2, 7, 0, 3};
static const uint8_t	F_offsets[] = {0, 2, 3, 6, 9, 12, 14};
static const uint8_t	F_num_fac[] = {2, 1, 3, 3, 3, 2, 2};

static const double		G_factors[] = {
	g35, -g35, 3.0 * g17, -g17, 6.0 * g5, -g5, -g5, 4.0 * g2, -3.0 * g2, -3.0 * g2,
	8.0, 3.0, 3.0, 6.0, -24.0, -24.0, 4.0 * g2, -3.0 * g2, -3.0 * g2, 3.0 * g5,
	-3.0 * g5, -0.5 * g5, 0.5 * g5,  g17,  -3.0 * g17, 0.25 * g35, 0.25 * g35, -1.50 * g35};
static const uint8_t    G_indices[] = {1, 6, 4, 11, 8, 1, 6, 13, 4, 11, 14, 0, 10, 3, 5, 12, 9, 2, 7, 5, 12, 0, 10, 2, 7, 0, 10, 3};
static const uint8_t	G_offsets[] = {0, 2, 4, 7, 10, 16, 19, 23, 25};
static const uint8_t	G_num_fac[] = {2, 2, 3, 3, 6, 3, 4, 2, 3};

#undef d3
#undef f5
#undef f15
#undef f3
#undef g35
#undef g17
#undef g5 
#undef g2 

static inline int spherical_momentum_num_components(int angl) {
	return angl * 2 + 1;
}

static inline int spherical_momentum_num_factors(int angl, int isph) {
	switch(angl) {
	case 0:
		ASSERT(isph < ARRAY_SIZE(S_num_fac));
		return S_num_fac[isph];
	case 1:
		ASSERT(isph < ARRAY_SIZE(P_num_fac));
		return P_num_fac[isph];
	case 2:
		ASSERT(isph < ARRAY_SIZE(D_num_fac));
		return D_num_fac[isph];
	case 3:
		ASSERT(isph < ARRAY_SIZE(F_num_fac));
		return F_num_fac[isph];
	case 4:
		ASSERT(isph < ARRAY_SIZE(G_num_fac));
		return G_num_fac[isph];
	default:
		ASSERT(false);
		return 0;
	}
}

static inline const double* spherical_momentum_factors(int angl, int isph) {
	switch(angl) {
	case 0:
		ASSERT(isph == 0);
		return S_factors;
	case 1:
		ASSERT(isph < ARRAY_SIZE(P_offsets));
		return P_factors + P_offsets[isph];
	case 2:
		ASSERT(isph < ARRAY_SIZE(D_offsets));
		return D_factors + D_offsets[isph];
	case 3:
		ASSERT(isph < ARRAY_SIZE(F_offsets));
		return F_factors + F_offsets[isph];
	case 4:
		ASSERT(isph < ARRAY_SIZE(G_offsets));
		return G_factors + G_offsets[isph];
	default:
		ASSERT(false);
		return NULL;
	}
}

static inline const uint8_t* spherical_momentum_indices(int angl, int isph) {
	switch(angl) {
	case 0:
		ASSERT(isph == 0);
		return S_indices;
	case 1:
		ASSERT(isph < ARRAY_SIZE(P_offsets));
		return P_indices + P_offsets[isph];
	case 2:
		ASSERT(isph < ARRAY_SIZE(D_offsets));
		return D_indices + D_offsets[isph];
	case 3:
		ASSERT(isph < ARRAY_SIZE(F_offsets));
		return F_indices + F_offsets[isph];
	case 4:
		ASSERT(isph < ARRAY_SIZE(G_offsets));
		return G_indices + G_offsets[isph];
	default:
		ASSERT(false);
		return NULL;
	}
}

typedef uint8_t lmn_t[3];

// S: 0
static const lmn_t S_lmn[1] = {{0,0,0}};
// P: x y z
static const lmn_t P_lmn[3] = {{1,0,0}, {0,1,0}, {0,0,1}};
// D: xx xy xz yy yz zz
static const lmn_t D_lmn[6] = {{2,0,0}, {1,1,0}, {1,0,1}, {0,2,0}, {0,1,1}, {0,0,2}};
// F: xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
static const lmn_t F_lmn[10] = {{3,0,0}, {2,1,0}, {2,0,1}, {1,2,0}, {1,1,1}, {1,0,2}, {0,3,0}, {0,2,1}, {0,1,2}, {0,0,3}};
// G: xxxx xxxy xxxz xxyy xxyz xxzz xyyy xyyz xyzz xzzz yyyy yyyz yyzz yzzz zzzz
static const lmn_t G_lmn[15] = {{4,0,0}, {3,1,0}, {3,0,1}, {2,2,0}, {2,1,1}, {2,0,2}, {1,3,0}, {1,2,1}, {1,1,2}, {1,0,3}, {0,4,0}, {0,3,1}, {0,2,2}, {0,1,3}, {0,0,4}};

static inline const lmn_t* cartesian_angular_momentum(int angl) {
	switch (angl) {
	case 0: return S_lmn;
	case 1: return P_lmn;
	case 2: return D_lmn;
	case 3: return F_lmn;
	case 4: return G_lmn;
	default: ASSERT(false); return NULL;
	}
}

typedef struct basis_func_t {
    int type;
    int count;
    double* exponents;
    double* normalization_coefficients;
} basis_func_t;

static inline basis_func_t get_basis_func(const basis_set_t* basis_set, int basis_func_idx) {
	basis_set_func_t func = basis_set->basis_func.data[basis_func_idx];
	return (basis_func_t) {
		.type = func.type,
		.count = func.param_count,
		.exponents = basis_set->param.exponents + func.param_offset,
		.normalization_coefficients = basis_set->param.coefficients + func.param_offset,
	};
}

static size_t basis_set_extract_atomic_basis_func_angl(basis_func_t* out_funcs, size_t cap_funcs, const basis_set_t* basis_set, int atomic_number, int angl) {
    size_t count = 0;

    basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_number);
    if (atom_basis) {
        int beg = atom_basis->basis_func_offset;
        int end = atom_basis->basis_func_offset + atom_basis->basis_func_count;
        for (int i = beg; i < end; ++i) {
            if (count == cap_funcs) break;
            if (basis_set->basis_func.data[i].type == angl) {
                out_funcs[count++] = get_basis_func(basis_set, i);
            }
        }
    }

    return count;
}

// =============================
// Public API Implementation
// =============================

md_molden_t* md_molden_create(md_allocator_i* backing) {
	ASSERT(backing);
	md_allocator_i* arena = md_arena_allocator_create(backing, MEGABYTES(1));
	ASSERT(arena);
	md_molden_t* molden = md_alloc(arena, sizeof(md_molden_t));
	if (!molden) {
		MD_LOG_ERROR("Failed to allocate memory for Molden object");
		return NULL;
	}
	MEMSET(molden, 0, sizeof(md_molden_t));
	molden->arena = arena;
	
	// Set defaults for basis function format
	molden->use_5d = true;
	molden->use_7f = true;
	molden->use_9g = true;

	return molden;
}

void md_molden_reset(md_molden_t* molden) {
	if (molden) {
		ASSERT(molden->arena);
		md_allocator_i* arena = molden->arena;
		md_arena_allocator_reset(arena);
		MEMSET(molden, 0, sizeof(md_molden_t));
		molden->arena = arena;
		molden->use_5d = true;
		molden->use_7f = true;
		molden->use_9g = true;
	}
}

void md_molden_destroy(md_molden_t* molden) {
	if (molden) {
		md_arena_allocator_destroy(molden->arena);
	} else {
		MD_LOG_DEBUG("Attempt to destroy NULL Molden object");
	}
}

bool md_molden_parse_file(md_molden_t* molden, str_t filename) {
	ASSERT(molden);
	if (!molden) {
		MD_LOG_ERROR("NULL Molden object provided");
		return false;
	}

	return parse_molden_file(molden, filename);
}

// Molecular properties
size_t md_molden_number_of_atoms(const md_molden_t* molden) {
	if (molden) return molden->number_of_atoms;
	return 0;
}

size_t md_molden_number_of_alpha_electrons(const md_molden_t* molden) {
	if (molden) return molden->number_of_alpha_electrons;
	return 0;
}

size_t md_molden_number_of_beta_electrons(const md_molden_t* molden) {
	if (molden) return molden->number_of_beta_electrons;
	return 0;
}

double md_molden_molecular_charge(const md_molden_t* molden) {
	if (molden) return molden->molecular_charge;
	return 0;
}

size_t md_molden_spin_multiplicity(const md_molden_t* molden) {
	if (molden) return molden->spin_multiplicity;
	return 0;
}

str_t md_molden_basis_set_ident(const md_molden_t* molden) {
	if (molden) return molden->basis_set_ident;
	return (str_t){0};
}

// Atomic data
const dvec3_t* md_molden_atom_coordinates(const md_molden_t* molden) {
	if (molden) return molden->atom_coordinates;
	return NULL;
}

const md_element_t* md_molden_atomic_numbers(const md_molden_t* molden) {
	if (molden) return molden->atomic_numbers;
	return NULL;
}

const int* md_molden_ao_to_atom_idx(const md_molden_t* molden) {
	if (molden) return molden->ao_to_atom_idx;
	return NULL;
}

// SCF data
md_molden_scf_type_t md_molden_scf_type(const md_molden_t* molden) {
	if (molden) return molden->scf.type;
	return MD_MOLDEN_SCF_TYPE_UNKNOWN;
}

size_t md_molden_scf_homo_idx(const md_molden_t* molden, md_molden_mo_type_t type) {
	if (molden) {
		if (type == MD_MOLDEN_MO_TYPE_ALPHA) {
			return molden->scf.homo_idx[0];
		} else if (type == MD_MOLDEN_MO_TYPE_BETA) {
			return molden->scf.homo_idx[1];
		}
	}
	return 0;
}

size_t md_molden_scf_lumo_idx(const md_molden_t* molden, md_molden_mo_type_t type) {
	if (molden) {
		if (type == MD_MOLDEN_MO_TYPE_ALPHA) {
			return molden->scf.lumo_idx[0];
		} else if (type == MD_MOLDEN_MO_TYPE_BETA) {
			return molden->scf.lumo_idx[1];
		}
	}
	return 0;
}

size_t md_molden_scf_number_of_atomic_orbitals(const md_molden_t* molden) {
	if (molden) {
		return molden->scf.alpha.coefficients.size[0];
	}
	return 0;
}

size_t md_molden_scf_number_of_molecular_orbitals(const md_molden_t* molden) {
	if (molden) {
		return molden->scf.alpha.coefficients.size[1];
	}
	return 0;
}

const double* md_molden_scf_mo_occupancy(const md_molden_t* molden, md_molden_mo_type_t type) {
	if (molden) {
		if (type == MD_MOLDEN_MO_TYPE_ALPHA) {
			return molden->scf.alpha.occupancy.data;
		} else if (type == MD_MOLDEN_MO_TYPE_BETA) {
			return molden->scf.beta.occupancy.data;
		}
	}
	return NULL;
}

const double* md_molden_scf_mo_energy(const md_molden_t* molden, md_molden_mo_type_t type) {
	if (molden) {
		if (type == MD_MOLDEN_MO_TYPE_ALPHA) {
			return molden->scf.alpha.energy.data;
		} else if (type == MD_MOLDEN_MO_TYPE_BETA) {
			return molden->scf.beta.energy.data;
		}
	}
	return NULL;
}

const str_t* md_molden_scf_mo_symmetry(const md_molden_t* molden, md_molden_mo_type_t type) {
	if (molden) {
		if (type == MD_MOLDEN_MO_TYPE_ALPHA) {
			return molden->scf.alpha.symmetry;
		} else if (type == MD_MOLDEN_MO_TYPE_BETA) {
			return molden->scf.beta.symmetry;
		}
	}
	return NULL;
}

const str_t* md_molden_scf_mo_spin(const md_molden_t* molden, md_molden_mo_type_t type) {
	if (molden) {
		if (type == MD_MOLDEN_MO_TYPE_ALPHA) {
			return molden->scf.alpha.spin;
		} else if (type == MD_MOLDEN_MO_TYPE_BETA) {
			return molden->scf.beta.spin;
		}
	}
	return NULL;
}

// Basis set and GTO extraction
bool md_molden_scf_extract_gto_data(md_gto_data_t* out_gto_data, const md_molden_t* molden, double cutoff_value, md_allocator_i* alloc) {
	if (!molden || !out_gto_data || !alloc) {
		MD_LOG_ERROR("Invalid arguments to md_molden_scf_extract_gto_data");
		return false;
	}

	if (molden->number_of_atoms == 0) {
		MD_LOG_ERROR("No atoms in Molden object");
		return false;
	}

	MEMSET(out_gto_data, 0, sizeof(md_gto_data_t));
	extract_gto_data(out_gto_data, molden->atom_coordinates, molden->atomic_numbers, molden->number_of_atoms, &molden->basis_set, alloc);

	(void)cutoff_value; // TODO: Apply cutoff if needed

	return true;
}

size_t md_molden_mo_gto_count(const md_molden_t* molden) {
	if (molden) {
		return molden->gto_data.num_pgtos;
	}
	return 0;
}

size_t md_molden_mo_gto_extract(md_gto_t gtos[], const md_molden_t* molden, size_t mo_idx, md_molden_mo_type_t type, double value_cutoff) {
	if (!gtos || !molden) {
		MD_LOG_ERROR("Invalid arguments to md_molden_mo_gto_extract");
		return 0;
	}

	const md_molden_orbital_t* orb = (type == MD_MOLDEN_MO_TYPE_ALPHA) ? &molden->scf.alpha : &molden->scf.beta;
	
	if (mo_idx >= orb->coefficients.size[1]) {
		MD_LOG_ERROR("MO index %zu out of range (max %zu)", mo_idx, orb->coefficients.size[1]);
		return 0;
	}

	size_t num_aos = orb->coefficients.size[0];
	
	// Extract molecular orbital coefficients for this MO
	double* mo_coeffs = md_alloc(md_get_temp_allocator(), num_aos * sizeof(double));
	if (!mo_coeffs) {
		MD_LOG_ERROR("Failed to allocate temporary memory for MO coefficients");
		return 0;
	}

	// Extract column mo_idx from coefficients matrix (stored as row-major AO x MO)
	for (size_t ao = 0; ao < num_aos; ao++) {
		mo_coeffs[ao] = orb->coefficients.data[ao * orb->coefficients.size[1] + mo_idx];
	}

	// Now extract PGTOs using the MO coefficients
	size_t count = 0;
	int natoms = (int)molden->number_of_atoms;
	int max_angl = compute_max_angular_momentum(&molden->basis_set, molden->atomic_numbers, molden->number_of_atoms);
	
	basis_func_t basis_funcs[128];
	size_t mo_coeff_idx = 0;

	// azimuthal quantum number: s,p,d,f,...
	for (int angl = 0; angl <= max_angl; angl++) {
		int nsph = spherical_momentum_num_components(angl);
		const lmn_t* lmn = cartesian_angular_momentum(angl);
		
		// magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...
		for (int isph = 0; isph < nsph; isph++) {
			// prepare Cartesian components
			int lx[8];
			int ly[8];
			int lz[8];
			int			      ncomp	= spherical_momentum_num_factors(angl, isph);
			const double*	fcarts  = spherical_momentum_factors(angl, isph);
			const uint8_t*	indices = spherical_momentum_indices(angl, isph);

			for (int icomp = 0; icomp < ncomp; icomp++) {
				int cartind = indices[icomp];
				lx[icomp] = lmn[cartind][0];
				ly[icomp] = lmn[cartind][1];
				lz[icomp] = lmn[cartind][2];
			}

			// go through atoms
			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				// process coordinates (convert from Ångström to Bohr)
				float x = (float)(molden->atom_coordinates[atomidx].x * ANGSTROM_TO_BOHR);
				float y = (float)(molden->atom_coordinates[atomidx].y * ANGSTROM_TO_BOHR);
				float z = (float)(molden->atom_coordinates[atomidx].z * ANGSTROM_TO_BOHR);

				int idelem = molden->atomic_numbers[atomidx];

				size_t num_basis_funcs = basis_set_extract_atomic_basis_func_angl(basis_funcs, ARRAY_SIZE(basis_funcs), &molden->basis_set, idelem, angl);

				// process atomic orbitals
				for (size_t funcidx = 0; funcidx < num_basis_funcs; funcidx++) {
					const double mo_coeff = (mo_coeff_idx < num_aos) ? mo_coeffs[mo_coeff_idx++] : 0.0;

					// process primitives
					basis_func_t basis_func = basis_funcs[funcidx];
					ASSERT(basis_func.type == angl);
					const int        nprims = basis_func.count;
					const double* exponents = basis_func.exponents;
					const double* normcoefs = basis_func.normalization_coefficients;

					for (int iprim = 0; iprim < nprims; iprim++) {
						double alpha = exponents[iprim];
						double coef1 = normcoefs[iprim];

						// transform from Cartesian to spherical harmonics
						for (int icomp = 0; icomp < ncomp; icomp++) {
							double fcart = fcarts[icomp];
							double coeff = coef1 * fcart * mo_coeff;

							// Apply cutoff
							if (value_cutoff > 0.0 && fabs(coeff) < value_cutoff) {
								continue;
							}

							gtos[count].x		= x;
							gtos[count].y		= y;
							gtos[count].z		= z;
							gtos[count].coeff	= (float)coeff;
							gtos[count].alpha	= (float)alpha; 
							gtos[count].cutoff	= FLT_MAX;
							gtos[count].i		= (uint8_t)lx[icomp];
							gtos[count].j		= (uint8_t)ly[icomp];
							gtos[count].k		= (uint8_t)lz[icomp];
							gtos[count].l		= (uint8_t)angl;

							count += 1;
						}
					}
				}
			}
		}
	}

	return count;
}

// System integration
bool md_molden_system_init(md_system_t* sys, const md_molden_t* molden, md_allocator_i* alloc) {
	ASSERT(sys);
	ASSERT(molden);
	ASSERT(alloc);

	if (molden->number_of_atoms == 0) {
		MD_LOG_ERROR("The Molden object contains no atoms");
		return false;
	}

	size_t capacity = ROUND_UP(molden->number_of_atoms, 16);

	MEMSET(sys, 0, sizeof(md_system_t));

	md_array_resize(sys->atom.x,        capacity, alloc);
	md_array_resize(sys->atom.y,        capacity, alloc);
	md_array_resize(sys->atom.z,        capacity, alloc);
	md_array_resize(sys->atom.type_idx, capacity, alloc);
	md_array_resize(sys->atom.flags,    capacity, alloc);

	MEMSET(sys->atom.x,        0, md_array_bytes(sys->atom.x));
	MEMSET(sys->atom.y,        0, md_array_bytes(sys->atom.y));
	MEMSET(sys->atom.z,        0, md_array_bytes(sys->atom.z));
	MEMSET(sys->atom.type_idx, 0, md_array_bytes(sys->atom.type_idx));
	MEMSET(sys->atom.flags,    0, md_array_bytes(sys->atom.flags));

	md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unknown"), 0, 0.0f, 0.0f, 0, alloc);

	for (size_t i = 0; i < molden->number_of_atoms; ++i) {
		sys->atom.x[i] = (float)molden->atom_coordinates[i].x;
		sys->atom.y[i] = (float)molden->atom_coordinates[i].y;
		sys->atom.z[i] = (float)molden->atom_coordinates[i].z;
		
		md_atomic_number_t z = molden->atomic_numbers[i];
		str_t sym  = md_atomic_number_symbol(z);
		float mass = md_atomic_number_mass(z);
		float radius = md_atomic_number_vdw_radius(z);

		md_atom_type_idx_t type_idx = md_atom_type_find_or_add(&sys->atom.type, sym, z, mass, radius, 0, alloc);
		sys->atom.type_idx[i] = type_idx;
	}

	sys->atom.count = molden->number_of_atoms;

	return true;
}

static bool molden_sys_init_from_str(md_system_t* mol, str_t str, const void* arg, md_allocator_i* alloc) {
	(void)mol;
	(void)str;
	(void)arg;
	(void)alloc;
	MD_LOG_ERROR("Molden string parsing is not implemented");
	return false;
}

static bool molden_sys_init_from_file(md_system_t* mol, str_t filename, const void* arg, md_allocator_i* alloc) {
	(void)arg;
	md_molden_t* molden = md_molden_create(md_get_heap_allocator());

	bool success = false;
	if (md_molden_parse_file(molden, filename)) {
		success = md_molden_system_init(mol, molden, alloc);
	}

	md_molden_destroy(molden);
	return success;
}

static md_system_loader_i molden_loader = {
	molden_sys_init_from_str,
	molden_sys_init_from_file
};

md_system_loader_i* md_molden_system_loader(void) {
	return &molden_loader;
}

// =============================
// Parsing Implementation
// =============================

static inline double compute_overlap(basis_func_t func, int i, int j) {
	const double fab  = 1.0 / (func.exponents[i] + func.exponents[j]);
	const double fab2 = fab * fab;
	const double ovl = func.normalization_coefficients[i] * func.normalization_coefficients[j] * pow(PI * fab, 1.5);

	switch (func.type) {
	case 0: return ovl;
	case 1: return 0.5 * fab * ovl;
	case 2: return 3.0 * fab2 * ovl;
	case 3: return 7.5 * fab2 * fab * ovl;
	case 4: return 420.0 * fab2 * fab2 * ovl;
	case 5: return 1890.0 * fab2 * fab2 * fab * ovl;
	case 6: return 41580.0 * fab2 * fab2 * fab2 * ovl;
	default:
		ASSERT(false);
		return 0;
	}
}

static void rescale_basis_func(basis_func_t func) {
	const double fpi = 2.0 / PI;

	for (int i = 0; i < func.count; i++) {
		func.normalization_coefficients[i] *= pow(func.exponents[i] * fpi, 0.75);
	}

	if (func.type < 0 || 6 < func.type) {
		MD_LOG_DEBUG("Invalid basis function type supplied in rescaling");
		return;
	}

	static const double f_table[] = {
		0,
		2.0,
		1.15470053837925152902, // 2.0 / sqrt(3.0)
		1.03279555898864450271, // 4.0 / sqrt(15.0)
		0.19518001458970663587, // 2.0 / sqrt(105.0)
		0.13012000972647109058, // 4.0 / sqrt(945.0)
		0.03923265908909997910, // 4.0 / sqrt(10395.0)
	};

	double f = f_table[func.type];
	double e = (double)func.type * 0.5;

	for (int i = 0; i < func.count; i++) {
		func.normalization_coefficients[i] *= pow(f * func.exponents[i], e);
	}
}

static void normalize_basis_set(basis_set_t* basis_set) {
	for (size_t func_idx = 0; func_idx < basis_set->basis_func.count; ++func_idx) {
		basis_func_t func = get_basis_func(basis_set, (int)func_idx);
		// uncontracted basis, set expansion coeficient to 1.0
		if (func.count == 1) func.normalization_coefficients[0] = 1.0;

		// normalize primitive GBFs
		rescale_basis_func(func);

		// compute overlap
		double ovl = 0.0;
		for (int i = 0; i < func.count; i++) {
			ovl += compute_overlap(func, i, i);
			for (int j = i + 1; j < func.count; j++) {
				ovl += 2.0 * compute_overlap(func, i, j);
			}
		}

		// renormalize primitive BFs
		ovl = 1.0 / sqrt(ovl);
		for (int i = 0; i < func.count; i++) {
			func.normalization_coefficients[i] *= ovl;
		}
	}
}

static size_t extract_ao_to_atom_idx(int* out_ao_to_atom, const md_element_t* atomic_numbers, size_t number_of_atoms, const basis_set_t* basis_set) {
	int natoms = (int)number_of_atoms;
	int max_angl = compute_max_angular_momentum(basis_set, atomic_numbers, number_of_atoms);

	size_t count = 0;

    basis_func_t basis_funcs[128];

	// azimuthal quantum number: s,p,d,f,...
	for (int angl = 0; angl <= max_angl; angl++) {
		int nsph = spherical_momentum_num_components(angl);
		// magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...
		for (int isph = 0; isph < nsph; isph++) {
			int	ncomp = spherical_momentum_num_factors(angl, isph);

			// go through atoms
			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				int idelem = atomic_numbers[atomidx];
				size_t num_ao = basis_set_extract_atomic_basis_func_angl(basis_funcs, ARRAY_SIZE(basis_funcs), basis_set, idelem, angl);

				for (size_t iao = 0; iao < num_ao; iao++) {
					out_ao_to_atom[count] = atomidx;
					count += 1;
				}
			}
		}
	}
	return count;
}

static void extract_gto_data(md_gto_data_t* out_data, const dvec3_t* atom_coordinates, const md_element_t* atomic_numbers, size_t number_of_atoms, const basis_set_t* basis_set, md_allocator_i* alloc) {
	int natoms = (int)number_of_atoms;
	int max_angl = compute_max_angular_momentum(basis_set, atomic_numbers, number_of_atoms);

	uint32_t coeff_idx = 0; // same as the cgto_idx

	// azimuthal quantum number: s,p,d,f,...
	for (int angl = 0; angl <= max_angl; angl++) {
		int nsph = spherical_momentum_num_components(angl);
		const lmn_t* lmn = cartesian_angular_momentum(angl);
		// magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...
		for (int isph = 0; isph < nsph; isph++) {
			// prepare Cartesian components (Maximum number of components should be 6 here for the currently supported basis sets)
			int lx[8];
			int ly[8];
			int lz[8];
			int			      ncomp	= spherical_momentum_num_factors(angl, isph);
			const double*	fcarts  = spherical_momentum_factors(angl, isph);
			const uint8_t*	indices = spherical_momentum_indices(angl, isph);

			for (int icomp = 0; icomp < ncomp; icomp++) {
				int cartind = indices[icomp];

				lx[icomp] = lmn[cartind][0];
				ly[icomp] = lmn[cartind][1];
				lz[icomp] = lmn[cartind][2];
			}

			// go through atoms
			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				// process coordinates
				// Conversion from Ångström to Bohr
				float x = (float)(atom_coordinates[atomidx].x * ANGSTROM_TO_BOHR);
				float y = (float)(atom_coordinates[atomidx].y * ANGSTROM_TO_BOHR);
				float z = (float)(atom_coordinates[atomidx].z * ANGSTROM_TO_BOHR);

				int idelem = atomic_numbers[atomidx];

				basis_func_t basis_funcs[128];
				size_t num_basis_funcs = basis_set_extract_atomic_basis_func_angl(basis_funcs, ARRAY_SIZE(basis_funcs), basis_set, idelem, angl);

				// process atomic orbitals
				for (size_t funcidx = 0; funcidx < num_basis_funcs; funcidx++) {
                    vec4_t   cgto_xyzr = {x, y, z, FLT_MAX};
                    uint32_t cgto_offset = (uint32_t)out_data->num_pgtos;

					// process primitives
					basis_func_t basis_func = basis_funcs[funcidx];
					ASSERT(basis_func.type == angl);
					const int        nprims = basis_func.count;
					const double* exponents = basis_func.exponents;
					const double* normcoefs = basis_func.normalization_coefficients;

					size_t new_pgto_count = out_data->num_pgtos + (size_t)nprims * (size_t)ncomp;
					md_array_ensure(out_data->pgto_alpha,  new_pgto_count, alloc);
					md_array_ensure(out_data->pgto_coeff,  new_pgto_count, alloc);
					md_array_ensure(out_data->pgto_radius, new_pgto_count, alloc);
					md_array_ensure(out_data->pgto_ijkl,   new_pgto_count, alloc);

					for (int iprim = 0; iprim < nprims; iprim++) {
						double alpha = exponents[iprim];
						double coef1 = normcoefs[iprim];

						// transform from Cartesian to spherical harmonics
						for (int icomp = 0; icomp < ncomp; icomp++) {
							uint32_t packed_ijkl = md_gto_pack_ijkl(lx[icomp], ly[icomp], lz[icomp], angl);
							md_array_push_no_grow(out_data->pgto_alpha, (float)alpha);
							md_array_push_no_grow(out_data->pgto_coeff, (float)(coef1 * fcarts[icomp]));
							md_array_push_no_grow(out_data->pgto_radius, FLT_MAX);
							md_array_push_no_grow(out_data->pgto_ijkl,  packed_ijkl);
							out_data->num_pgtos += 1;
						}
					}
					md_array_push(out_data->cgto_xyzr, cgto_xyzr, alloc);
					md_array_push(out_data->cgto_offset, cgto_offset, alloc);
					out_data->num_cgtos += 1;
				}
			}
		}
	}

	if (out_data->num_cgtos > 0) {
		md_array_push(out_data->cgto_offset, (uint32_t)out_data->num_pgtos, alloc);
	}
}

static bool parse_atoms_block(md_molden_t* molden, md_buffered_reader_t* reader, str_t* next_block_line) {
	md_allocator_i* alloc = molden->arena;
	str_t line;

	MD_LOG_DEBUG("Entered parse_atoms_block");

	// Read all atom lines into temporary storage
	str_t atom_lines[256];  // Max 256 atoms
	size_t atom_count = 0;
	
	while (md_buffered_reader_extract_line(&line, reader) && atom_count < 256) {
		MD_LOG_DEBUG("parse_atoms_block: read line: %.*s", (int)MIN(line.len, 50), line.ptr);
		line = str_trim(line);
		if (line.len == 0) continue;
		if (line.ptr[0] == '[') {
			// Hit next block - save it for main parser
			MD_LOG_DEBUG("parse_atoms_block: hit next block marker");
			*next_block_line = line;
			break;
		}
		atom_lines[atom_count++] = line;
	}

	MD_LOG_DEBUG("parse_atoms_block: read %zu atoms", atom_count);

	if (atom_count == 0) {
		MD_LOG_ERROR("No atoms found in [Atoms] block");
		return false;
	}

	// Allocate arrays
	molden->number_of_atoms = atom_count;
	molden->atom_coordinates = md_alloc(alloc, atom_count * sizeof(dvec3_t));
	molden->atomic_numbers = md_alloc(alloc, atom_count * sizeof(md_element_t));

	if (!molden->atom_coordinates || !molden->atomic_numbers) {
		MD_LOG_ERROR("Failed to allocate memory for atoms");
		return false;
	}

	// Parse atoms
	str_t tokens[10];
	for (size_t atom_idx = 0; atom_idx < atom_count; atom_idx++) {
		str_t atom_line = atom_lines[atom_idx];
		
		// Format: Element AtomNum AtomicNum X Y Z
		// Example: O     1    8    0.000000    0.000000    0.119262
		size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &atom_line);
		if (num_tokens < 6) {
			MD_LOG_ERROR("Invalid [Atoms] line: expected at least 6 tokens, got %zu", num_tokens);
			return false;
		}

		// Parse atomic number (third token, index 2)
		int atomic_num = (int)parse_int(tokens[2]);
		molden->atomic_numbers[atom_idx] = (md_element_t)atomic_num;

		// Parse coordinates (tokens 3, 4, 5)
		double x = parse_float(tokens[3]);
		double y = parse_float(tokens[4]);
		double z = parse_float(tokens[5]);

		// Convert to Ångström if needed
		if (molden->coord_unit == MD_MOLDEN_COORD_UNIT_AU) {
			x *= BOHR_TO_ANGSTROM;
			y *= BOHR_TO_ANGSTROM;
			z *= BOHR_TO_ANGSTROM;
		}

		molden->atom_coordinates[atom_idx].x = x;
		molden->atom_coordinates[atom_idx].y = y;
		molden->atom_coordinates[atom_idx].z = z;
	}

	return true;
}

static bool parse_gto_block(md_molden_t* molden, md_buffered_reader_t* reader, str_t* next_block_line) {
	md_allocator_i* alloc = molden->arena;
	str_t line;
	str_t tokens[10];

	// Initialize basis set structures
	basis_set_t* basis_set = &molden->basis_set;
	MEMSET(basis_set, 0, sizeof(basis_set_t));

	// Reserve space for parameters
	size_t max_params = 1024;
	basis_set->param.exponents = md_alloc(alloc, max_params * sizeof(double));
	basis_set->param.coefficients = md_alloc(alloc, max_params * sizeof(double));
	basis_set->param.count = 0;

	// Reserve space for basis functions
	size_t max_funcs = 256;
	basis_set->basis_func.data = md_alloc(alloc, max_funcs * sizeof(basis_set_func_t));
	basis_set->basis_func.count = 0;

	// Reserve space for atom basis (up to element 118)
	size_t max_elements = 119; // 0 is null, 1-118 are elements
	basis_set->atom_basis.data = md_alloc(alloc, max_elements * sizeof(basis_set_basis_t));
	basis_set->atom_basis.count = max_elements;
	MEMSET(basis_set->atom_basis.data, 0, max_elements * sizeof(basis_set_basis_t));

	int current_atom_idx = -1;
	int current_atomic_num = -1;

	while (md_buffered_reader_extract_line(&line, reader)) {
		line = str_trim(line);
		if (line.len == 0) continue;
		if (line.ptr[0] == '[') {
			// Hit next block
			*next_block_line = line; break;
		}

		str_t line_copy = line; size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line_copy);
		if (num_tokens < 1) continue;

		// Check if this is an atom index line (two numbers)
		if (num_tokens == 2) {
			int atom_idx, dummy;
			atom_idx = (int)parse_int(tokens[0]); dummy = (int)parse_int(tokens[1]); if (atom_idx != 0 && dummy != 0) {
				current_atom_idx = atom_idx;
				if (current_atom_idx >= 1 && current_atom_idx <= (int)molden->number_of_atoms) {
					current_atomic_num = molden->atomic_numbers[current_atom_idx - 1];
					
					// Initialize atom basis if not yet done
					if (current_atomic_num >= 0 && current_atomic_num < (int)basis_set->atom_basis.count) {
						basis_set_basis_t* atom_basis = &basis_set->atom_basis.data[current_atomic_num];
						if (atom_basis->basis_func_count == 0) {
							atom_basis->basis_func_offset = (uint16_t)basis_set->basis_func.count;
						}
					}
				}
				continue;
			}
		}

		// Parse basis function line (e.g., "s   3 1.00")
		if (num_tokens >= 2 && current_atomic_num >= 0) {
			int angl_type = char_to_angular_momentum_type(tokens[0].ptr[0]);
			if (angl_type < 0) continue;

			int num_primitives = 0;
			num_primitives = (int)parse_int(tokens[1]); if (num_primitives <= 0) {
				MD_LOG_ERROR("Invalid number of primitives in GTO block");
				return false;
			}

			// Record basis function
			if (basis_set->basis_func.count >= max_funcs) {
				MD_LOG_ERROR("Too many basis functions");
				return false;
			}

			basis_set_func_t* func = &basis_set->basis_func.data[basis_set->basis_func.count];
			func->type = (uint8_t)angl_type;
			func->param_count = (uint8_t)num_primitives;
			func->param_offset = (uint16_t)basis_set->param.count;
			basis_set->basis_func.count++;

			// Update atom basis
			if (current_atomic_num < (int)basis_set->atom_basis.count) {
				basis_set_basis_t* atom_basis = &basis_set->atom_basis.data[current_atomic_num];
				atom_basis->basis_func_count++;
				atom_basis->max_type = MAX(atom_basis->max_type, (uint8_t)angl_type);
			}

			// Parse primitives (exponent and coefficient pairs)
			for (int i = 0; i < num_primitives; i++) {
				if (!md_buffered_reader_extract_line(&line, reader)) {
					MD_LOG_ERROR("Unexpected end of file while reading GTO primitives");
					return false;
				}
				line = str_trim(line);
				if (line.len == 0) {
					i--;
					continue;
				}

				str_t line_copy2 = line;
				size_t prim_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line_copy2);
				if (prim_tokens < 2) {
					MD_LOG_ERROR("Invalid primitive line in GTO block");
					return false;
				}

				if (basis_set->param.count >= max_params) {
					MD_LOG_ERROR("Too many GTO parameters");
					return false;
				}

				double exponent = parse_float(tokens[0]);
				double coefficient = parse_float(tokens[1]);

				basis_set->param.exponents[basis_set->param.count] = exponent;
				basis_set->param.coefficients[basis_set->param.count] = coefficient;
				basis_set->param.count++;
			}
		}
	}

	// Normalize basis set
	normalize_basis_set(basis_set);

	return true;
}

static bool parse_mo_block(md_molden_t* molden, md_buffered_reader_t* reader, str_t* next_block_line) {
md_allocator_i* alloc = molden->arena;
str_t line;
str_t tokens[10];

// Read all MO lines into temporary storage
str_t mo_lines[2048];  // Max 2048 lines for MO block
size_t line_count = 0;

while (md_buffered_reader_extract_line(&line, reader) && line_count < 2048) {
line = str_trim(line);
if (line.len == 0) continue;
if (line.ptr[0] == '[') break;
mo_lines[line_count++] = line;
}

// First pass: count MOs and AOs
size_t num_mos = 0;
size_t num_aos = 0;

for (size_t i = 0; i < line_count; i++) {
line = mo_lines[i];

if (str_find_str(NULL, line, STR_LIT("Sym="))) {
num_mos++;
} else if (num_mos == 1 && num_aos == 0) {
// Count AO coefficients in first MO
str_t line_copy = line; 
size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line_copy);
if (num_tokens >= 2) {
int idx = (int)parse_int(tokens[0]);
num_aos = MAX(num_aos, (size_t)idx);
}
}
}

if (num_mos == 0 || num_aos == 0) {
MD_LOG_ERROR("No molecular orbitals found in [MO] block");
return false;
}

// Allocate MO data
molden->scf.alpha.coefficients.size[0] = num_aos;
molden->scf.alpha.coefficients.size[1] = num_mos;
molden->scf.alpha.coefficients.data = md_alloc(alloc, num_aos * num_mos * sizeof(double));
molden->scf.alpha.energy.size = num_mos;
molden->scf.alpha.energy.data = md_alloc(alloc, num_mos * sizeof(double));
molden->scf.alpha.occupancy.size = num_mos;
molden->scf.alpha.occupancy.data = md_alloc(alloc, num_mos * sizeof(double));
molden->scf.alpha.symmetry = md_alloc(alloc, num_mos * sizeof(str_t));
molden->scf.alpha.spin = md_alloc(alloc, num_mos * sizeof(str_t));

if (!molden->scf.alpha.coefficients.data || !molden->scf.alpha.energy.data || 
    !molden->scf.alpha.occupancy.data || !molden->scf.alpha.symmetry || !molden->scf.alpha.spin) {
MD_LOG_ERROR("Failed to allocate memory for MO data");
return false;
}

MEMSET(molden->scf.alpha.coefficients.data, 0, num_aos * num_mos * sizeof(double));
MEMSET(molden->scf.alpha.energy.data, 0, num_mos * sizeof(double));
MEMSET(molden->scf.alpha.occupancy.data, 0, num_mos * sizeof(double));
MEMSET(molden->scf.alpha.symmetry, 0, num_mos * sizeof(str_t));
MEMSET(molden->scf.alpha.spin, 0, num_mos * sizeof(str_t));

// Second pass: parse MOs
size_t mo_idx = 0;
str_t current_sym = {0};
double current_ene = 0.0;
str_t current_spin = {0};
double current_occup = 0.0;

for (size_t i = 0; i < line_count && mo_idx <= num_mos; i++) {
line = mo_lines[i];

// Parse MO properties
if (str_find_str(NULL, line, STR_LIT("Sym="))) {
// Start of new MO
if (mo_idx > 0) {
// Finalize previous MO
molden->scf.alpha.symmetry[mo_idx - 1] = current_sym;
molden->scf.alpha.spin[mo_idx - 1] = current_spin;
molden->scf.alpha.energy.data[mo_idx - 1] = current_ene;
molden->scf.alpha.occupancy.data[mo_idx - 1] = current_occup;
}
mo_idx++;

// Parse symmetry label
size_t eq_pos;
if (str_find_char(&eq_pos, line, '=')) {
current_sym = str_trim(str_substr(line, eq_pos + 1, line.len));
current_sym = str_copy(current_sym, alloc);
}

} else if (str_find_str(NULL, line, STR_LIT("Ene="))) {
size_t eq_pos;
if (str_find_char(&eq_pos, line, '=')) {
str_t value = str_trim(str_substr(line, eq_pos + 1, line.len));
current_ene = parse_float(value);
}

} else if (str_find_str(NULL, line, STR_LIT("Spin="))) {
size_t eq_pos;
if (str_find_char(&eq_pos, line, '=')) {
current_spin = str_trim(str_substr(line, eq_pos + 1, line.len));
current_spin = str_copy(current_spin, alloc);
}

} else if (str_find_str(NULL, line, STR_LIT("Occup="))) {
size_t eq_pos;
if (str_find_char(&eq_pos, line, '=')) {
str_t value = str_trim(str_substr(line, eq_pos + 1, line.len));
current_occup = parse_float(value);
}

} else {
// Parse coefficient
str_t line_copy = line; 
size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line_copy);
if (num_tokens >= 2 && mo_idx > 0) {
int idx = (int)parse_int(tokens[0]);
double coeff = parse_float(tokens[1]);
if (idx >= 1 && idx <= (int)num_aos && mo_idx - 1 < num_mos) {
molden->scf.alpha.coefficients.data[(idx - 1) * num_mos + (mo_idx - 1)] = coeff;
}
}
}
}

// Finalize last MO
if (mo_idx > 0 && mo_idx <= num_mos) {
molden->scf.alpha.symmetry[mo_idx - 1] = current_sym;
molden->scf.alpha.spin[mo_idx - 1] = current_spin;
molden->scf.alpha.energy.data[mo_idx - 1] = current_ene;
molden->scf.alpha.occupancy.data[mo_idx - 1] = current_occup;
}

// Calculate HOMO/LUMO indices and electron count
molden->scf.type = MD_MOLDEN_SCF_TYPE_RESTRICTED;
molden->scf.homo_idx[0] = 0;
molden->scf.lumo_idx[0] = 0;

size_t num_alpha_electrons = 0;
for (size_t i = 0; i < num_mos; i++) {
if (molden->scf.alpha.occupancy.data[i] > 0.5) {
molden->scf.homo_idx[0] = i;
num_alpha_electrons += (size_t)molden->scf.alpha.occupancy.data[i];
} else if (molden->scf.lumo_idx[0] == 0 || molden->scf.lumo_idx[0] == molden->scf.homo_idx[0]) {
molden->scf.lumo_idx[0] = i;
}
}

molden->number_of_alpha_electrons = num_alpha_electrons;
molden->number_of_beta_electrons = num_alpha_electrons;  // Restricted calc
molden->spin_multiplicity = 1;

return true;
}

static bool parse_molden_file(md_molden_t* molden, str_t filename) {
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (!file) {
		MD_LOG_ERROR("Failed to open Molden file: %.*s", (int)filename.len, filename.ptr);
		return false;
	}

	MD_LOG_DEBUG("Opened Molden file successfully");

	char* buf = md_alloc(md_get_temp_allocator(), MEGABYTES(4));
	if (!buf) {
		MD_LOG_ERROR("Failed to allocate buffer for reading Molden file");
		md_file_close(file);
		return false;
	}

	md_buffered_reader_t reader = md_buffered_reader_from_file(buf, MEGABYTES(4), file);
	str_t line;

	bool success = true;
	bool found_atoms = false;
	bool found_gto = false;
	bool found_mo = false;
	str_t next_block = {0};

	MD_LOG_DEBUG("Starting to read lines from file");
	
	int line_count = 0;
	while (md_buffered_reader_extract_line(&line, &reader)) {
process_line:
		line_count++;
		if (line_count < 10) {
			MD_LOG_DEBUG("Read line %d: %.*s", line_count, (int)MIN(line.len, 50), line.ptr);
		}
		line = str_trim(line);
		if (line.len == 0) continue;

		// Check for block headers
		if (str_eq_ignore_case(line, STR_LIT("[5D]"))) {
			molden->use_5d = true;
			MD_LOG_DEBUG("Found [5D] marker");
		} else if (str_eq_ignore_case(line, STR_LIT("[7F]"))) {
			molden->use_7f = true;
			MD_LOG_DEBUG("Found [7F] marker");
		} else if (str_eq_ignore_case(line, STR_LIT("[9G]"))) {
			molden->use_9g = true;
			MD_LOG_DEBUG("Found [9G] marker");
		} else if (str_begins_with(line, STR_LIT("[Atoms]")) || str_begins_with(line, STR_LIT("[atoms]")) || str_begins_with(line, STR_LIT("[ATOMS]"))) {
			MD_LOG_DEBUG("Parsing [Atoms] block");
			// Determine coordinate unit from the [Atoms] line
			if (str_find_str(NULL, line, STR_LIT("AU"))) {
				molden->coord_unit = MD_MOLDEN_COORD_UNIT_AU;
			} else {
				molden->coord_unit = MD_MOLDEN_COORD_UNIT_ANGS;
			}
			
			if (!parse_atoms_block(molden, &reader, &next_block)) {
				MD_LOG_ERROR("Failed to parse [Atoms] block");
				success = false;
				break;
			}
			MD_LOG_DEBUG("Successfully parsed [Atoms] block with %zu atoms", molden->number_of_atoms);
			found_atoms = true;
			if (next_block.len > 0) { line = next_block; next_block.len = 0; goto process_line; }
		} else if (str_eq_ignore_case(line, STR_LIT("[GTO]"))) {
			MD_LOG_DEBUG("Parsing [GTO] block");
			if (!parse_gto_block(molden, &reader, &next_block)) {
				MD_LOG_ERROR("Failed to parse [GTO] block");
				success = false;
				break;
			}
			MD_LOG_DEBUG("Successfully parsed [GTO] block");
			found_gto = true;
			if (next_block.len > 0) { line = next_block; next_block.len = 0; goto process_line; }
		} else if (str_eq_ignore_case(line, STR_LIT("[MO]"))) {
			MD_LOG_DEBUG("Parsing [MO] block");
			if (!parse_mo_block(molden, &reader, &next_block)) {
				MD_LOG_ERROR("Failed to parse [MO] block");
				success = false;
				break;
			}
			MD_LOG_DEBUG("Successfully parsed [MO] block");
			found_mo = true;
		}
	}

	md_file_close(file);

	MD_LOG_DEBUG("Parsing complete. Found: atoms=%d, gto=%d, mo=%d", found_atoms, found_gto, found_mo);

	if (success && found_atoms && found_gto) {
		// Extract GTO data and AO-to-atom mapping
		molden->ao_to_atom_idx = md_alloc(molden->arena, molden->scf.alpha.coefficients.size[0] * sizeof(int));
		if (molden->ao_to_atom_idx) {
			extract_ao_to_atom_idx(molden->ao_to_atom_idx, molden->atomic_numbers, molden->number_of_atoms, &molden->basis_set);
		}

		// Extract GTO data for visualization
		extract_gto_data(&molden->gto_data, molden->atom_coordinates, molden->atomic_numbers, molden->number_of_atoms, &molden->basis_set, molden->arena);
	}

	return success && found_atoms && found_gto;
}
