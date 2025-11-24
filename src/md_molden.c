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

	// TODO: Implement actual parsing
	MD_LOG_ERROR("Molden file parsing not yet implemented");
	(void)filename;
	return false;
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

	// TODO: Implement GTO extraction
	MD_LOG_ERROR("GTO data extraction not yet implemented");
	(void)cutoff_value;
	return false;
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

	// TODO: Implement MO GTO extraction
	MD_LOG_ERROR("MO GTO extraction not yet implemented");
	(void)mo_idx;
	(void)type;
	(void)value_cutoff;
	return 0;
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
