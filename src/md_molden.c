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
#include <string.h>
#include <ctype.h>

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.5291772109029999
#define HARTREE_TO_EV 27.2114079527

// Molden file sections
typedef enum {
    MOLDEN_SECTION_NONE = 0,
    MOLDEN_SECTION_ATOMS,
    MOLDEN_SECTION_GTO,
    MOLDEN_SECTION_MO,
    MOLDEN_SECTION_FREQ,
    MOLDEN_SECTION_GEOMETRIES,
} molden_section_t;

// Basis function data
typedef struct {
    uint8_t angular_momentum; // s=0, p=1, d=2, f=3, etc.
    size_t num_primitives;
    size_t primitive_offset; // Offset into exponents/coefficients arrays
} molden_basis_func_t;

// Basis set per atom
typedef struct {
    size_t atom_idx;
    size_t num_funcs;
    size_t func_offset; // Offset into basis_func array
} molden_atom_basis_t;

// Molecular orbital data
typedef struct {
    double energy;
    double occupancy;
    size_t coefficient_offset; // Offset into mo_coefficients array
} molden_mo_t;

// Main Molden data structure
struct md_molden_t {
    struct md_allocator_i* alloc;
    
    // Molecule structure
    struct {
        size_t count;
        dvec3_t* coord;    // Coordinates (in Angstroms)
        uint8_t* atomic_number;
    } atoms;
    
    // Basis set information
    struct {
        size_t num_primitives;
        double* exponents;
        double* coefficients;
        
        size_t num_funcs;
        molden_basis_func_t* funcs;
        
        size_t num_atom_basis;
        molden_atom_basis_t* atom_basis;
    } basis;
    
    // Molecular orbitals
    struct {
        size_t num_aos; // Number of atomic orbitals
        size_t num_mos; // Number of molecular orbitals
        
        molden_mo_t* alpha_mos;
        molden_mo_t* beta_mos;
        
        double* mo_coefficients; // All MO coefficients stored sequentially
        
        int* ao_to_atom; // Maps AO index to atom index
    } orbitals;
    
    // Molecular properties
    struct {
        double charge;
        size_t multiplicity;
        size_t num_alpha_electrons;
        size_t num_beta_electrons;
        md_molden_scf_type_t scf_type;
    } properties;
    
    // Flags for available data
    struct {
        bool has_basis : 1;
        bool has_mos : 1;
        bool has_beta_mos : 1;
        bool units_angs : 1; // true = Angstroms, false = Bohr
    } flags;
};

// Helper function to skip whitespace
static inline str_t skip_whitespace(str_t str) {
    while (str.len > 0 && isspace(str.ptr[0])) {
        str.ptr++;
        str.len--;
    }
    return str;
}

// Helper function to read a line
static inline str_t extract_line(str_t* str) {
    str_t line = {0};
    if (str->len == 0) return line;
    
    const char* end = str->ptr;
    const char* max = str->ptr + str->len;
    
    while (end < max && *end != '\n' && *end != '\r') {
        end++;
    }
    
    line.ptr = str->ptr;
    line.len = end - str->ptr;
    
    // Skip line ending
    if (end < max && *end == '\r') end++;
    if (end < max && *end == '\n') end++;
    
    str->ptr = end;
    str->len = max - end;
    
    return line;
}

// Helper to check if line starts with section header
static molden_section_t parse_section_header(str_t line) {
    line = skip_whitespace(line);
    if (line.len == 0 || line.ptr[0] != '[') return MOLDEN_SECTION_NONE;
    
    if (str_eq_n_ignore_case(line, STR_LIT("[Atoms]"), 7)) {
        return MOLDEN_SECTION_ATOMS;
    } else if (str_eq_n_ignore_case(line, STR_LIT("[GTO]"), 5)) {
        return MOLDEN_SECTION_GTO;
    } else if (str_eq_n_ignore_case(line, STR_LIT("[MO]"), 4)) {
        return MOLDEN_SECTION_MO;
    } else if (str_eq_n_ignore_case(line, STR_LIT("[FREQ]"), 6)) {
        return MOLDEN_SECTION_FREQ;
    } else if (str_eq_n_ignore_case(line, STR_LIT("[GEOMETRIES]"), 12)) {
        return MOLDEN_SECTION_GEOMETRIES;
    }
    
    return MOLDEN_SECTION_NONE;
}

// Parse [Atoms] section
static bool parse_atoms_section(md_molden_t* molden, str_t* content) {
    // First pass: count atoms
    str_t temp = *content;
    size_t atom_count = 0;
    
    str_t line = extract_line(&temp);
    bool units_angs = str_eq_n_ignore_case(line, STR_LIT("Angs"), 4) || 
                      str_eq_n_ignore_case(line, STR_LIT("(Angs)"), 6);
    molden->flags.units_angs = units_angs;
    
    while (temp.len > 0) {
        line = extract_line(&temp);
        line = skip_whitespace(line);
        if (line.len == 0) continue;
        if (line.ptr[0] == '[') break; // Next section
        atom_count++;
    }
    
    if (atom_count == 0) return false;
    
    // Allocate arrays
    molden->atoms.count = atom_count;
    molden->atoms.coord = md_alloc(molden->alloc, sizeof(dvec3_t) * atom_count);
    molden->atoms.atomic_number = md_alloc(molden->alloc, sizeof(uint8_t) * atom_count);
    
    if (!molden->atoms.coord || !molden->atoms.atomic_number) {
        MD_LOG_ERROR("Failed to allocate memory for atoms");
        return false;
    }
    
    // Second pass: parse atom data
    extract_line(content); // Skip header line
    
    for (size_t i = 0; i < atom_count; i++) {
        line = extract_line(content);
        
        char element[4] = {0};
        int dummy_idx, atomic_num;
        double x, y, z;
        
        // Format: Element Index AtomicNumber X Y Z
        if (sscanf(line.ptr, "%3s %d %d %lf %lf %lf", 
                   element, &dummy_idx, &atomic_num, &x, &y, &z) != 6) {
            MD_LOG_ERROR("Failed to parse atom line: %.*s", (int)line.len, line.ptr);
            return false;
        }
        
        molden->atoms.atomic_number[i] = (uint8_t)atomic_num;
        
        // Convert Bohr to Angstrom if needed
        if (!units_angs) {
            x *= BOHR_TO_ANGSTROM;
            y *= BOHR_TO_ANGSTROM;
            z *= BOHR_TO_ANGSTROM;
        }
        
        molden->atoms.coord[i].x = x;
        molden->atoms.coord[i].y = y;
        molden->atoms.coord[i].z = z;
    }
    
    return true;
}

// Parse [GTO] section
static bool parse_gto_section(md_molden_t* molden, str_t* content) {
    // Implementation for parsing Gaussian Type Orbitals
    // This is complex and would need multiple passes
    
    // For now, mark that we have basis set data
    molden->flags.has_basis = true;
    
    // Skip to next section for now
    while (content->len > 0) {
        str_t line = extract_line(content);
        line = skip_whitespace(line);
        if (line.len > 0 && line.ptr[0] == '[') {
            // Put line back
            content->ptr = line.ptr;
            content->len += line.len;
            break;
        }
    }
    
    return true;
}

// Parse [MO] section
static bool parse_mo_section(md_molden_t* molden, str_t* content) {
    // Implementation for parsing Molecular Orbitals
    // This section contains MO coefficients, energies, occupancies
    
    molden->flags.has_mos = true;
    
    // Skip to next section for now
    while (content->len > 0) {
        str_t line = extract_line(content);
        line = skip_whitespace(line);
        if (line.len > 0 && line.ptr[0] == '[') {
            // Put line back
            content->ptr = line.ptr;
            content->len += line.len;
            break;
        }
    }
    
    return true;
}

// Main parsing function
static bool molden_parse_str(md_molden_t* molden, str_t content) {
    if (!molden || content.len == 0) return false;
    
    // Initialize default properties
    molden->properties.charge = 0.0;
    molden->properties.multiplicity = 1;
    molden->properties.scf_type = MD_MOLDEN_SCF_TYPE_RESTRICTED;
    
    molden_section_t current_section = MOLDEN_SECTION_NONE;
    
    while (content.len > 0) {
        str_t line = extract_line(&content);
        line = skip_whitespace(line);
        
        if (line.len == 0) continue;
        
        // Check for section headers
        molden_section_t section = parse_section_header(line);
        if (section != MOLDEN_SECTION_NONE) {
            current_section = section;
            
            // Parse sections
            switch (current_section) {
                case MOLDEN_SECTION_ATOMS:
                    if (!parse_atoms_section(molden, &content)) {
                        MD_LOG_ERROR("Failed to parse [Atoms] section");
                        return false;
                    }
                    break;
                    
                case MOLDEN_SECTION_GTO:
                    if (!parse_gto_section(molden, &content)) {
                        MD_LOG_ERROR("Failed to parse [GTO] section");
                        return false;
                    }
                    break;
                    
                case MOLDEN_SECTION_MO:
                    if (!parse_mo_section(molden, &content)) {
                        MD_LOG_ERROR("Failed to parse [MO] section");
                        return false;
                    }
                    break;
                    
                default:
                    // Skip unknown sections
                    break;
            }
        }
    }
    
    return molden->atoms.count > 0;
}

// Public API Implementation

struct md_molden_t* md_molden_create(struct md_allocator_i* backing) {
    if (!backing) {
        MD_LOG_ERROR("Allocator is NULL");
        return NULL;
    }
    
    md_molden_t* molden = md_alloc(backing, sizeof(md_molden_t));
    if (!molden) {
        MD_LOG_ERROR("Failed to allocate md_molden_t");
        return NULL;
    }
    
    memset(molden, 0, sizeof(md_molden_t));
    molden->alloc = backing;
    
    return molden;
}

void md_molden_reset(struct md_molden_t* molden) {
    if (!molden) return;
    
    // Free all allocated data
    if (molden->atoms.coord) md_free(molden->alloc, molden->atoms.coord, molden->atoms.count * sizeof(dvec3_t));
    if (molden->atoms.atomic_number) md_free(molden->alloc, molden->atoms.atomic_number, molden->atoms.count * sizeof(uint8_t));
    
    if (molden->basis.exponents) md_free(molden->alloc, molden->basis.exponents, molden->basis.num_primitives * sizeof(double));
    if (molden->basis.coefficients) md_free(molden->alloc, molden->basis.coefficients, molden->basis.num_primitives * sizeof(double));
    if (molden->basis.funcs) md_free(molden->alloc, molden->basis.funcs, molden->basis.num_funcs * sizeof(molden_basis_func_t));
    if (molden->basis.atom_basis) md_free(molden->alloc, molden->basis.atom_basis, molden->basis.num_atom_basis * sizeof(molden_atom_basis_t));
    
    if (molden->orbitals.alpha_mos) md_free(molden->alloc, molden->orbitals.alpha_mos, molden->orbitals.num_mos * sizeof(molden_mo_t));
    if (molden->orbitals.beta_mos) md_free(molden->alloc, molden->orbitals.beta_mos, molden->orbitals.num_mos * sizeof(molden_mo_t));
    if (molden->orbitals.mo_coefficients) md_free(molden->alloc, molden->orbitals.mo_coefficients, molden->orbitals.num_mos * molden->orbitals.num_aos * sizeof(double));
    if (molden->orbitals.ao_to_atom) md_free(molden->alloc, molden->orbitals.ao_to_atom, molden->orbitals.num_aos * sizeof(int));
    
    // Reset structure
    memset(molden, 0, sizeof(md_molden_t));
    struct md_allocator_i* saved_alloc = molden->alloc;
    memset(molden, 0, sizeof(md_molden_t));
    molden->alloc = saved_alloc;
}

void md_molden_destroy(struct md_molden_t* molden) {
    if (!molden) return;
    
    struct md_allocator_i* alloc = molden->alloc;
    md_molden_reset(molden);
    md_free(alloc, molden, sizeof(md_molden_t));
}

bool md_molden_parse_file(struct md_molden_t* molden, str_t filename) {
    if (!molden) {
        MD_LOG_ERROR("Molden object is NULL");
        return false;
    }
    
    struct md_allocator_i* saved_alloc = molden->alloc;
    md_molden_reset(molden);
    molden->alloc = saved_alloc;
    
    bool result = false;
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    if (!file) {
        MD_LOG_ERROR("Failed to open file: %.*s", (int)filename.len, filename.ptr);
        return false;
    }
    
    // Get file size
    md_file_seek(file, 0, MD_FILE_END);
    int64_t file_size = md_file_tell(file);
    md_file_seek(file, 0, MD_FILE_BEG);
    
    if (file_size <= 0 || file_size > GIGABYTES(1)) {
        MD_LOG_ERROR("Invalid file size: %lld", (long long)file_size);
        md_file_close(file);
        return false;
    }
    
    // Allocate buffer and read file
    char* content = md_alloc(molden->alloc, (size_t)file_size + 1);
    if (!content) {
        MD_LOG_ERROR("Failed to allocate memory for file content");
        md_file_close(file);
        return false;
    }
    
    size_t bytes_read = md_file_read(file, content, (size_t)file_size);
    content[bytes_read] = '\0'; // Null terminate
    md_file_close(file);
    
    if (bytes_read != (size_t)file_size) {
        MD_LOG_ERROR("Failed to read complete file");
        md_free(molden->alloc, content, (size_t)file_size + 1);
        return false;
    }
    
    // Parse content
    str_t str = {content, bytes_read};
    result = molden_parse_str(molden, str);
    
    md_free(molden->alloc, content, (size_t)file_size + 1);
    
    return result;
}

// Molecule structure queries
size_t md_molden_number_of_atoms(const struct md_molden_t* molden) {
    return molden ? molden->atoms.count : 0;
}

size_t md_molden_number_of_alpha_electrons(const struct md_molden_t* molden) {
    return molden ? molden->properties.num_alpha_electrons : 0;
}

size_t md_molden_number_of_beta_electrons(const struct md_molden_t* molden) {
    return molden ? molden->properties.num_beta_electrons : 0;
}

double md_molden_molecular_charge(const struct md_molden_t* molden) {
    return molden ? molden->properties.charge : 0.0;
}

size_t md_molden_spin_multiplicity(const struct md_molden_t* molden) {
    return molden ? molden->properties.multiplicity : 1;
}

const dvec3_t* md_molden_atom_coordinates(const struct md_molden_t* molden) {
    return molden ? molden->atoms.coord : NULL;
}

const uint8_t* md_molden_atomic_numbers(const struct md_molden_t* molden) {
    return molden ? molden->atoms.atomic_number : NULL;
}

const int* md_molden_ao_to_atom_idx(const struct md_molden_t* molden) {
    return molden ? molden->orbitals.ao_to_atom : NULL;
}

// SCF and orbital information
md_molden_scf_type_t md_molden_scf_type(const struct md_molden_t* molden) {
    return molden ? molden->properties.scf_type : MD_MOLDEN_SCF_TYPE_UNKNOWN;
}

size_t md_molden_scf_homo_idx(const struct md_molden_t* molden, md_molden_mo_type_t type) {
    if (!molden) return 0;
    
    size_t num_electrons = (type == MD_MOLDEN_MO_TYPE_ALPHA) ? 
        molden->properties.num_alpha_electrons : 
        molden->properties.num_beta_electrons;
    
    // HOMO is the highest occupied orbital (0-indexed)
    return num_electrons > 0 ? num_electrons - 1 : 0;
}

size_t md_molden_scf_lumo_idx(const struct md_molden_t* molden, md_molden_mo_type_t type) {
    if (!molden) return 0;
    
    size_t num_electrons = (type == MD_MOLDEN_MO_TYPE_ALPHA) ? 
        molden->properties.num_alpha_electrons : 
        molden->properties.num_beta_electrons;
    
    // LUMO is the lowest unoccupied orbital
    return num_electrons;
}

size_t md_molden_scf_number_of_atomic_orbitals(const struct md_molden_t* molden) {
    return molden ? molden->orbitals.num_aos : 0;
}

size_t md_molden_scf_number_of_molecular_orbitals(const struct md_molden_t* molden) {
    return molden ? molden->orbitals.num_mos : 0;
}

const double* md_molden_scf_mo_occupancy(const struct md_molden_t* molden, md_molden_mo_type_t type) {
    if (!molden || !molden->flags.has_mos) return NULL;
    
    // Not implemented yet - would need to extract from MO data
    return NULL;
}

const double* md_molden_scf_mo_energy(const struct md_molden_t* molden, md_molden_mo_type_t type) {
    if (!molden || !molden->flags.has_mos) return NULL;
    
    // Not implemented yet - would need to extract from MO data
    return NULL;
}

// GTO extraction functions
bool md_molden_scf_extract_gto_data(md_gto_data_t* out_gto_data, const struct md_molden_t* molden, double cutoff_value, struct md_allocator_i* alloc) {
    if (!out_gto_data || !molden || !alloc) return false;
    
    // Not implemented yet - complex GTO extraction
    return false;
}

size_t md_molden_mo_gto_count(const struct md_molden_t* molden) {
    if (!molden || !molden->flags.has_basis) return 0;
    
    // Not implemented yet
    return 0;
}

size_t md_molden_mo_gto_extract(md_gto_t gtos[], const struct md_molden_t* molden, size_t mo_idx, md_molden_mo_type_t type, double value_cutoff) {
    if (!gtos || !molden || !molden->flags.has_mos) return 0;
    
    // Not implemented yet - complex GTO extraction
    return 0;
}

// Molecule system integration
bool md_molden_system_init(struct md_system_t* sys, const md_molden_t* molden, struct md_allocator_i* alloc) {
    if (!sys || !molden || !alloc) return false;
    
    if (molden->atoms.count == 0) return false;
    
    // Initialize basic molecule structure
    memset(sys, 0, sizeof(md_system_t));
    
    sys->atom.count = molden->atoms.count;
    sys->atom.x = md_alloc(alloc, sizeof(float) * molden->atoms.count);
    sys->atom.y = md_alloc(alloc, sizeof(float) * molden->atoms.count);
    sys->atom.z = md_alloc(alloc, sizeof(float) * molden->atoms.count);
    sys->atom.type_idx = md_alloc(alloc, sizeof(md_atom_type_idx_t) * molden->atoms.count);
    
    if (!sys->atom.x || !sys->atom.y || !sys->atom.z || !sys->atom.type_idx) {
        MD_LOG_ERROR("Failed to allocate memory for molecule system");
        return false;
    }
    
    // Setup atom types
    sys->atom.type.count = molden->atoms.count;
    sys->atom.type.z = md_alloc(alloc, sizeof(md_atomic_number_t) * molden->atoms.count);
    
    if (!sys->atom.type.z) {
        MD_LOG_ERROR("Failed to allocate memory for atom types");
        return false;
    }
    
    // Copy atom data
    for (size_t i = 0; i < molden->atoms.count; i++) {
        sys->atom.x[i] = (float)molden->atoms.coord[i].x;
        sys->atom.y[i] = (float)molden->atoms.coord[i].y;
        sys->atom.z[i] = (float)molden->atoms.coord[i].z;
        sys->atom.type_idx[i] = (md_atom_type_idx_t)i;
        sys->atom.type.z[i] = molden->atoms.atomic_number[i];
    }
    
    return true;
}

// System loader implementation
static bool molden_loader_init_from_str(struct md_system_t* sys, str_t string, const void* arg, struct md_allocator_i* alloc) {
    (void)arg;
    md_molden_t* molden = md_molden_create(alloc);
    if (!molden) return false;
    
    bool result = false;
    if (molden_parse_str(molden, string)) {
        result = md_molden_system_init(sys, molden, alloc);
    }
    
    md_molden_destroy(molden);
    return result;
}

static bool molden_loader_init_from_file(struct md_system_t* sys, str_t filename, const void* arg, struct md_allocator_i* alloc) {
    (void)arg;
    md_molden_t* molden = md_molden_create(alloc);
    if (!molden) return false;
    
    bool result = false;
    if (md_molden_parse_file(molden, filename)) {
        result = md_molden_system_init(sys, molden, alloc);
    }
    
    md_molden_destroy(molden);
    return result;
}

static struct md_system_loader_i molden_loader = {
    molden_loader_init_from_str,
    molden_loader_init_from_file,
};

struct md_system_loader_i* md_molden_system_loader(void) {
    return &molden_loader;
}
