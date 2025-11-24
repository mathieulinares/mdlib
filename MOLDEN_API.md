# Molden API Documentation

## Overview

The Molden API (`md_molden.h` and `md_molden.c`) provides functionality to parse and work with Molden file format data in MDLIB. This enables VIAMD to load and visualize molecules with quantum chemistry data including molecular orbitals, basis sets, and electronic structure information.

## Molden File Format

Molden is a widely-used file format in quantum chemistry for storing:
- Molecular structure (atom positions, elements)
- Basis set definitions (Gaussian Type Orbitals)
- Molecular orbitals and their properties
- SCF (Self-Consistent Field) calculation results
- Vibrational frequencies and normal modes
- Geometry optimization trajectories

## API Functions

### Object Lifecycle

#### `md_molden_t* md_molden_create(md_allocator_i* backing)`
Creates a new Molden parser object.
- **Parameters**: `backing` - Memory allocator to use
- **Returns**: Pointer to new md_molden_t object, or NULL on failure

#### `void md_molden_reset(md_molden_t* molden)`
Resets the parser state, freeing all parsed data.
- **Parameters**: `molden` - Molden object to reset

#### `void md_molden_destroy(md_molden_t* molden)`
Destroys the Molden object and frees all memory.
- **Parameters**: `molden` - Molden object to destroy

### File Parsing

#### `bool md_molden_parse_file(md_molden_t* molden, str_t filename)`
Parses a Molden format file.
- **Parameters**: 
  - `molden` - Molden parser object
  - `filename` - Path to Molden file
- **Returns**: true on success, false on failure
- **Supported sections**: 
  - `[Atoms]` - Molecular structure (Angs or Bohr units)
  - `[GTO]` - Basis set (partial support)
  - `[MO]` - Molecular orbitals (partial support)

### Molecule Structure Queries

#### `size_t md_molden_number_of_atoms(const md_molden_t* molden)`
Returns the number of atoms in the molecule.

#### `const dvec3_t* md_molden_atom_coordinates(const md_molden_t* molden)`
Returns array of atomic coordinates (in Angstroms).
- **Returns**: Array of length `md_molden_number_of_atoms()`

#### `const uint8_t* md_molden_atomic_numbers(const md_molden_t* molden)`
Returns array of atomic numbers.
- **Returns**: Array of length `md_molden_number_of_atoms()`

#### `double md_molden_molecular_charge(const md_molden_t* molden)`
Returns the molecular charge.

#### `size_t md_molden_spin_multiplicity(const md_molden_t* molden)`
Returns the spin multiplicity (2S+1).

### Electron Configuration

#### `size_t md_molden_number_of_alpha_electrons(const md_molden_t* molden)`
Returns the number of alpha electrons.

#### `size_t md_molden_number_of_beta_electrons(const md_molden_t* molden)`
Returns the number of beta electrons.

### SCF and Molecular Orbitals

#### `md_molden_scf_type_t md_molden_scf_type(const md_molden_t* molden)`
Returns the SCF calculation type.
- **Returns**: 
  - `MD_MOLDEN_SCF_TYPE_RESTRICTED` - RHF/RKS
  - `MD_MOLDEN_SCF_TYPE_UNRESTRICTED` - UHF/UKS
  - `MD_MOLDEN_SCF_TYPE_UNKNOWN`

#### `size_t md_molden_scf_homo_idx(const md_molden_t* molden, md_molden_mo_type_t type)`
Returns the index of the HOMO (Highest Occupied Molecular Orbital).
- **Parameters**: `type` - `MD_MOLDEN_MO_TYPE_ALPHA` or `MD_MOLDEN_MO_TYPE_BETA`

#### `size_t md_molden_scf_lumo_idx(const md_molden_t* molden, md_molden_mo_type_t type)`
Returns the index of the LUMO (Lowest Unoccupied Molecular Orbital).

#### `size_t md_molden_scf_number_of_atomic_orbitals(const md_molden_t* molden)`
Returns the number of atomic orbitals (basis functions).

#### `size_t md_molden_scf_number_of_molecular_orbitals(const md_molden_t* molden)`
Returns the number of molecular orbitals.

#### `const int* md_molden_ao_to_atom_idx(const md_molden_t* molden)`
Returns mapping from atomic orbital index to atom index.

### GTO (Gaussian Type Orbital) Extraction

#### `bool md_molden_scf_extract_gto_data(md_gto_data_t* out_gto_data, const md_molden_t* molden, double cutoff_value, md_allocator_i* alloc)`
Extracts GTO data for visualization.
- **Status**: API defined, implementation pending

#### `size_t md_molden_mo_gto_count(const md_molden_t* molden)`
Returns the number of GTOs for molecular orbital visualization.
- **Status**: API defined, implementation pending

#### `size_t md_molden_mo_gto_extract(md_gto_t gtos[], const md_molden_t* molden, size_t mo_idx, md_molden_mo_type_t type, double value_cutoff)`
Extracts GTOs for a specific molecular orbital.
- **Status**: API defined, implementation pending

### System Integration

#### `bool md_molden_system_init(md_system_t* sys, const md_molden_t* molden, md_allocator_i* alloc)`
Initializes an md_system_t from Molden data.
- **Parameters**:
  - `sys` - System structure to initialize
  - `molden` - Parsed Molden data
  - `alloc` - Allocator for system data
- **Returns**: true on success

#### `md_system_loader_i* md_molden_system_loader(void)`
Returns the system loader interface for Molden files.
- **Supported extensions**: `.molden`, `.input`, `.mold`

## Usage Example

```c
#include <md_molden.h>
#include <md_system.h>
#include <core/md_allocator.h>

// Create parser
md_allocator_i* alloc = md_get_heap_allocator();
md_molden_t* molden = md_molden_create(alloc);

// Parse file
if (md_molden_parse_file(molden, STR_LIT("molecule.molden"))) {
    // Query atom data
    size_t num_atoms = md_molden_number_of_atoms(molden);
    const dvec3_t* coords = md_molden_atom_coordinates(molden);
    const uint8_t* atomic_nums = md_molden_atomic_numbers(molden);
    
    // Initialize system for VIAMD
    md_system_t sys = {0};
    if (md_molden_system_init(&sys, molden, alloc)) {
        // Use system in VIAMD...
        md_system_free(&sys, alloc);
    }
}

md_molden_destroy(molden);
```

## Implementation Status

### ✅ Implemented
- Object lifecycle management
- File parsing infrastructure
- [Atoms] section parser with unit support (Angstrom/Bohr)
- Atom coordinate and atomic number queries
- Molecular property queries
- System integration (md_system_t)
- Unit tests (100% pass rate)

### 🚧 In Progress
- [GTO] section parser for basis sets
- [MO] section parser for molecular orbitals
- GTO extraction for visualization

### 📋 Planned
- [FREQ] section for vibrational frequencies
- [GEOMETRIES] section for optimization trajectories
- Complete orbital coefficient extraction
- Density matrix support
- Additional test coverage

## Testing

Current test suite (`test_molden.c`):
- ✅ `molden.create_destroy` - Object lifecycle
- ✅ `molden.parse_water_molecule` - Atom parsing
- ✅ `molden.system_init` - System integration
- ✅ `molden.reset` - State reset

Test data: `test_data/water.molden` (H2O molecule)

## API Compatibility

The Molden API follows the same design patterns as the VeloxChem API (`md_vlx.h`):
- Consistent naming conventions
- Similar function signatures
- Compatible with VIAMD loader requirements
- Seamless integration with md_system_t

## Performance Considerations

- File parsing uses buffered reading
- Memory allocated using md_array for efficient cleanup
- Coordinate conversions (Bohr ↔ Angstrom) performed during parsing
- GTO extraction will support cutoff values for optimization

## Error Handling

All functions validate inputs and return appropriate error codes:
- NULL pointers handled gracefully
- File I/O errors logged via MD_LOG_ERROR
- Parse errors include line information when available
- Memory allocation failures detected and reported

## Future Enhancements

1. **Basis Set Support**: Complete [GTO] section parsing
2. **Orbital Visualization**: Full MO coefficient extraction and GTO generation
3. **Frequency Analysis**: Parse vibrational modes
4. **Optimization Paths**: Support geometry optimization trajectories
5. **Performance**: Optimize large file parsing with streaming
6. **Validation**: Add file format validation and warnings
