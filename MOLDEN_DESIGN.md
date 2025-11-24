# Molden File Format Support - Design Document

## Overview

This document describes the design and implementation of Molden file format support in MDLIB, following the patterns established by the VeloxChem (VLX) implementation.

## Design Goals

1. **API Consistency**: Mirror the VLX API structure to maintain consistency across the codebase
2. **Minimal Dependencies**: Use text-based parsing (no HDF5 dependency like VLX)
3. **Memory Efficiency**: Use arena allocator pattern for efficient memory management
4. **Extensibility**: Design for future additions (vibrational modes, frequency analysis, etc.)

## Molden Format Overview

The Molden format is a text-based file format used in quantum chemistry to store:
- Molecular geometries
- Basis set information (Gaussian Type Orbitals)
- Molecular orbital coefficients
- Energies and occupancies
- Symmetry labels

### Key Format Blocks

| Block Name | Purpose | Status |
|------------|---------|--------|
| `[Molden Format]` | Version identifier | Planned |
| `[Atoms]` | Atomic coordinates (Angs or AU) | Planned |
| `[GTO]` | Gaussian basis set definition | Planned |
| `[MO]` | Molecular orbital data | Planned |
| `[5D]` | Use 5 d-orbitals instead of 6 | Planned |
| `[7F]` | Use 7 f-orbitals instead of 10 | Planned |
| `[9G]` | Use 9 g-orbitals instead of 15 | Planned |

## Data Structure Design

### Core Structure: `md_molden_t`

The main data structure mirrors `md_vlx_t` with the following components:

```c
typedef struct md_molden_t {
    basis_set_t basis_set;              // Basis set definition
    str_t basis_set_ident;              // Basis set name/identifier
    
    size_t number_of_atoms;             // Number of atoms
    size_t number_of_alpha_electrons;   // Alpha electrons
    size_t number_of_beta_electrons;    // Beta electrons
    
    double molecular_charge;            // Total charge
    size_t spin_multiplicity;           // 2S+1
    
    dvec3_t* atom_coordinates;          // Atomic positions (Ångström)
    md_element_t* atomic_numbers;       // Atomic numbers
    int* ao_to_atom_idx;               // AO to atom mapping
    
    md_molden_scf_t scf;               // SCF/MO data
    md_gto_data_t gto_data;            // GTO representation
    
    bool use_5d, use_7f, use_9g;       // Basis function formats
    
    md_allocator_i* arena;             // Memory allocator
} md_molden_t;
```

### Basis Set Representation

Following the VLX pattern exactly:

```c
typedef struct basis_set_func_t {
    uint8_t  type;          // Angular momentum (s=0, p=1, d=2, ...)
    uint8_t  param_count;   // Number of primitives
    uint16_t param_offset;  // Offset into parameter arrays
} basis_set_func_t;

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
    
    struct {
        size_t count;
        basis_set_basis_t* data;  // Indexed by atomic number
    } atom_basis;
} basis_set_t;
```

### Molecular Orbital Data

```c
typedef struct md_molden_orbital_t {
    md_molden_2d_data_t coefficients;  // [num_ao x num_mo]
    md_molden_1d_data_t energy;        // [num_mo]
    md_molden_1d_data_t occupancy;     // [num_mo]
    str_t* symmetry;                    // [num_mo] optional
    str_t* spin;                        // [num_mo] optional
} md_molden_orbital_t;

typedef struct md_molden_scf_t {
    md_molden_scf_type_t type;
    size_t homo_idx[2];  // [alpha, beta]
    size_t lumo_idx[2];  // [alpha, beta]
    
    md_molden_orbital_t alpha;
    md_molden_orbital_t beta;
} md_molden_scf_t;
```

## API Design

### Naming Convention

All public functions follow the pattern: `md_molden_<category>_<function>`

Examples:
- `md_molden_create()`
- `md_molden_parse_file()`
- `md_molden_number_of_atoms()`
- `md_molden_scf_mo_energy()`

This mirrors the VLX naming: `md_vlx_<category>_<function>`

### API Categories

1. **Lifecycle Management**
   - `md_molden_create()` - Allocate and initialize
   - `md_molden_reset()` - Clear and reset to initial state
   - `md_molden_destroy()` - Free all resources

2. **File I/O**
   - `md_molden_parse_file()` - Parse Molden file

3. **Molecular Properties**
   - `md_molden_number_of_atoms()`
   - `md_molden_number_of_alpha_electrons()`
   - `md_molden_number_of_beta_electrons()`
   - `md_molden_molecular_charge()`
   - `md_molden_spin_multiplicity()`
   - `md_molden_basis_set_ident()`

4. **Atomic Data**
   - `md_molden_atom_coordinates()`
   - `md_molden_atomic_numbers()`
   - `md_molden_ao_to_atom_idx()`

5. **SCF/MO Data**
   - `md_molden_scf_type()`
   - `md_molden_scf_homo_idx()`
   - `md_molden_scf_lumo_idx()`
   - `md_molden_scf_number_of_atomic_orbitals()`
   - `md_molden_scf_number_of_molecular_orbitals()`
   - `md_molden_scf_mo_occupancy()`
   - `md_molden_scf_mo_energy()`
   - `md_molden_scf_mo_symmetry()`
   - `md_molden_scf_mo_spin()`

6. **GTO Operations**
   - `md_molden_scf_extract_gto_data()`
   - `md_molden_mo_gto_count()`
   - `md_molden_mo_gto_extract()`

7. **System Integration**
   - `md_molden_system_init()`
   - `md_molden_system_loader()`

## Mapping Molden Format to VLX Structures

### Atoms Block Mapping

**Molden Format:**
```
[Atoms] Angs
C     1    6    0.000000    0.000000    0.000000
H     2    1    1.089000    0.000000    0.000000
...
```

**Maps to:**
- `atom_coordinates[]` - positions in Ångström
- `atomic_numbers[]` - element Z values
- `number_of_atoms` - count

**Coordinate Units:** Molden supports both `Angs` (Ångström) and `AU` (Bohr). All coordinates are internally stored in Ångström, with conversion applied during parsing if needed.

### GTO Block Mapping

**Molden Format:**
```
[GTO]
  1 0
 s   3 1.00
      71.6168370              0.15432897
      13.0450960              0.53532814
       3.5305122              0.44463454
 p   2 1.00
       2.9412494             -0.09996723
       0.6834831              0.39951283
```

**Maps to:**
- `basis_set.param.exponents[]` - Gaussian exponents
- `basis_set.param.coefficients[]` - Contraction coefficients
- `basis_set.basis_func[]` - Basis function definitions
- `basis_set.atom_basis[]` - Per-atom basis assignments

### MO Block Mapping

**Molden Format:**
```
[MO]
 Sym= A
 Ene= -11.2179
 Spin= Alpha
 Occup= 2.000000
  1   0.994510
  2  -0.025910
  ...
```

**Maps to:**
- `scf.alpha.energy.data[]` - Orbital energies
- `scf.alpha.occupancy.data[]` - Occupation numbers
- `scf.alpha.coefficients.data[]` - MO coefficients (2D: AO x MO)
- `scf.alpha.symmetry[]` - Symmetry labels (optional)
- `scf.alpha.spin[]` - Spin labels (optional)

## Key Differences from VLX

### 1. File Format
- **VLX**: Binary HDF5 format (requires HDF5 library)
- **Molden**: Text-based format (only standard C library needed)

### 2. Data Organization
- **VLX**: Multiple files (`.out`, `.h5`, `_NTO.h5`, etc.)
- **Molden**: Single file containing all data

### 3. Basis Function Conventions
- **VLX**: Uses VeloxChem internal conventions
- **Molden**: Supports both Cartesian and spherical harmonics
  - `[5D]` - Use 5 spherical d-orbitals (default)
  - `[7F]` - Use 7 spherical f-orbitals (default)
  - `[9G]` - Use 9 spherical g-orbitals (default)

### 4. Symmetry Information
- **VLX**: Limited symmetry labels
- **Molden**: Explicit symmetry labels for each MO (Sym= field)

### 5. Coordinate Units
- **VLX**: Always in Ångström
- **Molden**: Can be Ångström or Bohr (specified in [Atoms] block)

## Implementation Status

### Phase 1: Data Structure and API Design ✅
- [x] Define `md_molden_t` structure
- [x] Define API function signatures
- [x] Create header file (`md_molden.h`)
- [x] Create stub implementation (`md_molden.c`)
- [x] Document design decisions

### Phase 2: Parsing Implementation (Future)
- [ ] Implement `parse_molden_file()`
- [ ] Parse `[Atoms]` block
- [ ] Parse `[GTO]` block
- [ ] Parse `[MO]` block
- [ ] Handle coordinate unit conversion (Bohr ↔ Ångström)
- [ ] Handle basis function format flags (`[5D]`, `[7F]`, `[9G]`)

### Phase 3: GTO and Orbital Extraction (Future)
- [ ] Implement `md_molden_scf_extract_gto_data()`
- [ ] Implement `md_molden_mo_gto_extract()`
- [ ] Implement basis set normalization
- [ ] Implement spherical/Cartesian transformation

### Phase 4: Integration and Testing (Future)
- [ ] Add unit tests
- [ ] Test with real Molden files
- [ ] VIAMD integration
- [ ] Performance optimization

## Memory Management

Following VLX pattern:
- Uses **arena allocator** for efficient bulk allocation
- All data freed in single `md_molden_destroy()` call
- No manual free calls for individual allocations
- Memory layout optimized for cache efficiency

## Error Handling

Following VLX conventions:
- Boolean return values for operations that can fail
- `MD_LOG_ERROR()` for error reporting
- NULL checks on all pointer parameters
- Graceful degradation (return 0/NULL on error)

## Future Extensions

### Potential Additions
1. **Vibrational Analysis**: Parse `[FREQ]` and `[FR-COORD]` blocks
2. **CI/MCSCF Data**: Parse `[GEOCONV]` and other CI blocks
3. **Optimization Trajectories**: Support geometry optimization steps
4. **Property Surfaces**: Parse density and potential surface data

### VIAMD Integration
- Create `molden.h`/`molden.cpp` wrapper in VIAMD repository
- Implement visualization for:
  - Molecular orbitals
  - Electron density
  - Orbital energies and occupancies
- Follow `veloxchem.cpp` integration pattern

## References

### Molden Format Documentation
- Official specification: http://www.cmbi.ru.nl/molden/molden_format.html
- Basis set conventions: Various QM software (Gaussian, GAMESS, etc.)

### VeloxChem Implementation
- `md_vlx.h` - Header file with API definitions
- `md_vlx.c` - Implementation with HDF5 parsing and GTO extraction

### Related Standards
- Gaussian cube format (for density visualization)
- WFN/WFX formats (alternative wavefunction formats)

## Contributors

Design and initial implementation following VLX patterns established in MDLIB.

---

**Document Version**: 1.0  
**Last Updated**: 2025-11-24  
**Status**: Phase 1 Complete (Data Structures and API Design)
