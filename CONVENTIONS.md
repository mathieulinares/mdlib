# MDLIB Code Conventions and Patterns

This document captures important conventions and patterns discovered during Molden implementation.

## Build System

### Commands
```bash
cmake -B build -DMD_UNITTEST=ON -DMD_ENABLE_MOLDEN=ON
cmake --build build -j$(nproc)
```

### Optional Features
- Use CMake `option()` for features: `option(MD_ENABLE_FEATURE "Description" OFF)`
- Add conditional compilation: `#ifdef MD_FEATURE`
- Conditionally add files to SRC_FILES in CMakeLists.txt
- Pattern: VLX (lines 138-144), Molden (lines 146-149)

## File Format Support Pattern

### Required Functions
All file format parsers must implement:
1. `md_<format>_create()` - Allocate with backing allocator
2. `md_<format>_reset()` - Clear to initial state
3. `md_<format>_destroy()` - Free all resources
4. `md_<format>_parse_file()` - Parse file
5. `md_<format>_system_init()` - Convert to md_system_t
6. `md_<format>_system_loader()` - Return loader interface

### Example: VLX and Molden
- md_vlx.h (lines 40-43, 45, 177, 179)
- md_molden.h mirrors this exactly

## Memory Management

### Arena Allocator Pattern
**Always use arena allocator for file parsers:**
```c
md_allocator_i* arena = md_arena_allocator_create(backing, MEGABYTES(1));
md_<format>_t* obj = md_alloc(arena, sizeof(md_<format>_t));
```

**Benefits:**
- Single destroy() call frees everything
- No manual free() calls needed
- Cache-efficient memory layout

**References:**
- md_vlx.c: lines 2936-2947 (create), 2950-2957 (reset/destroy)
- md_molden.c: lines 219-236 (follows same pattern)

## Data Structure Organization

### Quantum Chemistry Data Hierarchy
```
md_<format>_t
├── basis_set_t
│   ├── param (exponents, coefficients)
│   ├── basis_func (per-function data)
│   └── atom_basis (per-atom assignments)
├── scf_t
│   ├── alpha (orbital_t)
│   │   ├── coefficients (2D: AO x MO)
│   │   ├── energy (1D: per MO)
│   │   └── occupancy (1D: per MO)
│   └── beta (orbital_t)
└── gto_data_t (atomic orbitals as GTOs)
```

### Pattern Source
- VLX structures: md_vlx.c lines 76-183
- Molden structures: md_molden.c lines 38-115 (mirrors VLX)

## Naming Conventions

### Public API
- Format: `md_<format>_<category>_<function>`
- Examples:
  - `md_molden_number_of_atoms()`
  - `md_molden_scf_mo_energy()`
  - `md_molden_system_init()`

### Internal Structures
- Use typedef for main structure: `typedef struct md_<format>_t md_<format>_t;`
- Internal basis/orbital types use _t suffix
- Keep internal structures in .c file (not exposed in .h)

## Error Handling

### Patterns
- Return `bool` for operations that can fail
- Return `NULL` or `0` for getters when object is invalid
- Use `MD_LOG_ERROR()` for error messages
- Always check parameters: `if (!param) { MD_LOG_ERROR(...); return false; }`

### Examples
- md_vlx.c: throughout (search for `MD_LOG_ERROR`)
- md_molden.c: lines 224-227, 245-249, etc.

## Constants and Units

### Coordinate Conversion
```c
#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.5291772109029999
```

**Internal storage:** Always Ångström  
**Conversion:** Apply during parsing if needed (Molden supports both AU and Angs)

### Source
- md_vlx.c: lines 18-19
- md_molden.c: lines 36-37

## Documentation

### Header Comments
- Each public function should have doc comment
- Document array lengths in return values
- Note optional vs required fields
- Reference similar VLX functions when applicable

### Example
```c
// Get MO energies (array of length number_of_molecular_orbitals)
const double* md_molden_scf_mo_energy(const struct md_molden_t* molden, md_molden_mo_type_t type);
```

## Testing

### Build Verification
Always test with feature:
- ON: `cmake -B build -DMD_ENABLE_FEATURE=ON`
- OFF: `cmake -B build -DMD_ENABLE_FEATURE=OFF`

Both must build successfully without warnings.

## Future Work Pattern

When implementing parsing (Phase 2):
1. Study VLX parsing in md_vlx.c
2. Use md_buffered_reader_t for file I/O
3. Use md_temp_allocator for temporary data
4. Apply same normalization/transformation as VLX
5. Test with real files from quantum chemistry codes

## References

### Key Files
- **VLX Implementation:** src/md_vlx.h, src/md_vlx.c
- **Molden Implementation:** src/md_molden.h, src/md_molden.c  
- **Design Doc:** MOLDEN_DESIGN.md
- **Build System:** CMakeLists.txt
- **Type Definitions:** src/md_types.h
- **GTO Utilities:** src/md_gto.h, src/md_gto.c

### Documentation
- Molden format: http://www.cmbi.ru.nl/molden/molden_format.html
- VeloxChem: https://github.com/VeloxChem/VeloxChem
