#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <md_types.h>
#include <md_gto.h>
#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_system_t;
struct md_system_loader_i;

/*
 * Molden File Format Support for MDLIB
 * 
 * This module provides support for reading Molden format files, which are commonly used
 * in quantum chemistry for storing molecular orbitals, basis sets, and related data.
 * 
 * Design follows the VeloxChem (md_vlx) implementation patterns to maintain consistency
 * across the MDLIB codebase.
 * 
 * Molden Format Blocks Supported:
 * - [Molden Format]: Version information
 * - [Atoms]: Atomic coordinates and elements (Angs or AU)
 * - [GTO]: Gaussian Type Orbital basis set information
 * - [MO]: Molecular orbital coefficients, energies, occupancies, symmetry labels
 * - [5D], [7F], [9G]: Basis function format specifiers
 * 
 * Key Differences from VLX:
 * - Molden files are text-based (unlike VLX's HDF5 format)
 * - Molden uses different conventions for spherical/Cartesian basis functions
 * - Molden may include symmetry labels for MOs
 * - No separate excited state or vibrational analysis files (single file format)
 */

// Molecular orbital type (analogous to md_vlx_mo_type_t)
typedef enum {
	MD_MOLDEN_MO_TYPE_ALPHA = 0,
	MD_MOLDEN_MO_TYPE_BETA = 1,
} md_molden_mo_type_t;

// Coordinate unit type (Molden supports both Angstrom and Bohr)
typedef enum {
	MD_MOLDEN_COORD_UNIT_ANGS = 0,
	MD_MOLDEN_COORD_UNIT_AU = 1,
} md_molden_coord_unit_t;

// SCF type classification
typedef enum {
	MD_MOLDEN_SCF_TYPE_UNKNOWN = 0,
	MD_MOLDEN_SCF_TYPE_RESTRICTED,
	MD_MOLDEN_SCF_TYPE_UNRESTRICTED,
} md_molden_scf_type_t;

typedef struct md_molden_t md_molden_t;

// =============================
// Basic operations
// =============================

// Create a new Molden data structure
struct md_molden_t* md_molden_create(struct md_allocator_i* backing);

// Reset the Molden data structure to initial state
void md_molden_reset(struct md_molden_t* molden);

// Destroy and free the Molden data structure
void md_molden_destroy(struct md_molden_t* molden);

// Parse Molden data from a file
// Supports standard Molden format (.molden, .mold, .mol)
bool md_molden_parse_file(struct md_molden_t* molden, str_t filename);

// =============================
// Molecular properties
// =============================

// Get the number of atoms in the molecule
size_t md_molden_number_of_atoms(const struct md_molden_t* molden);

// Get the number of alpha electrons
size_t md_molden_number_of_alpha_electrons(const struct md_molden_t* molden);

// Get the number of beta electrons
size_t md_molden_number_of_beta_electrons(const struct md_molden_t* molden);

// Get molecular charge
double md_molden_molecular_charge(const struct md_molden_t* molden);

// Get spin multiplicity (2S+1)
size_t md_molden_spin_multiplicity(const struct md_molden_t* molden);

// Get basis set identifier/name
str_t md_molden_basis_set_ident(const struct md_molden_t* molden);

// =============================
// Atomic data
// =============================

// Get atom coordinates (array of length number_of_atoms)
// Returns coordinates in Ångström
const dvec3_t* md_molden_atom_coordinates(const struct md_molden_t* molden);

// Get atomic numbers (array of length number_of_atoms)
const md_element_t* md_molden_atomic_numbers(const struct md_molden_t* molden);

// Maps AO index to atom index, length is number of atomic orbitals
const int* md_molden_ao_to_atom_idx(const struct md_molden_t* molden);

// =============================
// SCF / Molecular Orbital data
// =============================

// Get the SCF type (restricted, unrestricted)
md_molden_scf_type_t md_molden_scf_type(const struct md_molden_t* molden);

// Get HOMO index for the specified spin type
size_t md_molden_scf_homo_idx(const struct md_molden_t* molden, md_molden_mo_type_t type);

// Get LUMO index for the specified spin type
size_t md_molden_scf_lumo_idx(const struct md_molden_t* molden, md_molden_mo_type_t type);

// Get number of atomic orbitals
size_t md_molden_scf_number_of_atomic_orbitals(const struct md_molden_t* molden);

// Get number of molecular orbitals
size_t md_molden_scf_number_of_molecular_orbitals(const struct md_molden_t* molden);

// Get MO occupancies (array of length number_of_molecular_orbitals)
const double* md_molden_scf_mo_occupancy(const struct md_molden_t* molden, md_molden_mo_type_t type);

// Get MO energies (array of length number_of_molecular_orbitals)
const double* md_molden_scf_mo_energy(const struct md_molden_t* molden, md_molden_mo_type_t type);

// Get MO symmetry labels (array of length number_of_molecular_orbitals)
// Returns NULL if symmetry labels are not available
const str_t* md_molden_scf_mo_symmetry(const struct md_molden_t* molden, md_molden_mo_type_t type);

// Get MO spin labels (array of length number_of_molecular_orbitals)
// Returns NULL if spin labels are not available (restricted calculation)
const str_t* md_molden_scf_mo_spin(const struct md_molden_t* molden, md_molden_mo_type_t type);

// =============================
// Basis set and GTO data
// =============================

// Extract GTO data for atomic orbitals
// Similar to md_vlx_scf_extract_gto_data
bool md_molden_scf_extract_gto_data(md_gto_data_t* out_gto_data, const struct md_molden_t* molden, double cutoff_value, struct md_allocator_i* alloc);

// Get the maximum count of GTOs needed to represent molecular orbitals
size_t md_molden_mo_gto_count(const struct md_molden_t* molden);

// Extract GTOs for a specific molecular orbital
// gtos: array to hold the extracted gtos (must be pre-allocated with size from md_molden_mo_gto_count)
// molden: pointer to valid Molden object
// mo_idx: Molecular Orbital Index (0-based)
// type: Molecular orbital type (Alpha / Beta)
// value_cutoff: Cutoff value for GTO radius of influence (0 == no cutoff)
// Returns the number of GTOs written
size_t md_molden_mo_gto_extract(md_gto_t gtos[], const struct md_molden_t* molden, size_t mo_idx, md_molden_mo_type_t type, double value_cutoff);

// =============================
// System integration
// =============================

// Initialize md_system_t from Molden data
bool md_molden_system_init(struct md_system_t* sys, const struct md_molden_t* molden, struct md_allocator_i* alloc);

// Get system loader interface for Molden files
struct md_system_loader_i* md_molden_system_loader(void);

#ifdef __cplusplus
}
#endif
