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

typedef enum {
    MD_MOLDEN_MO_TYPE_ALPHA = 0,
    MD_MOLDEN_MO_TYPE_BETA = 1,
} md_molden_mo_type_t;

typedef enum {
    MD_MOLDEN_SCF_TYPE_UNKNOWN = 0,
    MD_MOLDEN_SCF_TYPE_RESTRICTED,
    MD_MOLDEN_SCF_TYPE_UNRESTRICTED,
} md_molden_scf_type_t;

typedef struct md_molden_t md_molden_t;

// Basic operations

struct md_molden_t* md_molden_create(struct md_allocator_i* backing);
void md_molden_reset(struct md_molden_t* molden);
void md_molden_destroy(struct md_molden_t* molden);

// Parse Molden format files (.molden, .input, etc.)
bool md_molden_parse_file(struct md_molden_t* molden, str_t filename);

// Molecule structure queries
size_t md_molden_number_of_atoms(const struct md_molden_t* molden);
size_t md_molden_number_of_alpha_electrons(const struct md_molden_t* molden);
size_t md_molden_number_of_beta_electrons(const struct md_molden_t* molden);

double md_molden_molecular_charge(const struct md_molden_t* molden);
size_t md_molden_spin_multiplicity(const struct md_molden_t* molden);

const dvec3_t* md_molden_atom_coordinates(const struct md_molden_t* molden);
const uint8_t* md_molden_atomic_numbers(const struct md_molden_t* molden);

// Maps AO index to atom index, length is number of atomic orbitals
const int* md_molden_ao_to_atom_idx(const struct md_molden_t* molden);

// SCF type and molecular orbital information
md_molden_scf_type_t md_molden_scf_type(const struct md_molden_t* molden);
size_t md_molden_scf_homo_idx(const struct md_molden_t* molden, md_molden_mo_type_t type);
size_t md_molden_scf_lumo_idx(const struct md_molden_t* molden, md_molden_mo_type_t type);

size_t md_molden_scf_number_of_atomic_orbitals(const struct md_molden_t* molden);
size_t md_molden_scf_number_of_molecular_orbitals(const struct md_molden_t* molden);

const double* md_molden_scf_mo_occupancy(const struct md_molden_t* molden, md_molden_mo_type_t type);
const double* md_molden_scf_mo_energy(const struct md_molden_t* molden, md_molden_mo_type_t type);

// Atomic orbital GTO data
bool md_molden_scf_extract_gto_data(md_gto_data_t* out_gto_data, const struct md_molden_t* molden, double cutoff_value, struct md_allocator_i* alloc);

// Extract Molecular Orbital (MO) GTOs
// Provides an upper limit to the number of GTOs
size_t md_molden_mo_gto_count(const struct md_molden_t* molden);

// Extract GTOs for a specific molecular orbital
// Returns the number of GTOs written
// gtos: array to hold the extracted GTOs
// molden: pointer to valid Molden object
// mo_idx: Molecular Orbital Index
// type: Molecular orbital type: Alpha / Beta
// value_cutoff: cutoff value for GTO radius of influence (0 == no cutoff)
size_t md_molden_mo_gto_extract(md_gto_t gtos[], const struct md_molden_t* molden, size_t mo_idx, md_molden_mo_type_t type, double value_cutoff);

// MOLECULE SYSTEM INTEGRATION
bool md_molden_system_init(struct md_system_t* sys, const md_molden_t* molden, struct md_allocator_i* alloc);

struct md_system_loader_i* md_molden_system_loader(void);

#ifdef __cplusplus
}
#endif
