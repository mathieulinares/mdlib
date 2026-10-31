#pragma once

#include <stdbool.h>
#include <stdint.h>

#include <core/md_str.h>

struct md_system_t;
struct md_system_state_t;

#ifdef __cplusplus
extern "C" {
#endif

// True only for supported classic NetCDF files carrying the Amber trajectory convention.
bool md_amber_netcdf_file_is_trajectory(str_t filename);

// Coordinates and cell from frame zero of AMBER, or the single AMBERRESTART state.
// Checks the atom count before modifying the state. NetCDF-4 remains unsupported.
bool md_amber_netcdf_state_init_from_file(struct md_system_state_t* state, str_t filename, size_t num_atoms);

// Publishes an Amber NetCDF trajectory as a streamed run. Coordinates and cells are Angstrom;
// time is picoseconds when present, and otherwise the run uses unitless frame ordinals.
bool md_amber_netcdf_system_publish_run(struct md_system_t* sys, str_t filename, str_t run, uint32_t flags);

#ifdef __cplusplus
}
#endif
