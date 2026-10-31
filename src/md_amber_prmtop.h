#pragma once

#include <stdbool.h>
#include <core/md_str.h>

struct md_system_t;
struct md_system_state_t;

#ifdef __cplusplus
extern "C" {
#endif

bool md_amber_prmtop_file_is_topology(str_t filename);

// A topology has no coordinates. coordinates must name a matching ASCII rst7/inpcrd, a
// classic NetCDF AMBERRESTART, or an Amber trajectory (whose first frame initializes the state).
// NetCDF input requires MD_AMBER_NETCDF. Atom order is the caller's responsibility.
bool md_amber_prmtop_system_init_from_file(struct md_system_t* sys, struct md_system_state_t* state,
                                         str_t filename, str_t coordinates);

// Applies topology metadata to a loaded system whose atom and residue order matches the topology.
bool md_amber_prmtop_system_supplement_from_file(struct md_system_t* sys, str_t filename);

#ifdef __cplusplus
}
#endif
