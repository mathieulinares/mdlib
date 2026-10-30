#pragma once

#include <stdbool.h>

#include <core/md_str.h>

struct md_system_t;

#ifdef __cplusplus
extern "C" {
#endif

// Recognizes an Amber mdout file by its energy-record markers. Used for ambiguous .out files
// and the conventional extensionless filename "mdout".
bool md_amber_energy_file_is_mdout(str_t filename);

// Recognizes an Amber mden file (written with ntwe > 0) by its "L0  Nsteps  time(ps)  Etot" header.
bool md_amber_energy_file_is_mden(str_t filename);

// Either of the two above.
bool md_amber_energy_file_is_energy(str_t filename);

// Publishes the raw samples in an Amber mdout or mden file under "<run>/amber_energy", on their own
// time axis, and only when that time axis covers every frame of the run. The format is told from the
// contents; terms mden shares with mdout are published under their mdout names.
bool md_amber_energy_system_supplement_from_file(struct md_system_t* sys, str_t filename, str_t run);

#ifdef __cplusplus
}
#endif
