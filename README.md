[![Windows (MSVC)](https://github.com/scanberg/mold/actions/workflows/windows.yml/badge.svg)](https://github.com/scanberg/mold/actions/workflows/windows.yml)
[![Ubuntu 22.04 (GCC 11)](https://github.com/scanberg/mold/actions/workflows/ubuntu22.yml/badge.svg)](https://github.com/scanberg/mold/actions/workflows/ubuntu22.yml)
[![Ubuntu 24.04 (GCC 13)](https://github.com/scanberg/mold/actions/workflows/ubuntu24.yml/badge.svg)](https://github.com/scanberg/mold/actions/workflows/ubuntu24.yml)
[![MacOS (Clang)](https://github.com/scanberg/mold/actions/workflows/macos.yml/badge.svg)](https://github.com/scanberg/mold/actions/workflows/macos.yml)

# mdlib
Molecular Dynamics Library

Collection of utilities and procedures for processing molecular dynamics data.

## script
Scripting language for defining selections, computing properties and evaulating expressions.

## gl
OpenGL draw functionality for drawing moleculecules under different representations.

## util
Utilities for computing secondary structures, com, rmsd, internal reference frame (rotation), covalent bonds, etc.

## Amber files
The Amber topology reader loads `.prmtop`/`.parm7` with ASCII restart coordinates, or
supplements a matching structure with topology metadata and bonds. Separate readers publish
`mdout` and `mden` energy output, including restarted segments.

Enable `-DMD_ENABLE_AMBER_NETCDF=ON` to build the optional NetCDF reader. It requires the
NetCDF-C CMake package and supports convention-checked classic and 64-bit-offset trajectories,
as well as initial coordinates from a classic Amber NetCDF restart or the first trajectory
frame. Trajectory atom order must match the loaded structure or topology. NetCDF-4, legacy
text trajectories, CHAMBER topologies and velocities/forces are not supported.
