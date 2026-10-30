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

## Amber NetCDF trajectories
Enable `-DMD_ENABLE_AMBER_NETCDF=ON` to build the optional Amber trajectory reader. It requires the
NetCDF-C CMake package and supports convention-checked classic NetCDF and 64-bit-offset files; the
molecular structure and atom order must come from a separately loaded structure file. This
trajectory reader does not support NetCDF-4, legacy Amber text trajectories, or topology/restart
files. A separate reader supports Amber energy output, both the `mdout` text log and the `mden`
energy file (`ntwe`).
