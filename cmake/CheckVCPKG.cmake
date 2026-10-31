# We do not want to have a hard dependency on VCPKG
# As it is only used for big external libs
# So if VCPKG_ROOT is defined and no toolchain was supplied, use vcpkg's toolchain
if (NOT DEFINED CMAKE_TOOLCHAIN_FILE
    AND (NOT DEFINED ENV{CMAKE_TOOLCHAIN_FILE} OR "$ENV{CMAKE_TOOLCHAIN_FILE}" STREQUAL "")
    AND DEFINED ENV{VCPKG_ROOT})
    if(WIN32)
        set(VCPKG_TARGET_TRIPLET "x64-windows-static" CACHE STRING "VCPKG Target Triplet to use")
    endif()

    # VCPKG_ROOT uses backslashes on Windows. Convert to a CMake path, otherwise the path gets
    # written verbatim into CMakeSystem.cmake as include("C:\...") and fails to parse as an invalid escape.
    file(TO_CMAKE_PATH "$ENV{VCPKG_ROOT}" _md_vcpkg_root)
    set(CMAKE_TOOLCHAIN_FILE "${_md_vcpkg_root}/scripts/buildsystems/vcpkg.cmake")
    message(STATUS "VCPKG_ROOT is defined, using vcpkg-toolchain: ${CMAKE_TOOLCHAIN_FILE}")
endif()

# This file is included before project(), whether the toolchain was selected above or explicitly
# supplied by the caller. Manifest features must be set before vcpkg's toolchain is loaded.
set(_md_vcpkg_toolchain "")
if (DEFINED CMAKE_TOOLCHAIN_FILE)
    set(_md_vcpkg_toolchain "${CMAKE_TOOLCHAIN_FILE}")
elseif (DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
    set(_md_vcpkg_toolchain "$ENV{CMAKE_TOOLCHAIN_FILE}")
endif()
file(TO_CMAKE_PATH "${_md_vcpkg_toolchain}" _md_vcpkg_toolchain)
get_filename_component(_md_vcpkg_toolchain_name "${_md_vcpkg_toolchain}" NAME)

if (_md_vcpkg_toolchain_name STREQUAL "vcpkg.cmake"
    AND (NOT DEFINED VCPKG_MANIFEST_MODE OR VCPKG_MANIFEST_MODE)
    AND (MD_ENABLE_HDF5 OR MD_ENABLE_AMBER_NETCDF OR VIAMD_ENABLE_QM OR VIAMD_ENABLE_VELOXCHEM))
    # When mdlib is a subdirectory, its manifest is not next to the top-level CMakeLists.txt.
    get_filename_component(MD_VCPKG_MANIFEST_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
    if (NOT DEFINED VCPKG_MANIFEST_DIR
        AND NOT EXISTS "${CMAKE_SOURCE_DIR}/vcpkg.json"
        AND EXISTS "${MD_VCPKG_MANIFEST_DIR}/vcpkg.json")
        set(VCPKG_MANIFEST_DIR "${MD_VCPKG_MANIFEST_DIR}")
        message(STATUS "Using mdlib vcpkg manifest: ${VCPKG_MANIFEST_DIR}/vcpkg.json")
    endif()

    if (DEFINED VCPKG_MANIFEST_DIR)
        get_filename_component(_md_vcpkg_selected_manifest "${VCPKG_MANIFEST_DIR}" ABSOLUTE)
    else()
        set(_md_vcpkg_selected_manifest "${CMAKE_SOURCE_DIR}")
    endif()
    if (MD_ENABLE_AMBER_NETCDF AND _md_vcpkg_selected_manifest STREQUAL MD_VCPKG_MANIFEST_DIR)
        list(APPEND VCPKG_MANIFEST_FEATURES amber-netcdf)
        list(REMOVE_DUPLICATES VCPKG_MANIFEST_FEATURES)
    endif()
endif()