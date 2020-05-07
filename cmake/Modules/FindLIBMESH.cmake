# Finds the libMesh installation using CMake's PkgConfig
# module and creates a libmesh imported target

find_path(LIBMESH_PC NAMES libmesh-opt.pc
          HINTS ${LIBMESH_DIR} $ENV{LIBMESH_ROOT}
          PATHS ENV LD_LIBRARY_PATH
          PATH_SUFFIXES lib/pkgconfig pkgconfig
          NO_DEFAULT_PATH)

set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH OFF)
set(PKG_CONFIG_USE_CMAKE_ENVIRONMENT_PATH OFF)
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${LIBMESH_PC}")
message(STATUS "PKG_CONFIG_PATH: $ENV{PKG_CONFIG_PATH}")
find_package(PkgConfig REQUIRED)

pkg_check_modules(LIBMESH REQUIRED IMPORTED_TARGET libmesh)
pkg_check_modules(LIBMESH_DBG REQUIRED IMPORTED_TARGET libmesh-dbg)
pkg_check_modules(LIBMESH_OPT REQUIRED IMPORTED_TARGET libmesh-opt)

# older versions of CMake won't check the system libraries
# if '-L' appears in the package config string. This
# block makes sure the system locations are checked
# if any of the libraries aren't found.
get_cmake_property(ALL_VARIABLES VARIABLES)
foreach(VARIABLE ${ALL_VARIABLES})
  if (var MATCHES "^pkgcfg_lib_LIBMESH" AND ${var} MATCHES "NOTFOUND$")
    # get the library name
    string(REGEX REPLACE "^(.*[_])" "" LIBNAME ${var})
    find_library(${var} NAMES ${LIBNAME})
  endif()
endforeach()


message(STATUS "Found LIBMESH in ${LIBMESH_PC}")
