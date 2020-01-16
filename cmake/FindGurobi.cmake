# Based on https://gitlab.inf.unibe.ch/CGG-public/cmake-library

# Once done this will define
#  GUROBI_FOUND - System has Gurobi
#  Targets:
#       Gurobi::GurobiC   - only the C interface
#       Gurobi::GurobiCXX - C and C++ interface
#       Gurobi::Gurobi    - same as GurobiCXX

find_path(GUROBI_HOME
          NAMES include/gurobi_c++.h
          PATHS
          $ENV{GUROBI_HOME}
          "/opt/gurobi/linux64/"
          )

find_path(GUROBI_INCLUDE_DIR
    NAMES gurobi_c++.h
    HINTS
    "${GUROBI_HOME}/include"
    )

set(GUROBI_BIN_DIR "${GUROBI_HOME}/bin")
set(GUROBI_LIB_DIR "${GUROBI_HOME}/lib")

file(GLOB GUROBI_LIBRARY_LIST
    RELATIVE ${GUROBI_LIB_DIR}
    ${GUROBI_LIB_DIR}/libgurobi*.so
    ${GUROBI_LIB_DIR}/gurobi*.lib
    )

# Ignore libgurobiXY_light.so, libgurobi.so (without version):
string(REGEX MATCHALL
    "gurobi([0-9]+)\\..*"
    GUROBI_LIBRARY_LIST
    "${GUROBI_LIBRARY_LIST}"
    )

string(REGEX REPLACE
    ".*gurobi([0-9]+)\\..*"
    "\\1"
    GUROBI_LIBRARY_VERSIONS
    "${GUROBI_LIBRARY_LIST}")
list(LENGTH GUROBI_LIBRARY_VERSIONS GUROBI_NUMVER)

if (GUROBI_NUMVER EQUAL 0)
elseif (GUROBI_NUMVER EQUAL 1)
    list(GET GUROBI_LIBRARY_VERSIONS 0 GUROBI_LIBRARY_VERSION)
else()
    # none or more than one versioned library -let's try without suffix,
    # maybe the user added a symlink to the desired library
    message(STATUS "Found more than one Gurobi library version (${GUROBI_LIBRARY_VERSIONS}), trying without suffix. Set GUROBI_LIBRARY if you want to pick a certain one.")
    set(GUROBI_LIBRARY_VERSION "")
endif()

find_library(GUROBI_LIBRARY
    NAMES "gurobi${GUROBI_LIBRARY_VERSION}"
    PATHS
    ${GUROBI_LIB_DIR}
)

if (NOT GUROBI_FOUND)
	message(STATUS "Found no Gurobi library. If you have Gurobi installed, make sure GUROBI_HOME is set properly.")
	set(GUROBI_LIBRARY "")
	set(GUROBI_INCLUDE_DIR "")
endif()
