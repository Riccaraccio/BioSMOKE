# --------------------------------------------------------------------
# Preliminary Settings
# --------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16)

project(
  BioSMOKEpp
  DESCRIPTION ""
  LANGUAGES C CXX Fortran
  VERSION 0.22.0)

# Set module path in order to use custom CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# ----------------------------------------------------------------------------
# Operating system
# ----------------------------------------------------------------------------
# Detect operating system and define SYSTEM macros
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(SYSTEM_DEFINITIONS SYSTEM_LINUX)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(SYSTEM_DEFINITIONS SYSTEM_DARWIN)
elseif(CMAKE_SYSTEM_NAME STREQUAL "AIX")
  set(SYSTEM_DEFINITIONS SYSTEM_AIX)
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
  set(SYSTEM_DEFINITIONS SYSTEM_WINDOWS)
endif()

# Apply system-specific definitions to the target
if(SYSTEM_DEFINITIONS)
  message(STATUS "System definition: ${SYSTEM_DEFINITIONS}")
endif()

# ----------------------------------------------------------------------------
# Compiler options
# ----------------------------------------------------------------------------
# include(${CMAKE_MODULE_PATH}/CompilerOptions.cmake)

# ----------------------------------------------------------------------------
# Options
# ----------------------------------------------------------------------------
include(${CMAKE_MODULE_PATH}/Options.cmake)

# ----------------------------------------------------------------------------
# Compulsory libraries
# ----------------------------------------------------------------------------
# Eigen https://eigen.tuxfamily.org/index.php?title=Main_Page
find_package(Eigen3 REQUIRED NO_MODULE)

# Boost https://www.boost.org/
find_package(
  Boost REQUIRED
  COMPONENTS date_time
             filesystem
             program_options
             system
             regex
             timer
             chrono)


# OpenSMOKE++ Send an email to alberto.cuoci@polimi.it
find_package(OpenSMOKEpp REQUIRED)

# --------------------------------------------------------------------
# MKL / OpenBLAS
# --------------------------------------------------------------------
set(LINEAR_ALGEBRA_INCLUDE)
set(LINEAR_ALGEBRA_LIB)
if(OPENSMOKE_USE_MKL)
  set(BLA_VENDOR Intel10_64lp_seq)
  find_package(BLAS REQUIRED COMPONENTS mkl_core mkl_sequential)
  find_package(LAPACK REQUIRED)
endif()

if(OPENSMOKE_USE_OPENBLAS)
  # For OpenBLAS we are not using the cmake standard defined find blas stuff
  # because it doesnt provide correctly the include directory so in order to do
  # that we are using the one automatically generated when building from source
  # the OpenBLAS library usually located in ".../lib/cmake/openblas"
  find_package(OpenBLAS REQUIRED)
  message(STATUS "Found OpenBLAS: version ${OpenBLAS_VERSION}")
  list(APPEND LINEAR_ALGEBRA_INCLUDE ${OpenBLAS_INCLUDE_DIRS})
  # The library is handled in this way because we need to link against the
  # static version maybe AC can confirm if this is the right thing
  list(APPEND LINEAR_ALGEBRA_LIB ${OpenBLAS_INCLUDE_DIRS}/../lib/libopenblas.a)
endif()

# --------------------------------------------------------------------
# Project Tree directory
# --------------------------------------------------------------------
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/source)
add_subdirectory(${PROJECT_SOURCE_DIR})
