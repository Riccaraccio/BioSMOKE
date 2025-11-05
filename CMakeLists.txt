# ----------------------------------------------------------------------------
# Preliminary Settings
# ----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16)

project(
  BioSMOKEpp
  DESCRIPTION ""
  LANGUAGES C CXX Fortran
  VERSION 0.22.0)

# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# Apply system-specific definitions to the target
if(SYSTEM_DEFINITIONS)
  message(STATUS "System definition: ${SYSTEM_DEFINITIONS}")
endif()

# ----------------------------------------------------------------------------
# Compiler options
# ----------------------------------------------------------------------------
include(${CMAKE_MODULE_PATH}/CompilerOptions.cmake)

# ----------------------------------------------------------------------------
# Options
# ----------------------------------------------------------------------------
include(${CMAKE_MODULE_PATH}/Options.cmake)

# ----------------------------------------------------------------------------
# Compulsory libraries
# ----------------------------------------------------------------------------
# Eigen3: https://eigen.tuxfamily.org/index.php?title=Main_Page
find_package(Eigen3 REQUIRED NO_MODULE)

# ----------------------------------------------------------------------------
# Boost: https://www.boost.org/
find_package(
  Boost REQUIRED
  COMPONENTS date_time
             filesystem
             program_options
             system
             regex
             timer
             chrono)

# ----------------------------------------------------------------------------
# OpenSMOKE++: Send an email to <alberto.cuoci@polimi.it>
find_package(OpenSMOKEpp REQUIRED)

# ----------------------------------------------------------------------------
# Linear Algebra acceleration
# ----------------------------------------------------------------------------
set(LINEAR_ALGEBRA_INCLUDE)
set(LINEAR_ALGEBRA_LIB)
set(OPENBLAS_VALUE 0)
set(MKL_VALUE 0)

# ----------------------------------------------------------------------------
# MKL
if (OPENSMOKE_USE_MKL)
  set (BLA_VENDOR Intel10_64lp)
  set (BLA_STATIC ON)
  find_package (BLAS)
  
  # For MKL, BLAS includes LAPACK functionality
  set(LAPACK_LIBRARIES ${BLAS_LIBRARIES})
  set(LAPACK_FOUND TRUE)
 
  message (STATUS "MKL_LIBRARIES = ${BLAS_LIBRARIES}")
  message (STATUS "MKL_DIR       = ${BLAS_DIR}")
  message (STATUS "LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}")
  list(APPEND LINEAR_ALGEBRA_INCLUDE ${BLAS_DIR})
  list(APPEND LINEAR_ALGEBRA_LIB ${BLAS_LIBRARIES})
  add_compile_definitions(OPENSMOKE_USE_MKL=1)
endif()

# ----------------------------------------------------------------------------
# OpenBLAS: http://www.openmathlib.org/OpenBLAS/
if(BIOSMOKE_USE_OPENBLAS)
  # For OpenBLAS we are not using the cmake standard defined find blas stuff
  # because it doesnt provide correctly the include directory so in order to do
  # that we are using the one automatically generated when building from source
  # the OpenBLAS library usually located in ".../lib/"
  find_package(OpenBLAS REQUIRED)
  message(STATUS "Found OpenBLAS VERSION: ${OpenBLAS_VERSION}")
  message(STATUS "Found OpenBLAS: ${OpenBLAS_INCLUDE_DIRS}")
  set(OPENBLAS_VALUE 1)
  set(MKL_VALUE 0)
  list(APPEND LINEAR_ALGEBRA_INCLUDE ${OpenBLAS_INCLUDE_DIRS})
  list(APPEND LINEAR_ALGEBRA_LIB ${OpenBLAS_INCLUDE_DIRS}/../lib/libopenblas.a)
endif()

# ----------------------------------------------------------------------------
# Include install definitions
# ----------------------------------------------------------------------------
include(${CMAKE_MODULE_PATH}/InstallationOptions.cmake)

# ----------------------------------------------------------------------------
# Project Tree directory
# ----------------------------------------------------------------------------
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/source)
add_subdirectory(${PROJECT_SOURCE_DIR})
