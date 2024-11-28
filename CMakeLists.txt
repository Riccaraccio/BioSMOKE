# --------------------------------------------------------------------
# Preliminary Settings
# --------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16)

project(
  BioSMOKE
  DESCRIPTION ""
  LANGUAGES CXX
  VERSION 0.1.0) 

# Set module path in order to use custom CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/source)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Building project in mode: ${CMAKE_BUILD_TYPE}")

# --------------------------------------------------------------------
# Options
# --------------------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Gw")  # -> whole program global optimization
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /GS-") # -> buffer security check: no 
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /GL")  # -> whole program optimization: enable link-time code generation (disables Zi)
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /GF")  # -> enable string pooling)
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj")  # -> enable string pooling)
endif()
# --------------------------------------------------------------------
# Compulsory libraries
# --------------------------------------------------------------------

# Boost (Compulsory)
find_package(
  Boost REQUIRED
  COMPONENTS date_time
             filesystem
             program_options
             system
             regex
             timer
             chrono)

# Eigen (Compulsory)
find_package(Eigen3 REQUIRED NO_MODULE)

# OpenSmoke++ (Compulsory)
find_package(OpenSMOKEpp REQUIRED)

set(COMPULSORY_LIBRARIES ${Boost_LIBRARIES})

set(COMPULSORY_INCLUDE ${Boost_INCLUDE_DIRS}
                       ${OPENSMOKE_INCLUDE_DIR} ${Config_INCLUDE_DIR})

# --------------------------------------------------------------------
# Semi-Compulsory libraries (compulsory just for some solvers)
# --------------------------------------------------------------------

set(LINEAR_ALGEBRA_INCLUDE)
set(LINEAR_ALGEBRA_LIB)

if(OPENSMOKE_USE_MKL AND OPENSMOKE_USE_OPENBLAS)
  message(
    FATAL_ERROR
      "Solvers can be compiled with just one BLAS distribution at a time! Choose between MKL and OpenBLAS"
  )
endif()

# Intel MKL (Optional but strongly rcommended: choose one between MKL=BLAS and
# OpenBLAS)
option(OPENSMOKE_USE_MKL "Activate MKL/BLAS support for OpenSMOKEpp" FALSE)
if(OPENSMOKE_USE_MKL)
  set(BLA_VENDOR Intel10_64lp_seq)
  find_package(BLAS REQUIRED COMPONENTS mkl_core mkl_sequential)
  message(STATUS $ENV{BLAS_ROOT}/include)
  list(APPEND LINEAR_ALGEBRA_INCLUDE ${BLAS_ROOT}/include)
  list(APPEND LINEAR_ALGEBRA_LIB ${BLAS_LIBRARIES})
endif()

# OpenBLAS (Optional but strongly Recommended: choose one between MKL and
# OpenBLAS)
option(OPENSMOKE_USE_OPENBLAS "Activate OpenBLAS support for OpenSMOKEpp" FALSE)
if(OPENSMOKE_USE_OPENBLAS)
  set(BLA_VENDOR OpenBLAS)
  set(BLA_STATIC ON)
  find_package(BLAS REQUIRED)
  list(APPEND LINEAR_ALGEBRA_INCLUDE $ENV{BLAS_ROOT}/include)
  list(APPEND LINEAR_ALGEBRA_LIB ${BLAS_LIBRARIES})
endif()

# --------------------------------------------------------------------
# Compilation flags and options
# --------------------------------------------------------------------

option(OPENSMOKE_USE_OPENMP "Activate support for OpenMP" FALSE)
if(OPENSMOKE_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  list(APPEND OPTIONAL_LIBRARIES ${OpenMP_CXX_LIBRARIES})
  # add_compile_definitions(OPENSMOKE_USE_OPENMP=1)
endif()
add_compile_options(-fPIC -m64 -O3)
add_link_options(-fPIC -m64 -O3)

# --------------------------------------------------------------------
# Project Tree directory
# --------------------------------------------------------------------
message(STATUS "${PROJECT_SOURCE_DIR}")
add_subdirectory(${PROJECT_SOURCE_DIR})
