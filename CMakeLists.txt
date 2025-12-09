# --------------------------------------------------------------------
# Preliminary Settings
# --------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16)

project(
  BioSMOKE
  DESCRIPTION ""
  LANGUAGES CXX
  VERSION 0.1.0) 

function(apply_compiler_options TARGET_NAME)
  # Set C++ standard for the target
  target_compile_features(${TARGET_NAME} PUBLIC cxx_std_17)

  # Compiler-specific options
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Release build options
    target_compile_options(${TARGET_NAME} PRIVATE
      $<$<CONFIG:Release>:-O3>
      $<$<CONFIG:Release>:-fPIC>
      # The following are dangerous here but they are annoying
      $<$<CONFIG:Release>:-Wno-unused-value>
      $<$<CONFIG:Release>:-Wno-switch>
      $<$<CONFIG:Release>:-Wno-return-type-c-linkage>
      $<$<CONFIG:Release>:-Wno-self-assign-field>
      $<$<CONFIG:Release>:-Wno-dangling-else>
      $<$<CONFIG:Release>:-Wno-mismatched-new-delete>
      $<$<CONFIG:Release>:-Wno-delete-abstract-non-virtual-dtor>
      $<$<CONFIG:Release>:-Wno-tautological-undefined-compare>
    )

    # Debug build options
    target_compile_options(${TARGET_NAME} PRIVATE
      $<$<CONFIG:Debug>:-O0>
      $<$<CONFIG:Debug>:-g>
    )

    # Debug mode warnings (only for project targets, not dependencies)
    target_compile_options(${TARGET_NAME} PRIVATE
      $<$<CONFIG:Debug>:-Wall>
      $<$<CONFIG:Debug>:-Wextra>
      $<$<CONFIG:Debug>:-Wunused>
      $<$<CONFIG:Debug>:-Wswitch-default>
      $<$<CONFIG:Debug>:-Wuninitialized>
      $<$<CONFIG:Debug>:-Woverloaded-virtual>
      $<$<CONFIG:Debug>:-Wwrite-strings>
      $<$<CONFIG:Debug>:-Wignored-qualifiers>
      $<$<CONFIG:Debug>:-Wmissing-braces>
      $<$<CONFIG:Debug>:-Wreturn-type>
      $<$<CONFIG:Debug>:-Wswitch>
      $<$<CONFIG:Debug>:-Wmissing-field-initializers>
      $<$<CONFIG:Debug>:-Wno-unknown-pragmas>
      $<$<CONFIG:Debug>:-Wsign-compare>
    )

  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # MSVC-specific options
    target_compile_options(${TARGET_NAME} PRIVATE
      $<$<CONFIG:Release>:/O2>
      $<$<CONFIG:Debug>:/Od>
      $<$<CONFIG:Debug>:/W4>
      $<$<CONFIG:Debug>:/Zi>
    )
  endif()
endfunction()

# Set module path in order to use custom CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/source)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Building project in mode: ${CMAKE_BUILD_TYPE}")

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
