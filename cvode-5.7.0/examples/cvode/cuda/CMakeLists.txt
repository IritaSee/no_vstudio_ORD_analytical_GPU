# -----------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody Balos @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
# CMakeLists.txt for CVODE CUDA examples.
#
# This file is generated from a template using variables
# set at configuration time. It can be used as a template for
# other user CMakeLists configuration files.
# -----------------------------------------------------------------

# Set the minimum required cmake version
cmake_minimum_required(VERSION 3.10)

# Set cache variables for compilers and flags
set(CMAKE_C_COMPILER
  C:/Program Files (x86)/Microsoft Visual Studio/2019/Community/VC/Tools/MSVC/14.28.29333/bin/Hostx64/x64/cl.exe
  CACHE FILEPATH "C compiler")

set(CMAKE_C_FLAGS
  "/DWIN32 /D_WINDOWS /W3"
  CACHE STRING "C compiler flags")

set(CMAKE_CXX_COMPILER
  C:/Program Files (x86)/Microsoft Visual Studio/2019/Community/VC/Tools/MSVC/14.28.29333/bin/Hostx64/x64/cl.exe
  CACHE FILEPATH "C++ compiler")

set(CMAKE_CXX_FLAGS
  "/DWIN32 /D_WINDOWS /W3 /GR /EHsc"
  CACHE STRING "C++ compiler flags")

set(CMAKE_CUDA_COMPILER
  C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.1/bin/nvcc.exe
  CACHE FILEPATH "CUDA compiler")

set(CMAKE_CUDA_HOST_COMPILER
  C:/Program Files (x86)/Microsoft Visual Studio/2019/Community/VC/Tools/MSVC/14.28.29333/bin/Hostx64/x64/cl.exe
  CACHE FILEPATH "CUDA host compiler")

# Specify project name and languages
project(CVODE_cuda_examples C CXX CUDA)

# Need this as long as CUDA libraries like cuSOLVER are not
# available through some other way.
find_package(CUDA REQUIRED)

# Enable testing
include(CTest)

# ------------------------------------------------------------------------------

# Specify the path to SUNDIALS header files
set(SUNDIALS_INCLUDE_DIR
  C:/Users/Marcell-CML/Documents/Visual Studio 2019/Libraries/cvode-5.7.0/include
  CACHE PATH "Location of SUNDIALS header files")

# Specify the path to SUNDIALS libraries
set(SUNDIALS_LIBRARY_DIR
  C:/Users/Marcell-CML/Documents/Visual Studio 2019/Libraries/cvode-5.7.0/lib
  CACHE PATH "Location of SUNDIALS libraries")

# Find the SUNDIALS libraries
find_library(SUNDIALS_SOLVER_LIB
  sundials_cvode ${SUNDIALS_LIBRARY_DIR}
  DOC "CVODE library")

find_library(SUNDIALS_NVEC_LIB
  sundials_nveccuda ${SUNDIALS_LIBRARY_DIR}
  DOC "NVECTOR_CUDA library")

find_library(SUNDIALS_CUSPARSEMAT_LIB
  sundials_sunmatrixcusparse ${SUNDIALS_LIBRARY_DIR}
  DOC "SUNMATRIX_CUSPARSE library")

find_library(SUNDIALS_CUSOLVERSP_LIB
  sundials_sunlinsolcusolversp
  ${SUNDIALS_LIBRARY_DIR}
  DOC "SUNLINSOL_CUSOLVERSP linear solver library")

# For SUNDIALS module examples the solver library is not needed
if(NOT SUNDIALS_SOLVER_LIB)
  set(SUNDIALS_SOLVER_LIB "")
endif()

# For some examples the SUNDIALS_CUSPARSEMAT_LIB library is not needed
if(NOT SUNDIALS_CUSPARSEMAT_LIB)
  set(SUNDIALS_CUSPARSEMAT_LIB "")
endif()

# For some examples the SUNLINSOL_CUSOLVERSP library is not needed
if(NOT SUNDIALS_CUSOLVERSP_LIB)
  set(SUNDIALS_CUSOLVERSP_LIB "")
endif()

# Set additional libraries
set(SUNDIALS_EXTRA_LIBS -L${SUNDIALS_LIBRARY_DIR}  CACHE STRING "Additional libraries")

# List of SUNDIALS libraries
set(SUNDIALS_LIBRARIES
  -L${SUNDIALS_LIBRARY_DIR}
  ${SUNDIALS_SOLVER_LIB}
  ${SUNDIALS_NVEC_LIB}
  ${SUNDIALS_CUSPARSEMAT_LIB}
  ${SUNDIALS_EXTRA_LIBS})

# ------------------------------------------------------------------------------

# Set the names of the examples to be built and their dependencies
set(examples  cvAdvDiff_kry_cuda cvAdvDiff_kry_cuda_managed cvAdvDiff_diag_cuda)
set(examples_dependencies )
if(examples)
  list(REMOVE_DUPLICATES examples)
endif()

# Create targets for each example
foreach(example ${examples})

  # example source files
  add_executable(${example} ${example}.cu ${examples_dependencies})

  # set c++11
  target_compile_features(${example} PRIVATE cxx_std_11)

  # directories to include
  target_include_directories(${example} PRIVATE ${SUNDIALS_INCLUDE_DIR})

  # libraries to link against
  target_link_libraries(${example} ${SUNDIALS_LIBRARIES})

  # add the example to ctest
  add_test(NAME ${example} COMMAND ${example})

endforeach()

# ------------------------------------------------------------------------------

# Set the names of the examples to be built and their dependencies
set(examples_cusolver )
if(examples_cusolver)
  list(REMOVE_DUPLICATES examples_cusolver)
endif()

# List of SUNDIALS libraries
set(SUNDIALS_LIBRARIES
  ${SUNDIALS_SOLVER_LIB}
  ${SUNDIALS_NVEC_LIB}
  ${SUNDIALS_CUSPARSEMAT_LIB}
  ${SUNDIALS_CUSOLVERSP_LIB}
  ${SUNDIALS_EXTRA_LIBS})

# Create targets for each example
foreach(example ${examples_cusolver})

  # example source files
  add_executable(${example} ${example}.cu ${examples_dependencies})

  # set c++11
  target_compile_features(${example} PRIVATE cxx_std_11)

  # directories to include
  target_include_directories(${example} PRIVATE ${SUNDIALS_INCLUDE_DIR})

  # libraries to link against
  target_link_libraries(${example} ${SUNDIALS_LIBRARIES})
  target_link_libraries(${example} ${CUDA_cusolver_LIBRARY})
  target_link_libraries(${example} ${CUDA_cusparse_LIBRARY})

  # add the example to ctest
  add_test(NAME ${example} COMMAND ${example})

endforeach()
