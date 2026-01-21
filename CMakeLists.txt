cmake_minimum_required(VERSION 3.16)
project(petsc_hypre_ams_gmres C)

# Set C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find PETSc
# Users should set PETSC_DIR and PETSC_ARCH environment variables
if(NOT DEFINED ENV{PETSC_DIR})
    message(FATAL_ERROR "PETSC_DIR environment variable is not set. Please set it to your PETSc installation directory.")
endif()

set(PETSC_DIR $ENV{PETSC_DIR})
set(PETSC_ARCH $ENV{PETSC_ARCH})

# Include PETSc configuration
if(EXISTS "${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables")
    # PETSc with PETSC_ARCH
    set(PETSC_CONF_DIR "${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf")
elseif(EXISTS "${PETSC_DIR}/lib/petsc/conf/petscvariables")
    # PETSc without PETSC_ARCH (prefix install)
    set(PETSC_CONF_DIR "${PETSC_DIR}/lib/petsc/conf")
else()
    message(FATAL_ERROR "Cannot find PETSc configuration files. Check PETSC_DIR and PETSC_ARCH.")
endif()

# Parse PETSc configuration
file(READ "${PETSC_CONF_DIR}/petscvariables" PETSC_VARS_CONTENT)

# Extract compiler and flags from petscvariables
string(REGEX MATCH "PETSC_CC_INCLUDES = [^\n]+" PETSC_CC_INCLUDES_LINE "${PETSC_VARS_CONTENT}")
string(REGEX REPLACE "PETSC_CC_INCLUDES = " "" PETSC_CC_INCLUDES "${PETSC_CC_INCLUDES_LINE}")

# Find PETSc include directories
if(PETSC_ARCH)
    set(PETSC_INCLUDE_DIRS
        "${PETSC_DIR}/include"
        "${PETSC_DIR}/${PETSC_ARCH}/include"
    )
    set(PETSC_LIB_DIR "${PETSC_DIR}/${PETSC_ARCH}/lib")
else()
    set(PETSC_INCLUDE_DIRS "${PETSC_DIR}/include")
    set(PETSC_LIB_DIR "${PETSC_DIR}/lib")
endif()

# Find MPI
find_package(MPI REQUIRED)

# Create executable
add_executable(petsc_hypre_ams_gmres petsc_hypre_ams_gmres_full.c)

# Include directories
target_include_directories(petsc_hypre_ams_gmres PRIVATE
    ${PETSC_INCLUDE_DIRS}
    ${MPI_C_INCLUDE_PATH}
)

# Link directories
target_link_directories(petsc_hypre_ams_gmres PRIVATE ${PETSC_LIB_DIR})

# Link libraries
target_link_libraries(petsc_hypre_ams_gmres
    petsc
    ${MPI_C_LIBRARIES}
    m
)


# MPI compile options
if(MPI_C_COMPILE_FLAGS)
    target_compile_options(petsc_hypre_ams_gmres PRIVATE ${MPI_C_COMPILE_FLAGS})
endif()

# Platform-specific settings
if(WIN32)
    # Windows-specific settings
    target_compile_definitions(petsc_hypre_ams_gmres PRIVATE _USE_MATH_DEFINES)
elseif(UNIX)
    # Unix-specific settings
    target_link_libraries(petsc_hypre_ams_gmres dl)
endif()

# Print configuration summary
message(STATUS "=== PETSc AMS GMRES Configuration ===")
message(STATUS "PETSC_DIR: ${PETSC_DIR}")
message(STATUS "PETSC_ARCH: ${PETSC_ARCH}")
message(STATUS "PETSC_INCLUDE_DIRS: ${PETSC_INCLUDE_DIRS}")
message(STATUS "PETSC_LIB_DIR: ${PETSC_LIB_DIR}")
message(STATUS "MPI_C_INCLUDE_PATH: ${MPI_C_INCLUDE_PATH}")
message(STATUS "MPI_C_LIBRARIES: ${MPI_C_LIBRARIES}")
message(STATUS "=====================================")
