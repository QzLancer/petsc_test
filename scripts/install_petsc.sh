#!/bin/bash
# install_petsc.sh - Script to install PETSc with Hypre support on Linux
#
# This script downloads and builds PETSc with the configuration needed
# for the petsc_hypre_ams_gmres example.
#
# Usage:
#   chmod +x install_petsc.sh
#   ./install_petsc.sh [install_dir]
#
# Default install directory: $HOME/petsc

set -e

# Installation directory
INSTALL_DIR="${1:-$HOME/petsc}"

echo "=========================================="
echo "PETSc Installation Script"
echo "Installation directory: $INSTALL_DIR"
echo "=========================================="

# Check for required dependencies
echo "Checking dependencies..."

check_command() {
    if ! command -v $1 &> /dev/null; then
        echo "ERROR: $1 is not installed."
        echo "Please install it first:"
        echo "  Ubuntu/Debian: sudo apt-get install $2"
        echo "  CentOS/RHEL:   sudo yum install $3"
        exit 1
    fi
    echo "  Found: $1"
}

check_command gcc "build-essential" "gcc"
check_command g++ "build-essential" "gcc-c++"
check_command gfortran "gfortran" "gcc-gfortran"
check_command python3 "python3" "python3"
check_command make "make" "make"
check_command git "git" "git"

# Check for MPI
if command -v mpicc &> /dev/null; then
    echo "  Found: MPI (mpicc)"
    USE_DOWNLOAD_MPICH=""
else
    echo "  MPI not found, will download MPICH"
    USE_DOWNLOAD_MPICH="--download-mpich"
fi

# Check for BLAS/LAPACK using pkg-config or ldconfig
if pkg-config --exists blas 2>/dev/null || ldconfig -p 2>/dev/null | grep -q libblas; then
    echo "  Found: BLAS"
    USE_DOWNLOAD_FBLASLAPACK=""
else
    echo "  BLAS not found, will download fblaslapack"
    USE_DOWNLOAD_FBLASLAPACK="--download-fblaslapack"
fi

echo ""
echo "Downloading PETSc..."

# Check if installation directory exists and is not empty
if [ -d "$INSTALL_DIR" ] && [ "$(ls -A "$INSTALL_DIR" 2>/dev/null)" ]; then
    echo "WARNING: Directory $INSTALL_DIR already exists and is not empty."
    read -p "Do you want to remove it and proceed? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Installation aborted."
        exit 1
    fi
    rm -rf "$INSTALL_DIR"
fi

# Create directory and clone
mkdir -p "$(dirname "$INSTALL_DIR")"
git clone -b release https://gitlab.com/petsc/petsc.git "$INSTALL_DIR"
cd "$INSTALL_DIR"

echo ""
echo "Configuring PETSc with complex scalar type and Hypre..."

# Configure PETSc
./configure \
    --with-cc=gcc \
    --with-cxx=g++ \
    --with-fc=gfortran \
    --with-scalar-type=complex \
    --download-hypre \
    $USE_DOWNLOAD_MPICH \
    $USE_DOWNLOAD_FBLASLAPACK \
    --with-debugging=0 \
    COPTFLAGS='-O3' \
    CXXOPTFLAGS='-O3' \
    FOPTFLAGS='-O3'

echo ""
echo "Building PETSc..."

# Build
make all

echo ""
echo "Testing PETSc installation..."

# Test
make check

# Get the PETSC_ARCH
PETSC_ARCH=$(ls -d arch-* | head -n1)

echo ""
echo "=========================================="
echo "PETSc installation completed successfully!"
echo "=========================================="
echo ""
echo "Add the following to your ~/.bashrc or ~/.profile:"
echo ""
echo "  export PETSC_DIR=$INSTALL_DIR"
echo "  export PETSC_ARCH=$PETSC_ARCH"
echo ""
echo "Then reload your shell or run:"
echo "  source ~/.bashrc"
echo ""
echo "To build the example, run:"
echo "  cd /path/to/petsc_test"
echo "  make"
echo "=========================================="
