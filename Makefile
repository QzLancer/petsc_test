# Makefile for PETSc Hypre AMS GMRES Example
#
# Prerequisites:
#   - PETSc installed with Hypre support
#   - PETSC_DIR and PETSC_ARCH environment variables set
#
# Usage:
#   make          - Build the executable
#   make run      - Run with default options
#   make clean    - Clean build artifacts
#   make help     - Show this help

# Check if PETSC_DIR is set
ifndef PETSC_DIR
    $(error PETSC_DIR is not set. Please set it to your PETSc installation directory)
endif

# Include PETSc configuration
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# Target executable
TARGET = petsc_hypre_ams_gmres

# Source files
SOURCES = petsc_hypre_ams_gmres_full.c

# Object files
OBJECTS = $(SOURCES:.c=.o)

# Default target
all: $(TARGET)

# Link the executable using PETSc's standard method
$(TARGET): $(OBJECTS)
	${CLINKER} -o $@ $^ ${PETSC_KSP_LIB}
	${RM} $(OBJECTS)

# Run the executable
run: $(TARGET)
	${MPIEXEC} -n 1 ./$(TARGET)

# Run with multiple processes
run-mpi: $(TARGET)
	${MPIEXEC} -n 4 ./$(TARGET)

# Clean build artifacts
clean::
	${RM} $(TARGET) $(OBJECTS) *.o

# Help target
help:
	@echo "PETSc Hypre AMS GMRES Example - Build System"
	@echo ""
	@echo "Prerequisites:"
	@echo "  - PETSc installed with Hypre support"
	@echo "  - PETSC_DIR environment variable set"
	@echo "  - PETSC_ARCH environment variable set (if applicable)"
	@echo ""
	@echo "Available targets:"
	@echo "  make          - Build the executable"
	@echo "  make run      - Run with single process"
	@echo "  make run-mpi  - Run with 4 MPI processes"
	@echo "  make clean    - Clean build artifacts"
	@echo "  make help     - Show this help"
	@echo ""
	@echo "Example usage:"
	@echo "  export PETSC_DIR=/path/to/petsc"
	@echo "  export PETSC_ARCH=arch-linux-c-debug"
	@echo "  make"
	@echo "  make run"

.PHONY: all run run-mpi clean help
