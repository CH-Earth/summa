#========================================================================
# PART 0: Define directory paths
#========================================================================

# Define core directory below which everything resides. This is the
# parent directory of the 'build' directory
F_MASTER = /home/mclark/summa

# Core directory that contains source code
F_KORE_DIR = $(F_MASTER)/build/source

# Location of the compiled modules
MOD_PATH = $(F_MASTER)/build/source/lapack

# Define the directory for the executables
EXE_PATH = $(F_MASTER)/build/source/lapack

#========================================================================
# PART 1: Assemble all of the sub-routines
#========================================================================

# Define directories
NUMREC_DIR = $(F_KORE_DIR)/numrec
LAPACK_DIR = $(F_KORE_DIR)/lapack

#
# Numerical Recipes utilities
SUMMA_NRUTIL= \
        nrtype.f90
NRUTIL = $(patsubst %, $(NUMREC_DIR)/%, $(SUMMA_NRUTIL))

# numerical recipes routines
SUMMA_NUMREC= \
		luSolv_numrec.f90
NUMREC = $(patsubst %, $(LAPACK_DIR)/%, $(SUMMA_NUMREC))

# ... stitch together all programs
LUTEST = $(NRUTIL) $(NUMREC)

#========================================================================
# PART 2: Define the libraries, driver programs, and executables
#========================================================================

# Define the Fortran Compiler
FC  = ifort

# Define the lapack library
LAPK_PATH = /usr
LIBLAPACK = -L$(LAPK_PATH)/lib -llapack

DRIVER = test_lusolve.f90

# Define the executable
DRIVER__EX = test_lusolve.exe


#========================================================================
# PART 3: Compile the puppy
#========================================================================

# Define flags
FLAGS = -debug -warn all -check all -FR -O0 -auto -WB -traceback -g -fltconsistency -fpe0

all: compile link clean

compile:
	$(FC) $(FLAGS) -c $(LUTEST) $(DRIVER)

link:
	$(FC) *.o $(LIBLAPACK) -o $(DRIVER__EX)

clean:
	rm -f *.o
	rm -f *.mod
