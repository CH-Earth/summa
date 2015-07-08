#========================================================================
# Makefile to compile SUMMA
#========================================================================
#
# To troubleshoot your paths and setup, type 'make check'
#

#========================================================================
# PART 0: Define directory paths
#========================================================================

# Define core directory below which everything resides. This is the
# parent directory of the 'build' directory
F_MASTER =

# Core directory that contains source code
F_KORE_DIR = $(F_MASTER)/build/source

# Location of the compiled modules
MOD_PATH = $(F_MASTER)/build

# Define the directory for the executables
EXE_PATH = $(F_MASTER)/bin

#========================================================================
# PART 1: Assemble all of the SUMMA sub-routines
#========================================================================

# Define directories
DRIVER_DIR = $(F_KORE_DIR)/driver
HOOKUP_DIR = $(F_KORE_DIR)/hookup
NETCDF_DIR = $(F_KORE_DIR)/netcdf
DSHARE_DIR = $(F_KORE_DIR)/dshare
NUMREC_DIR = $(F_KORE_DIR)/numrec
NOAHMP_DIR = $(F_KORE_DIR)/noah-mp
ENGINE_DIR = $(F_KORE_DIR)/engine
#
# Numerical Recipes utilities
SUMMA_NRUTIL= \
        nrtype.f90 \
		nr_utility.f90
NRUTIL = $(patsubst %, $(NUMREC_DIR)/%, $(SUMMA_NRUTIL))

#
# Numerical recipes procedures
# NOTE: all numerical recipes procedures are now replaced with free versions
SUMMA_NRPROC= \
		expIntegral.f90 \
		spline_int.f90
NRPROC = $(patsubst %, $(NUMREC_DIR)/%, $(SUMMA_NRPROC))

# Hook-up modules (set files and directory paths)
SUMMA_HOOKUP= \
        summaFileManager.f90
HOOKUP = $(patsubst %, $(HOOKUP_DIR)/%, $(SUMMA_HOOKUP))


# Data modules
SUMMA_DATAMS= \
		var_lookup.f90 \
		data_struc.f90 \
		multiconst.f90
DATAMS = $(patsubst %, $(DSHARE_DIR)/%, $(SUMMA_DATAMS))

# utility modules
SUMMA_UTILMS= \
		time_utils.f90
UTILMS = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_UTILMS))

# Model guts
SUMMA_MODGUT= \
MODGUT = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_MODGUT))

# Solver
SUMMA_SOLVER= \
		vegPhenlgy.f90 \
		diagn_evar.f90 \
		groundwatr.f90 \
		vegSWavRad.f90 \
		vegNrgFlux.f90 \
		ssdNrgFlux.f90 \
		vegLiqFlux.f90 \
		snowLiqFlx.f90 \
		soilLiqFlx.f90 \
		systemSolv.f90 \
		coupled_em.f90
SOLVER = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_SOLVER))

# Define routines for SUMMA preliminaries
SUMMA_PRELIM= \
		conv_funcs.f90 \
		ascii_util.f90 \
		sunGeomtry.f90 \
		get_ixname.f90 \
		mDecisions.f90 \
		snow_utils.f90 \
		soil_utils.f90 \
		updatState.f90 \
		convE2Temp.f90 \
		allocspace.f90 \
		read_metad.f90 \
		ffile_info.f90 \
		read_attrb.f90 \
		read_pinit.f90 \
		pOverwrite.f90 \
		read_param.f90 \
		paramCheck.f90 \
		read_icond.f90
PRELIM = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_PRELIM))

SUMMA_NOAHMP= \
		module_model_constants.F \
		module_sf_noahutl.F \
		module_sf_myjsfc.F \
		module_sf_sfclay.F \
		module_sf_noahlsm.F \
		module_sf_noahmplsm.F
NOAHMP = $(patsubst %, $(NOAHMP_DIR)/%, $(SUMMA_NOAHMP))

# Define routines for the SUMMA model runs
SUMMA_MODRUN = \
		var_derive.f90 \
		read_force.f90 \
		derivforce.f90 \
		snowAlbedo.f90 \
		canopySnow.f90 \
		tempAdjust.f90 \
		snwCompact.f90 \
		layerMerge.f90 \
		layerDivide.f90 \
		volicePack.f90 \
		qTimeDelay.f90
MODRUN = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_MODRUN))

# Define routines for the solver
SUMMA_MSOLVE = \

# Define NetCDF routines
SUMMA_NETCDF = \
		def_output.f90 \
		modelwrite.f90
NETCDF = $(patsubst %, $(NETCDF_DIR)/%, $(SUMMA_NETCDF))

# ... stitch together common programs
COMM_ALL = $(NRUTIL) $(NRPROC) $(HOOKUP) $(DATAMS) $(UTILMS)

# ... stitch together SUMMA programs
SUMMA_ALL = $(PRELIM) $(MODRUN) $(SOLVER) $(NETCDF)

#========================================================================
# PART 2: Define the libraries, driver programs, and executables
#========================================================================

# Define the Fortran Compiler
FC  = gfortran-mp-4.8

# Define the NetCDF libraries and path to include files
ifeq "$(FC)" "gfortran"
 NCDF_PATH = /usr/local/netcdf-4.3.0/gcc-4.7.2
 LAPK_PATH = /usr
 # define the lapack libraries
 LIBLAPACK = -L$(LAPK_PATH)/lib -llapack -lblas -latlas
endif
ifeq "$(FC)" "gfortran-mp-4.8"
 NCDF_PATH = /opt/local
 LAPK_PATH = /opt/local
 # define the lapack libraries
 LIBLAPACK = -L$(LAPK_PATH)/lib -llapack -lblas -latlas
endif
ifeq "$(FC)" "ifort"
 NCDF_PATH = /opt/netcdf-4.3.0+ifort-12.1
 LAPK_PATH = /usr
 # define the lapack libraries
 LIBLAPACK = -L$(LAPK_PATH)/lib -llapack
endif
ifeq "$(FC)" "pgf90"
 NCDF_PATH = /usr/local/netcdf-4.3.0/pgi-13.8
 LAPK_PATH = /usr
 # define the lapack libraries
 LIBLAPACK = -L$(LAPK_PATH)/lib -llapack
endif
LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff
INCNETCDF = -I$(NCDF_PATH)/include

# Define the driver routine
SUMMA_DRIVER= \
	multi_driver.f90
DRIVER = $(patsubst %, $(DRIVER_DIR)/%, $(SUMMA_DRIVER))

# Define the executable
DRIVER__EX = summa.exe

#========================================================================
# PART 3: Compile the puppy
#========================================================================

# Define flags
ifeq "$(FC)" "gfortran"
 FLAGS_NOAH = -ffree-form -fdefault-real-8 -ffree-line-length-none -fmax-errors=0 -g -fbacktrace
 FLAGS_COMM = -g -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace
 FLAGS_SUMMA = -g -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace
endif
ifeq "$(FC)" "gfortran-mp-4.8"
 FLAGS_NOAH = -ffree-form -fdefault-real-8 -ffree-line-length-none -fmax-errors=0 -g -fbacktrace
 FLAGS_COMM = -g -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace
 FLAGS_SUMMA = -g -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace
endif
ifeq "$(FC)" "ifort"
 FLAGS_NOAH = -warn nounused -autodouble -noerror_limit -FR -O0 -auto -WB -traceback -g -fltconsistency
 FLAGS_COMM = -debug -warn all -check all -FR -O0 -auto -WB -traceback -g -fltconsistency -fpe0
 FLAGS_SUMMA = -debug -warn all -check all -FR -O0 -auto -WB -traceback -g -fltconsistency -fpe0
endif
ifeq "$(FC)" "pgf90"
 FLAGS1 = -Mbackslash -g -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=divz,inv -Mprof=lines,time
 FLAGS_NOAH = -g -Mr8 -Mfree
 FLAGS_COMM = -Mbackslash -g -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=fp -traceback
 FLAGS_SUMMA = -Mbackslash -g -Mbounds -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=divz,inv,ovf -traceback
endif

#FLAGS = -O3 -W -v

#.SUFFIXES: .f .o .f90

# make sure that the paths are defined. These are just some high level checks
ifndef F_MASTER
 $(error F_MASTER is undefined)
endif
ifndef FC
 $(error FC is undefined: Specify your compiler)
endif
ifndef FLAGS_SUMMA
 $(error Specify flags for your compiler: $(FC))
endif
ifndef NCDF_PATH
 $(error NCDF_PATH is undefined)
endif

# Compile
all: compile_noah compile_comm compile_summa link clean install

check:
	$(info)
	$(info Displaying make variables:)
	$(info   F_MASTER   : $(F_MASTER))
	$(info   EXE_PATH   : $(EXE_PATH))
	$(info   FC         : $(FC))
	$(info   NCDF_PATH  : $(NCDF_PATH))
	$(info   LIBNETCDF  : $(LIBNETCDF))
	$(info   INCNETCDF  : $(INCNETCDF))
	$(info   LAPK_PATH  : $(LAPK_PATH))
	$(info   LIBLAPACK  : $(LIBLAPACK))
	$(info   FLAGS_NOAH : $(FLAGS_NOAH))
	$(info   FLAGS_COMM : $(FLAGS_COMM))
	$(info   FLAGS_SUMMA: $(FLAGS_SUMMA))
	$(info Add $(NCDF_PATH)/lib to your LD_LIBRARY_PATH)
	$(info)

# compile Noah-MP routines
compile_noah:
	$(FC) $(FLAGS_NOAH) -c $(NOAHMP)

# compile common routines
compile_comm:
	$(FC) $(FLAGS_COMM) -c $(COMM_ALL)

# compile SUMMA routines
compile_summa:
	$(FC) $(FLAGS_SUMMA) -c $(SUMMA_ALL) $(DRIVER) \
	$(INCNETCDF)

# link routines
link:
	$(FC) *.o $(LIBNETCDF) $(LIBLAPACK) -o $(DRIVER__EX)

# Remove object files
clean:
	rm -f *.o
	rm -f *.mod
	rm -f soil_veg_gen_parm__genmod.f90

# Copy the executable to the bin directory
install:
	@mkdir -p $(EXE_PATH)
	@mv $(DRIVER__EX) $(EXE_PATH)
	$(info Add $(NCDF_PATH)/lib to your LD_LIBRARY_PATH)
	$(info $(DRIVER__EX) successfully installed in $(EXE_PATH))
