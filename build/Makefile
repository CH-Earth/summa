#========================================================================
# Makefile to compile the FUSE-style snow model 
#========================================================================
#
#
#
#========================================================================
# PART 0: Define directory paths
#========================================================================

# Define core directory below which everything resides
F_MASTER = /d1/mclark/FUSE_SNOW/

# Core directory that contains source code
F_KORE_DIR = $(F_MASTER)build/source/

# Location of the compiled modules
MOD_PATH = $(F_MASTER)build/

# Define the directory for the executables
EXE_PATH = $(F_MASTER)bin
#
#========================================================================
# PART 1: Assemble all of the FUSE sub-routines
#========================================================================

# Define directories
DRIVER_DIR = $(F_KORE_DIR)driver
HOOKUP_DIR = $(F_KORE_DIR)hookup
NETCDF_DIR = $(F_KORE_DIR)netcdf
DSHARE_DIR = $(F_KORE_DIR)dshare
NUMREC_DIR = $(F_KORE_DIR)numrec
ENGINE_DIR = $(F_KORE_DIR)engine
#
# Numerical Recipes utilities
FUSE_NRUTIL= \
        nrtype.f90 \
		nr_utility.f90
NRUTIL = $(patsubst %, $(NUMREC_DIR)/%, $(FUSE_NRUTIL))

#
# Numerical recipes procedures
FUSE_NRPROC= \
		spline_int.f90
NRPROC = $(patsubst %, $(NUMREC_DIR)/%, $(FUSE_NRPROC))

# Hook-up modules (set files and directory paths)
FUSE_HOOKUP= \
        snow_fileManager.f90
HOOKUP = $(patsubst %, $(HOOKUP_DIR)/%, $(FUSE_HOOKUP))


# Data modules
FUSE_DATAMS= \
		var_lookup.f90 \
		data_struc.f90 \
		multiconst.f90
DATAMS = $(patsubst %, $(DSHARE_DIR)/%, $(FUSE_DATAMS))

# Model guts
FUSE_MODGUT= \
MODGUT = $(patsubst %, $(ENGINE_DIR)/%, $(FUSE_MODGUT))

# Solver
FUSE_SOLVER= \
		tridagSolv.f90 \
		diagn_evar.f90 \
		phseChange.f90 \
		heatTransf.f90 \
		snowHydrol.f90 \
		soilHydrol.f90 \
		groundwatr.f90 \
        picardSolv.f90 \
        coupled_em.f90
SOLVER = $(patsubst %, $(ENGINE_DIR)/%, $(FUSE_SOLVER))

# Define routines for FUSE preliminaries
FUSE_PRELIM= \
		conv_funcs.f90 \
		ascii_util.f90 \
		get_ixname.f90 \
		mDecisions.f90 \
		snow_utils.f90 \
		soil_utils.f90 \
		ConvE2Temp.f90 \
		time_utils.f90 \
		allocspace.f90 \
		read_metad.f90 \
		ffile_info.f90 \
		read_pinit.f90 \
		read_param.f90 \
		read_icond.f90
PRELIM = $(patsubst %, $(ENGINE_DIR)/%, $(FUSE_PRELIM))

# Define routines for the FUSE model runs
FUSE_MODRUN = \
		var_derive.f90 \
		read_force.f90 \
		derivforce.f90 \
		newsnwfall.f90 \
		layerMerge.f90 \
		surfAlbedo.f90 \
		snwDensify.f90
MODRUN = $(patsubst %, $(ENGINE_DIR)/%, $(FUSE_MODRUN))

# Define routines for the solver
FUSE_MSOLVE = \
		

# Define NetCDF routines
FUSE_NETCDF = \
		def_output.f90 \
		modelwrite.f90
NETCDF = $(patsubst %, $(NETCDF_DIR)/%, $(FUSE_NETCDF))

# ... and stitch it all together...
FUSE_ALL = $(UTILMS) $(NRUTIL) $(NRPROC) $(HOOKUP) $(DATAMS) \
		$(MODGUT) $(PRELIM) $(MODRUN) $(SOLVER) \
        $(NETCDF)

#========================================================================
# PART 2: Define the libraries, driver programs, and executables
#========================================================================

# Define the Fortran Compiler
FC  = ifort

# Define the NetCDF libraries and path to include files
ifeq "$(FC)" "gfortran"
 LIBNETCDF = -L/opt/netcdf-4.1.3/gcc-4.4.5/lib -lnetcdff -lnetcdf -L/opt/hdf5-1.8.7/lib -lhdf5_hl -lhdf5 -L/opt/lib -lsz -lz
 INCNETCDF = -I/opt/netcdf-4.1.3/gcc-4.4.5/include
endif
ifeq "$(FC)" "ifort"
 HDF_PATH  = /usr/local/hdf5-1.8.7
 NCDF_PATH = /usr/local/netcdf-4.1.3/ifort-12.0.5
 LOCL_PATH = /usr/local
 LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf -L$(HDF_PATH)/lib -lhdf5_hl -lhdf5 -L$(LOCL_PATH)/lib -lsz -lz
 INCNETCDF = -I$(NCDF_PATH)/include -I/$(HDF_PATH)/include
endif
ifeq "$(FC)" "pgf90"
 LIBNETCDF = -L/usr/local/netcdf4-pgi/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz
 INCNETCDF = -I/usr/local/netcdf4-pgi/include
endif

# Define the driver routine
FUSE_DRIVER= \
	multi_driver.f90
DRIVER = $(patsubst %, $(DRIVER_DIR)/%, $(FUSE_DRIVER))

# Define the executable
DRIVER__EX = fuse_snow.exe

#========================================================================
# PART 3: Compile the puppy
#========================================================================

# Define flags
ifeq "$(FC)" "gfortran"
 FLAGS = -static -Wall -g -ffree-line-length-none
endif
ifeq "$(FC)" "ifort"
 FLAGS1 = -warn all -noerror_limit -FR -O0 -auto -WB -traceback -g -fltconsistency
 FLAGS = -static -debug -warn all -check all -FR -O0 -auto -WB -traceback -g -fltconsistency -fpe0
endif
ifeq "$(FC)" "pgf90"
 FLAGS1 = -Bstatic -Mbackslash -g -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=divz,inv -Mprof=lines,time
 FLAGS = -Bstatic -Mbackslash -g -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=divz,inv
endif

#FLAGS = -O3 -W -v

#.SUFFIXES: .f .o .f90

# Compile
all: check compile install clean

check:
	echo test
	echo $(FLAGS)
	echo $(FC)

# compile FUSE into a static library
compile:
	$(FC) $(FLAGS) $(FUSE_ALL) $(DRIVER) \
	$(LIBNETCDF) $(INCNETCDF) -o $(DRIVER__EX)

# Remove object files
clean:
						rm -f *.o
						rm -f *.mod
#
# Copy the executable to the bin directory
install:
		mv $(DRIVER__EX) $(EXE_PATH)
