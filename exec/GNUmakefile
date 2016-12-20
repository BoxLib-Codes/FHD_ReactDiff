# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= ../../BoxLib
FHDSRC_HOME ?= ../src

# "not debug" - "t" builds an optimized code
NDEBUG    := t
# build an MPI executable
MPI       := t
# OpenMP not fully supported yet - only some routines are threaded
OMP       :=
# Fortran and C compilers
COMP      := gfortran
CCOMP     := gcc
# need this to compile bl_rng.f90 (module to store random number engines and distributions)
CXX11     = t

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

VPATH_LOCATIONS += ..

include $(FHDSRC_HOME)/src_reactDiff/GPackage.mak
VPATH_LOCATIONS += $(FHDSRC_HOME)/src_reactDiff

include $(FHDSRC_HOME)/src_common/GPackage.mak
VPATH_LOCATIONS += $(FHDSRC_HOME)/src_common

include $(FHDSRC_HOME)/src_random/GPackage.mak
VPATH_LOCATIONS += $(FHDSRC_HOME)/src_random

include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
