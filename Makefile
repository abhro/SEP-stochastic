
SHELL = /bin/bash

# Fortran compiler name
FC = gfortran

# options for the Fortran compiler
FFLAGS =

#FFLAGS += -O
FFLAGS += -fopenmp -ggdb3
FFLAGS += -cpp
FFLAGS += -fcheck=all
#FFLAGS += -Wall -Wno-unused-variable
FFLAGS += -fbacktrace
FFLAGS += -std=f2008
FFLAGS += -g -O0 -fvar-tracking -fopenmp
#FFLAGS += -g -O0 -ffpe-trap=zero,overflow,underflow -fopenmp

#LDFLAGS += -O
#LDFLAGS += -g -O0 -fvar-tracking
#LDFLAGS += -g -O0 -ffpe-trap=zero,overflow,underflow
#LDFLAGS += -O0 -g

# Options for compiling and linking with GNU Scientific Library (GSL) Fortran version
#FGSL_OPTS_COMP =
#FGSL_OPTS_LINK = -L/usr/local/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-9.4.0/gsl-2.7-x2znrtoned4nhaqpluk6oh4jma6tjixe/lib \
#            -lfgsl -lgsl -lgslcblas -lm

FFLAGS += $(FGSL_OPTS_COMP)
LDFLAGS += $(FGSL_OPTS_LINK)

PROGEXT = # .exe for windows
OBJEXT = .o

COMPILE.f90 = $(FC) $(FFLAGS) $(TARGET_ARCH) -c
LINK.f90 = $(FC) $(FFLAGS) $(LDFLAGS) $(TARGET_ARCH)

MODULE_DIR = modules
PROG_DIR = programs

# All of the Fortran "module" files, they live under the $(MODULE_DIR) directory
MODULE_SRCS      = $(wildcard $(MODULE_DIR)/*.f90)
MODULE_BASENAMES = $(patsubst $(MODULE_DIR)/%, %, $(patsubst %.f90, %, $(MODULE_SRCS)))
MODULE_OBJS      = $(addsuffix $(OBJEXT), $(MODULE_BASENAMES))
# module name must be same as file name
MODULE_MODS      = $(addsuffix .mod,$(MODULE_BASENAMES))

PROG_SRCS = $(wildcard $(PROG_DIR)/*.f90)
PROG_BASENAMES = $(patsubst $(PROG_DIR)/%, %, $(patsubst %.f90, %, $(PROG_SRCS)))
PROG_OBJS = $(addsuffix $(OBJEXT),                   $(PROG_BASENAMES))
PROG_BINS = $(addsuffix $(PROGEXT),                  $(PROG_BASENAMES))

# List of all source files
SRC       = $(MODULE_SRCS) $(PROG_SRCS)
BASENAMES = $(MODULE_BASENAMES) $(PROG_BASENAMES)
OBJS      = $(addsuffix $(OBJEXT), $(BASENAMES))

.SUFFIXES: .o .for .f .f90 .f95 .f03
.PHONY: modules executables all clean clean-logs clean-objs printvars
.DEFAULT_GOAL := all

$(PROG_BINS): $(OBJS)
	$(FC) $(OUTFLAG)$(PROG_BIN) $(FFLAGS) $(SRC)

modules: $(MODULE_OBJS) $(MODULE_MODS)

subroutines: $(SUBR_OBJS)


COMPILE.f90 = $(FC) $(FFLAGS) $(TARGET_ARCH) -c
LINK.f90 = $(FC) $(FFLAGS) $(LDFLAGS) $(TARGET_ARCH)

# how to compile any Fortran source into an object file
#
# the $(OUTPUT_OPTION) is copies to mimic make's internal
# database, even though it's never set anywhere in this file
.f90$(OBJEXT):
	$(COMPILE.f90) $(OUTPUT_OPTION) $<
%$(OBJEXT): %.f90
	$(COMPILE.f90) $< $(OUTPUT_OPTION)
%.mod %$(OBJEXT): modules/%.f90
	$(COMPILE.f90) $<
# For an object file, it needs to depend on its module file, because.... reasons
# https://stackoverflow.com/a/67591729/3174398
%$(OBJEXT): %.mod

# how to compile (link) any object file to make an executable, fortran specific
%: %$(OBJEXT)
	@#$(LINK.f90) $^ $(MODULE_OBJS) $(LDLIBS) -o $@
	@#$(LINK.f90) $^ $(MODULE_MODS) $(LDLIBS) -o $@
	$(LINK.f90) $^ $(LDLIBS) -o $@

all: $(MODULE_OBJS) executables

# these targets create executables
executables: $(MODULE_OBJS) $(PROG_BINS)

clean: clean-objs clean-modules
	$(RM) -vf $(PROG_BINS)

clean-logs:
	$(RM) -vf logs/*.out.txt logs/*.err.txt

clean-modules:
	$(RM) -vf $(MODULE_MODS)

clean-objs:
	$(RM) -vf $(OBJS)

# manual dependencies because automatic dependency generation isn't a thing
mtrx.o: param.o epv.o dxx.o dmumu.o file_op.o
dxx.o: dmumu.o
cme_cross.o: sim3d_utils.o datetime_utils.o
sim3d_utils.o: fb.o epv.o mtrx.o
fb.o: param.o

# this depedency is only relevant when making object files
# because they need the .mod files for function signatures
$(PROG_OBJS): $(MODULE_MODS)

# link together the object file with all the modules.
# this should ideally be done with automatic dependency
# generation, since not all programs use all modules.
# alas, we don't have that in fortran
$(PROG_BINS): $(MODULE_OBJS)

# for debugging the Makefile
printvars:
	@$(foreach V, $(sort $(.VARIABLES)), \
		$(if $(filter-out environment% default automatic, $(origin $V)), \
		$(info $V = $($V) ($(value $V)))))

