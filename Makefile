
SHELL = /bin/bash

ifeq ($(toolset),pgi)
  #MPIBIN = /usr/local/mpich/1.2.5/ip/up/pgi/ssh/bin/
  # Requires the pgi module
  FC = $(MPIBIN)pgf95

  #FFLAGS += -ansi -O -w4
  #FFLAGS += -mp
  FFLAGS += -g -mp -mcmodel=medium

  #LDFLAGS += -Mmpi=mpich -g -O0
  #LDFLAGS += -O0 -g -mp

# default assume gcc
else
  MPIBIN =
  FC = $(MPIBIN)gfortran
  #FC = f95

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

endif

FGSL_OPTS_COMP =

FFLAGS += $(FGSL_OPTS_COMP)

FGSL_OPTS_LINK = -L/usr/local/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-9.4.0/gsl-2.7-x2znrtoned4nhaqpluk6oh4jma6tjixe/lib \
            -lfgsl -lgsl -lgslcblas -lm

LDFLAGS += $(FGSL_OPTS_LINK)

PROGEXT = # .exe for windows

COMPILE.f90 = $(FC) $(FFLAGS) $(TARGET_ARCH) -c
LINK.f90 = $(FC) $(FFLAGS) $(LDFLAGS) $(TARGET_ARCH)

MODULE_SRC_BASENAMES = param file_op epv fb dmumu dxx random datetime_utils \
                       mtrx loadptcl sim3d_utils cme_cross rksolvers
MODULE_SRCS = $(addprefix modules/, $(addsuffix .f90, $(MODULE_SRC_BASENAMES)))
MODULE_OBJS = $(addsuffix .o,                         $(MODULE_SRC_BASENAMES))
MODULE_MODS = $(addsuffix .mod,                       $(MODULE_SRC_BASENAMES))

PROG_SRC_BASENAMES = 2damrhistss maggrid mapb2s combiner seedgen sim3d shockfront
#PROG_SRC_BASENAMES += sim3d_em
PROG_SRCS = $(addprefix programs/, $(addsuffix .f90, $(PROG_SRC_BASENAMES)))
PROG_OBJS = $(addsuffix .o,                          $(PROG_SRC_BASENAMES))
PROG_BINS = $(addsuffix $(PROGEXT),                  $(PROG_SRC_BASENAMES))

.SUFFIXES: .o .for .f .f90 .f95 .f03
.PHONY: modules executables all clean clean-logs clean-objs
.DEFAULT_GOAL := all

# how to compile any fortran source into an object file
#
# the $(OUTPUT_OPTION) is copies to mimic make's internal
# database, even though it's never set anywhere in this file
.f90.o:
	$(COMPILE.f90) $(OUTPUT_OPTION) $<
%.o: %.f90
	$(COMPILE.f90) $< $(OUTPUT_OPTION)
%.mod: %.f90
	$(COMPILE.f90) $<

# how to compile (link) any object file to make an executable, fortran specific
%: %.o
	@#$(LINK.f90) $^ $(MODULE_OBJS) $(LDLIBS) -o $@
	@#$(LINK.f90) $^ $(MODULE_MODS) $(LDLIBS) -o $@
	$(LINK.f90) $^ $(LDLIBS) -o $@

all: $(MODULE_OBJS) executables

# these targets create executables
executables: $(MODULE_OBJS) $(PROG_BINS)

clean: clean-objs clean-modules
	rm -vf $(PROG_BINS)

clean-logs:
	rm -vf logs/*.out.txt logs/*.err.txt

clean-modules:
	rm -vf *.mod

clean-objs:
	rm -vf *.o

# these targets create modules: object files (.o) and module files (.mod)
modules: $(MODULE_OBJS) $(MODULE_MODS)

mtrx.o: param.o epv.o dxx.o dmumu.o sim3d_utils.o
dxx.o: dmumu.o

# this depedency is only relevant when making object files
# because they need the .mod files for function signatures
$(PROG_OBJS): $(MODULE_MODS)

# link together the object file with all the modules.
# this should ideally be done with automatic dependency
# generation, since not all programs use all modules.
# alas, we don't have that in fortran
$(PROG_BINS): $(MODULE_OBJS)
