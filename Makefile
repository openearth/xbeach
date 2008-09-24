.SUFFIXES:
# 
# external variables that control the working of this makefile:
# 
# USEMPI:  if defined, make an MPI version of the program: xbeach.mpi
#          else make a serial version: xbeach
# USEXLF:  if defined: use xlf compiler (IBM)
#          else use gfortran
# TESTING: if defined: compile with -O3
#          else compile with -O0
# 
ifdef USEMPE
USEMPI=yes
endif

LINKFLAGS =  #-lefence

ifdef USEMPI
  program = xbeach.mpi
  ifdef USEMPE
    LINKFLAGS += -lmpe_f2cmpi -llmpe -lmpe #-mpe=mpilog
    F90FLAGS += -DUSEMPE
  endif
else
  program = xbeach
endif
all: $(program)

OBJS:= boundaryconditions.o \
	constants.o \
	flow_timestep.o \
	initialize.o \
	interp.o \
	math_tools.o \
	morphevolution.o \
	params.o \
	readkey.o \
	readtide.o \
	roelvink.o \
	spaceparams.o \
	waveparams.o \
	wave_stationary.o \
	wave_timestep.o \
	xmpi.o \
	general_mpi.o \
	varoutput.o \
	timestep.o \
	xbeach.o \
        mnemonic.o

SRCS = $(patsubst %.o,%.F90,$(OBJS))
#	wave_dist.o \

$(program): $(OBJS)
	$(F90) $(F90FLAGS) -o $@ $(OBJS) $(LINKFLAGS) 
	@echo Do not forget to make install

testgenmoduleobjs = testgenmodule.o xmpi.o general_mpi.o 

testdetsubobjs = testdetsub.o general_mpi.o

demoobjs = demo.o spaceparams.o xmpi.o mnemonic.o readkey.o

# mpe_mpilog
# mpe_mpitrace
# mpe_mpianim
# mpe_mpicheck
testgenmodule: $(testgenmoduleobjs)
	$(F90) $(F90FLAGS) -o $@ $(testgenmoduleobjs) $(LINKFLAGS)

testdetsub: $(testdetsubobjs)
	$(F90) $(F90FLAGS) -o $@ $(testdetsubobjs) $(LINKFLAGS)

demo: $(demoobjs)
	$(F90) $(F90FLAGS) -o $@ $(demoobjs) $(LINKFLAGS)

install: all
	mkdir -p ../bin
	cp $(program) ../bin

testje:
	@echo OBJS: $(OBJS) 
	@echo USEMPI: $(USEMPI)

%.o: %.F90
	@echo compiling $<
	$(F90) -c $(F90FLAGS) $<

%.gen: spaceparams.tmpl makeincludes
	@echo $@ | ./makeincludes

%.mod: 
	@echo compiling $< to create $@
	@$(F90) -c $(F90FLAGS) $<

#
# separate rule for math_tools and varoutput: 
# -O2 and -Wall gives spurious warnings with gfortran
#
math_tools.o: math_tools.F90
	@echo compiling $<
	$(F90) -c $(F90FLAGSNOWALL) math_tools.F90

math_tools.mod:
	@echo compiling $< to create $@
	@$(F90) -c $(F90FLAGSNOWALL) math_tools.F90

varoutput.o: varoutput.F90
	@echo compiling $<
	$(F90) -c $(F90FLAGSNOWALL) varoutput.F90

outputmod.mod:
	@echo compiling $< to create $@
	@$(F90) -c $(F90FLAGSNOWALL) varoutput.F90

DEPENDENCIES dep:
	awk -f ./makedepo $(SRCS) testgenmodule.F90 demo.F90 > DEPENDENCIES

s.ind: space_ind.gen

s.inp: space_inp.gen

include DEPENDENCIES

makeincludes: makeincludes.F90
	$(F90_NO_MPI) $(F90FLAGS) -o $@ makeincludes.F90
	
# some mpi.mod files have a definition of MPI_Wtime,
# some don't. Comment next line out if your mpi.mod 
# lacks a definition

F90FLAGS += -DHAVE_MPI_WTIME

ifdef TESTING
  OPT=-O0
else
  OPT=-O3
endif
F90FLAGS+=-g $(OPT) -I.
ifdef USEXLF
  F90=xlf_r
else
  F90:=gfortran
endif
F90_NO_MPI := $(F90)
ifeq ($(F90),gfortran)
  F90FLAGS += -Wall
endif
ifdef USEMPI
  ifdef USEXLF
    F90:=mpfort
  else
    F90:=mpif90
  endif
  F90FLAGS +=  -DUSEMPI
endif

ifeq "$(filter -O3,$(F90FLAGS))" "-O3"
  F90FLAGSNOWALL=$(filter-out -Wall,$(F90FLAGS))
  else
  F90FLAGSNOWALL=$(F90FLAGS)
endif

ifdef USEXLF
  comma=,
  F90FLAGS := $(subst -D,-WF$(comma)-D,$(F90FLAGSNOWALL)) -qfree=f90
endif

clean:
	rm -f *.o *.mod core testgenmodule demo a.out

realclean: clean
	rm -f DEPENDENCIES *.gen makeincludes xbeach xbeach.mpi tags

show:
	@echo F90: $(F90)
	@echo F90FLAGS: $(F90FLAGS)
	@echo USEMPI: $(USEMPI)
	@echo F90FLAGSNOWALL :$(F90FLAGSNOWALL)

