.SUFFIXES:
ifdef USEMPE
USEMPI=yes
endif

LINKFLAGS = #-lefence

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
	xbeach.o 

SRCS = $(patsubst %.o,%.F90,$(OBJS))
#	wave_dist.o \

$(program): $(OBJS)
	$(F90) $(F90FLAGS) -o $@ $(OBJS) $(LINKFLAGS) 
	@echo Do not forget to make install

testgenmoduleobjs = testgenmodule.o xmpi.o general_mpi.o 

testdetsubobjs = testdetsub.o general_mpi.o

# mpe_mpilog
# mpe_mpitrace
# mpe_mpianim
# mpe_mpicheck
testgenmodule: $(testgenmoduleobjs)
	$(F90) $(F90FLAGS) -o $@ $(testgenmoduleobjs) $(LINKFLAGS)

testdetsub: $(testdetsubobjs)
	$(F90) $(F90FLAGS) -o $@ $(testdetsubobjs) $(LINKFLAGS)

install: all
	mkdir -p ../bin
	cp $(program) ../bin


testje:
	@echo OBJS: $(OBJS) 
	@echo USEMPI: $(USEMPI)

%.o: %.F90
	@echo compiling $<
	$(F90) -c $(F90FLAGS) $<

%.mod: 
	@echo compiling $< to create $@
	@$(F90) -c $(F90FLAGS) $<

#
# separate rule for math_tools: 
# -O2 and -Wall gives spurious warnings with gfortran
#
math_tools.o: math_tools.F90
	@echo compiling $<
	$(F90) -c $(filter-out -Wall,$(F90FLAGS)) math_tools.F90

math_tools.mod:
	@echo compiling $< to create $@
	@$(F90) -c $(F90FLAGS) math_tools.F90

DEPENDENCIES dep:
	./makedepo $(SRCS) testgenmodule.F90 > DEPENDENCIES

include DEPENDENCIES

F90FLAGS+=-g -O2 -I. # -pg # -fprofile-arcs -ftest-coverage
F90:=gfortran
ifeq ($(F90),gfortran)
  F90FLAGS += -Wall
endif
ifdef USEMPI
  F90:=mpif90
ifdef USETAU
  F90:=tau_f90.sh
endif
  F90FLAGS +=  -DUSEMPI
endif

clean:
	rm -f *.o *.mod $(program) core testgenmodule

depclean: clean
	rm -f DEPENDENCIES
show:
	@echo F90: $(F90)
	@echo F90FLAGS: $(F90FLAGS)
	@echo USEMPI: $(USEMPI)
