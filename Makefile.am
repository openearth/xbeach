bin_PROGRAMS=makeincludes xbeach


ACLOCAL_AMFLAGS = -I config -I m4
makeincludes_SOURCES= \
	makeincludes.F90
makeincludes_LDFLAGS=

xbeach_SOURCES= \
	s.ind \
	s.inp \
	xmpi.F90 \
	general_mpi.F90 \
	readkey.F90 \
	mnemonic.F90 \
	interp.F90 \
	math_tools.F90 \
	params.F90 \
	spaceparams.F90 \
	waveparams.F90 \
	constants.F90 \
	roelvink.F90 \
	wave_timestep.F90 \
	wave_stationary.F90 \
        boundaryconditions.F90 \
	morphevolution.F90 \
	groundwater.F90 \
	drifters.F90 \
	readtide.F90 \
	readwind.F90 \
	initialize.F90 \
	timestep.F90 \
	flow_timestep.F90 \
	varianceupdate.F90 \
	varoutput.F90 \
	$(extra_sources) \
	xbeach.F90 



BUILT_SOURCES= \
	spacedecl.gen	\
	mnemonic.gen \
	indextos.gen \
	space_alloc_scalars.gen \
	space_alloc_arrays.gen \
	space_ind.gen \
	space_inp.gen \
	chartoindex.gen \
	s.ind \
	s.inp 

s.ind: space_ind.gen
s.inp: space_inp.gen
# this is not a proper way
%.gen: spaceparams.tmpl makeincludes
	echo $@ | ./makeincludes

check_SCRIPTS= \
	./testbed/tools/testT01.sh \
	./testbed/tools/testrun01.sh
TESTS=./testbed/tools/testT01.sh ./testbed/tools/testrun01.sh

if USEMPI
xbeach_FCFLAGS=-DUSEMPI -DHAVE_MPI_WTIME
endif


# EXTRA_DIST=s.ind s.inp
EXTRA_DIST=spaceparams.tmpl
SUBDIRS=doc


