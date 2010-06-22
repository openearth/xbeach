bin_PROGRAMS=makeincludes xbeach


ACLOCAL_AMFLAGS = -I config -I m4
makeincludes_SOURCES= \
	makeincludes.F90
makeincludes_LDFLAGS=

if USENETCDF
extra_xbeach_sources=ncoutput.F90
endif

xbeach_SOURCES= \
	s.ind \
	s.inp \
	xmpi.F90 \
	general_mpi.F90 \
	logging.F90 \
	filefunctions.F90 \
	readkey.F90 \
	mnemonic.F90 \
	interp.F90 \
	math_tools.F90 \
	general_fileio.F90 \
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
	flow_secondorder.F90 \
	solver.F90 \
	nonh.F90 \
	flow_timestep.F90 \
	postprocess.F90 \
	varianceupdate.F90 \
	varoutput.F90 \
	${extra_xbeach_sources} \
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
	s.inp \
	version.dat

version.dat: 
	echo "Build_Revision = '`svnversion`'" > version.dat
	echo "Build_URL = '`svn info | grep URL`'" >> version.dat
	echo "Build_Date = '`date`'" >> version.dat

s.ind: space_ind.gen
s.inp: space_inp.gen


spacedecl.gen: makeincludes$(EXEEXT)
	echo $@ | ./makeincludes
mnemonic.gen: makeincludes$(EXEEXT)
	echo $@ | ./makeincludes
indextos.gen: makeincludes$(EXEEXT)
	echo $@ | ./makeincludes
space_alloc_scalars.gen: makeincludes$(EXEEXT)
	echo $@ | ./makeincludes
space_alloc_arrays.gen: makeincludes$(EXEEXT)
	echo $@ | ./makeincludes
space_ind.gen: makeincludes$(EXEEXT)
	echo $@ | ./makeincludes
space_inp.gen: makeincludes$(EXEEXT)
	echo $@ | ./makeincludes
chartoindex.gen: makeincludes$(EXEEXT)
	echo $@ | ./makeincludes

# for testing we need the testbed
check: testbed

# Link the testbed to the current directory, autotools will not look below it's root. So we link it in. 
# Look for it in the two directories above (trunk/.. and branches/foo/../..
testbed:
	if [ -d ../testbed ]; \
	 then ln -s ../testbed testbed; \
	elif [ -d ../../testbed ]; \
	 then ln -s ../../testbed testbed; \
	else echo "Could not find testbed in .. or ../.. relative to `pwd`. The testbed can be found at https://repos.deltares.nl/repos/XBeach/testbed"; \
	fi; \
	for i in `ls testbed/input`; \
	do \
	 echo "adding $$i"; \
	 ln -sf ../test.sh  testbed/tools/autotests/autotest_$$i ; \
	 ln -sf ../test.sh  testbed/tools/autotests/mpi_autotest_$$i ; \
	done
testdir=./testbed/tools/autotests
if USEMPI
check_SCRIPTS= \
	${testdir}/mpi_autotest_*
TESTS=${testdir}/mpi_autotest_*
else
check_SCRIPTS= \
	${testdir}/autotest_*
TESTS=${testdir}/autotest_*
endif

# Set flags off by default
xbeach_FCFLAGS=

# if mpi, then add some flags
if USEMPI
xbeach_FCFLAGS+=-DUSEMPI -DHAVE_MPI_WTIME
endif
if USENETCDF
xbeach_FCFLAGS+=-DUSENETCDF ${NETCDF_CFLAGS}
# Manualy add the netcdff (fortran dll)
LIBS+=${NETCDF_LIBS}  -lnetcdff
endif


clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	rm -rf testbed; \
	rm -rf run

# EXTRA_DIST=s.ind s.inp
EXTRA_DIST=spaceparams.tmpl s.ind s.inp nh_pars.inc version.def VS2008 xbeach.dsp xbeach.dsw
SUBDIRS=doc
CLEANFILES=*.mod *.MOD *.gen version.dat 


