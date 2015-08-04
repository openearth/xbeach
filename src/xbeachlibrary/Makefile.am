ACLOCAL_AMFLAGS = -I m4 -I config

lib_LTLIBRARIES=libxbeach.la
libxbeach_la_SOURCES= \
	sleeper.F90 \
	xmpi.F90 \
	general_mpi.F90 \
	typesandkinds.F90 \
	logging.F90 \
	paramsconst.F90 \
	filefunctions.F90 \
	readkey.F90 \
	mnemonic.F90 \
	interp.F90 \
	math_tools.F90 \
	params.F90 \
	getkey.F90 \
	spaceparamsdef.F90 \
	indextos.F90 \
	spaceparams.F90 \
	wetcells.F90 \
	vegetation.F90 \
	wave_functions.F90 \
	waveparams.F90 \
	waveparamsnew.F90 \
	constants.F90 \
	roelvink.F90 \
	wave_stationary.F90 \
	wave_instationary.F90 \
	wave_directions.F90 \
	wave_timestep.F90 \
	solver.F90 \
	flow_secondorder.F90 \
	nonh.F90 \
	boundaryconditions.F90 \
	morphevolution.F90 \
	beachwizard.F90 \
	groundwater.F90 \
	drifters.F90 \
	readtide.F90 \
	readwind.F90 \
	initialize.F90 \
	timestep.F90 \
	bedroughness.F90 \
	flow_timestep.F90 \
	postprocess.F90 \
	varianceupdate.F90 \
	varoutput.F90 \
	ncoutput.F90 \
	output.F90 \
        ship.F90 \
	mnemoniciso.F90 \
	loopcounters.F90 \
	libxbeach.F90 \
	libxbeach_dynamic.F90 \
	debugging.F90 \
	introspection.F90
BUILT_SOURCES= \
	version.dat  \
	spacedecl.gen	\
	mnemonic.gen \
	indextos.gen \
	space_alloc_arrays.gen \
	space_ind.gen \
	space_inp.gen \
	chartoindex.gen \
	getkey.gen \
	index_allocate.gen \
	index_deallocate.gen \
	index_reallocate.gen \
	index_allocated.gen

version.dat: 
	echo "Build_Revision = '`svnversion`'" > version.dat
	echo "Build_URL = '`svn info | grep -m 1 URL`'" >> version.dat
	echo "Build_Date = '`date`'" >> version.dat

getkey.gen: params.F90


# We should also specify the dependencies of source files and gen files. If we want the source to be rebuild when gen files get updated
#
# wwvv: Here it is:
# (the .gen files are automatically build because of the dependencies in
# deps)

mi=$(builddir)/../makeincludes/makeincludes$(EXEEXT)
spacedecl.gen: $(mi) spaceparams.tmpl
	echo $@ | $(mi)
mnemonic.gen: $(mi) spaceparams.tmpl
	echo $@ | $(mi)
indextos.gen: $(mi) spaceparams.tmpl
	echo $@ | $(mi)
index_allocate.gen: $(mi) spaceparams.tmpl
	echo $@ | $(mi)
index_deallocate.gen: $(mi) spaceparams.tmpl
	echo $@ | $(mi)
index_reallocate.gen: $(mi) spaceparams.tmpl
	echo $@ | $(mi)
index_allocated.gen: $(mi) spaceparams.tmpl
	echo $@ | $(mi)
space_alloc_arrays.gen: $(mi) spaceparams.tmpl
	echo $@ | $(mi)
space_ind.gen: $(mi) spaceparams.tmpl
	echo $@ | $(mi)
space_inp.gen: $(mi) spaceparams.tmpl
	echo $@ | $(mi)
chartoindex.gen: $(mi) spaceparams.tmpl
	echo $@ | $(mi)
getkey.gen: $(mi) params.F90
	echo $@ | $(mi)
parameters.inc: $(mi) params.F90
	echo $@ | $(mi)

# Set flags off by default
# include link to top_srcdir so it can find config.h, or should this be builddir
libxbeach_la_FCFLAGS=-I$(top_builddir)
# if mpi, then add some flags
if USEMPI
libxbeach_la_FCFLAGS+=-DUSEMPI -DHAVE_MPI_WTIME
endif
if USENETCDF
libxbeach_la_FCFLAGS+=-DUSENETCDF ${NETCDF_CFLAGS} ${NETCDF_FORTRAN_CFLAGS}
# Manualy add the netcdff (fortran dll)
LIBS+=${NETCDF_LIBS} ${NETCDF_FORTRAN_LIBS} -lnetcdff
endif

# EXTRA_DIST=s.ind s.inp
EXTRA_DIST=spaceparams.tmpl s.ind s.inp nh_pars.inc version.def
CLEANFILES=*.mod *.MOD *.gen version.dat parameters.inc

depclean: clean
	rm -f $(deps)

dep $(deps):
	$(MAKEDEPO) -p libxbeach_la- -s .lo $(libxbeach_la_SOURCES) > $(deps)

include $(deps)
