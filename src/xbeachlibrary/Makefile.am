ACLOCAL_AMFLAGS = -I m4 -I config

lib_LTLIBRARIES=libxbeach.la
if USENETCDF
extra_xbeach_sources=ncoutput.F90
endif
libxbeach_la_SOURCES= \
	s.ind \
	s.inp \
	xmpi.F90 \
	general_mpi.F90 \
	typesandkinds.F90 \
	logging.F90 \
	filefunctions.F90 \
	readkey.F90 \
	mnemonic.F90 \
	interp.F90 \
	math_tools.F90 \
	params.F90 \
	getkey.F90 \
	spaceparams.F90 \
	wave_functions.F90 \
	waveparams.F90 \
	waveparamsnew.F90 \
	constants.F90 \
	roelvink.F90 \
	wave_stationary.F90 \
	wave_instationary.F90 \
	wave_timestep.F90 \
        boundaryconditions.F90 \
	morphevolution.F90 \
	beachwizard.F90 \
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
	output.F90 \
        ship.F90 \
	mnemoniciso.F90 \
	libxbeach.F90 \
	libxbeach_dynamic.F90 \
	introspection.F90
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

mi=$(builddir)/../makeincludes/makeincludes$(EXEEXT)
spacedecl.gen: $(mi)
	echo $@ | $(mi)
mnemonic.gen: $(mi)
	echo $@ | $(mi)
indextos.gen: $(mi)
	echo $@ | $(mi)
space_alloc_scalars.gen: $(mi)
	echo $@ | $(mi)
space_alloc_arrays.gen: $(mi)
	echo $@ | $(mi)
space_ind.gen: $(mi)
	echo $@ | $(mi)
space_inp.gen: $(mi)
	echo $@ | $(mi)
chartoindex.gen: $(mi)
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
CLEANFILES=*.mod *.MOD *.gen version.dat 



