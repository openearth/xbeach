ACLOCAL_AMFLAGS = -I m4 -I config

lib_LTLIBRARIES=libxbeach.la
libxbeach_la_SOURCES= \
	s.ind \
	s.inp \
	sleeper.F90 \
	xmpi.F90 \
	general_mpi.F90 \
	typesandkinds.F90 \
	iso_c_utils.f90 \
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
	spaceparams.F90 \
	compute_tide_zs0.F90 \
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
	vsm_u_XB.f90 \
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
	debugging.F90 \
	introspection.F90 \
	xbeach_bmi.f90

BUILT_SOURCES= \
	version.dat  \
	spacedecl.inc \
	chartoindex.inc \
	get_var.inc \
	get_var_shape.inc \
	getkey.inc \
	index_allocate.inc \
	index_allocated.inc \
	index_deallocate.inc \
	index_reallocate.inc \
	indextos.inc \
	mnemonic.inc \
	set_var.inc \
	space_alloc_arrays.inc \
	space_alloc_scalars.inc \
	space_distribute.inc \
	space_ind.inc \
	space_inp.inc \
	writelog.inc \
	writeloginterface.inc

version.dat: 
	echo "Build_Revision = '`svnversion $(top_srcdir)`'" > version.dat
	echo "Build_URL = '`svn info $(top_srcdir) | grep -m 1 URL`'" >> version.dat
	echo "Build_Date = '`date`'" >> version.dat

logging.F90: version.dat writelog.inc writeloginterface.inc

spacedecl.inc:
	$(top_srcdir)/scripts/generate.py

# Set flags off by default
# include link to top_srcdir so it can find config.h, or should this be builddir
# wwvv changed following
#libxbeach_la_FCFLAGS=-I$(top_builddir) -ffree-line-length-none -fPIC
# into
libxbeach_la_FCFLAGS=-I$(top_builddir) -fPIC
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
EXTRA_DIST=s.ind s.inp nh_pars.inc version.def
CLEANFILES=*.mod *.MOD *.inc version.dat

#depclean: clean
#	rm -f $(deps)

#dep $(deps):
#	$(MAKEDEPO) -p libxbeach_la- -s .lo $(libxbeach_la_SOURCES) > $(deps)

#include $(deps)
