ACLOCAL_AMFLAGS = -I m4 -I config
bin_PROGRAMS=xbeach

xbeach_SOURCES=\
	input.F90\
	xbeach.F90
xbeach_LDADD=\
	$(builddir)/../xbeachlibrary/libxbeach.la
BUILT_SOURCES= \
	$(builddir)/../xbeachlibrary/version.def\
	$(builddir)/../xbeachlibrary/version.dat

# Set flags off by default
xbeach_FCFLAGS=-I$(builddir)/../xbeachlibrary
if USEMPI
xbeach_FCFLAGS+=-DUSEMPI -DHAVE_MPI_WTIME 
endif
if USENETCDF
# Why don't we use NETCDF_FCFLAGS?? or FFLAGS?
xbeach_FCFLAGS+=-DUSENETCDF ${NETCDF_CFLAGS}
# Manualy add the netcdff (fortran dll)
LIBS+=${NETCDF_LIBS}  -lnetcdff
endif


CLEANFILES=*.mod *.MOD *.gen version.dat 


