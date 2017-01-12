#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Revision$

AC_PREREQ(2.59)
AC_INIT(xbeach, [v18], [d.roelvink@unesco-ihe.org])

dnl if you're in the source directory you should see xbeach/xbeach.F90
AC_CONFIG_SRCDIR([src/xbeach/xbeach.F90])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
LT_INIT
AM_INIT_AUTOMAKE
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_HEADERS([config.h])

# dependency stuff:

AC_SUBST(deps,[deps])
AC_SUBST(MAKEDEPO,['$(top_srcdir)/makedepo.py])'

# Checks for programs.
# use doxygen
dnl DX_INIT_DOXYGEN(XBeach, Doxyfile, [doc])
dnl DX_DOXYGEN_FEATURE(ON)
dnl DX_DOT_FEATURE(ON)
dnl DX_HTML_FEATURE(ON)
dnl DX_CHM_FEATURE(OFF)
dnl DX_CHI_FEATURE(OFF)
dnl DX_MAN_FEATURE(OFF)
dnl DX_RTF_FEATURE(OFF)
dnl DX_XML_FEATURE(OFF)
dnl DX_PDF_FEATURE(ON)
dnl DX_PS_FEATURE(OFF)

# prefer compilers we tested
AC_PROG_FC([gfortran ifort ifc g95 xlf95 f95 fort  efc pgf95 lf95 ftn xlf90 f90 pgf90 pghpf epcf90 g77 xlf f77 frt pgf77 cf77 fort77 fl32 af77])
AC_PROG_INSTALL
# Checks for libraries.

AC_SUBST(FC_NOMPI, $FC)
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],[Build with mpi support @<:@default=no@:>@])], [],[with_mpi=no])
AS_IF([test "x$with_mpi" != xno],
            [FCFLAGS_save=$FCFLAGS
	    AC_LANG_PUSH(Fortran)
            ACX_MPI
            AC_LANG_POP
            FC=$MPIFC
	    FCFLAGS=$FCFLAGS_save],
            [])
# define usempi (same as HAVE_MPI)
AM_CONDITIONAL(USEMPI, test "x$with_mpi" != "xno")

# check for functions
AC_LANG_PUSH([Fortran])
AC_FC_FREEFORM
AC_FC_LINE_LENGTH([unlimited])
ACX_FORTRAN_ISNAN
ACX_FORTRAN_IEEE_ARITHMETIC

AC_LANG_POP
#AC_SUBST(with_netcdf)

dnl not prefered, but added for compatibility
AC_PATH_PROG(PKG_CONFIG,pkg-config)
dnl prefer to use this one:
dnl PKG_PROG_PKG_CONFIG

AC_ARG_WITH(netcdf,[AS_HELP_STRING([--with-netcdf],[Build with netcdf support @<:@default=no@:>@])], [],[with_netcdf=no])
AS_IF([test "x$with_netcdf" != "xno"],
        [
        # For netcdf >= 4.2 (separate fortran)
        # Check for netcdf 4 if not found...
        PKG_CHECK_MODULES(NETCDF_FORTRAN, netcdf-fortran,[],[PKG_CHECK_MODULES(NETCDF, netcdf >= 4)])
        AC_DEFINE(HAVE_NETCDF,1,[Defined if you have NETCDF library.])
        ],
        [])
AM_CONDITIONAL(USENETCDF, test "x$with_netcdf" != "xno")


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.


AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/xbeachlibrary/Makefile
                 src/xbeach/Makefile])
AC_OUTPUT
