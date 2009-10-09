#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Revision$

AC_PREREQ(2.59)
AC_INIT(xbeach, [v18], [d.roelvink@unesco-ihe.org])

AC_CONFIG_SRCDIR([boundaryconditions.F90])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(m4)
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([config.h])


# Checks for programs.
AC_PROG_FC
AC_PROG_F77
AC_PROG_INSTALL
# Checks for libraries.

AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],[Build with mpi support @<:@default=check@:>@])], [],[])
AS_IF([test "x$with_mpi" != xno],
	    [AC_LANG_PUSH(Fortran)
	    ACX_MPI
	    AC_LANG_POP
	    FC=$MPIFC],
	    [])
# define usempi (same as HAVE_MPI)
AM_CONDITIONAL(USEMPI, test "x$with_mpi" != "xno")


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.


AC_CONFIG_FILES([Makefile
		 doc/Makefile])
AC_OUTPUT
