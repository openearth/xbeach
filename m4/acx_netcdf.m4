AC_DEFUN([ACX_NETCDF], [
acx_netcdf_ok=no
FCFLAGS_NETCDF=""
LIBS_NETCDF=""

dnl Check if the library was given in the command line
AC_ARG_WITH(netcdf, [AS_HELP_STRING([--with-netcdf=DIR], [http://www.unidata.ucar.edu/packages/netcdf/])])
case $with_netcdf in
  yes | "") ;;
  no ) acx_netcdf_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) LIBS_NETCDF="$with_netcdf" ;;
  *) LIBS_NETCDF="-l$with_netcdf" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_netcdf_save_LIBS="$LIBS"
acx_netcdf_save_FCFLAGS="$FCFLAGS"

dnl The tests
if test $acx_netcdf_ok = no; then
  AC_MSG_CHECKING([for netcdf])  
  # If LIBS_NETCDF has been passed with --with-netcdf just test this
  if test "$LIBS_NETCDF"; then
    netcdf_fcflags=""; netcdf_libs="$LIBS_NETCDF"
FCFLAGS="$netcdf_fcflags $acx_netcdf_save_FCFLAGS"
LIBS="$netcdf_libs $acx_netcdf_save_LIBS"
AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
use netcdf
integer :: ncid
integer :: status
status=nf90_close(ncid)
]), [acx_netcdf_ok=yes; FCFLAGS_NETCDF="$netcdf_fcflags"; LIBS_NETCDF="$netcdf_libs"], [])
  else
    for netcdf_fcflags in "" -I/usr/include/netcdf-3; do
      for netcdf_libsL in ""; do
        for netcdf_libsl in "" -lnetcdf "-lnetcdff -lnetcdf"; do
	  if test "$netcdf_libsL" -a "$netcdf_libsl"; then
	    netcdf_libs="$netcdf_libsL $netcdf_libsl"
	  else
	    netcdf_libs="$netcdf_libsL$netcdf_libsl"
	  fi
FCFLAGS="$netcdf_fcflags $acx_netcdf_save_FCFLAGS"
LIBS="$netcdf_libs $acx_netcdf_save_LIBS"
AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
use netcdf
integer :: ncid
integer :: status
status=nf90_close(ncid)
]), [acx_netcdf_ok=yes; FCFLAGS_NETCDF="$netcdf_fcflags"; LIBS_NETCDF="$netcdf_libs"], [])
          if test $acx_netcdf_ok != no; then break; fi
        done
        if test $acx_netcdf_ok != no; then break; fi
      done
      if test $acx_netcdf_ok != no; then break; fi
    done
  fi
  if test $acx_netcdf_ok = no -o -z "$FCFLAGS_NETCDF$LIBS_NETCDF"; then
    AC_MSG_RESULT([$acx_netcdf_ok])
  else
    AC_MSG_RESULT([$acx_netcdf_ok ($FCFLAGS_NETCDF $LIBS_NETCDF)])
  fi
fi

AC_SUBST(FCFLAGS_NETCDF)
AC_SUBST(LIBS_NETCDF)
FCFLAGS="$acx_netcdf_save_FCFLAGS"
LIBS="$acx_netcdf_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_netcdf_ok" = xyes; then
  AC_DEFINE(HAVE_NETCDF,1,[Defined if you have NETCDF library.])
  $1
else
  if test $acx_netcdf_ok != disable; then
    AC_MSG_WARN([Could not find netcdf library. 
                *** Will compile without netcdf support])
  fi
  $2
fi
])dnl ACX_NETCDF
