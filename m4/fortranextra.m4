# AC_FC_LINE_LENGTH([LENGTH], [ACTION-IF-SUCCESS],
#                   [ACTION-IF-FAILURE = FAILURE])
# ------------------------------------------------
# Look for a compiler flag to make the Fortran (FC) compiler accept long lines
# in the current (free- or fixed-format) source code, and adds it to FCFLAGS.
# The optional LENGTH may be 80, 132 (default), or `unlimited' for longer
# lines.  Note that line lengths above 254 columns are not portable, and some
# compilers (hello ifort) do not accept more than 132 columns at least for
# fixed format.  Call ACTION-IF-SUCCESS (defaults to nothing) if successful
# (i.e. can compile code using new extension) and ACTION-IF-FAILURE (defaults
# to failing with an error message) if not.  (Defined via DEFUN_ONCE to
# prevent flag from being added to FCFLAGS multiple times.)
# You should call AC_FC_FREEFORM or AC_FC_FIXEDFORM to set the desired format
# prior to using this macro.
#
# The known flags are:
# -f{free,fixed}-line-length-N with N 72, 80, 132, or 0 or none for none.
# -ffree-line-length-none: GNU gfortran
#       -qfixed=132 80 72: IBM compiler (xlf)
#                -Mextend: Cray
#            -132 -80 -72: Intel compiler (ifort)
#                          Needs to come before -extend_source because ifort
#                          accepts that as well with an optional parameter and
#                          doesn't fail but only warns about unknown arguments.
#          -extend_source: SGI compiler
#     -W NN (132, 80, 72): Absoft Fortran
#          +extend_source: HP Fortran (254 in either form, default is 72 fixed,
#                          132 free)
#                   -wide: Lahey/Fujitsu Fortran (255 cols in fixed form)
#                      -e: Sun Fortran compiler (132 characters)
AC_DEFUN([AC_FC_LINE_LENGTH],
[AC_LANG_PUSH([Fortran])dnl
m4_case(m4_default([$1], [132]),
  [unlimited], [ac_fc_line_len_string=unlimited
                       ac_fc_line_len=0
                       ac_fc_line_length_test='
      subroutine longer_than_132(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,'\
'arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19)'],
  [132],            [ac_fc_line_len=132
                       ac_fc_line_length_test='
      subroutine longer_than_80(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,'\
'arg10)'],
  [80],             [ac_fc_line_len=80
                       ac_fc_line_length_test='
      subroutine longer_than_72(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)'],
  [m4_warning([Invalid length argument `$1'])])
: ${ac_fc_line_len_string=$ac_fc_line_len}
AC_CACHE_CHECK(
[for Fortran flag needed to accept $ac_fc_line_len_string column source lines],
               [ac_cv_fc_line_length],
[ac_cv_fc_line_length=unknown
ac_fc_line_length_FCFLAGS_save=$FCFLAGS
for ac_flag in none \
               -ffree-line-length-none -ffixed-line-length-none \
               -ffree-line-length-$ac_fc_line_len \
               -ffixed-line-length-$ac_fc_line_len \
               -qfixed=$ac_fc_line_len -Mextend \
               -$ac_fc_line_len -extend_source \
               "-W $ac_fc_line_len" +extend_source -wide -e
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_line_length_FCFLAGS_save $ac_flag"
  AC_COMPILE_IFELSE([[$ac_fc_line_length_test
      end subroutine]],
                    [ac_cv_fc_line_length=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fc_line_length_FCFLAGS_save
])
if test "x$ac_cv_fc_line_length" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Fortran does not accept long source lines], 77)])
else
  if test "x$ac_cv_fc_line_length" != xnone; then
    FCFLAGS="$FCFLAGS $ac_cv_fc_line_length"
  fi
  $2
fi
AC_LANG_POP([Fortran])dnl
])# AC_FC_LINE_LENGTH

# AC_FC_PREPROCESS([ACTION-IF-SUCCESS],
#                   [ACTION-IF-FAILURE = FAILURE])
# Look for a compiler flag to make the Fortran (FC) compiler preprocess source
# in the current (free- or fixed-format) source code, and adds it to FCFLAGS.
# Call ACTION-IF-SUCCESS (defaults to nothing) if successful
# (i.e. can compile code using new extension) and ACTION-IF-FAILURE (defaults
# to failing with an error message) if not.  (Defined via DEFUN_ONCE to
# prevent flag from being added to FCFLAGS multiple times.)
# You should call AC_FC_FREEFORM or AC_FC_FIXEDFORM to set the desired format
# prior to using this macro.
#
# The known flags are:
# -fpp ifort
# -cpp gfortran
dnl this AC_FC_PREPROCESS macro checks for a preprocess flags.
AC_DEFUN([AC_FC_PREPROCESS],
[AC_LANG_PUSH([Fortran])dnl
AC_CACHE_CHECK(
[for Fortran flag needed to preprocess source],
               [ac_cv_fc_preprocess],
[ac_cv_fc_preprocess=unknown
ac_fc_preprocess_FCFLAGS_save=$FCFLAGS
dnl a fortran program with some preprocess statement
ac_fc_preprocess_test='program a
#if defined (A)
#endif
end program'

for ac_flag in -fpp \
               -cpp \
               "-x f95-cpp-input"
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_preprocess_FCFLAGS_save -Werror $ac_flag"
  AC_COMPILE_IFELSE([[$ac_fc_preprocess_test]],
                    [ac_cv_fc_preprocess=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fc_preprocess_FCFLAGS_save
])
if test "x$ac_cv_fc_preprocess" = xunknown; then
  m4_default([$2],
             [AC_MSG_ERROR([Fortran does do preprocessing])])
else
  if test "x$ac_cv_fc_preprocess" != xnone; then
    FCFLAGS="$FCFLAGS $ac_cv_fc_preprocess"
  fi
  $1
fi
AC_LANG_POP([Fortran])dnl
])# AC_FC_PREPROCESS



# AC_F77_LINE_LENGTH([LENGTH], [ACTION-IF-SUCCESS],
#                   [ACTION-IF-FAILURE = FAILURE])
# ------------------------------------------------
# Look for a compiler flag to make the Fortran (F77) compiler accept long lines
# in the current (free- or fixed-format) source code, and adds it to FFLAGS.
# The optional LENGTH may be 80, 132 (default), or `unlimited' for longer
# lines.  Note that line lengths above 254 columns are not portable, and some
# compilers (hello ifort) do not accept more than 132 columns at least for
# fixed format.  Call ACTION-IF-SUCCESS (defaults to nothing) if successful
# (i.e. can compile code using new extension) and ACTION-IF-FAILURE (defaults
# to failing with an error message) if not.  (Defined via DEFUN_ONCE to
# prevent flag from being added to FFLAGS multiple times.)
# You should call AC_F77_FREEFORM or AC_F77_FIXEDFORM to set the desired format
# prior to using this macro.
#
# The known flags are:
# -f{free,fixed}-line-length-N with N 72, 80, 132, or 0 or none for none.
# -ffree-line-length-none: GNU gfortran
#       -qfixed=132 80 72: IBM compiler (xlf)
#                -Mextend: Cray
#            -132 -80 -72: Intel compiler (ifort)
#                          Needs to come before -extend_source because ifort
#                          accepts that as well with an optional parameter and
#                          doesn't fail but only warns about unknown arguments.
#          -extend_source: SGI compiler
#     -W NN (132, 80, 72): Absoft Fortran
#          +extend_source: HP Fortran (254 in either form, default is 72 fixed,
#                          132 free)
#                   -wide: Lahey/Fujitsu Fortran (255 cols in fixed form)
#                      -e: Sun Fortran compiler (132 characters)
AC_DEFUN([AC_F77_LINE_LENGTH],
[AC_LANG_PUSH([Fortran 77])dnl
m4_case(m4_default([$1], [132]),
  [unlimited], [ac_f77_line_len_string=unlimited
                       ac_f77_line_len=0
                       ac_f77_line_length_test='
      subroutine longer_than_132(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,'\
'arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19)'],
  [132],            [ac_f77_line_len=132
                       ac_f77_line_length_test='
      subroutine longer_than_80(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,'\
'arg10)'],
  [80],             [ac_f77_line_len=80
                       ac_f77_line_length_test='
      subroutine longer_than_72(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)'],
  [m4_warning([Invalid length argument `$1'])])
: ${ac_f77_line_len_string=$ac_f77_line_len}
AC_CACHE_CHECK(
[for Fortran flag needed to accept $ac_f77_line_len_string column source lines],
               [ac_cv_f77_line_length],
[ac_cv_f77_line_length=unknown
ac_f77_line_length_FFLAGS_save=$FFLAGS
for ac_flag in none \
               -ffree-line-length-none -ffixed-line-length-none \
               -ffree-line-length-$ac_f77_line_len \
               -ffixed-line-length-$ac_f77_line_len \
               -qfixed=$ac_f77_line_len -Mextend \
               -$ac_f77_line_len -extend_source \
               "-W $ac_f77_line_len" +extend_source -wide -e
do
  test "x$ac_flag" != xnone && FFLAGS="$ac_f77_line_length_FFLAGS_save $ac_flag"
  AC_COMPILE_IFELSE([[$ac_f77_line_length_test
      end subroutine]],
                    [ac_cv_f77_line_length=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FFLAGS=$ac_f77_line_length_FFLAGS_save
])
if test "x$ac_cv_f77_line_length" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Fortran does not accept long source lines], 77)])
else
  if test "x$ac_cv_f77_line_length" != xnone; then
    FFLAGS="$FFLAGS $ac_cv_f77_line_length"
  fi
  $2
fi
AC_LANG_POP([Fortran 77])dnl
])# AC_F77_LINE_LENGTH



# AX_FC_HAS_MODULE([MODULE])
# Check if a module is usable from fortran.
# Sets the FC_HAS_MODULE variable
# ------------------------------------------------
AC_DEFUN([AX_FC_HAS_MODULE],
[AC_CACHE_CHECK([whether $FC has the module $1],
[ax_cv_fc_has_module],
[AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
AC_COMPILE_IFELSE([
program moduletest
use $1
end program moduletest
],
[ax_cv_fc_has_module="yes"],
[ax_cv_fc_has_module="no"])
])
if test x"$ax_cv_fc_has_module" = xyes; then 
AC_DEFINE([HAVE_$1],1,[Define if Fortran supports module $1.])
fi
AM_CONDITIONAL([HAVE_$1], [test x"$HAVE_$1" = x1])
])
