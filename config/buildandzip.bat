echo off

echo solution dir	: %1
echo solution file	: %2
echo VS version		: %3

IF %3==VS2010 GOTO VS2010

set devenv_path="c:\Program Files\Microsoft Visual Studio 9.0\Common7\IDE"
IF EXIST "C:\Program Files (x86)\" set devenv_path="c:\Program Files (x86)\Microsoft Visual Studio 9.0\Common7\IDE"
GOTO COMPILE

:VS2010
set devenv_path="c:\Program Files\Microsoft Visual Studio 10.0\Common7\IDE"
IF EXIST "C:\Program Files (x86)\" set devenv_path="c:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE"

:COMPILE
echo devenv		: %devenv_path%

IF EXIST build.log del build.log
echo on
echo clean build
%devenv_path%\devenv.exe "%1\%2" /Clean "Release|Win32" /Project xbeach /Out clean_release.log
type clean_release.log
%devenv_path%\devenv.exe "%1\%2" /Clean "Release|Win32" /Project xbeach /Out clean_netcdf.log
type clean_netcdf.log
%devenv_path%\devenv.exe "%1\%2" /Clean "Release|Win32" /Project xbeach /Out clean_mpi.log
type clean_mpi.log

echo off
IF %ERRORLEVEL% == 0 GOTO BUILD
echo ##teamcity[buildStatus status='FAILURE' text='{build.status.text} in cleanup, check build log']
EXIT /B 1

:BUILD
echo on
echo build xbeach test project
%devenv_path%\devenv.exe "%1\%2" /Build "Release|Win32" /Project xbeach /Out build_exe.log
type build_exe.log
%devenv_path%\devenv.exe "%1\%2" /Build "netcdf_Release|Win32" /Project xbeach /Out build_netcdf.log
type build_netcdf.log
%devenv_path%\devenv.exe "%1\%2" /Build "MPI_Release|Win32" /Project xbeach /Out build_mpi.log
type build_mpi.log
echo off

IF %ERRORLEVEL% == 0 GOTO END
echo ##teamcity[buildStatus status='FAILURE' text='{build.status.text} in compilation, check build log']
EXIT /B 1

:END