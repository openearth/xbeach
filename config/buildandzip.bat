echo off

echo solution dir	: %1
echo solution file	: %2
echo VS version		: %3

echo ##teamcity[progressMessage 'Removing old builds']
IF EXIST "%1\dist\win32" (
rmdir "%1\dist\win32" /s /q
mkdir "%1\dist\win32"
)
IF EXIST "%1\dist\x64" (
rmdir "%1\dist\x64" /s /q
mkdir "%1\dist\x64"
)

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
echo ##teamcity[progressMessage 'Cleaning Release||Win32']
%devenv_path%\devenv.exe "%1\%2" /Clean "Release|Win32" /Project xbeach /Out clean_release.log
type clean_release.log
echo ##teamcity[progressMessage 'Cleaning netcdf_Release||Win32']
%devenv_path%\devenv.exe "%1\%2" /Clean "netcdf_Release|Win32" /Project xbeach /Out clean_netcdf.log
type clean_netcdf.log
echo ##teamcity[progressMessage 'Cleaning MPI_Release||Win32']
%devenv_path%\devenv.exe "%1\%2" /Clean "MPI_Release|Win32" /Project xbeach /Out clean_mpi.log
type clean_mpi.log
echo ##teamcity[progressMessage 'Cleaning Release||x64']
%devenv_path%\devenv.exe "%1\%2" /Clean "Release|x64" /Project xbeach /Out clean_release_x64.log
type clean_release_x64.log
echo ##teamcity[progressMessage 'Cleaning MPI_Release||x64']
%devenv_path%\devenv.exe "%1\%2" /Clean "MPI_Release|x64" /Project xbeach /Out clean_mpi.log
type clean_mpi.log

echo off
IF %ERRORLEVEL% == 0 GOTO BUILD
echo ##teamcity[buildStatus status='FAILURE' text='{build.status.text} in cleanup, check build log']
EXIT /B 1

:BUILD
echo on
echo build and zip xbeach executables
echo ##teamcity[progressMessage 'Building Release||Win32']
%devenv_path%\devenv.exe "%1\%2" /Build "Release|Win32" /Project xbeach /Out build_exe.log
type build_exe.log
echo ##teamcity[progressMessage 'Building netcdf_Release||Win32']
%devenv_path%\devenv.exe "%1\%2" /Build "netcdf_Release|Win32" /Project xbeach /Out build_netcdf.log
type build_netcdf.log
echo ##teamcity[progressMessage 'Building MPI_Release||Win32']
%devenv_path%\devenv.exe "%1\%2" /Build "MPI_Release|Win32" /Project xbeach /Out build_mpi.log
type build_mpi.log
echo ##teamcity[progressMessage 'Building Release||x64']
%devenv_path%\devenv.exe "%1\%2" /Build "Release|x64" /Project xbeach /Out build_release_x64.log
type build_release_x64.log
echo ##teamcity[progressMessage 'Building MPI_Release||x64']
%devenv_path%\devenv.exe "%1\%2" /Build "MPI_Release|x64" /Project xbeach /Out build_mpi_x64.log
type build_mpi_x64.log
echo off

IF %ERRORLEVEL% == 0 GOTO END
echo ##teamcity[buildStatus status='FAILURE' text='{build.status.text} in compilation, check build log']
EXIT /B 1

:END