echo off

echo solution dir	: %1
echo VS version		: %2
echo Configuration	: %3
echo Platform (Win32/x64): %4

echo Selecting devenv to compile with (VS2008 or VS2010)

IF %2==VS2010 GOTO VS2010

set devenv_path="c:\Program Files\Microsoft Visual Studio 9.0\Common7\IDE"
IF EXIST "C:\Program Files (x86)\" set devenv_path="c:\Program Files (x86)\Microsoft Visual Studio 9.0\Common7\IDE"
GOTO COMPILE

:VS2010
set devenv_path="c:\Program Files\Microsoft Visual Studio 10.0\Common7\IDE"
IF EXIST "C:\Program Files (x86)\" set devenv_path="c:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE"

:COMPILE
echo devenv		: %devenv_path%

set solutionpath="%1/XBeach_%2.sln"
echo solution path	: %solutionpath%

IF EXIST build.log del build.log
echo clean build
echo ##teamcity[progressMessage 'Cleaning %3 - %4']
echo on
%devenv_path%\devenv.exe %solutionpath% /Clean "%3|%4" /Project xbeach /Out clean.log
type clean.log
echo off

IF %ERRORLEVEL% == 0 GOTO BUILD
echo ##teamcity[buildStatus status='FAILURE' text='{build.status.text} in cleanup, check build log']
EXIT /B 1

:BUILD
echo build and zip xbeach executables
echo ##teamcity[progressMessage 'Building %3 - %4']
echo on
%devenv_path%\devenv.exe %solutionpath% /Build "%3|%4" /Project xbeach /Out build.log
type build.log
%devenv_path%\devenv.exe %solutionpath% /Build "%3|%4" /Project xbeachlibrary_dynamic /Out build_library.log
type build_library.log
echo off

IF %ERRORLEVEL% == 0 GOTO END
echo ##teamcity[buildStatus status='FAILURE' text='{build.status.text} in compilation, check build log']
EXIT /B 1

:END
echo ##teamcity[progressMessage 'Build finished %3 - %4']