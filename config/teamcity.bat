echo off

echo solution dir	: %1
echo solution file	: %2
echo Configuration name	: %3
echo VS version		: %4
echo Project name	: xbeachlibrary_test

IF %4==VS2010 GOTO VS2010

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
%devenv_path%\devenv.exe "%1\%2" /Clean "%3|Win32" /Project xbeachlibrary_test /Out build.log
type build.log

echo off
IF %ERRORLEVEL% == 0 GOTO BUILD
echo ##teamcity[buildStatus status='FAILURE' text='{build.status.text} in cleanup, check build log']
EXIT /B 1

:BUILD
IF EXIST build.log del build.log
echo on
echo build xbeach test project
%devenv_path%\devenv.exe "%1\%2" /Build "%3|Win32" /Project xbeachlibrary_test /Out build.log
type build.log
echo off

IF %ERRORLEVEL% == 0 GOTO NEXT
echo ##teamcity[buildStatus status='FAILURE' text='{build.status.text} in compilation, check build log']
EXIT /B 1

:NEXT
type build.log 

echo on
echo run tests
%1\test\xbeachlibrary_test\bin\%3\xbeachlibrary_test.exe