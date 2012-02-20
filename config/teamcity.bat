echo off

echo solution dir	: %1
echo solution file	: %2
echo Configuration name	: %3
echo Project name	: xbeachlibrary_test

set devenv_path="c:\Program Files\Microsoft Visual Studio 10.0\Common7\IDE"
IF EXIST "C:\Program Files (x86)\" set devenv_path="c:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE"

echo devenv		: %devenv_path%

echo on
echo clean build
%devenv_path%\devenv.exe "%1\%2" /Clean "%3|Win32" /Project xbeachlibrary_test

echo build xbech test project
%devenv_path%\devenv.exe "%1\%2" /Build "%3|Win32" /Project xbeachlibrary_test

echo run tests
%1\test\xbeachlibrary_test\bin\%3\xbeachlibrary_test.exe