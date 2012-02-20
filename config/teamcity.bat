echo off
rem %1 path to devenv exe

echo solution dir	: %2
echo solution file	: %3
echo Configuration name	: %4
echo Project name	: xbeachlibrary_test

echo clean build
%1\devenv.exe "%2\%3" /Clean %4 /Project xbeachlibrary_test

echo build xbech test project
%1\devenv.exe "%2\%3" /Build %4 /Project xbeachlibrary_test

echo run tests
%2\test\xbeachlibrary_test\bin\%4\xbeachlibrary_test.exe