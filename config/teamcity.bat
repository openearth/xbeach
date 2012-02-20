echo off

echo solution dir	: %1
echo solution file	: %2
echo Configuration name	: %3
echo Project name	: xbeachlibrary_test

echo on
echo clean build
"c:\Program Files\Microsoft Visual Studio 10.0\Common7\IDE\devenv.exe" "%1\%2" /Clean %3 /Project xbeachlibrary_test

echo build xbech test project
"c:\Program Files\Microsoft Visual Studio 10.0\Common7\IDE\devenv.exe" "%1\%2" /Build %3 /Project xbeachlibrary_test

echo run tests
%1\test\xbeachlibrary_test\bin\%3\xbeachlibrary_test.exe