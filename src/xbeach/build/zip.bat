echo off

rem 1: SolutionDir
rem 2: ConfigurationName
rem 3: platform (win32 / x64)

echo on

rem copy zip command
copy build\7za.exe %1\dist\%3\%2\

cd %1\dist\%3\%2\

7za.exe a xbeach_%2_%3.zip *
7za.exe d xbeach_%2_%3.zip 7za.exe

cd %1\src\xbeach\

rem del zip command
del %1\dist\%3\%2\7za.exe
