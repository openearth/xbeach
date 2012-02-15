echo off

rem 1: SolutionDir
rem 2: ConfigurationName

echo on

rem copy zip command
copy build\7za.exe %1\dist\win32\%2\

cd %1\dist\win32\%2\

7za.exe a xbeach_%2.zip *
7za.exe d xbeach_%2.zip 7za.exe

cd %1\src\xbeach\

rem del zip command
del %1\dist\win32\%2\7za.exe