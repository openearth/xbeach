echo off
rem 1: SolutionDir
rem 2: ConfigurationName
rem 3: platform (win32 / x64)

set type=%2
set type=%type:netcdf_=%
set type=%type:MPI_=%
echo on

cd %1src\makeincludes\bin\%3

rem copy input files
copy %1src\xbeachlibrary\params.F90 %1src\makeincludes\bin\%3
copy %1src\xbeachlibrary\*.tmpl %1src\makeincludes\bin\%3

rem del intermediate output
del %1src\makeincludes\bin\%3\parameters.inc

rem create gen files
%1src\makeincludes\bin\%3\%type%\makeincludes

rem delete input
del %1src\makeincludes\bin\%3\params.F90

rem copy gen files to right location
copy %1src\makeincludes\bin\%3\*.gen %1src\xbeachlibrary\includes\genfiles\

cd %1src\xbeachlibrary\
