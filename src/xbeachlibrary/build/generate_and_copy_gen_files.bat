echo off
rem 1: SolutionDir
rem 2: ConfigurationName

set type=%2
set type=%type:netcdf_=%
set type=%type:MPI_=%
echo on

cd %1src\makeincludes\bin

rem copy input files
copy %1src\xbeachlibrary\params.F90 %1src\makeincludes\bin
copy %1src\xbeachlibrary\*.tmpl %1src\makeincludes\bin

rem del intermediate output
del %1src\makeincludes\bin\parameters.inc

rem create gen files
%1src\makeincludes\bin\%type%\makeincludes

rem delete input
del %1src\makeincludes\bin\params.F90

rem copy gen files to right location
copy %1src\makeincludes\bin\*.gen %1src\xbeachlibrary\includes\genfiles\

cd %1src\xbeachlibrary\