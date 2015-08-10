echo off

rem read command line parameters
set SolutionDir=%1
set ConfigurationName=%2
set Platform=%3

rem replace quotes
set SolutionDir=%SolutionDir:"=%
set ConfigurationName=%ConfigurationName:"=%
set Platform=%Platform:"=%

set type=%ConfigurationName%
set type=%type:netcdf_=%
set type=%type:MPI_=%

echo on

set currentDir=%cd%

cd "%SolutionDir%scripts
dist\generate.exe
cd "%currentDir%"