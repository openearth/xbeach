echo off

rem read command line parameters
set SolutionDir=%1
set TargetDir=%2
set Platform=%3
set Configuration=%4

rem replace quotes
set SolutionDir=%SolutionDir:"=%
set TargetDir=%TargetDir:"=%
set Platform=%Platform:"=%
set Configuration=%Configuration:"=%

echo on

copy "%SolutionDir%\lib\%Platform%\all\*.dll" "%TargetDir%"
copy "%SolutionDir%\lib\%Platform%\netcdf\%Configuration%\*.dll" "%TargetDir%"