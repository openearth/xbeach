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

rem copy "%SolutionDir%\lib\win32\vs\pg\*.dll" "%TargetDir%"
rem copy "%SolutionDir%\lib\win32\all\hdf5\dll\*.dll" "%TargetDir%"
rem copy "%SolutionDir%\lib\win32\all\szip\dll\*.dll" "%TargetDir%"
rem copy "%SolutionDir%\lib\win32\all\zlib\dll\*.dll" "%TargetDir%"