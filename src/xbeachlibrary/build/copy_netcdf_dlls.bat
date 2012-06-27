echo off

rem read command line parameters
set SolutionDir=%1
set TargetDir=%2

rem replace quotes
set SolutionDir=%SolutionDir:"=%
set TargetDir=%TargetDir:"=%

echo on

copy "%SolutionDir%\lib\win32\vs\pg\*.dll" "%TargetDir%"
copy "%SolutionDir%\lib\win32\all\hdf5\dll\*.dll" "%TargetDir%"
copy "%SolutionDir%\lib\win32\all\szip\dll\*.dll" "%TargetDir%"
copy "%SolutionDir%\lib\win32\all\zlib\dll\*.dll" "%TargetDir%"