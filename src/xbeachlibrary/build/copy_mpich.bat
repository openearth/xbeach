echo off
echo SolutionDir = %1
echo TargetDir = %2
echo platform (win32 / x64) = %3

echo remove old mpich dir in target dir
if exist "%2\mpich\" rmdir %2\mpich\ /s /q

if exist "C:\Program Files (x86)\" if %3==win32 (
echo 64-bit machine detected, while compiling 32 bit
if exist "c:\Program Files (x86)\MPICH2\" (
echo copy locally installed 32 bit files from "C:\Program Files (x86)\MPICH2\*.*" to %2\mpich\
xcopy "C:\Program Files (x86)\MPICH2\*.*" %2\mpich\ /s
GOTO:EOF
)
) else (
rem Not 32 bit on 64 bit machine, try default install dir
if exist "C:\Program Files\MPICH2\" (
echo copy locally installed files from "C:\Program Files\MPICH2\*.*" to %2\mpich\
xcopy "C:\Program Files\MPICH2\*.*" %2\mpich\ /s
GOTO:EOF
)
)

rem copy files in tree
if %3==x64 (
echo copy included files from "%1\lib\win64\mpich\*.*" to %2\mpich\
xcopy "%1\lib\win64\mpich\*.*" %2\mpich\ /s
GOTO:EOF
)

echo copy included files from "%1\lib\win32\mpich\*.*" to %2\mpich\
xcopy "%1\lib\win32\mpich\*.*" %2\mpich\ /s
