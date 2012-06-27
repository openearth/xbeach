echo off

rem read command line parameters
set SolutionDir=%1
set ConfigurationName=%2

rem replace quotes
set SolutionDir=%SolutionDir:"=%
set ConfigurationName=%ConfigurationName:"=%

echo on

cd "%SolutionDir%test\xbeachlibrary_test\bin\%ConfigurationName%"

rem run tests and create html
start xbeachlibrary_test.exe

:loopbusy

tasklist /FI "IMAGENAME eq xbeachlibrary_test.exe" 2>NUL | find /I /N "xbeachlibrary_test.exe">NUL
if "%ERRORLEVEL%"=="0" goto loopbusy

ftnunit.html

cd "%SolutionDir%test\xbeachlibrary_test"