echo off
rem %1 = SolutionDir
rem %2 = ConfigurationName
echo on

cd %1test\xbeachlibrary_test\bin\%2

rem run tests and create html
start xbeachlibrary_test.exe

:loopbusy

tasklist /FI "IMAGENAME eq xbeachlibrary_test.exe" 2>NUL | find /I /N "xbeachlibrary_test.exe">NUL
if "%ERRORLEVEL%"=="0" goto loopbusy

ftnunit.html

cd %1test\xbeachlibrary_test