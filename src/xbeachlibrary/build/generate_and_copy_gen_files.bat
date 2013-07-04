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

cd "%SolutionDir%src\makeincludes\bin\%Platform%"

if not EXIST "%SolutionDir%src\makeincludes\bin\%Platform%\params.F90" goto creategen
if not EXIST "%SolutionDir%src\makeincludes\bin\%Platform%\spaceparams.tmpl" goto creategen

fc "%SolutionDir%src\xbeachlibrary\params.F90" "%SolutionDir%src\makeincludes\bin\%Platform%\params.F90" > nul
if errorlevel 1 goto creategen
fc "%SolutionDir%src\xbeachlibrary\spaceparams.tmpl" "%SolutionDir%src\makeincludes\bin\%Platform%\spaceparams.tmpl" > nul
if errorlevel 1 goto creategen

goto end

:creategen
rem copy input files
xcopy "%SolutionDir%src\xbeachlibrary\params.F90" "%SolutionDir%src\makeincludes\bin\%Platform%" /D /Y
xcopy "%SolutionDir%src\xbeachlibrary\spaceparams.tmpl" "%SolutionDir%src\makeincludes\bin\%Platform%" /D /Y

rem del intermediate output
del "%SolutionDir%src\makeincludes\bin\%Platform%\parameters.inc"

rem create gen files
"%SolutionDir%src\makeincludes\bin\%Platform%\%type%\makeincludes"

rem copy gen files to right location
xcopy "%SolutionDir%src\makeincludes\bin\%Platform%\*.gen" "%SolutionDir%src\xbeachlibrary\includes\genfiles\" /D /Y

cd "%SolutionDir%src\xbeachlibrary\"

:end
cd "%currentDir%"