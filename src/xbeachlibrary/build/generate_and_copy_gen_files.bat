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

cd "%SolutionDir%src\makeincludes\bin\%Platform%"

rem copy input files
copy "%SolutionDir%src\xbeachlibrary\params.F90" "%SolutionDir%src\makeincludes\bin\%Platform%"
copy "%SolutionDir%src\xbeachlibrary\*.tmpl" "%SolutionDir%src\makeincludes\bin\%Platform%"

rem del intermediate output
del "%SolutionDir%src\makeincludes\bin\%Platform%\parameters.inc"

rem create gen files
"%SolutionDir%src\makeincludes\bin\%Platform%\%type%\makeincludes"

rem delete input
del "%SolutionDir%src\makeincludes\bin\%Platform%\params.F90"

rem copy gen files to right location
copy "%SolutionDir%src\makeincludes\bin\%Platform%\*.gen" "%SolutionDir%src\xbeachlibrary\includes\genfiles\"

cd "%SolutionDir%src\xbeachlibrary\"
