echo off

rem read command line parameters
set SolutionDir=%1
set ConfigurationName=%2
set TargetDir=%3
set Platform=%4

rem replace quotes
set SolutionDir=%SolutionDir:"=%
set ConfigurationName=%ConfigurationName:"=%
set TargetDir=%TargetDir:"=%
set Platform=%Platform:"=%

echo on

rem clear output dir
if exist "%SolutionDir%\dist\%Platform%\%ConfigurationName%\XBeach.Library.dll" del "%SolutionDir%\dist\%Platform%\%ConfigurationName%\XBeach.Library.dll" /q /f
if not exist "%SolutionDir%\dist\%Platform%\%ConfigurationName%\" mkdir "%SolutionDir%\dist\%Platform%\%ConfigurationName%\"

rem copy additional dll files from static to bin folder
copy "%SolutionDir%\src\xbeachlibrary\bin\dynamic\%Platform%\%ConfigurationName%\xbeachlibrary_dynamic.dll" "%SolutionDir%\dist\%Platform%\%ConfigurationName%\XBeach.Library.dll"
copy "%SolutionDir%\src\xbeachlibrary\bin\bmi\%Platform%\%ConfigurationName%\xbeachlibrary_bmi.dll" "%SolutionDir%\dist\%Platform%\%ConfigurationName%\XBeach.BMI.dll"
