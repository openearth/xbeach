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
if exist "%SolutionDir%\dist\%Platform%\%ConfigurationName%\" rmdir "%SolutionDir%\dist\%Platform%\%ConfigurationName%\" /q /s
mkdir "%SolutionDir%\dist\%Platform%\%ConfigurationName%\"

rem copy exe
copy "%TargetDir%\*.exe" "%SolutionDir%\dist\%Platform%\%ConfigurationName%\"
copy "%TargetDir%\*.dll" "%SolutionDir%\dist\%Platform%\%ConfigurationName%\"

rem copy additional dll files from static output folder
copy "%SolutionDir%\src\xbeachlibrary\bin\static\%Platform%\%ConfigurationName%\*.dll" "%SolutionDir%\dist\%Platform%\%ConfigurationName%\"

rem copy additional dll files from static to bin folder
copy "%SolutionDir%\src\xbeachlibrary\bin\static\%Platform%\%ConfigurationName%\*.dll" "%TargetDir%\"

rem copy current manual
copy "%SolutionDir%\doc\manual\xbeach_manual.pdf" "%SolutionDir%\dist\%Platform%\%ConfigurationName%\"

rem copy license
copy "%SolutionDir%\LICENSE" "%SolutionDir%\dist\%Platform%\%ConfigurationName%\"