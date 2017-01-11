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

rem determine paths
set DstDir="%SolutionDir%\dist\%Platform%\%ConfigurationName%"
set SrcDir1="%SolutionDir%\src\xbeachlibrary\bin\dynamic\%Platform%\%ConfigurationName%"
set SrcDir2="%SolutionDir%\src\xbeachlibrary\bin\bmi\%Platform%\%ConfigurationName%"

echo on

rem clear output dir
if exist "%DstDir%\XBeach.Library.dll" del "%DstDir%\XBeach.Library.dll" /q /f
if exist "%DstDir%\XBeach.BMI.dll" del "%DstDir%\XBeach.BMI.dll" /q /f
if not exist "%DstDir%" mkdir "%DstDir%"

rem copy dll files
if exist "%SrcDir1%\xbeachlibrary_dynamic.dll" copy "%SrcDir1%\xbeachlibrary_dynamic.dll" "%DstDir%\XBeach.Library.dll"
if exist "%SrcDir2%\xbeachlibrary_bmi.dll" copy "%SrcDir2%\xbeachlibrary_bmi.dll" "%DstDir%\XBeach.BMI.dll"
