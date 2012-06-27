echo off

rem read command line parameters
set SolutionDir=%1
set ConfigurationName=%2
set Platform=%3

rem replace quotes
set SolutionDir=%SolutionDir:"=%
set ConfigurationName=%ConfigurationName:"=%
set Platform=%Platform:"=%

echo on

rem copy zip command
copy build\7za.exe "%SolutionDir%\dist\%Platform%\%ConfigurationName%\"

cd "%SolutionDir%\dist\%Platform%\%ConfigurationName%\"

7za.exe a "xbeach_%ConfigurationName%_%Platform%.zip" *
7za.exe d "xbeach_%ConfigurationName%_%Platform%.zip" 7za.exe

cd "%SolutionDir%\src\xbeach\"

rem del zip command
del "%SolutionDir%\dist\%Platform%\%ConfigurationName%\7za.exe"