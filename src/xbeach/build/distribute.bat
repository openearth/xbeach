echo off

rem 1: SolutionDir
rem 2: ConfigurationName
rem 3: TargetDir
rem 4: platform (win32 \ x64)

echo on

rem clear output dir
if exist %1\dist\%4\%2\ rmdir %1\dist\%4\%2\ /q /s
mkdir %1\dist\%4\%2\

rem copy exe
copy %3\*.exe %1\dist\%4\%2\
copy %3\*.dll %1\dist\%4\%2\

rem copy additional dll files from static output folder
copy %1\src\xbeachlibrary\bin\static\%4\%2\*.dll %1\dist\%4\%2\

rem copy additional dll files from static to bin folder
copy %1\src\xbeachlibrary\bin\static\%4\%2\*.dll %3\

rem copy current manual
copy %1\doc\manual\xbeach_manual.pdf %1\dist\%4\%2\

rem copy license
copy %1\LICENSE %1\dist\%4\%2\
