echo off

echo solution dir	: %1

echo ##teamcity[progressMessage 'Removing old builds']
IF EXIST "%1\dist\win32" (
rmdir "%1\dist\win32" /s /q
mkdir "%1\dist\win32"
)
IF EXIST "%1\dist\x64" (
rmdir "%1\dist\x64" /s /q
mkdir "%1\dist\x64"
)