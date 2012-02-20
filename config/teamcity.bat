echo off
rem %1 path to devenv exe

echo solution file	: %2
echo Configuration name	: %3
echo Project name	: %4

echo clean build
%1\devenv.exe %2 /Clean %3

echo xbeach build solution
%1\devenv.exe %2 /Build %3 /Project %4
