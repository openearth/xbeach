@echo off

SET exe=C:\Program Files\TortoiseSVN\bin\SubWCRev.exe
SET ok=0

IF EXIST "%exe%" (
  "%exe%" "%cd%" "%cd%\includes\version\version.tpl" "%cd%\includes\version\version.dat"
  IF ERRORLEVEL 1 (SET ok=0) ELSE (SET ok=1)
)

IF %ok%==0 (
  copy "%cd%\includes\version\version.000" "%cd%\includes\version\version.dat"
)