@echo off

echo ##teamcity[progressStart 'Build started']

echo Build_Revision = '%3' > ../version.dat
echo Build_Date     = '%date% %time%' >> ../version.dat
echo Build_URL      = '%4' >> ../version.dat

"%VS90COMNTOOLS%..\IDE\devenv.exe" /rebuild "%1|%2" XBeach.sln

IF ERRORLEVEL 1 (
	echo ##teamcity[buildStatus status='FAILED']
) ELSE (
	echo ##teamcity[buildStatus status='SUCCES']
)

echo ##teamcity[progressStart 'Build finished']

exit %ERRORLEVEL%