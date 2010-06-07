echo Build_Revision = '%3' > ../version.dat
echo Build_Date     = '%date% %time%' >> ../version.dat
echo Build_URL      = '%4' >> ../version.dat

"%VS90COMNTOOLS%..\IDE\devenv.exe" /rebuild "%1|%2" XBeach.sln