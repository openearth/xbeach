echo Build_Revision = '%1' > ../version.dat
"%VS90COMNTOOLS%..\IDE\devenv.exe" /rebuild "%2|%3" XBeach.sln