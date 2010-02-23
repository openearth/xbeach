echo Build_Revision = '%1' > ../version.dat
"%VS90COMNTOOLS%..\IDE\devenv.exe" /rebuild "Release|Win32" XBeach.sln