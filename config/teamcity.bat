echo Build_Revision = '%env.BUILD_VCS_NUMBER%' > ../version.dat
"%VS90COMNTOOLS%..\IDE\devenv.exe" /rebuild "Release|Win32" XBeach.sln