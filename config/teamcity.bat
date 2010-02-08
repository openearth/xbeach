echo Build_Revision = '%env.BUILD_VCS_NUMBER%' //  ' Mixed revisions' // 'modifications  M' > ../version.dat
"%VS90COMNTOOLS%..\IDE\devenv.exe" /rebuild "Release|Win32" XBeach.sln