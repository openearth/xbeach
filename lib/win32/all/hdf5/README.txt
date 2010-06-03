This directory contains the binary (release) distribution of 
HDF5 1.8.1 that was compiled on an Intel PC running Windows XP,
using Visual Studio 2005/Intel Fortran 9.1.  It was built with
the following options: 
    -- C/C++/Fortran libraries, both static and shared
    -- SZIP (encoder enabled) and ZLIB
    -- Static and shared HDF5 tools

The contents of this directory are:

    COMPILE.txt         - Instructions for using binaries
    COPYING.txt         - Copyright notice
    INSTALL_Windows.txt - Install instructions from the source code.
                          Section IV discusses how to compile an
                          application
    README.txt          - This file
    RELEASE.txt         - Detailed information regarding this release
    bin\                - HDF5 static Utilities
    bindll\             - HDF5 shared Utilities
    dll\                - HDF5 dlls 
    include\            - HDF5 include files
    lib\                - HDF5 libraries 
    modsdll\            - HDF5 Fortran modules 
  

The binaries were built with ZLIB compression enabled (zlib 1.2.3). Therefore 
you MUST link with the GNU ZLIB library when compiling with these binaries. We 
provide the pre-compiled binary distribution for Windows for ZLIB 1.2.3 from 
our ftp server at:

    ftp://ftp.hdfgroup.org/lib-external/zlib/1.2/bin/windows

These binaries were also built with SZIP 2.1 compression (encoder ENABLED).  
When compiling with these binaries, you must link with the SZIP library.  
Information regarding SZIP can be found at:

   http://www.hdfgroup.org/doc_resource/SZIP/

You can obtain the SZIP source code and binaries from our ftp server at:
  
   ftp://ftp.hdfgroup.org/lib-external/szip/2.1/bin/windows/

Source code can be found on the THG ftp server in:

   ftp://ftp.hdfgroup.org/HDF5/current/src

Please send questions, comments, and suggestions to:

    help@hdfgroup.org

