# Microsoft Developer Studio Project File - Name="xbeach" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=xbeach - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "xbeach.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "xbeach.mak" CFG="xbeach - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "xbeach - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "xbeach - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE "xbeach - Win32 Release_NoMPI" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "xbeach - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /browser /check:bounds /check:power /check:overflow /compile_only /define:"BARRIER" /define:"USEMPI" /define:"HAVE_MPI_WTIME" /fpp /include:"c:\Program Files\MPICH2\include\\" /nologo /real_size:64 /traceback /warn:nofileopt
# SUBTRACT F90 /check:underflow
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x413 /d "NDEBUG"
# ADD RSC /l 0x413 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib fmpich2s.lib mpi.lib /nologo /stack:0x5f5e100 /subsystem:console /machine:I386 /libpath:"c:\Program Files\MPICH2\lib\\"
# SUBTRACT LINK32 /pdb:none

!ELSEIF  "$(CFG)" == "xbeach - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /browser /check:bounds /check:power /check:overflow /compile_only /dbglibs /debug:full /fltconsistency /fpconstant /fpe:0 /fpp /nologo /real_size:64 /traceback /warn:argument_checking /warn:declarations /warn:nofileopt /warn:unused
# SUBTRACT F90 /check:underflow
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x413 /d "_DEBUG"
# ADD RSC /l 0x413 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /stack:0x3d0900 /subsystem:console /debug /machine:I386 /pdbtype:sept
# SUBTRACT LINK32 /pdb:none

!ELSEIF  "$(CFG)" == "xbeach - Win32 Release_NoMPI"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "xbeach___Win32_Release_NoMPI"
# PROP BASE Intermediate_Dir "xbeach___Win32_Release_NoMPI"
# PROP BASE Ignore_Export_Lib 0
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release_NoMPI"
# PROP Intermediate_Dir "Release_NoMPI"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /browser /check:bounds /check:power /check:overflow /compile_only /define:"BARRIER" /define:"USEMPI" /define:"HAVE_MPI_WTIME" /fpp /include:"c:\Program Files\MPICH2\include\\" /nologo /real_size:64 /traceback /warn:nofileopt
# SUBTRACT BASE F90 /check:underflow
# ADD F90 /browser /check:bounds /check:power /check:overflow /compile_only /fpp /list /nologo /real_size:64 /show:include /traceback /warn:nofileopt
# SUBTRACT F90 /check:underflow
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x413 /d "NDEBUG"
# ADD RSC /l 0x413 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib fmpich2s.lib mpi.lib /nologo /stack:0x5f5e100 /subsystem:console /machine:I386 /libpath:"c:\Program Files\MPICH2\lib\\"
# SUBTRACT BASE LINK32 /pdb:none
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /stack:0x5f5e100 /subsystem:console /machine:I386
# SUBTRACT LINK32 /pdb:none

!ENDIF 

# Begin Target

# Name "xbeach - Win32 Release"
# Name "xbeach - Win32 Debug"
# Name "xbeach - Win32 Release_NoMPI"
# Begin Source File

SOURCE=.\boundaryconditions.f90
DEP_F90_BOUND=\
	".\nh_pars.inc"\
	".\s.ind"\
	".\s.inp"\
	".\space_ind.gen"\
	".\space_inp.gen"\
	
NODEP_F90_BOUND=\
	".\Release\interp.mod"\
	".\Release\params.mod"\
	".\Release\readkey_module.mod"\
	".\Release\spaceparams.mod"\
	".\Release\wave_timestep_module.mod"\
	".\Release\waveparams.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\constants.f90
# End Source File
# Begin Source File

SOURCE=.\drifters.f90
DEP_F90_DRIFT=\
	".\s.ind"\
	".\s.inp"\
	".\space_ind.gen"\
	".\space_inp.gen"\
	
NODEP_F90_DRIFT=\
	".\Release\mnemmodule.mod"\
	".\Release\params.mod"\
	".\Release\readkey_module.mod"\
	".\Release\spaceparams.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\flow_secondorder.F90
DEP_F90_FLOW_=\
	".\nh_pars.inc"\
	
NODEP_F90_FLOW_=\
	".\Release\params.mod"\
	".\Release\spaceparams.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\flow_timestep.f90
DEP_F90_FLOW_T=\
	".\s.ind"\
	".\s.inp"\
	".\space_ind.gen"\
	".\space_inp.gen"\
	
NODEP_F90_FLOW_T=\
	".\Release\boundaryconditions.mod"\
	".\Release\flow_secondorder_module.mod"\
	".\Release\nonh_module.mod"\
	".\Release\params.mod"\
	".\Release\spaceparams.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\general_fileio.F90
NODEP_F90_GENER=\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\general_mpi.F90
DEP_F90_GENERA=\
	"c:\Program Files\MPICH2\include\mpi.mod"\
	
NODEP_F90_GENERA=\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\groundwater.f90
NODEP_F90_GROUN=\
	".\Release\params.mod"\
	".\Release\readkey_module.mod"\
	".\Release\spaceparams.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\initialize.f90
NODEP_F90_INITI=\
	".\Release\interp.mod"\
	".\Release\params.mod"\
	".\Release\readkey_module.mod"\
	".\Release\spaceparams.mod"\
	".\Release\wave_timestep_module.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\interp.f90
# End Source File
# Begin Source File

SOURCE=.\math_tools.f90
NODEP_F90_MATH_=\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mnemonic.F90
DEP_F90_MNEMO=\
	".\chartoindex.gen"\
	".\mnemonic.gen"\
	
# End Source File
# Begin Source File

SOURCE=.\morphevolution.f90
DEP_F90_MORPH=\
	".\s.ind"\
	".\s.inp"\
	".\space_ind.gen"\
	".\space_inp.gen"\
	
NODEP_F90_MORPH=\
	".\Release\interp.mod"\
	".\Release\params.mod"\
	".\Release\readkey_module.mod"\
	".\Release\spaceparams.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\nonh.F90
DEP_F90_NONH_=\
	".\nh_pars.inc"\
	
NODEP_F90_NONH_=\
	".\Release\flow_secondorder_module.mod"\
	".\Release\params.mod"\
	".\Release\solver_module.mod"\
	".\Release\spaceparams.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\nonhydrostatic.F90
NODEP_F90_NONHY=\
	".\Release\params.mod"\
	".\Release\spaceparams.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\params.f90
NODEP_F90_PARAM=\
	".\Release\general_fileio.mod"\
	".\Release\readkey_module.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\readkey.f90
NODEP_F90_READK=\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\readtide.f90
DEP_F90_READT=\
	".\s.ind"\
	".\s.inp"\
	".\space_ind.gen"\
	".\space_inp.gen"\
	
NODEP_F90_READT=\
	".\Release\params.mod"\
	".\Release\readkey_module.mod"\
	".\Release\spaceparams.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\readwind.F90
DEP_F90_READW=\
	".\s.ind"\
	".\s.inp"\
	".\space_ind.gen"\
	".\space_inp.gen"\
	
NODEP_F90_READW=\
	".\Release\params.mod"\
	".\Release\readkey_module.mod"\
	".\Release\spaceparams.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\roelvink.f90
NODEP_F90_ROELV=\
	".\Release\params.mod"\
	".\Release\spaceparams.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\solver.F90
DEP_F90_SOLVE=\
	".\nh_pars.inc"\
	
NODEP_F90_SOLVE=\
	".\Release\params.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\spaceparams.F90
DEP_F90_SPACE=\
	".\indextos.gen"\
	".\space_alloc_arrays.gen"\
	".\space_alloc_scalars.gen"\
	".\spacedecl.gen"\
	
NODEP_F90_SPACE=\
	".\Release\general_mpi_module.mod"\
	".\Release\mnemmodule.mod"\
	".\Release\params.mod"\
	".\Release\readkey_module.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\timestep.f90
NODEP_F90_TIMES=\
	".\Release\params.mod"\
	".\Release\spaceparams.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\varianceupdate.F90
NODEP_F90_VARIA=\
	".\Release\mnemmodule.mod"\
	".\Release\params.mod"\
	".\Release\spaceparams.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\varoutput.f90
NODEP_F90_VAROU=\
	".\Release\mnemmodule.mod"\
	".\Release\params.mod"\
	".\Release\readkey_module.mod"\
	".\Release\spaceparams.mod"\
	".\Release\varianceupdate_module.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\wave_stationary.f90
DEP_F90_WAVE_=\
	".\s.ind"\
	".\s.inp"\
	".\space_ind.gen"\
	".\space_inp.gen"\
	
NODEP_F90_WAVE_=\
	".\Release\params.mod"\
	".\Release\roelvink_module.mod"\
	".\Release\spaceparams.mod"\
	".\Release\wave_timestep_module.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\wave_timestep.f90
DEP_F90_WAVE_T=\
	".\s.ind"\
	".\s.inp"\
	".\space_ind.gen"\
	".\space_inp.gen"\
	
NODEP_F90_WAVE_T=\
	".\Release\mnemmodule.mod"\
	".\Release\params.mod"\
	".\Release\roelvink_module.mod"\
	".\Release\spaceparams.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\waveparams.f90
NODEP_F90_WAVEP=\
	".\Release\interp.mod"\
	".\Release\math_tools.mod"\
	".\Release\params.mod"\
	".\Release\readkey_module.mod"\
	".\Release\spaceparams.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\xbeach.f90
DEP_F90_XBEAC=\
	".\version.dat"\
	".\version.def"\
	
NODEP_F90_XBEAC=\
	".\Release\boundaryconditions.mod"\
	".\Release\config.h"\
	".\Release\drifter_module.mod"\
	".\Release\flow_timestep_module.mod"\
	".\Release\groundwaterflow.mod"\
	".\Release\initialize.mod"\
	".\Release\morphevolution.mod"\
	".\Release\outputmod.mod"\
	".\Release\params.mod"\
	".\Release\readkey_module.mod"\
	".\Release\readtide_module.mod"\
	".\Release\readwind_module.mod"\
	".\Release\spaceparams.mod"\
	".\Release\timestep_module.mod"\
	".\Release\wave_stationary_module.mod"\
	".\Release\wave_timestep_module.mod"\
	".\Release\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\xmpi.F90
DEP_F90_XMPI_=\
	"c:\Program Files\MPICH2\include\mpi.mod"\
	
# End Source File
# End Target
# End Project
