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
# ADD F90 /browser /check:bounds /check:power /check:overflow /compile_only /fpp /nologo /real_size:64 /traceback /warn:nofileopt
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
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /stack:0x5f5e100 /subsystem:console /machine:I386
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
# ADD F90 /browser /check:bounds /check:power /check:overflow /compile_only /dbglibs /debug:full /fltconsistency /fpconstant /fpe:0 /fpp /nologo /real_size:64 /traceback /warn:argument_checking /warn:nofileopt
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

!ENDIF 

# Begin Target

# Name "xbeach - Win32 Release"
# Name "xbeach - Win32 Debug"
# Begin Source File

SOURCE=.\boundaryconditions.f90
DEP_F90_BOUND=\
	".\interp.mod"\
	".\params.mod"\
	".\s.ind"\
	".\s.inp"\
	".\spaceparams.mod"\
	".\wave_timestep_module.mod"\
	".\waveparams.mod"\
	".\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\constants.f90
# End Source File
# Begin Source File

SOURCE=.\flow_timestep.f90
DEP_F90_FLOW_=\
	".\params.mod"\
	".\s.ind"\
	".\s.inp"\
	".\spaceparams.mod"\
	".\xmpi_module.mod"\
	
NODEP_F90_FLOW_=\
	".\Release\mpi.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\general_mpi.F90
DEP_F90_GENER=\
	".\xmpi_module.mod"\
	
NODEP_F90_GENER=\
	".\Release\mpi.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\initialize.f90
DEP_F90_INITI=\
	".\interp.mod"\
	".\params.mod"\
	".\readkey_module.mod"\
	".\spaceparams.mod"\
	".\wave_timestep_module.mod"\
	".\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\interp.f90
# End Source File
# Begin Source File

SOURCE=.\math_tools.f90
DEP_F90_MATH_=\
	".\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\morphevolution.f90
DEP_F90_MORPH=\
	".\params.mod"\
	".\s.ind"\
	".\s.inp"\
	".\spaceparams.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\params.f90
DEP_F90_PARAM=\
	".\readkey_module.mod"\
	".\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\readkey.f90
DEP_F90_READK=\
	".\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\readtide.f90
DEP_F90_READT=\
	".\params.mod"\
	".\readkey_module.mod"\
	".\s.ind"\
	".\s.inp"\
	".\spaceparams.mod"\
	".\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\roelvink.f90
DEP_F90_ROELV=\
	".\params.mod"\
	".\spaceparams.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\spaceparams.F90
DEP_F90_SPACE=\
	".\general_mpi_module.mod"\
	".\params.mod"\
	".\readkey_module.mod"\
	".\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\timestep.f90
DEP_F90_TIMES=\
	".\params.mod"\
	".\spaceparams.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\varoutput.f90
DEP_F90_VAROU=\
	".\params.mod"\
	".\readkey_module.mod"\
	".\spaceparams.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\wave_stationary.f90
DEP_F90_WAVE_=\
	".\params.mod"\
	".\roelvink_module.mod"\
	".\s.ind"\
	".\s.inp"\
	".\spaceparams.mod"\
	".\wave_timestep_module.mod"\
	".\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\wave_timestep.f90
DEP_F90_WAVE_T=\
	".\params.mod"\
	".\roelvink_module.mod"\
	".\s.ind"\
	".\s.inp"\
	".\spaceparams.mod"\
	".\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\waveparams.f90
DEP_F90_WAVEP=\
	".\interp.mod"\
	".\math_tools.mod"\
	".\params.mod"\
	".\readkey_module.mod"\
	".\spaceparams.mod"\
	".\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\xbeach.f90
DEP_F90_XBEAC=\
	".\boundaryconditions.mod"\
	".\flow_timestep_module.mod"\
	".\initialize.mod"\
	".\morphevolution.mod"\
	".\outputmod.mod"\
	".\params.mod"\
	".\readtide_module.mod"\
	".\spaceparams.mod"\
	".\timestep_module.mod"\
	".\wave_stationary_module.mod"\
	".\wave_timestep_module.mod"\
	".\xmpi_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\xmpi.F90
NODEP_F90_XMPI_=\
	".\Release\mpi.mod"\
	
# End Source File
# End Target
# End Project
