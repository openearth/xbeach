!==============================================================================
!                               NH_PARAMETERS                        
!==============================================================================
  integer, parameter   :: rKind        = 8       !Default precision for reals
  integer, parameter   :: iKind        = 4       !Default precision for integers
  integer, parameter   :: iPrint       = 1000    !Unit number print file
  integer, parameter   :: iScreen      = 6       !Unit number screen
  
  integer, parameter   :: iUnit_U      = 1001    !Unit number u boundary file
  integer, parameter   :: iUnit_ZsFile   = 1002  !Unit number for reading initial zs
  integer, parameter   :: iFileNameLen = 80      !Length of filenames
  
  character(len=iFilenameLen), parameter :: cParamsFile = 'params.txt'
  character(len=iFilenameLen), parameter :: cZsFile     = 'zsfile.ini'