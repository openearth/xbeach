module paramsconst
   implicit none
   save
   integer, parameter :: TURB_NONE                   =  0
   integer, parameter :: TURB_BORE_AVERAGED          =  1
   integer, parameter :: TURB_WAVE_AVERAGED          =  2

   integer, parameter :: BREAK_ROELVINK1             =  1
   integer, parameter :: BREAK_BALDOCK               =  2
   integer, parameter :: BREAK_ROELVINK2             =  3
   integer, parameter :: BREAK_ROELVINK_DALY         =  4
   integer, parameter :: BREAK_JANSSEN               =  5

   integer, parameter :: LATERALWAVE_NEUMANN         =  1
   integer, parameter :: LATERALWAVE_WAVECREST       =  2
   integer, parameter :: LATERALWAVE_CYCLIC          =  3

   integer, parameter :: LEFTWAVE_NEUMANN            =  1
   integer, parameter :: LEFTWAVE_WAVECREST          =  2

   integer, parameter :: RIGHTWAVE_NEUMANN           =  1
   integer, parameter :: RIGHTWAVE_WAVECREST         =  2

   integer, parameter :: INSTAT_STAT                 =  0
   integer, parameter :: INSTAT_BICHROM              =  1
   integer, parameter :: INSTAT_TS_1                 =  2
   integer, parameter :: INSTAT_TS_2                 =  3
   integer, parameter :: INSTAT_JONS                 =  4
   integer, parameter :: INSTAT_SWAN                 =  5
   integer, parameter :: INSTAT_VARDENS              =  6
   integer, parameter :: INSTAT_REUSE                =  7
   integer, parameter :: INSTAT_TS_NONH              =  8
   integer, parameter :: INSTAT_OFF                  =  9
   integer, parameter :: INSTAT_STAT_TABLE           =  10
   integer, parameter :: INSTAT_JONS_TABLE           =  11

   integer, parameter :: GRIDFORM_XBEACH             =  1
   integer, parameter :: GRIDFORM_DELFT3D            =  2

   integer, parameter :: FRONT_ABS_1D                =  0
   integer, parameter :: FRONT_ABS_2D                =  1
   integer, parameter :: FRONT_WALL                  =  2
   integer, parameter :: FRONT_WLEVEL                =  3
   integer, parameter :: FRONT_NONH_1D               =  4
   integer, parameter :: FRONT_WAVEFLUME             =  5

   integer, parameter :: LR_NEUMANN                  =  0
   integer, parameter :: LR_WALL                     =  1
   integer, parameter :: LR_NO_ADVEC                 =  2
   integer, parameter :: LR_NEUMANN_V                =  3

   integer, parameter :: BACK_WALL                   =  0
   integer, parameter :: BACK_ABS_1D                 =  1
   integer, parameter :: BACK_ABS_2D                 =  2
   integer, parameter :: BACK_WLEVEL                 =  3

   integer, parameter :: TIDETYPE_INSTANT            =  0
   integer, parameter :: TIDETYPE_VELOCITY           =  1

   integer, parameter :: PAULREVERE_LAND             =  0
   integer, parameter :: PAULREVERE_SEA              =  1

   integer, parameter :: BEDFRICTION_CHEZY           =  0
   integer, parameter :: BEDFRICTION_CF              =  1
   integer, parameter :: BEDFRICTION_WHITE_COLEBROOK =  2
   integer, parameter :: BEDFRICTION_MANNING         =  3
   integer, parameter :: BEDFRICTION_WHITE_COLEBROOK_GRAINSIZE =  4

   integer, parameter :: GWSCHEME_LAMINAR            =  0
   integer, parameter :: GWSCHEME_TURBULENT          =  1

   integer, parameter :: GWHEADMODEL_PARABOLIC       =  0
   integer, parameter :: GWHEADMODEL_EXPONENTIAL     =  1

   integer, parameter :: SOLVER_SIPP                 =  0
   integer, parameter :: SOLVER_TRIDIAGG             =  1

   integer, parameter :: FORM_SOULSBY_VANRIJN        =  0
   integer, parameter :: FORM_VANTHIEL_VANRIJN       =  1

   integer, parameter :: WAVEFORM_RUESSINK_VANRIJN   =  0
   integer, parameter :: WAVEFORM_VANTHIEL           =  1

   integer, parameter :: OUTPUTFORMAT_FORTRAN        =  0
   integer, parameter :: OUTPUTFORMAT_NETCDF         =  1
   integer, parameter :: OUTPUTFORMAT_DEBUG          =  2

   integer, parameter :: SCHEME_UPWIND_1             =  0
   integer, parameter :: SCHEME_LAX_WENDROFF         =  1
   integer, parameter :: SCHEME_UPWIND_2             =  2

   integer, parameter :: TURBADV_NONE                =  0
   integer, parameter :: TURBADV_LAGRANGIAN          =  1
   integer, parameter :: TURBADV_EULERIAN            =  2

   integer, parameter :: BDSLPEFFMAG_NONE            =  0
   integer, parameter :: BDSLPEFFMAG_ROELV_TOTAL     =  1
   integer, parameter :: BDSLPEFFMAG_ROELV_BED       =  2
   integer, parameter :: BDSLPEFFMAG_SOULS_TOTAL     =  3
   integer, parameter :: BDSLPEFFMAG_SOULS_BED       =  4

   integer, parameter :: BDSLPEFFINI_NONE            =  0
   integer, parameter :: BDSLPEFFINI_TOTAL           =  1
   integer, parameter :: BDSLPEFFINI_BED             =  2

   integer, parameter :: BDSLPEFFDIR_NONE            =  0
   integer, parameter :: BDSLPEFFDIR_TALMON          =  1
   
   integer, parameter :: OUTPUTPRECISION_SINGLE      =  0
   integer, parameter :: OUTPUTPRECISION_DOUBLE      =  1

end module paramsconst
