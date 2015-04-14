!****h* General/Constants
!
!  SYNOPSIS
!     Module defining common constants
!
!  COPYRIGHT
!     (C) 2004 by WL|Delft Hydraulics
!
!  AUTHOR
!     Arjen Markus
!
!  DESCRIPTION
!  This module contains integer and real parameters:
!  - kind parameters for single (sp) and double (dp) precision
!
! SOURCE
!
! wwvv never used
!      this would be a nice place to hold things as pi, twopi
!      and so on

module constants
  implicit none
  save

  integer, parameter     :: sp = kind(1.0)
  integer, parameter     :: dp = kind(1.0d00)


end module constants
!****
