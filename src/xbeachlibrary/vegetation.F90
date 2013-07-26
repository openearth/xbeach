!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Copyright (C) 2011 UNESCO-IHE, WL|Delft Hydraulics and Delft University !
! Dano Roelvink, Ap van Dongeren, Ad Reniers, Jamie Lescinski,            !
! Jaap van Thiel de Vries, Robert McCall                                  !       
!                                                                         !
! d.roelvink@unesco-ihe.org                                               !
! UNESCO-IHE Institute for Water Education                                !
! P.O. Box 3015                                                           !
! 2601 DA Delft                                                           !
! The Netherlands                                                         !
!                                                                         !
! This library is free software; you can redistribute it and/or           !
! modify it under the terms of the GNU Lesser General Public              !
! License as published by the Free Software Foundation; either            !
! version 2.1 of the License, or (at your option) any later version.      !
!                                                                         !
! This library is distributed in the hope that it will be useful,         !
! but WITHOUT ANY WARRANTY; without even the implied warranty of          !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        !
! Lesser General Public License for more details.                         !
!                                                                         !
! You should have received a copy of the GNU Lesser General Public        !
! License along with this library; if not, write to the Free Software     !
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307     !
! USA                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vegetation_module
implicit none
private
type vegie
   character*256                      :: name
   integer ,                  pointer :: nsec        ! Number of vegetation stands per unit horizontal area [m-2]
   real*8 , dimension(:)    , pointer :: ah          ! Vegetation height [m] 
   real*8 , dimension(:)    , pointer :: Cd          ! Vertically integrated drag coefficient [-]
   real*8 , dimension(:)    , pointer :: bv          ! Width of vegetation stands
   real*8 , dimension(:)    , pointer :: Dragterm    ! Number of vegetation stands per unit horizontal area [m-2]
   integer , dimension(:)   , pointer :: N           ! Number of vegetation stands per unit horizontal area [m-2]
   integer , dimension(:,:) , pointer :: vegtype     ! spatial mapping of vegetation types [-]
end type

type(vegie), dimension(:), allocatable, save            :: veg

public vegie_init
public swvegatt
public lwvegatt

contains

subroutine vegie_init(s,par)
    use params
    use spaceparams
    use readkey_module
    use filefunctions
    use logging_module
    use interp

    IMPLICIT NONE

    type(parameters)                            :: par
    type(spacepars)                             :: s

    character(1)                                :: ch
    integer                                     :: i,j,fid,ier
    integer, save                               :: nsec           ! Number of vegetation sections in the vertical per specie [-]
    !integer, dimension(s%nx+1,s%ny+1)           :: vegtype
    
    include 's.ind'
    include 's.inp'
    
    ! Read files with vegetation properties:
    ! file 1: list of species
    ! file 2: vegetation properties per specie (could be multiple files)
    ! file 3: distribution of species oevr space
    
    fid=create_new_fid() ! see filefunctions.F90
    call check_file_exist(par%vegiefile)
    open(fid,file=par%vegiefile)     
    ier=0
    i=0
    do while (ier==0)
       read(fid,'(a)',iostat=ier)ch
       if (ier==0)i=i+1
    enddo
    par%nveg=i                        
    rewind(fid)   
      
    allocate(veg(par%nveg))
    do i=1,par%nveg
       read(fid,'(a)')veg(i)%name      ! set vegetation name 
    enddo
    close(fid)
    
    ! allocate and read specie specific vegetation properties
    do i=1,par%nveg  ! for each specie
       call check_file_exist(veg(i)%name)
       allocate (veg(i)%nsec) 
       veg(i)%nsec = readkey_int(veg(i)%name,'nsec',  1,        1,      10  )
       allocate (veg(i)%ah(veg(i)%nsec)) 
       allocate (veg(i)%Cd(veg(i)%nsec))
       allocate (veg(i)%bv(veg(i)%nsec))
       allocate (veg(i)%N(veg(i)%nsec))
       allocate (veg(i)%Dragterm(veg(i)%nsec))
       veg(i)%ah   = readkey_dblvec(veg(i)%name,'ah',nsec,size(veg(i)%ah),0.1d0,0.05d0,20.d0)       ! Think about default values...
       veg(i)%Cd   = readkey_dblvec(veg(i)%name,'Cd',nsec,size(veg(i)%Cd),0.1d0,0.05d0,20.d0)       ! Think about default values... 
       veg(i)%bv   = readkey_dblvec(veg(i)%name,'bv',nsec,size(veg(i)%bv),0.1d0,0.05d0,20.d0)       ! Think about default values...
       veg(i)%N    = nint(readkey_dblvec(veg(i)%name,'N',nsec,size(veg(i)%bv),0.1d0,0.05d0,20.d0))  ! Jaap: quick&dirt transform real into integer
       veg(i)%Dragterm = 0.5d0/par%g/par%px*veg(i)%Cd*veg(i)%bv*veg(i)%N ! Drag coefficient based on first part equation 6.5 Suzuki, 2011
    enddo
    
    ! read spatial distribution of species:
    ! vegtype = 1 corresponds to first vegetation specified in vegiefile
    allocate (veg(1)%vegtype(s%nx+1,s%ny+1))
    fid=create_new_fid() ! see filefunctions.F90
    call check_file_exist(par%vegiemapfile)
    open(fid,file=par%vegiemapfile)
    do j=1,s%ny+1    ! Is this the right way to do it in a module
       read(fid,*,iostat=ier)(veg(1)%vegtype,i=1,s%nx+1)
       if (ier .ne. 0) then
          !Jaap doesn't work
          !call report_file_read_error(par%vegiemapfile)
       endif
    end do
    close(fid)
    
    !distrubue vegtype over species
    !allocate (veg(i)%vegtype(s%nx+1,s%ny+1))
    !do i=1,par%nveg
    !   where (vegtype == i)
    !      veg(i)%vegtype = i
    !   endwhere
    !enddo   
  
end subroutine vegie_init  

subroutine swvegatt(s,par)
    use params
    use spaceparams
    use readkey_module
    use filefunctions
    use interp
    
    type(parameters)                            :: par
    type(spacepars)                             :: s

    integer                                     :: i,j,m,ind  ! indices of actual x,y point
    real*8                                      :: aht,hterm,htermold,Dvgt
    real*8, dimension(s%nx+1,s%ny+1)            :: Dvg,kmr
    integer, dimension(s%nx+1,s%ny+1)           :: vegtype    
    
    include 's.ind'
    include 's.inp'
       
    kmr = min(max(s%k, 0.01d0), 100.d0)      
    
    !vegtype = 0
    !do i = 1,par%nveg
    !   vegtype = vegtype + veg(i)%vegtype;
    !enddo
    
    ! Ser dissipation in vegetation to zero everywhere for a start
    Dvg = 0.d0
    do j=1,ny+1
       do i=1,nx+1
          ind = vegtype(i,j)
          htermold = 0.d0
          do m=1,veg(ind)%nsec
             ! restrict vegetation height to water depth
             aht = min(veg(ind)%ah(m),hh(i,j)) 
             ! compute dissipation based on aht
             hterm = sinh(kmr(i,j)*aht)**3+3*sinh(kmr(i,j)*aht)/(3.d0*kmr(i,j)*cosh(kmr(i,j)*hh(i,j))**3)
             ! correct for lower elevated dissipation layers 
             Dvgt = veg(ind)%Dragterm(m)*(kmr(i,j)*par%g/s%sigm(i,j)/2.d0)**3*(hterm-htermold)*s%H(i,j)**3
             ! update htermold
             htermold = hterm
             ! set dissipation per layer
             Dvg(i,j) = Dvg(i,j) + Dvgt
          enddo
       enddo
    enddo
    ! store dissipation due to vegetation in s%
    s%Dveg = Dvg
    
end subroutine swvegatt

subroutine lwvegatt(s,par)
    use params
    use spaceparams
    use readkey_module
    use filefunctions
    use interp
    
    type(parameters)                            :: par
    type(spacepars)                             :: s
    type(vegie), dimension(:), pointer          :: veg

    !integer                                     :: 
    !real*8                                      :: 
    
    include 's.ind'
    include 's.inp'
    
end subroutine lwvegatt
    
    
 
end module vegetation_module