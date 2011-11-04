!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
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
module ship_module
implicit none
type ship
   character*256                   :: name
   real*8                          :: dx
   real*8                          :: dy
   integer                         :: nx
   integer                         :: ny
   character*256                   :: shipgeom
   real*8, dimension(:,:), pointer :: depth
   character*256                   :: shiptrack
   integer                         :: track_nt
   real*8, dimension(:)  , pointer :: track_t
   real*8, dimension(:)  , pointer :: track_x
   real*8, dimension(:)  , pointer :: track_y
   real*8, dimension(:)  , pointer :: track_dir
end type

contains



  subroutine shipwave(s,par)
    use params
    use xmpi_module
    use spaceparams
    use readkey_module
    use filefunctions
    use interp

    IMPLICIT NONE

    type(parameters)                            :: par
    type(spacepars)                             :: s
    type(ship), dimension(:), pointer, save     :: sh

    integer                                     :: i,j,fid,ix,iy,ier,it,shp_indx,i1,j1
    logical, save                               :: firstship=.true.
    integer, save                               :: nships
    character(1)                                :: ch
    real*8                                      :: xship,yship,dirship,radius,cosdir,sindir
    real*8                                      :: x1,y1,xrel,yrel

    include 's.ind'
    include 's.inp'
    
    if (firstship) then
    
    ! Read ship names (== filenames with ship geometry and track data)

      fid=create_new_fid()
      open(fid,file=par%shipfile)
      ier=0
      i=0
      do while (ier==0)
         read(fid,'(a)',iostat=ier)ch
         if (ier==0)i=i+1
      enddo
      nships=i
      rewind(fid)
      
      allocate(sh(nships))
      do i=1,nships
        read(fid,'(a)')sh(i)%name 
      enddo
      close(fid)
      
      do i=1,nships
      ! Read ship geometry
        sh(i)%dx  = readkey_dbl(sh(i)%name,'dx',  5.d0,   0.d0,      100.d0)
        sh(i)%dy  = readkey_dbl(sh(i)%name,'dy',  5.d0,   0.d0,      100.d0)
        sh(i)%nx  = readkey_int(sh(i)%name,'nx',  20,        1,      1000  )
        sh(i)%ny  = readkey_int(sh(i)%name,'ny',  20,        1,      1000  )
        sh(i)%shipgeom = readkey_name(sh(i)%name,'shipgeom',required=.true.)
        sh(i)%shiptrack = readkey_name(sh(i)%name,'shiptrack',required=.true.)

        allocate (sh(i)%depth(sh(i)%nx+1,sh(i)%ny+1))
        fid=create_new_fid()
        open(fid,file=sh(i)%shipgeom)
        do iy=1,sh(i)%ny+1
           read(fid,*)(sh(i)%depth(ix,iy),ix=1,sh(i)%nx+1)
        end do
        close(fid)

      ! Read t,x,y of ship position

        fid=create_new_fid()
        open(fid,file=sh(i)%shiptrack)
        ier=0
        it=0
        do while (ier==0)
           read(fid,'(a)',iostat=ier)ch
           if (ier==0)it=it+1
        enddo
        sh(i)%track_nt=it
        rewind(fid)
        allocate(sh(i)%track_t(sh(i)%track_nt))
        allocate(sh(i)%track_x(sh(i)%track_nt))
        allocate(sh(i)%track_y(sh(i)%track_nt))
        allocate(sh(i)%track_dir(sh(i)%track_nt))
        do it=1,sh(i)%track_nt
           read(fid,*)sh(i)%track_t(it),sh(i)%track_x(it),sh(i)%track_y(it)
        enddo
        close(fid)
        
     !  Compute ship direction

        sh(i)%track_dir(1)=atan2(sh(i)%track_y(2)-sh(i)%track_y(1),sh(i)%track_x(2)-sh(i)%track_x(1))
        do it=2,sh(i)%track_nt-1
           sh(i)%track_dir(it)=atan2(sh(i)%track_y(it+1)-sh(i)%track_y(it-1),sh(i)%track_x(it+1)-sh(i)%track_x(it-1))
        enddo
        it=sh(i)%track_nt
        sh(i)%track_dir(it)=atan2(sh(i)%track_y(it)-sh(i)%track_y(it-1),sh(i)%track_x(it)-sh(i)%track_x(it-1))

      enddo ! loop over ships

    endif !  firstship

    s%ph = 0.d0
    do i=1,nships

     !  Compute pressure head (m) on XBeach grid

        call linear_interp(sh(i)%track_t,sh(i)%track_x,sh(i)%track_nt,par%t,xship,shp_indx)
        call linear_interp(sh(i)%track_t,sh(i)%track_y,sh(i)%track_nt,par%t,yship,shp_indx)
        call linear_interp(sh(i)%track_t,sh(i)%track_dir,sh(i)%track_nt,par%t,dirship,shp_indx)
        radius=max(sh(i)%nx*sh(i)%dx,sh(i)%ny*sh(i)%dy)/2
        cosdir=cos(dirship)
        sindir=sin(dirship)
        do iy=1,s%ny+1
           do ix=1,s%nx+1
              x1 =  (xz(ix,iy)-xship)*cosdir + (yz(ix,iy)-yship)*sindir
              y1 = -(xz(ix,iy)-xship)*sindir + (yz(ix,iy)-yship)*cosdir
              xrel=x1/sh(i)%dx+sh(i)%nx/2
              yrel=y1/sh(i)%dy+sh(i)%ny/2
              i1=floor(xrel)
              j1=floor(yrel)
              if (i1>=0 .and. i1<SH(i)%nx .and. j1>=0 .and. j1<sh(i)%ny) then
                 s%ph(ix,iy)=s%ph(ix,iy)+(1.d0-(xrel-float(i1)))*(1.d0-(yrel-float(j1)))*sh(i)%depth(i1+1,j1+1)  &
                                        +(      xrel-float(i1) )*(1.d0-(yrel-float(j1)))*sh(i)%depth(i1+2,j1+1  )  &
                                        +(1.d0-(xrel-float(i1)))*(      yrel-float(j1) )*sh(i)%depth(i1+1,j1+2)  &
                                        +(      xrel-float(i1) )*(      yrel-float(j1) )*sh(i)%depth(i1+2,j1+2)
              endif
           enddo
        enddo
    enddo
     
    if (firstship) then

      ! apply initial condition

        s%zs=s%zs-s%ph
    endif
     
    firstship=.false.

  end subroutine shipwave

end module