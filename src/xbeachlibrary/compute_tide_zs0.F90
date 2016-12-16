!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (C) 2007 UNESCO-IHE, WL|Delft Hydraulics and Delft University !
! Robert McCall, Dano Roelvink, Ap van Dongeren                           !
!                                                                         !
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
module compute_tide_module


   implicit none
   save

   real*8,dimension(4),private  :: ndistcorners, sdistcorners
   
   public tide_init         ! routine called on global xmaster for full grid to initialise tide 
   public tide_boundary_timestep  ! routine called by all processes to update tide boundary during the run
   private timeinterp_tide  ! routine called by any processes to interpolate tide time series to current timestep
   private boundaryinterp_tide ! routine called by any processes to interpolate current tide level at model corners to current
                              ! level along the model boundary
   private boundaryinterp_tide_complex ! routine called by any processes to compute tide levels on complex boundaries 
                                       ! (interpolation,  splitting, MPI, etc.)
   private fill_tide_grid     ! routine called only by xmaster process to fill an initial field with zs0 from zs0 time series
                              ! and values generated at the boundary. Copies method of boundary_tide_complex to initialise in the 
                              ! centre of the domain

contains

subroutine tide_init(s,par)
   ! This routine is called by xmaster on full grid
      use params
      use spaceparams
      use paramsconst
      use filefunctions
      use logging_module

      implicit none
      save

      type(parameters),intent(in)                 :: par
      type(spacepars),target                      :: s
      integer                                     :: i,j,fid,ier
      logical                                     :: exists
      
      ! Set global domain corner distances
      ndistcorners(1) = s%ndist(1,1)
      ndistcorners(4) = s%ndist(s%nx+1,1)
      ndistcorners(3) = s%ndist(s%nx+1,s%ny+1)
      ndistcorners(2) = s%ndist(1,s%ny+1)
      sdistcorners(1) = s%sdist(1,1)
      sdistcorners(4) = s%sdist(s%nx+1,1)
      sdistcorners(3) = s%sdist(s%nx+1,s%ny+1)
      sdistcorners(2) = s%sdist(1,s%ny+1)
      
      ! if there is a file with initial zs0, then use that
      inquire(file=par%zsinitfile,exist=exists)
      if (exists) then
         fid = create_new_fid()
         open(fid,file=par%zsinitfile)
         do j=1,s%ny+1
            read(fid,*,iostat=ier)(s%zs0(i,j),i=1,s%nx+1)
            if (ier .ne. 0) then
               call report_file_read_error(par%zsinitfile)
            endif
         enddo
         close(fid)   
      else  ! no initial zs0 file exists       
         ! Call tidal water level at four corners of the model at current time (returns zs01-zs04)
         call timeinterp_tide(s,par,0.d0)
      
         ! call tide levels at the boundaries
         call boundaryinterp_tide(s,par,globalonly=.true.,correctforwetz=.false.) ! in this case we want all to be carried out on 
                                                                                  ! global grid, not local grids, and not correct
                                                                                  ! for wetz, as this is not yet initialised
         
         ! call zs0 water levels at intermediate points
         call fill_tide_grid(s)
      endif  
      
end subroutine tide_init

subroutine tide_boundary_timestep(s,par)
      use params
      use spaceparams
#ifdef USEMPI
      use xmpi_module
#endif
      
      implicit none
      save

      type(parameters),intent(in)                 :: par
      type(spacepars),target                      :: s
      logical,save                                :: initialised = .false.
      integer                                     :: i
      
      ! First we need common data from xmaster regarding corner points (initialisation done only on xmaster)
#ifdef USEMPI
      if (.not. initialised) then
         do i=1,4
            call xmpi_bcast(ndistcorners(i))
            call xmpi_bcast(sdistcorners(i))
         enddo
         initialised = .true.
      endif
#endif            

      ! Call tidal water level at four corners of the model at current time (returns zs01-zs04)
      call timeinterp_tide(s,par,par%t)
      
      ! call tide levels at the boundaries
      call boundaryinterp_tide(s,par) ! in this case we want all to be carried out on all local grids
      
      if(par%tidetype==TIDETYPE_INSTANT) then
         call fill_tide_grid(s)
      endif
        

end subroutine tide_boundary_timestep

subroutine timeinterp_tide(s,par,timenow)
      ! This routine may be called by all processes and interpolates tide time series to current time step at all four corners of
      ! the global model grid
      use params
      use spaceparams
      use paramsconst
      use interp

      implicit none
      save

      type(parameters),intent(in)                 :: par
      type(spacepars),target                      :: s
      real*8, intent(in)                          :: timenow
      integer                                     :: indt


      select case (par%tideloc)
         case (0)
            s%zs01 = par%zs0
            s%zs02 = par%zs0
            s%zs03 = par%zs0
            s%zs04 = par%zs0
         case (1)
            ! Need to interpolate to the correct moment in time. First point in tidal
            ! record not necessarily == 0.0
            call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,1),s%tidelen,timenow, s%zs01, indt)
            s%zs02 = s%zs01
            s%zs03 = s%zs01
            s%zs04 = s%zs01
         case (2)
            select case (par%paulrevere)
               case(PAULREVERE_LAND)
                  ! One time series on bay and one on sea
                  call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,1),s%tidelen,timenow, s%zs01, indt)
                  s%zs02 = s%zs01
                  call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,2),s%tidelen,timenow, s%zs03, indt)
                  s%zs04 = s%zs03
               case(PAULREVERE_SEA)
                  ! One time series for each sea corner (copied also to bay corners)
                  call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,1),s%tidelen,timenow, s%zs01, indt)
                  s%zs03 = s%zs01
                  call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,2),s%tidelen,timenow, s%zs02, indt)
                  s%zs04 = s%zs02
            end select
         case (4)
            ! One time series per corner
            call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,1),s%tidelen,timenow, s%zs01, indt)
            call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,2),s%tidelen,timenow, s%zs02, indt)
            call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,3),s%tidelen,timenow, s%zs03, indt)
            call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,4),s%tidelen,timenow, s%zs04, indt)
      end select
         
end subroutine timeinterp_tide

subroutine boundaryinterp_tide(s,par,globalonly,correctforwetz)
      ! This routine interpolates in space the current tide at the four corners of the model to the (local MPI subdomain) grid 
      ! boundary points
      use params
      use spaceparams
      use paramsconst
      use interp
      use xmpi_module

      implicit none
      save

      type(parameters),intent(in)                 :: par
      type(spacepars),target                      :: s
      logical,intent(in),optional                 :: globalonly,correctforwetz
      logical                                     :: globalonly_local,correctforwetz_local
      integer                                     :: i,j
      
      if(present(globalonly)) then
         globalonly_local = globalonly
      else
         globalonly_local = .false.
      endif
      
      if(present(correctforwetz)) then
         correctforwetz_local = correctforwetz
      else
         correctforwetz_local = .true.
      endif
            
      select case (par%tideloc)
         case (0,1)
            ! No spatially-varying boundary conditions
            if (xmpi_istop .or. globalonly_local) s%zs0(1:2,:) = s%zs01  ! carry out also if running in global grid
            if (xmpi_isbot .or. globalonly_local) s%zs0(s%nx:s%nx+1,:) = s%zs01
            if (xmpi_isleft .or. globalonly_local) s%zs0(:,1) = s%zs01
            if (par%ny>0) then ! stop unnecessary additional work
               if (xmpi_isright .or. globalonly_local) s%zs0(:,s%ny+1) = s%zs01
            endif
         case (2)
            select case (par%paulrevere)
               case(PAULREVERE_LAND)
                  ! One time series on bay and one on sea, easy for onshore and offshore boundaries
                  if (xmpi_istop .or. globalonly_local) s%zs0(1:2,:) = s%zs01
                  if (xmpi_isbot .or. globalonly_local) s%zs0(s%nx:s%nx+1,:) = s%zs03
                  call boundaryinterp_tide_complex(s,globalonly_local,3) ! left
                  if (par%ny>0) then
                     call boundaryinterp_tide_complex(s,globalonly_local,4) ! right
                  endif
               case(PAULREVERE_SEA)
                  ! One time series for each of the two sea points (copied to land points)
                  if (xmpi_isleft .or. globalonly_local) s%zs0(:,1) = s%zs01
                  if (par%ny>0) then
                     if (xmpi_isright .or. globalonly_local) s%zs0(:,s%ny+1) = s%zs02
                  endif
                  call boundaryinterp_tide_complex(s,globalonly_local,1) ! front
                  call boundaryinterp_tide_complex(s,globalonly_local,2) ! back
                  if (xmpi_istop .or. globalonly_local) s%zs0(2,:) = s%zs0(1,:)
                  if (xmpi_isbot .or. globalonly_local) s%zs0(s%nx,:) = s%zs0(s%nx+1,:)
            end select
         case(4)
            ! All four corners can differ
            call boundaryinterp_tide_complex(s,globalonly_local,1) ! front
            call boundaryinterp_tide_complex(s,globalonly_local,2) ! back
            call boundaryinterp_tide_complex(s,globalonly_local,3) ! left
            if (par%ny>0) then
               call boundaryinterp_tide_complex(s,globalonly_local,4) ! right
            endif
            if (xmpi_istop .or. globalonly_local) s%zs0(2,:) = s%zs0(1,:)
            if (xmpi_isbot .or. globalonly_local) s%zs0(s%nx,:) = s%zs0(s%nx+1,:)
         end select
         
         ! now set zs0 to zb level in all cells where the boundary neighbour is dry
         !
         ! wetz not initialised until after tide initialised, so only apply this step after the first time step
         if (correctforwetz_local) then
            do j=1,s%ny+1
               ! front
               if(s%wetz(2,j)==0) then 
                  s%zs0(1,j) = min(s%zb(1,j),s%zb(2,j))
               endif
               ! back
               if(s%wetz(s%nx,j)==0) then
                  s%zs0(s%nx+1,j) = min(s%zb(s%nx,j),s%zb(s%nx+1,j))
               endif
            enddo
            do i = 1,s%nx+1
               if (par%ny>0) then
                  if (s%wetz(i,2)==0) then
                     s%zs0(i,1) = min(s%zb(i,1),s%zb(i,2))
                  endif
                  if (s%wetz(i,s%ny)==0) then
                     s%zs0(i,s%ny+1) = min(s%zb(i,s%ny+1),s%zb(i,s%ny))
                  endif
               else
                  where(s%wetz==0)
                     s%zs0 = s%zb
                  endwhere
               endif
            enddo
         endif ! correctforwetz_local
      
end subroutine boundaryinterp_tide

subroutine boundaryinterp_tide_complex(s,globalonly_local,boundaryID)
      use spaceparams
      use interp
      use xmpi_module

      implicit none
      save

      type(spacepars),target                      :: s
      logical,intent(in)                          :: globalonly_local
      integer,intent(in)                          :: boundaryID
      logical                                     :: mpidiscretisation
      integer                                     :: i,j
      real*8                                      :: maxzbg,maxzbl
      real*8                                      :: distmaxzbg,distmaxzbl
      integer                                     :: relevantxmpisize
      integer                                     :: imin,imax,jmin,jmax
      logical                                     :: relevantxmpidomain
      real*8                                      :: tide1,tide2,disttide1,disttide2
      real*8                                      :: tidegradient1,tidegradient2
      real*8                                      :: fac,ddist
      
      
      ! Based on boundaryID (1 = front, 2 = back, 3 = left, 4 = right), we need to set some common switches
#ifdef USEMPI
      select case (boundaryID)
         case (1,2)
            relevantxmpisize = xmpi_n
         case (3,4)
            relevantxmpisize = xmpi_m
      end select
#else
      relevantxmpisize = 0
#endif
      select case (boundaryID)
         case (1) ! front
            imin = 1
            imax = 1
            jmin = 1
            jmax = s%ny+1
            relevantxmpidomain = xmpi_istop
            tide1 = s%zs01
            tide2 = s%zs02
            disttide1 = ndistcorners(1)
            disttide2 = ndistcorners(2)
            ! we don't use tidal gradient on front and back boundaries
            tidegradient1 = 0.d0 ! (s%zs04-s%zs01)/(sdistcorners(4)-sdistcorners(1))
            tidegradient2 = 0.d0 !(s%zs03-s%zs02)/(sdistcorners(3)-sdistcorners(2))
         case (2) ! back
            imin = s%nx+1
            imax = s%nx+1
            jmin = 1
            jmax = s%ny+1
            relevantxmpidomain = xmpi_isbot
            tide1 = s%zs04
            tide2 = s%zs03
            disttide1 = ndistcorners(4)
            disttide2 = ndistcorners(3)
            tidegradient1 = 0.d0 ! (s%zs04-s%zs01)/(sdistcorners(4)-sdistcorners(1))
            tidegradient2 = 0.d0 ! (s%zs03-s%zs02)/(sdistcorners(3)-sdistcorners(2))
         case (3) ! left (flow = right mpi)
            imin = 1
            imax = s%nx+1
            jmin = s%ny+1
            jmax = s%ny+1
            relevantxmpidomain = xmpi_isright
            tide1 = s%zs02
            tide2 = s%zs03
            disttide1 = sdistcorners(2)
            disttide2 = sdistcorners(3)
            tidegradient1 = (s%zs02-s%zs01)/(ndistcorners(2)-ndistcorners(1))
            tidegradient2 = (s%zs03-s%zs04)/(ndistcorners(3)-ndistcorners(4))
         case (4) ! right (flow = left mpi)
            imin = 1
            imax = s%nx+1
            jmin = 1
            jmax = 1
            relevantxmpidomain = xmpi_isleft
            tide1 = s%zs01
            tide2 = s%zs04
            disttide1 = sdistcorners(1)
            disttide2 = sdistcorners(4)
            tidegradient1 = (s%zs02-s%zs01)/(ndistcorners(2)-ndistcorners(1))
            tidegradient2 = (s%zs03-s%zs04)/(ndistcorners(3)-ndistcorners(4))
      end select
      !        
      ! We decide if there is MPI subdomain discretization along the boundary. This is not the case if run
      ! without MPI, or if MPI is passing global field (optional input "globalonly"), or if the MPI division is 
      ! such that there is discretisation along the boundary (i.e., xmpi_m = 1, or xmpi_n = 1)
      !
#ifdef USEMPI                     
      if (globalonly_local .or. relevantxmpisize==1) then ! only runnnig on global s, or domain has no MPI discretization
                                                          ! in cross-shore direction                                                               
         mpidiscretisation = .false.
      else
         mpidiscretisation = .true.
      endif
#else
      mpidiscretisation = .false.
#endif
      !
      ! find maximum bed level on boundary of each (or single) domain
      maxzbl = maxval(s%zb(imin:imax,jmin:jmax))
      ! find sdist at which maximum bed level occurs
      select case (boundaryID)
         case(1,2)
            distmaxzbl = s%sdist(imin,maxval(maxloc(s%zb(imin,:))))
         case(3,4)
            distmaxzbl = s%sdist(maxval(maxloc(s%zb(:,jmin))),jmin)
      end select
      !                
      ! if not running in single domain, then we need to communicate maximum bed levels and location thereof
      if (mpidiscretisation) then
#ifdef USEMPI 
         ! remove data from domains not on the boundary
         if(.not. relevantxmpidomain) maxzbl = -huge(0.d0)  
         ! store results in one array to be communicated
         maxzbg = maxzbl
         ! find the maximum of local maxima (in precompile statement to allow to compile)            
         call xmpi_allreduce(maxzbg,MPI_MAX)
         ! now find the corresponding distance value if the maximum value is in your domain
         if(abs(maxzbg-maxzbl)<1.d-8) then
            distmaxzbg = distmaxzbl
         else
            distmaxzbg = - huge(0.d0)
         endif
         ! and communicate location of global maximum       
         call xmpi_allreduce(distmaxzbg,MPI_MAX)
#endif                     
      else
         maxzbg = maxzbl
         distmaxzbg = distmaxzbl
      endif
                  
      ! Now 4 cases for imposing water levels on boundary:
      ! 1) Water level offshore (OWL) and bay (BWL) both lower than maximum bed level (zbmax), then OWB from offshore
      !    to smaxzbg, and BWL from smaxzbg to bay
      ! 2) OWL and BWL greater than maxzb, then linear interpolation along boundary
      ! 3) OWL greater than maxzb, but BWL not, then OWL from offshore to smaxzbg and interpolation from smaxzbg to
      !    bay
      ! 4) BWL greater than maxzb, but OWL not, then BWL from bay to smaxzbg and interpolation from smaxzbg to 
      !    offshore
      !
      
      ! Only carry this out on the relevant mpi subdomains or if running everything on the global grid (xmaster is not 
      ! necessarily also top, bottom, left and right) 
      if (relevantxmpidomain .or. globalonly_local) then    
         ! (1)
         if (tide1<=maxzbg .and. tide2<=maxzbg) then
            select case (boundaryID)
               case(1,2) ! top, bottom
                  where(s%sdist(imin,:)<=distmaxzbg)
                     s%zs0(imin,:) = tide1
                  elsewhere
                     s%zs0(imin,:) = tide2
                  endwhere
               case(3,4) ! left, right
                  where(s%ndist(:,jmin)<=distmaxzbg)
                     s%zs0(:,jmin) = tide1
                     s%dzs0dn(:,jmin) = tidegradient1*s%wetz(:,jmin)
                  elsewhere
                     s%zs0(:,jmin) = tide2
                     s%dzs0dn(:,jmin) = tidegradient2*s%wetz(:,jmin)
                  endwhere
            end select
         ! (2)      
         elseif (tide1>maxzbg .and. tide2>maxzbg) then
            select case (boundaryID)
               case(1,2) ! top, bottom
                  ! don't want to call this using linear_interp function, because that requires having the whole array of ndist
                  ddist = 1.d0/(disttide2-disttide1)
                  do j = jmin,jmax
                     fac = (s%ndist(imin,j)-disttide1)*ddist
                     s%zs0(imin,j) = (1.d0-fac)*tide1+fac*tide2
                  enddo
               case(3,4) ! left, right
                  ddist = 1.d0/(disttide2-disttide1)
                  do i = imin,imax
                     fac = (s%sdist(i,jmin)-disttide1)*ddist
                     s%zs0(i,jmin) = (1.d0-fac)*tide1+fac*tide2
                     s%dzs0dn(i,jmin) = ((1.d0-fac)*tidegradient1+fac*tidegradient2)*s%wetz(i,jmin)
                  enddo
            end select
         ! (3)
         elseif (tide1>maxzbg .and. tide2<=maxzbg) then
            select case (boundaryID)
               case(1,2) ! top, bottom
                  ! don't want to call this using linear_interp function, because that requires having the whole array of ndist
                  ddist = 1.d0/(disttide2-distmaxzbg)
                  do j = jmin,jmax
                     if (s%ndist(imin,j)<=distmaxzbg) then
                        s%zs0(imin,j) = tide1
                     else
                        fac = (s%ndist(imin,j)-distmaxzbg)*ddist
                        s%zs0(imin,j) = (1.d0-fac)*tide1+fac*tide2
                     endif
                  enddo
               case (3,4) ! left,right
                  ddist = 1.d0/(disttide2-distmaxzbg)
                  do i = imin,imax
                     if (s%sdist(imin,jmin)<=distmaxzbg) then
                        s%zs0(i,jmin) = tide1
                        s%dzs0dn(i,jmin) = tidegradient1*s%wetz(i,jmin)
                     else
                        fac = (s%sdist(i,jmin)-distmaxzbg)*ddist
                        s%zs0(i,jmin) = (1.d0-fac)*tide1+fac*tide2
                        s%dzs0dn(i,jmin) = ((1.d0-fac)*tidegradient1+fac*tidegradient2)*s%wetz(i,jmin)
                     endif
                  enddo
            end select
         ! (4)
         elseif(tide1<=maxzbg .and. tide2>maxzbg) then
            select case (boundaryID)
               case(1,2) ! top, bottom
                  ! don't want to call this using linear_interp function, because that requires having the whole array of ndist
                  ddist = 1.d0/(distmaxzbg-disttide1)
                  do j = jmin,jmax
                     if (s%ndist(imin,j)<=distmaxzbg) then
                        fac = (s%ndist(imin,j)-disttide1)*ddist
                        s%zs0(imin,j) = (1.d0-fac)*tide1+fac*tide2
                     else
                        s%zs0(imin,j) = tide2
                     endif
                  enddo
               case (3,4) ! left,right
                  ddist = 1.d0/(distmaxzbg-disttide1)
                  do i = imin,imax
                     if (s%sdist(imin,jmin)<=distmaxzbg) then
                        fac = (s%sdist(i,jmin)-disttide1)*ddist
                        s%zs0(i,jmin) = (1.d0-fac)*tide1+fac*tide2
                        s%dzs0dn(i,jmin) = ((1.d0-fac)*tidegradient1+fac*tidegradient2)*s%wetz(i,jmin)
                     else
                        s%zs0(i,jmin) = tide2
                        s%dzs0dn(i,jmin) = tidegradient1*s%wetz(i,jmin)
                     endif
                  enddo
            end select 
         endif ! selection of tide combinations
         !
         ! make sure that zs0 is not below zb
         select case (boundaryID)
            case(1,2) ! top, bottom
               s%zs0(imin,:) = max(s%zs0(imin,:),s%zb(imin,:))
            case(3,4) ! left, right
               s%zs0(:,jmin) = max(s%zs0(:,jmin),s%zb(:,jmin))
         end select
      endif ! relevant mpi subdomain
   
end subroutine boundaryinterp_tide_complex

subroutine fill_tide_grid(s)
      ! This subroutine fills an initial zs0 field, done by xmaster. Not called during the remainder of the simulation.
      use spaceparams
      use interp

      implicit none
      save

      type(spacepars),target                      :: s
      integer                                     :: i,j
      real*8                                      :: maxzbg
      real*8                                      :: distmaxzbg
      real*8                                      :: tide1,tide2
      real*8                                      :: fac,ddist
      
      ! Run down each cross-shore grid line. Note that this will automatically be skipped if ny=0. This is fine, zs0 at ny=1
      ! is already set by boundary condition call
      do j=2,s%ny
         ! find maximum bed level in each grid line
         maxzbg = maxval(s%zb(:,j))
         ! find sdist at which maximum bed level occurs
         distmaxzbg = s%sdist(maxloc(s%zb(:,j),dim=1),j)
         ! find offshore (tide1) and bay (tide2) tide level
         tide1 = s%zs0(1,j)
         tide2 = s%zs0(s%nx+1,j)
                           
         ! Now 4 cases for imposing water levels on boundary:
         ! 1) Water level offshore (OWL) and bay (BWL) both lower than maximum bed level (zbmax), then OWB from offshore
         !    to smaxzbg, and BWL from smaxzbg to bay
         ! 2) OWL and BWL greater than maxzb, then linear interpolation along boundary
         ! 3) OWL greater than maxzb, but BWL not, then OWL from offshore to smaxzbg and interpolation from smaxzbg to
         !    bay
         ! 4) BWL greater than maxzb, but OWL not, then BWL from bay to smaxzbg and interpolation from smaxzbg to 
         !    offshore
         !
         !
         ! Note that we only fill in i = 3:nx-1, for some reason i=2 is currently always copied from i=1 and i=nx from i=nx+1.
         ! Best not disturb that. 
         ! (1)
         if (tide1<=maxzbg .and. tide2<=maxzbg) then
            where(s%sdist(3:s%nx-1,j)<=distmaxzbg)
               s%zs0(3:s%nx-1,j) = tide1
            elsewhere
               s%zs0(3:s%nx-1,j) = tide2
            endwhere
         ! (2)      
         elseif (tide1>maxzbg .and. tide2>maxzbg) then
            ddist = 1.d0/(s%sdist(s%nx+1,j)-s%sdist(1,j))
            do i = 3,s%nx-1
               fac = (s%sdist(i,j)-s%sdist(1,j))*ddist
               s%zs0(i,j) = (1.d0-fac)*tide1+fac*tide2
            enddo
         ! (3)
         elseif (tide1>maxzbg .and. tide2<=maxzbg) then          
            ddist = 1.d0/(s%sdist(s%nx+1,j)-distmaxzbg)
            do i = 3,s%nx-1
               if (s%sdist(i,j)<=distmaxzbg) then
                  s%zs0(i,j) = tide1
               else
                  fac = (s%sdist(i,j)-distmaxzbg)*ddist
                  s%zs0(i,j) = (1.d0-fac)*tide1+fac*tide2
               endif
            enddo
         ! (4)
         elseif(tide1<=maxzbg .and. tide2>maxzbg) then
            ddist = 1.d0/(distmaxzbg-s%sdist(1,j))
            do i = 3,s%nx-1
               if (s%sdist(i,j)<=distmaxzbg) then
                  fac = (s%sdist(i,j)-s%sdist(1,j))*ddist
                  s%zs0(i,j) = (1.d0-fac)*tide1+fac*tide2
               else
                  s%zs0(i,j) = tide2
               endif
            enddo         
         endif ! selection of tide combinations
         s%zs0(3:s%nx-1,j) = max(s%zs0(3:s%nx-1,j),s%zb(3:s%nx-1,j))
      enddo ! j=2,ny
      
end subroutine fill_tide_grid

end module
