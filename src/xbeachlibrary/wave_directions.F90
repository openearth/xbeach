module wave_directions_module
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Copyright (C) 2007 UNESCO-IHE, WL|Delft Hydraulics and Delft University !
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
   implicit none
   save   
      
contains
   
   subroutine wave_directions(s,par)
      use params
      use spaceparams
      use roelvink_module
      use wave_functions_module
      use xmpi_module
      use logging_module
      use paramsconst

      ! wwvv in my testcase, this routine was not called, so it is not
      ! tested. Nevertheless, I put in code for the parallel version.

      IMPLICIT NONE

      type(spacepars), target     :: s
      type(parameters)            :: par

      integer                     :: i,imax,i1
      integer                     :: j
      integer                     :: itheta,iter
      real*8,dimension(:),allocatable,save        :: dist,factor,e01
      real*8 , dimension(:,:)  ,allocatable,save  :: dhdx,dhdy,dudx,dudy,dvdx,dvdy
      real*8 , dimension(:,:)  ,allocatable,save  :: uorb
      real*8 , dimension(:,:)  ,allocatable,save  :: sinh2kh ! ,wm
      real*8 , dimension(:,:,:),allocatable,save  :: xadvec,yadvec,thetaadvec,dd
            real*8 , dimension(:),allocatable,save      :: Hprev
      real*8                                      :: Herr,dtw
      logical                                     :: stopiterate

      if (.not. allocated(e01)) then
         allocate(e01(1:s%ntheta_s))
         allocate(dist(1:s%ntheta_s))
         allocate(factor(1:s%ntheta_s))
         allocate(xadvec    (s%nx+1,s%ny+1,s%ntheta_s))
         allocate(yadvec    (s%nx+1,s%ny+1,s%ntheta_s))
         allocate(thetaadvec(s%nx+1,s%ny+1,s%ntheta_s))
         allocate(dd        (s%nx+1,s%ny+1,s%ntheta_s))

        
         allocate(dhdx(s%nx+1,s%ny+1))
         allocate(dhdy(s%nx+1,s%ny+1))
         allocate(dudx(s%nx+1,s%ny+1))
         allocate(dudy(s%nx+1,s%ny+1))
         allocate(dvdx(s%nx+1,s%ny+1))
         allocate(dvdy(s%nx+1,s%ny+1))
         allocate(uorb(s%nx+1,s%ny+1))
         allocate(sinh2kh(s%nx+1,s%ny+1))
         allocate(Hprev(s%ny+1))

      endif

   
      xadvec      = 0.0d0
      yadvec      = 0.0d0
      thetaadvec  = 0.0d0
      dd          = 0.0d0
      dhdx        = 0.0d0
      dhdy        = 0.0d0
      dudx        = 0.0d0
      dudy        = 0.0d0
      dvdx        = 0.0d0
      dvdy        = 0.0d0
      sinh2kh     = 0.0d0

      ! Slopes of water depth
      call slope2D(s%hhws,s%nx,s%ny,s%dsu,s%dnv,dhdx,dhdy,s%wete)
      ! Dano limit slopes used in refraction to avoid unrealistic refraction speeds
      dhdx=sign(1.d0,dhdx)*min(abs(dhdx),0.1d0)
      dhdy=sign(1.d0,dhdy)*min(abs(dhdy),0.1d0)
      if (par%wci==1) then
         call slope2D(s%uws,s%nx,s%ny,s%dsu,s%dnv,dudx,dudy,s%wete)
         call slope2D(s%vws,s%nx,s%ny,s%dsu,s%dnv,dvdx,dvdy,s%wete)
      else
         dudx = 0.d0
         dudy = 0.d0
         dvdx = 0.d0
         dvdy = 0.d0
      endif
      ! wwvv these slope routines are in wave_timestep, and are
      !   MPI-aware
      !
      ! Calculate once sinh(2kh)
      
      where(s%wete==1 .and. 2*s%hhws*s%k<=3000.d0)
         sinh2kh=sinh(min(2*s%k*s%hhws,10.0d0))
      elsewhere
         sinh2kh = 3000.d0
      endwhere
            
      forall (i=1:s%nx+1,j=1:s%ny+1,s%wete(i,j)==1)
         s%thetamean(i,j)=(sum(s%ee_s(i,j,:)*s%thet_s(i,j,:))/s%ntheta_s)/(max(sum(s%ee_s(i,j,:)),0.00001d0)/s%ntheta_s)
      endforall
      
      !
      ! split wave velocities in wave grid directions theta
      call compute_wave_direction_velocities(s,par,2,dhdx,dhdy,dudx,dudy,dvdx,dvdy,sinh2kh)

      !dtw=.9*minval(xz(2:nx+1)-xz(1:nx))/maxval(cgx_s)

      s%E(1,:)=max(sum(s%ee_s(1,:,:),2)*s%dtheta_s,0.d0)
      s%H(1,:)=sqrt(s%E(1,:)/par%rhog8)

      imax=s%nx
      !Dano  This is ok, since we will set mpiboundary to y in stationary mode
      !
      ! write to screen that waves are updated
      if(xmaster) call writelog('ls','(a,f0.2,a)','Computing wave directions at t = ',par%t,' s')
      do i=2,imax
         if(par%wci==0) then
            dtw=.5*minval(s%dsu(i:i+1,jmin_ee:jmax_ee))/maxval(s%cgx_s(i-1:i+1,jmin_ee:jmax_ee,:))
            dtw=min(dtw,.5*minval(s%dnv(i,jmin_ee:jmax_ee))/maxval(abs(s%cgy_s(i,jmin_ee:jmax_ee,:))))
            dtw=min(dtw,.5*s%dtheta/max(1.0d-6,maxval(abs(s%ctheta_s(i,jmin_ee:jmax_ee,:)))))
            !Dano: need to make sure all processes use the same dtw, min of all processes
#ifdef USEMPI
            call xmpi_allreduce(dtw,MPI_MIN)
#endif
         else
            dtw = par%dt
         endif
      
         Herr=1.
         iter=0
         stopiterate=.false.
         do while (stopiterate .eqv. .false.)
            iter=iter+1
            Hprev=s%H(i,:)
            !
            ! transform to wave action
            !
            i1=max(i-2,1)
            do itheta=1,s%ntheta_s
               where(s%wete(i1:i+1,:)==1)
                  s%ee_s(i1:i+1,:,itheta) = s%ee_s(i1:i+1,:,itheta)/s%sigm(i1:i+1,:)
               endwhere
            enddo
            !
            ! Upwind Euler timestep propagation
            !
            ! Robert: can we also use WARMBEAM here?
            if  (i>2.and. (par%scheme==SCHEME_UPWIND_2 .or. par%scheme==SCHEME_WARMBEAM)) then
               call advecxho(s%ee_s(i-2:i+1,:,:),s%cgx_s(i-2:i+1,:,:),xadvec(i-2:i+1,:,:),    &
               3,s%ny,s%ntheta_s,s%dnu(i-2:i+1,:),s%dsu(i-2:i+1,:),s%dsdnzi(i-2:i+1,:),par%scheme, &
               s%wete(i-2:i+1,:),par%dt,s%dsz)
            else
               call advecxho(s%ee_s(i-1:i+1,:,:),s%cgx_s(i-1:i+1,:,:),xadvec(i-1:i+1,:,:),    &
               2,s%ny,s%ntheta_s,s%dnu(i-1:i+1,:),s%dsu(i-1:i+1,:),s%dsdnzi(i-1:i+1,:),SCHEME_UPWIND_1,&
               s%wete(i-1:i+1,:),par%dt,s%dsz)
            endif
            call advecyho(s%ee_s(i,:,:),s%cgy_s(i,:,:),yadvec(i,:,:),                                  &
            0,s%ny,s%ntheta_s,s%dsv(i,:),s%dnv(i,:),s%dsdnzi(i,:),SCHEME_UPWIND_1,s%wete(i,:),par%dt,s%dnz)
            call advecthetaho(s%ee_s(i,:,:),s%ctheta_s(i,:,:),thetaadvec(i,:,:),0,s%ny,s%ntheta_s,s%dtheta,par%scheme,s%wete(i,:))

            s%ee_s(i,:,:)=s%ee_s(i,:,:)-dtw*(xadvec(i,:,:) + yadvec(i,:,:) + thetaadvec(i,:,:))
#ifdef USEMPI
            call xmpi_shift(s%ee_s(i-1:i,:,:),SHIFT_Y_R,1,2)
            call xmpi_shift(s%ee_s(i-1:i,:,:),SHIFT_Y_L,3,4)
#endif
            !
            ! transform back to wave energy
            !
            do itheta=1,s%ntheta_s
               where(s%wete(i1:i+1,:)==1)
                  s%ee_s(i1:i+1,:,itheta) = s%ee_s(i1:i+1,:,itheta)*s%sigm(i1:i+1,:)
                  s%ee_s(i1:i+1,:,itheta)=max(s%ee_s(i1:i+1,:,itheta),0.0d0)
               elsewhere
                  s%ee_s(i1:i+1,:,itheta) = 0.d0
               endwhere
            enddo
            !
            ! Energy integrated over wave directions,Hrms
            !
            forall(j=1:s%ny+1,s%wete(i,j)==1)
               s%E(i,j)=sum(s%ee_s(i,j,:))*s%dtheta
            endforall
            where(s%wete(i,:)==1)
               s%H(i,:)=sqrt(s%E(i,:)/par%rhog8)
            endwhere
            do itheta=1,s%ntheta_s
               where(s%wete(i,:)==1)
                  s%ee_s(i,:,itheta)=s%ee_s(i,:,itheta)/max(1.0d0,(s%H(i,:)/(par%gammax*s%hhws(i,:)))**2)
               endwhere
            enddo
            where(s%wete(i,:)==1)
               s%H(i,:)=min(s%H(i,:),par%gammax*s%hhws(i,:))
               s%E(i,:)=par%rhog8*s%H(i,:)**2
            endwhere
            if (par%snells==0) then !Dano not for SNellius
               where(s%wete(i,:)==1)
                  s%thetamean(i,:) = (sum(s%ee_s(i,:,:)*s%thet_s(i,:,:),2)/s%ntheta_s)/ &
                                       (max(sum(s%ee_s(i,:,:),2),0.000010d0)/s%ntheta_s)
               endwhere
            endif
            !
            ! Total dissipation

            select case(par%break)
               case(BREAK_ROELVINK1,BREAK_ROELVINK2)
               call roelvink       (par,s,i)
               case(BREAK_BALDOCK)
               call baldock        (par,s,i)
               case(BREAK_JANSSEN)
               call janssen_battjes(par,s,i)
            end select


            ! Dissipation by bed friction
            where(s%fw(i,:)>0.d0 .and. s%wete(i,:)==1 .and. s%hhws(i,:)>par%fwcutoff)
               uorb(i,:)=par%px*s%H(i,:)/par%Trep/sinh(min(max(s%k(i,:),0.01d0)*s%hhws(i,:),10.0d0))
               s%Df(i,:)=0.6666666d0/par%px*par%rho*s%fw(i,:)*uorb(i,:)**3
            elsewhere
               s%Df(i,:) = 0.d0
            end where
            !
            ! Distribution of dissipation over directions and frequencies
            !
            do itheta=1,s%ntheta_s
               where(s%wete(i,:)==1)
                  dd(i,:,itheta)=s%ee_s(i,:,itheta)*(s%D(i,:)+s%Df(i,:))/max(s%E(i,:),0.00001d0) ! Robert: missing Dveg here!
               elsewhere
                  dd(i,:,itheta)=0.d0
               endwhere
            end do
            !
            ! Euler step dissipation
            !
            do j=jmin_ee,jmax_ee
               do itheta=1,s%ntheta_s
                  if (dtw*dd(i,j,itheta)>s%ee_s(i,j,itheta) .and. s%wete(i,j)==1) then
                     dtw=min(dtw,.5*s%ee_s(i,j,itheta)/dd(i,j,itheta))
                  endif
               enddo
            enddo
            !Dano: need to make sure all processes use the same dtw, min of all processes
#ifdef USEMPI
            call xmpi_allreduce(dtw,MPI_MIN)
#endif
            do j=1,s%ny+1
               do itheta=1,s%ntheta_s
                  if(s%wete(i,j)==1) then
                     s%ee_s(i,j,itheta)=s%ee_s(i,j,itheta)-dtw*dd(i,j,itheta)
                     s%ee_s(i,j,itheta)=max(s%ee_s(i,j,itheta),0.0d0)
                  else
                     s%ee_s(i,j,itheta)=0.0d0
                  end if
               end do
            end do
            ! Lateral boundary condition
            if (xmpi_isleft .and. s%ny>0) then
               do itheta=1,s%ntheta_s
                  if (s%sinth_s(i,1,itheta)>=0.) then
                     s%ee_s(i,1,itheta)=s%ee_s(i,2,itheta)
                  endif
               enddo
               s%k(:,1)=s%k(:,2)
               s%sigm(:,1)=s%sigm(:,2)
            endif
            if (xmpi_isright .and. s%ny>0) then
               do itheta=1,s%ntheta_s
                  if (s%sinth_s(i,s%ny+1,itheta)<=0.) then
                     s%ee_s(i,s%ny+1,itheta)=s%ee_s(i,s%ny,itheta)
                  endif
               end do
               s%k(:,s%ny+1)=s%k(:,s%ny)
               s%sigm(:,s%ny+1)=s%sigm(:,s%ny)
            endif
            !
            ! Compute mean wave direction
            !
            if (par%snells==0) then
               where(s%wete(i,:)==1)
                  s%thetamean(i,:)=(sum(s%ee_s(i,:,:)*s%thet_s(i,:,:),2)/size(s%ee_s(i,:,:),2)) &
                  /(max(sum(s%ee_s(i,:,:),2),0.000010d0) /size(s%ee_s(i,:,:),2))
               elsewhere 
                  ! forward-copy wave directions on dry cells in case they flood before the next wave update
                  s%thetamean(i,:) = s%thetamean(i-1,:)
               endwhere
            endif
            !
            ! Energy integrated over wave directions,Hrms
            !
            where(s%wete(i,:)==1)
               s%E(i,:)=sum(s%ee_s(i,:,:),2)*s%dtheta
               s%H(i,:)=sqrt(s%E(i,:)/par%rhog8)
            elsewhere
               s%E(i,:)=0.d0
               s%H(i,:)=0.d0
            endwhere
            Herr=maxval(abs(Hprev(jmin_ee:jmax_ee)-s%H(i,jmin_ee:jmax_ee)),mask=s%wete(i,:)==1)

#ifdef USEMPI
            call xmpi_allreduce(Herr,MPI_MAX)
#endif
            ! Stopping criteria
            if (iter<par%maxiter) then
               if (Herr<par%maxerror) then
                  stopiterate=.true.
                  !if(xmaster) call writelog('ls','(a,i4,a,i4)','Wave propagation row ',i,', iteration ',iter)
               endif
            else
               stopiterate=.true.
               if(xmaster) call writelog('ls','(a,i4,a,i4,a,f5.4)','Wave propagation row ',i,', iteration ',iter,', error: ',Herr)
            endif
         enddo ! End while loop
      enddo ! End do i=2:s%nx loop

      s%ee_s(s%nx+1,:,:) = s%ee_s(s%nx,:,:)
      s%E(s%nx+1,:)    = s%E(s%nx,:)
      s%H(s%nx+1,:)    = s%H(s%nx,:)
      s%k(s%nx+1,:)   = s%k(s%nx,:)
      s%sigm(s%nx+1,:) = s%sigm(s%nx,:)
      s%cg(s%nx+1,:)   = s%cg(s%nx,:)
      s%c(s%nx+1,:)    = s%c(s%nx,:)
      s%thet_s(s%nx+1,:,:) = s%thet_s(s%nx,:,:)
      
   end subroutine wave_directions
end module wave_directions_module
