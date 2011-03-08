module drifter_module
contains
subroutine drifter(s,par)
    use params
    use spaceparams
    use readkey_module
    use xmpi_module
    use mnemmodule

    IMPLICIT NONE

    type(spacepars), target     :: s
    type(parameters)            :: par

    integer                     :: i
    integer                     :: reclen,wordsize
    integer, save               :: it_drifter=0
    integer, save               :: ndrifter

    real*8 , dimension(:)  ,allocatable,save  :: xdrift,ydrift,releasetime,retrievaltime
    logical, save                             :: first_drifter=.true.
    character*14                              :: fname
    character(256)                            :: drifterfile
    integer, dimension(2)                     :: ijz
    integer                                   :: iz,jz,id,jd
    real*8                                    :: alfad,alfaq1,alfaq2,alfaq3,alfaq4
    real*8                                    :: dxu1,dyu1,dxu2,dyu2,dxv1,dyv1,dxv2,dyv2,dxt,dyt
    real*8                                    :: ux1,uy1,ux2,uy2,vx1,vy1,vx2,vy2
    real*8                                    :: fxd,fyd

    include 's.ind'
    include 's.inp'

if (first_drifter) then                                             ! bas: should this not be in initialize ??
   first_drifter=.false.
   ! check how many drifters are present
   ! if (xmaster) then                               !!! Do not broadcast in xmaster mode!  Robert
      ndrifter   = readkey_int     ('params.txt','ndrifter',    0,         0,        50)
   ! endif
   if (ndrifter>0) then
      allocate(xdrift(ndrifter))
      allocate(ydrift(ndrifter))
      allocate(releasetime(ndrifter))
      allocate(retrievaltime(ndrifter))
      if (xmaster) then
         drifterfile = readkey_name('params.txt','drifterfile',bcast=.false.)
         open(10,file=drifterfile)
         do i=1,ndrifter
            read(10,*)xdrift(i),ydrift(i),releasetime(i),retrievaltime(i)
         enddo
         close(10)
      endif
#ifdef USEMPI
      call xmpi_bcast(xdrift)
      call xmpi_bcast(ydrift)
      call xmpi_bcast(releasetime)
      call xmpi_bcast(retrievaltime)
#endif
      if (xmaster) then
         inquire(iolength=wordsize) 1.d0
         reclen=wordsize*3
         do i=1,ndrifter
            write(fname(7:10),'(i4)')i+1000
            fname(1:7)='drifter'
            fname(11:14)='.dat'
            open(700+i,file=fname,form='unformatted',access='direct',recl=reclen)
         enddo
      endif
   endif
else
   if (ndrifter>0) then
      do i=1,ndrifter
         if (par%t>releasetime(i).and.par%t<retrievaltime(i)) then
         
            ! step1: snap to z point
            ijz = minloc(sqrt((xz-xdrift(i))**2+(yz-ydrift(i))**2))
            iz = ijz(1)
            jz = ijz(2)
            
            ! only continue if position is within domain                                                        bas: drifters are discarded once they reach the border of the domain, even in case of walls
            if (   (iz>1 .and. iz<nx .and. jz>1 .and. jz<ny)                                 .or. &
                 & (iz==1 .and. xdrift(i)>xz(iz,jz)) .or. (iz==nx .and. xdrift(i)<xz(iz,jz)) .or. &
                 & (jz==1 .and. ydrift(i)>yz(iz,jz)) .or. (jz==ny .and. ydrift(i)<yz(iz,jz))        ) then 
               
               alfad = 2*par%px - atan2(ydrift(i)-yz(iz,jz),xdrift(i)-xz(iz,jz))
               if (alfad < 0)        alfad = alfad + 2*par%px
               if (alfad > 2*par%px) alfad = alfad - 2*par%px
               
               ! step2: determine quadrant around z point where the drifter is in
               alfaq1 = alfau(iz,jz) + alfad
               alfaq2 = alfav(iz,jz) + alfad
               
               if (   alfaq1 > alfaq2 .or. sign(1.d0, alfaq1) /= sign(1.d0, alfaq2)) then
                  id = iz
                  jd = jz
               endif
               if (iz > 1) then
                  alfaq3 = alfau(iz-1,jz) + alfad - par%px
                  if (alfaq2 > alfaq3 .or. sign(1.d0, alfaq2) /= sign(1.d0, alfaq3)) then
                     id = iz-1
                     jd = jz
                  endif
               endif
               if (jz > 1) then
                  alfaq4 = alfav(iz,jz-1) + alfad - par%px
                  if (alfaq4 > alfaq1 .or. sign(1.d0, alfaq4) /= sign(1.d0, alfaq1)) then
                     id = iz
                     jd = jz-1
                  endif
               endif
               if (iz > 1 .and. jz > 1) then
                  if (alfaq3 > alfaq4 .or. sign(1.d0, alfaq3) /= sign(1.d0, alfaq4)) then
                     id = iz-1
                     jd = jz-1 
                  endif
               endif
               
               ! step3: compute distances of drifter to surrounding u and v points
               dxu1 = max(abs(xu(id,jd)   - xdrift(i)), tiny(0.d0))
               dyu1 = max(abs(yu(id,jd)   - ydrift(i)), tiny(0.d0))
               dxu2 = max(abs(xu(id+1,jd) - xdrift(i)), tiny(0.d0))
               dyu2 = max(abs(yu(id+1,jd) - ydrift(i)), tiny(0.d0))
               dxv1 = max(abs(xv(id,jd)   - xdrift(i)), tiny(0.d0))
               dyv1 = max(abs(yv(id,jd)   - ydrift(i)), tiny(0.d0))
               dxv2 = max(abs(xv(id,jd+1) - xdrift(i)), tiny(0.d0))
               dyv2 = max(abs(yv(id,jd+1) - ydrift(i)), tiny(0.d0))
               
               dxt  = 1/dxu1+1/dxu2+1/dxv1+1/dxv2
               dyt  = 1/dyu1+1/dyu2+1/dyv1+1/dyv2
               
               ! step4: convert flow velocity vectors from s and n to x and y coordinates
               ux1  = uu(id,jd)   * cos(alfau(id,jd))   - vu(id,jd)   * sin(alfau(id,jd))
               uy1  = uu(id,jd)   * sin(alfau(id,jd))   + vu(id,jd)   * cos(alfau(id,jd))
               ux2  = uu(id+1,jd) * cos(alfau(id+1,jd)) - vu(id+1,jd) * sin(alfau(id+1,jd))
               uy2  = uu(id+1,jd) * sin(alfau(id+1,jd)) + vu(id+1,jd) * cos(alfau(id+1,jd))
               vx1  = uv(id,jd)   * cos(alfav(id,jd))   - vv(id,jd)   * sin(alfav(id,jd))
               vy1  = uv(id,jd)   * sin(alfav(id,jd))   + vv(id,jd)   * cos(alfav(id,jd))
               vx2  = uv(id,jd+1) * cos(alfav(id,jd+1)) - vv(id,jd+1) * sin(alfav(id,jd+1))
               vy2  = uv(id,jd+1) * sin(alfav(id,jd+1)) + vv(id,jd+1) * cos(alfav(id,jd+1))
               
               ! step5: compute weighed average of flow velocities at drifter location in x and y direction
               fxd  = (ux1/dxu1 + ux2/dxu2 + vx1/dxv1 + vx2/dxv2)/dxt
               fyd  = (uy1/dyu1 + uy2/dyu2 + vy1/dyv1 + vy2/dyv2)/dyt
               
               ! step6: update drifter location based on flow velocities and timestep
               xdrift(i) = xdrift(i) + fxd*par%dt
               ydrift(i) = ydrift(i) + fyd*par%dt
#ifdef USEMPI
            else
            
               ! in case of MPI set drifter coordinates to huge if outside domain
               ! this allows the right coordinates to be communicated by allreduce
               
               xdrift(i)=huge(0.0d0)
               ydrift(i)=huge(0.0d0)  
#endif
            endif
            
#ifdef USEMPI
            call xmpi_allreduce(xdrift(i),MPI_MIN)
            call xmpi_allreduce(ydrift(i),MPI_MIN)
#endif

         endif
      enddo
   endif

   ! write drifter location to file   
   if (xmaster) then
      if (abs(mod(par%t,par%tintp))<1.d-6) then
         it_drifter=it_drifter+1                                        ! bas: should this not be in varoutput and it_drifter == itp ??
         do i=1,ndrifter
         
            ! set dummy value if release time has not passed yet
            if (par%t>releasetime(i).and.par%t<retrievaltime(i)) then
                write(700+i,rec=it_drifter) xdrift(i),ydrift(i),par%t
            else
                write(700+i,rec=it_drifter) -999d0,-999d0,par%t
            endif
         enddo
      endif 
   endif
   
endif


end subroutine drifter

end module drifter_module
