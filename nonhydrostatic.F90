!*************************************************************************************************
!                                  dynamic pressure correction module
!*************************************************************************************************

module nonhydrostatic

use params
use spaceparams
use logging_module

implicit none

!Declare everything private, only public when specifially stated
private ::
    !extra physical variables
	real*8,allocatable,dimension(:,:)       :: p0   		!previous dynamic pressure
	real*8,allocatable,dimension(:,:)       :: p1	    	!current dynamic pressure 
	real*8,allocatable,dimension(:,:)       :: wb0		    !previous vert. vel. bottom
	real*8,allocatable,dimension(:,:)       :: wb1  		!current vert. vel. bottom
	real*8,allocatable,dimension(:,:)       :: ws0	    	!previous vert. vel. surface
	real*8,allocatable,dimension(:,:)       :: ws1		    !current vert. vel. surface
	
	!Advection
    real*8,allocatable,dimension(:,:)       :: udwsdx        !Advection of ws in x-dir
    real*8,allocatable,dimension(:,:)       :: vdwsdy        !Advection of ws in y-dir
    real*8,allocatable,dimension(:,:)       :: udwbdx        !Advection of wb in x-dir
    real*8,allocatable,dimension(:,:)       :: vdwbdy        !Advection of wb in y-dir
    
    !Cont. eq. coeficients for velocities
	real*8,allocatable,dimension(:,:)       :: huplus       !Coeficient for u(i,j) in cont.eq.
	real*8,allocatable,dimension(:,:)       :: humin        !Coeficient for u(i-1,j) in cont.eq.
	real*8,allocatable,dimension(:,:)       :: hvplus       !Coeficient for v(i,j) in cont.eq.
	real*8,allocatable,dimension(:,:)       :: hvmin        !Coeficient for u(i,j-1) in cont.eq.

    !Grid variables
	real*8,allocatable,dimension(:,:)       :: dx           !distance between consequtive u-points
	real*8,allocatable,dimension(:,:)       :: dy           !distance between consequtive v-points
	real*8,allocatable,dimension(:,:)       :: dxm          !distance between consequtive s-points in i-direction
	real*8,allocatable,dimension(:,:)       :: dym          !distance between consequtive s-points in j-direction
    real*8,allocatable,dimension(:,:)       :: zbu          !Bottom level in u-point
	real*8,allocatable,dimension(:,:)       :: zbv          !Bottom level in v-point

    
    !Pressure coefficient matrices/vectors
	real*8,allocatable,dimension(:,:,:)     :: au           !Pressure coeficients for u(i,j)
	real*8,allocatable,dimension(:,:)       :: ur           !rhs of au(i,j)
	real*8,allocatable,dimension(:,:,:)     :: av           !Pressure coeficients for v(i,j)
	real*8,allocatable,dimension(:,:)       :: vr           !rhs of av(i,j)
	real*8,allocatable,dimension(:,:)       :: aw           !Pressure coeficients for w(i,j)
	real*8,allocatable,dimension(:,:,:,:)   :: a            !Pressure coeficients for the cont. eq.
	real*8,allocatable,dimension(:,:)       :: ar           !Rhs of the cont. eq.

    integer                                 :: nt = 0       !Counter for total number of steps
    
    !Declare subroutine(s) public
    public                                  :: nhcorrection

contains
	subroutine nh_init(s,par)
        !Initialisation routine for the nonhydrostatic subroutine
		type(spacepars),target                  :: s
		type(parameters)                        :: par
		integer i,j
		integer									:: iw,jw,ie,je
		
		call writelog('ls','', 'Setting up dynamic pressure correction')

		open(124,file='ws.dat',form='binary')       !Output surface velocity    
		open(128,file='p1.dat',form='binary')       !Output pressure
		open(129,file='wb.dat',form='binary')		!Output bottom velocity
		
		!Declare arrays, all arrays include the boundaries but they are never used
		allocate(p1(s%nx+1,s%ny+1))
		allocate(p0(s%nx+1,s%ny+1))
		allocate(ws0(s%nx+1,s%ny+1))
		allocate(ws1(s%nx+1,s%ny+1))
		allocate(wb0(s%nx+1,s%ny+1))
		allocate(wb1(s%nx+1,s%ny+1))
		allocate(huplus(s%nx+1,s%ny+1))
		allocate(hvplus(s%nx+1,s%ny+1))
		allocate(humin(s%nx+1,s%ny+1))
		allocate(hvmin(s%nx+1,s%ny+1))
		allocate(dx(s%nx+1,s%ny+1))
		allocate(dxm(s%nx+1,s%ny+1))
		allocate(dy(s%nx+1,s%ny+1))
		allocate(dym(s%nx+1,s%ny+1))

		allocate(au(s%nx+1,s%ny+1,0:1))
		allocate(av(s%nx+1,s%ny+1,0:1))
		allocate(aw(s%nx+1,s%ny+1))
		allocate(a(s%nx+1,s%ny+1,-1:1,-1:1))
		allocate(ar(s%nx+1,s%ny+1))
		allocate(ur(s%nx+1,s%ny+1))
		allocate(vr(s%nx+1,s%ny+1))
		allocate(zbu(s%nx+1,s%ny+1))
		allocate(zbv(s%nx+1,s%ny+1))
		allocate(udwsdx(s%nx+1,s%ny+1))
		allocate(vdwsdy(s%nx+1,s%ny+1))
		allocate(udwbdx(s%nx+1,s%ny+1))
		allocate(vdwbdy(s%nx+1,s%ny+1))
			
		wb0 = 0.0d0
		wb1 = 0.0d0
		ws0 = 0.0d0
		ws1 = 0.0d0
		p1	= 0.0d0
		p0  = 0.0d0
		huplus = 0.0d0
		hvplus = 0.0d0
		humin = 0.0d0
		hvmin = 0.0d0
        a   = 0.0d0
        ar  = 0.0d0
        au  = 0.0
        av  = 0.0
        ur  = 0.0
        vr  = 0.0
        udwsdx = 0.0
        vdwsdy = 0.0
        udwbdx = 0.0
        vdwbdy = 0.0
        zbu   = 0.0
        zbv   = 0.0
        
        !Setup array for dx for computational efficiency
        do j=1,s%ny+1        
		    do i=1,s%nx
		        dx(i,j)  = s%xz(i+1)-s%xz(i)
		        dxm(i,j) = s%x(i+1,j) - s%x(i,j)
		    enddo
		enddo
		dx(s%nx+1,:)  = dx(s%nx,:)    
		dxm(s%nx+1,:) = dxm(s%nx,:)
		
		!Setup array for dy for computational efficiency
		do j=1,s%ny
		    do i=1,s%nx+1
                dy(i,j)  = s%yz(j+1)-s%yz(j)
				dym(i,j) = s%y(i,j+1) - s%y(i,j)
		    enddo
		enddo    
		dy(:,s%ny+1)  = dy(:,s%ny)    
		dym(:,s%ny+1) = dym(:,s%ny)
	end subroutine nh_init

!*************************************************************************************************************
!*******************************    pressure correction routine **********************************************
!*************************************************************************************************************

	subroutine nhcorrection(s,par)
		type(spacepars),target                  :: s                    !Spacepar type from the XBeach model
		type(parameters)                        :: par
		integer									:: i,j,iw,jw,ie,je,m
		real*8                                  :: vol
	    
		if (.not. allocated(p1) ) then
			call nh_init(s,par)
		endif
			
!Bottom gradient
        do j=1,s%ny+1
			do i=1,s%nx+1
				ie = min(s%nx,i+1)
				je = min(s%ny,j+1)
				zbu(i,j)    = max(s%zb(i,j),s%zb(ie,j))
				zbv(i,j)    = max(s%zb(i,j),s%zb(i,je))
			enddo
		enddo	
	
!PCC
        !Pressure coefficients used in first fractional step
		do i=1,s%nx
		    do j=1,s%ny
		        iw = max(2,i-1)
		        jw = max(2,j-1)
		        ie = min(s%nx,i+1)
		        je = min(s%ny,j+1)
            !Pressure coefficients for U

                vol      =  dxm(i,j)*(s%hh(i,j)+s%hh(ie,j))    
                au(i,j,0)= -par%dt*(s%zs(i,j)-s%zb(ie,j))/vol !Dependence of u(i,j) on p(i-1,j)
                au(i,j,1)=  par%dt*(s%zs(ie,j)-s%zb(i,j))/vol !Dependence of u(i,j) on p(i,j)

            !Pressure coefficients for V				
                vol      =  dym(i,j)*(s%hh(i,j)+s%hh(i,je))    
                av(i,j,0)= -par%dt*(s%zs(i,j)-s%zb(i,je))/vol !Dependence of u(i,j) on p(i-1,j)
                av(i,j,1)=  par%dt*(s%zs(i,je)-s%zb(i,j))/vol !Dependence of u(i,j) on p(i,j)

			!Pressure coefficients for W				
			    aw(i,j)	 = 2.0d0*par%dt/(s%hh(i,j))           !Dependence of w(i,j) on p(i,j)
		    enddo
		enddo    	
        
!HS 
       !Coefficients for the continuity equation used in first fractional step			
	    do j=2,s%ny+1
    		do i=2,s%nx+1
			    iw = max(2,i-1)
			    jw = max(2,j-1)
				
			    humin(i,j)   =  (- s%hh(i,j) - 2*(s%zb(i,j) - zbu(i-1,j)))/dx(i,j)
				huplus(i,j)  =  (  s%hh(i,j) - 2*(zbu(i,j)  - s%zb(i,j)))/dx(i,j)
				hvmin(i,j)   =  (- s%hh(i,j) - 2*(s%zb(i,j) - zbv(i,j-1)))/dy(i,j)
				hvplus(i,j)  =  (  s%hh(i,j) - 2*(zbv(i,j)  - s%zb(i,j)))/dy(i,j)
		    enddo
	    enddo
              
!ADVECTION OF WS,WB
if (.false.) then
		do i=2,s%nx
            do j=2,s%ny
                iw = max(2,i-1)
			    jw = max(2,j-1)
			    ie = min(s%nx,i+1)
			    je = min(s%ny,j+1)
            
            !Calculate advection terms for w
                !In i-direction    
                udwsdx(i,j) = max(0.0,ur(i-1,j))*(ws0(i,j)-ws0(iw,j))/dxm(i-1,j) &
                            + min(0.0,ur(i,j))*(ws0(ie,j)-ws0(i,j))/dxm(i,j) 
                !In j-direction           
                vdwsdy(i,j) = max(0.0,vr(i,j-1))*(ws0(i,j)-ws0(i,jw))/dym(i,j-1) &
                            + min(0.0,vr(i,j))*(ws0(i,je)-ws0(i,j))/dym(i,j)
                
                udwbdx(i,j) = max(0.0,ur(i-1,j))*(wb0(i,j)-wb0(iw,j))/dxm(i-1,j) &
                            + min(0.0,ur(i,j))*(wb0(ie,j)-wb0(i,j))/dxm(i,j) 
                !In j-direction           
                vdwbdy(i,j) = max(0.0,vr(i,j-1))*(wb0(i,j)-wb0(i,jw))/dym(i,j-1) &
                            + min(0.0,vr(i,j))*(wb0(i,je)-wb0(i,j))/dym(i,j)                                        
            enddo
        enddo    
        !Update vertical velocity with explicit advection terms
         ws0   = ws0 - par%dt*udwsdx - par%dt*vdwsdy
         wb0   = wb0 - par%dt*udwbdx - par%dt*vdwbdy
endif
        
!FLOODING AND DRYING        
	do i=1,s%nx+1
	    do j=1,s%ny+1
	        iw = max(2,i-1)
		    jw = max(2,j-1)
		    ie = min(s%nx,i+1)
		    je = min(s%ny,j+1)
	        if (s%wetu(i,j) == 0) then
	            !All dependencies on u(i,j) = 0
	            huplus(i,j)  = 0.0
	            humin(ie,j)  = 0.0
	            ur(i,j)      = 0.0
	            au(i,j,:)    = 0.0   
	        endif
	        if (s%wetv(i,j) == 0) then
	            !All dependencies on v(i,j) = 0
                hvplus(i,j)  = 0.0
	            hvmin(i,je)  = 0.0
	            vr(i,j)      = 0.0
	            av(i,j,:)    = 0.0		        
	        endif
	        if(s%wetz(i,j) == 0) then
	            !All depdencies on p(i,j) = 0
	            ws0(i,j) = 0.0
	            wb0(i,j) = 0.0
	            aw(i,j)  = 1.0 !To keep the matrix from becoming undetirmined
	        endif 
	    enddo
	enddo    
    
!Set right hand side of u(i,j)
	    ur = s%uu
		
!************************************ FIRST FRACTIONAL STEP ******************************************       
    !Pressure in boundary not calculated, so main diagonal set to 1, rest 0
        a(1,:,:,:)      = 0.0 
	    a(1,:,0,0)      = 1.0
	    ar(1,:)	        = 0.0		       		
		do j=2,s%ny
	    !Front boundary
	    	if (trim(par%front) == 'nonh_1d') then
		   	    humin(2,j)  = 0.0    !Closed boundary, only used in basin tests
		    else
                au(1,j,:)   = 0.0 !dp/dx = 0.0
		    endif
			    
		!Back boundary
		    if (trim(par%back) == 'abs_1d') then
           	    huplus(s%nx,j)  = 0.0 !Closed boundary
           	    au(s%nx,j,:)    = 0.0
		    else
                au(s%nx,j,:)    = 0.0 !dp/dx = 0.0
		    endif
		!Inner matrix
			do i=2,s%nx
				    a(i,j,-1,0)  =  s%hh(i,j)*(-humin(i,j)  * au(i-1,j,0))                                    !Left diagonal         
				    a(i,j,1,0)   =  s%hh(i,j)*(-huplus(i,j) * au(i,j,1))                                      !Main diagonal
				    a(i,j,0,0)   =  s%hh(i,j)*(aw(i,j)-humin(i,j)*au(i-1,j,1)-huplus(i,j)*au(i,j,0))          !Right diagonal
                    ar(i,j)      =  s%hh(i,j)*(-ws0(i,j)-wb0(i,j) -humin(i,j)*ur(i-1,j)-huplus(i,j)*ur(i,j) &	  !RHS
				                   -hvmin(i,j)*vr(i,j-1)-hvplus(i,j)*vr(i,j))
			end do		
		enddo
	        
	!Pressure in boundary not calculated, so main diagonal set to 1, rest 0		
		a(s%nx+1,:,:,:) = 0.0
		a(s%nx+1,:,0,0) = 1.0
		ar(s%nx+1,:)    = 0.0
	!Solve matrix eq., store result in p1
		call sweep(a,ar,p1,s%nx,s%ny,1)
	!Pressure in boundary, equal to neighbouring cells to ensure dp/dx == 0 and dp/dy == 0	
		p1(1,:) = p1(2,:)                   
		p1(s%nx+1,:) = p1(s%nx,:)
		p1(:,1) = p1(:,2)
		p1(:,s%ny+1) = p1(:,s%ny)
		
	!Correct velocity
		do i=2,s%nx
            do j=2,s%ny  
				s%uu(i,j)=s%uu(i,j)-au(i,j,1)*p1(i+1,j)-au(i,j,0)*p1(i,j)
			end do   
		enddo
		
	!Update vertical velocity from u-equations			
	    do i=2,s%nx
		    do j=2,s%ny
		        iw = max(2,i-1)
		        jw = max(2,j-1)
		        wb1(i,j) = (s%uu(i,j)*(zbu(i,j)-s%zb(i,j)) + s%uu(iw,j)*(s%zb(i,j)-zbu(iw,j))) / dx(i,j) &
		                  +(s%vv(i,j)*(zbv(i,j)-s%zb(i,j)) + s%vv(i,jw)*(s%zb(i,j)-zbv(i,jw))) / dy(i,j)
				ws1(i,j) = aw(i,j)*(p1(i,j)) + wb0(i,j) - wb1(i,j) + ws0(i,j)
			enddo
		enddo
 
		ws0 = ws1
		wb0 = wb1
		ur = s%uu
		vr = s%vv 
        p0 = p1     !set p0 to intermediate pressure for output
        
!************************************ SECOND FRACTIONAL STEP ******************************************        
		!Implicit in y direction
		    a(:,1,:,:)	  = 0.0
			a(:,1,0,0)	  = 1.0
			ar(:,1) 	  = 0.0

		    do i=2,s%nx
		        if (trim(par%left) == 'wall') then
            	    hvmin(i,2)  = 0.0 !Closed boundary            	    
			    else
                    av(i,1,:)     = 0.0 !dp/dy = 0.0
			    endif
			    if (trim(par%right) == 'wall') then
            	    hvplus(i,s%ny)  = 0.0 !Closed boundary
            	    hvplus(i,s%ny)  = 0.0
            	else
            	    av(i,s%ny,:)    = 0.0 !dp/dy = 0.0
            	endif
				!Setup interiour pressure matrix in j direction
				do j=2,s%ny
					    a(i,j,0,-1) = s%hh(i,j)*(-hvmin(i,j)*av(i,j-1,0))                                      !left diagonal         
					    a(i,j,0,1)  = s%hh(i,j)*(-hvplus(i,j)*av(i,j,1))                                       !Right diagonal
					    a(i,j,0,0)  = s%hh(i,j)*(aw(i,j) -hvmin(i,j)*av(i,j-1,1)-hvplus(i,j)*av(i,j,0))		   !Main diagonal     
						ar(i,j)     = s%hh(i,j)*(-ws0(i,j)-wb0(i,j) -humin(i,j)*ur(i-1,j)-huplus(i,j)*ur(i,j) &	
					                  -hvmin(i,j)*vr(i,j-1)-hvplus(i,j)*vr(i,j))		
				end do
		    enddo

		    a(:,s%ny+1,:,:)  =  0.0
			a(:,s%ny+1,0,0)  =  1.0
			ar(:,s%ny+1)     =  0.0
		    
			call sweep(a,ar,p1,s%nx,s%ny,2)
        !Define pressures at boundaries to ensure that dp/dx=0 or dp/dy=0 	
		p1(:,1) = p1(:,2)
		p1(:,s%ny+1) = p1(:,s%ny)
		p1(1,:) = p1(2,:)           
		p1(s%nx+1,:) = p1(s%nx,:)

		!IMPLICIT STEP V-Direction
		do i=2,s%nx
            do j=2,s%ny  
				s%vv(i,j)=s%vv(i,j)-av(i,j,1)*p1(i,j+1)-av(i,j,0)*p1(i,j)
			end do   
		enddo
		
		!Update the vertical velocities
		    do i=2,s%nx
		        do j=2,s%ny
		            iw = max(2,i-1)
		            jw = max(2,j-1)
		            wb1(i,j) = (s%uu(i,j)*(zbu(i,j)-s%zb(i,j)) + s%uu(iw,j)*(s%zb(i,j)-zbu(iw,j))) / dx(i,j) &
		                     + (s%vv(i,j)*(zbv(i,j)-s%zb(i,j)) + s%vv(i,jw)*(s%zb(i,j)-zbv(i,jw))) / dy(i,j)
					ws1(i,j) = aw(i,j)*(p1(i,j)) + wb0(i,j) - wb1(i,j) + ws0(i,j)
			    enddo
		    enddo
            
		ws0 = ws1
		wb0 = wb1
		vr = s%vv 	

	    if (modulo(par%t,par%tint)==0) then
	        if(par%t>=par%tstart) then
	            write(124)ws1
			    write(128)p1
                write(129)wb1
		    endif    
        endif
		p0 = p1
	end subroutine nhcorrection

!*********************************************************************************************************************
!****************************************	SWEEP					**************************************************
!*********************************************************************************************************************

    
    subroutine sweep(a,ar,p,m,n,dir)
    !Subroutine to sweep several m by n tri-diagonal matrices for use in an adi/directional splitting structure
	!note, boundaries aren't solved!
        !VARS
		integer                          :: m,n           !Maximum number of entries in x,y direction
        real*8,dimension(:,:,:,:)        :: a             !pressure coeficient matrix
        real*8,dimension(:,:)            :: ar            !rhs
        real*8,dimension(:,:)            :: p             !Target vector for result
        real*8,dimension(m,n)            :: factor        !Scaling factor, needs to be eliminated for speed
        integer                          :: i,j           !Loop variables in x,y direction
        integer                          :: dir           !Adi direction, 1: X-direction,2: Y-direction
        
        if (dir==1) then
        !SWEEP in x-direction, y is explicit
            !scaling
                do i=1,m+1
                    a(i,2:n,1,2)        = a(i,2:n,1,2)/a(i,2:n,2,2)
                    a(i,2:n,3,2)        = a(i,2:n,3,2)/a(i,2:n,2,2)
                    ar(i,2:n)           = ar(i,2:n)/a(i,2:n,2,2)
                    a(i,2:n,2,2)        = 1.0d0
                end do
            !Forward sweep
                do i=1,m
                    factor(i,2:n)       = a(i+1,2:n,1,2)/a(i,2:n,2,2)
                    a(i+1,2:n,2,2)      = a(i+1,2:n,2,2) - factor(i,2:n)*a(i,2:n,3,2)
                    ar(i+1,2:n)         = ar(i+1,2:n)    - factor(i,2:n)*ar(i,2:n)
                end do           
            !Backward substitution
                p(m,2:n) = ar(m,2:n)/a(m,2:n,2,2)
                do i=m-1,1,-1
		    		p(i,2:n) = (ar(i,2:n)-a(i,2:n,3,2)*p(i+1,2:n))/a(i,2:n,2,2)
                end do
        else
        !SWEEP in y-direction, x is explicit
            !scaling
                do j=1,n+1
                    a(2:m,j,2,1)        = a(2:m,j,2,1)/a(2:m,j,2,2)
                    a(2:m,j,2,3)        = a(2:m,j,2,3)/a(2:m,j,2,2)
                    ar(2:m,j)           = ar(2:m,j)/a(2:m,j,2,2)
                    a(2:m,j,2,2)        = 1.0d0
                end do
            !Forward sweep
                do j=1,n
                    factor(2:m,j)       = a(2:m,j+1,2,1)/a(2:m,j,2,2)
                    a(2:m,j+1,2,2)      = a(2:m,j+1,2,2) - factor(2:m,j)*a(2:m,j,2,3)
                    ar(2:m,j+1)         = ar(2:m,j+1)    - factor(2:m,j)*ar(2:m,j)
                end do           
            !Backward substitution
                p(2:m,n) = ar(2:m,n)/a(2:m,n,2,2)
                do j=n-1,1,-1
    				p(2:m,j) = (ar(2:m,j)-a(2:m,j,2,3)*p(2:m,j+1))/a(2:m,j,2,2)
                end do
       endif    
    end subroutine sweep
end module nonhydrostatic

