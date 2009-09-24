      program enve4
      parameter (Nmax=8192)
      real eta(Nmax),etahi(Nmax),E(Nmax),etabi(Nmax,2),etalo(Nmax),
     &     dEdt(Nmax),ans(2*Nmax),cv(-Nmax/2:Nmax/2)
      real k(Nmax)

      character*78 filnam

c
c     Read input     
c-----------------------------------------------------------------------
      call mess ('Reading input.............................')
      open(1,file='envelo.inp')
      read(1,*)filnam
      read(1,*)Nt
      read(1,*)dtf
      read(1,*)fspl
      read(1,*)flo
      read(1,*)fNyq
      read(1,*)rho
      read(1,*)g
      read(1,*)h
      read(1,*)dx 
      close(1)
       
      rhog=rho*g 
      Nsamp=nint(.5/fNyq/dtf)
      dt=Nsamp*dtf 
      Nt=Nt/Nsamp
c
c     Read eta(t)      
c-----------------------------------------------------------------------
      call mess ('Reading eta(t).............................')
      call geteta (filnam,eta,Nt,Nsamp,0)
      
      call demean (eta,etabar,Nt)

      do 500 ix=2,1,-1

      x=(ix-1)*dx  

      call mess ('Computing components......................')
      call bound(eta,etahi,E,etabi(1,ix),etalo,
     &             Nmax,Nt,dt,fspl,k,x,rho,g,h)
      WRITE(*,'(6e12.5)')
     &    (eta(i),etahi(i),E(i),etabi(i,ix),etalo(i),k(i),i=1,15)
      write(*,*)Nmax,Nt,dt,fspl,x,rho,g,h

500   continue

      N=1
25    if (N.lt.Nt) then
         N=2*N
         goto 25
      endif       

     
      taumax=3.*dx/sqrt(g*h)
      kmax=nint(taumax/dt)
      call correl(etabi(1,1),etabi(1,2),N,ans)

         ic=-N/2
         do 26 i=N/2+1,N
            cv(ic)=ans(i)/Nt
            ic=ic+1
26       continue
         do 27 i=1,N/2
            cv(ic)=ans(i)/Nt
            ic=ic+1
27       continue 

         cvmin=1.e10
         cvmax=-1.e10 
         do 501 ik=-kmax,kmax
            if (cv(ik).gt.cvmax) then
               cvmax=cv(ik)
               ik0=ik
            endif
501      continue
         ikm=ik0-1
         ikp=ik0+1
         aa= (cv(ikm)-2.*cv(ik0)+cv(ikp))/(2.*dt*dt)
         bb= (cv(ikp)-cv(ikm))/(2.*dt)
         cc= cv(ik0)
         tmax=-.5*bb/aa
         cvmax=aa*tmax*tmax+bb*tmax+cc
         tmax=ik0*dt+tmax
         Cg=-dx/tmax
         write(*,*)' tmax  = ',tmax
         write(*,*)' ik0   = ',ik0
         write(*,*)' Cvmax = ',cvmax         
         write(*,*)' Cg    = ',Cg

         f1=2.*fspl
         cg1=cgfun(f1,h)
         f2=2.2*fspl
         cg2=cgfun(f2,h)

502      f2=f1+(Cg-cg1)/(cg2-cg1)*(f2-f1)
         cg2=cgfun(f2,h)
         write(*,*)f2,cg2
         if (abs(cg2-cg).gt..001) goto 502 


      E0=E(1)       
      Ebar=avrg(E,Nt)
      Hrms=sqrt(8.*Ebar/rhog)
      call deriv (E,dt,Nt,dEdt)

      call mess ('Writing output............................')
 
      do 10 ic=1,75
         if (filnam(ic:ic).eq.'.') then
            filnam(ic+1:ic+3)='ezs'
            goto 11
         endif
10    continue
11    continue
      open(2,file=filnam)
      write(2,'(a1,a9,4a10)')'*','t (s)','eta LF(m)','E (J/m2)',
     &                       'eta BI(m)','eta F(m)'
      write(2,'(a4)')'BL01'
      write(2,*)Nt,5
      write (2,'(f10.2,f10.5,f10.3,2f10.5)')
     &  ((i-1)*dt,etalo(i),E(i),etabi(i,1),etalo(i)-etabi(i,1),
     &    i=1,Nt)
      close (2)

      do 20 ic=1,75
         if (filnam(ic:ic).eq.'.') then
            filnam(ic+1:ic+3)='env'
            goto 21
         endif
20    continue
21    continue
      open(2,file=filnam)
      write(2,'(a4)')'BL01'
      write(2,*)Nt,3
      write (2,'(f10.2,2f10.5)')
     &     ((i-1)*dt,etahi(i),sqrt(2.*max(E(i),1.e-6)/rhog),
     &       i=1,Nt)
      close(2)

      do 30 ic=1,75
         if (filnam(ic:ic).eq.'.') then
            filnam(ic+1:ic+3)='det'
            goto 31
         endif
30    continue
31    continue
      open(2,file=filnam)
      write(2,'(f10.5)')(dEdt(i),i=1,Nt)
      close(2) 


      do 40 ic=1,75
         if (filnam(ic:ic).eq.'.') then
            filnam(ic+1:ic+3)='gen'
            goto 41
         endif
40    continue
41    continue
      open(2,file=filnam)
      write(2,'(i12  ,10x,a)')Nt    ,'Nt'
      write(2,'(f12.5,10x,a)')dt    ,'dt'
      write(2,'(f12.5,10x,a)')etabar,'etabar'
      write(2,'(f12.5,10x,a)')Hrms  ,'Hrms'
      write(2,'(f12.5,10x,a)')Ebar  ,'Ebar'
      write(2,'(f12.5,10x,a)')E0    ,'E0'
      write(2,'(f12.5,10x,a)')f2    ,'frep'
      close(2)

      end


      subroutine bound(eta,etahi,E,etabi,etalo,Nmax,Nt,
     &                 dt,fspl,k,x,rho,g,h)
      real eta(Nmax),etahi(Nmax),E(Nmax),etabi(Nmax),etalo(Nmax)
      real k(Nmax)

      pi   = 4.*atan(1.)
      rhog = rho*g
      gh   = g*h
      g2   = 2.*g
      gp5  = .5*g

      N=1
1        if (N.lt.Nt) then
         N=N+N
         goto 1
      endif
      N=N/2
      if (2*N.gt.Nmax) stop ' 2*N>Nmax!'

      df=1./(2*N*dt)
      dom= 2.*pi*df
      kspl=nint(fspl/df)  
      write(*,*)' N    = ',N
      write(*,*)' df   = ',df
      write(*,*)' fspl = ',fspl
      write(*,*)' df   = ',df
      write(*,*)' kspl = ',kspl

      do 5 i=kspl,N
         t=1./(i*df)
         call wavenr (h,t,k(i))
5        continue

      do 10 i=1,2*N
         E(i)=0.
         etalo(i)=0.
         etabi(i)=0.
10       etahi(i)=eta(i)
*
*        Compute Fourier transform of etahi
*----------------------------------------------------------------------
      call realft(etahi,N,1)   

      if (abs(x).lt.1.e-3) then

      do 100 j=0,kspl-1
         write(*,*)j,kspl-1
         SA    = 0.
         SB    = 0.
         Salfa = 0.
         Sbeta = 0.
         do 110 i=kspl,N-j-1
            an   = etahi(2*  i  +1)
            bn   = etahi(2*  i  +2)
            anpj = etahi(2*(i+j)+1)
            bnpj = etahi(2*(i+j)+2)
            A    = an*anpj + bn*bnpj
            B    = an*bnpj - anpj*bn
            SA   = SA+A
            SB   = SB+B
            if (j.gt.0) then
            Cg   = j*dom/(k(i+j)-k(i))
            C    = (2*i+j)*dom/(k(i)+k(i+j))
            R    = (g2*Cg-gp5*C)/(C*Cg*Cg-gh*C)
            alfa = A*R
            beta = B*R
            Salfa= Salfa+alfa
            Sbeta= Sbeta+beta 
            endif
110      continue
         E(2*j+1) = SA/N
         E(2*j+2) = SB/N
         etabi(2*j+1) = Salfa/N
         etabi(2*j+2) = Sbeta/N 
         etalo(2*j+1) = etahi(2*j+1)
         etalo(2*j+2) = etahi(2*j+2)
         etahi(2*j+1) = 0.
         etahi(2*j+2) = 0.
100      continue
      write(*,'(6e12.5)')(E(ii),ii=1,12) 
      etabi(1) = 0.
      etabi(2) = 0.

      else

      do 200 j=0,kspl-1
         write(*,*)j,kspl-1
         SA    = 0.
         SB    = 0.
         Salfa = 0.
         Sbeta = 0.
         do 210 i=kspl,N-j-1
            an   = etahi(2*  i  +1)
            bn   = etahi(2*  i  +2)
            anpj = etahi(2*(i+j)+1)
            bnpj = etahi(2*(i+j)+2)
            dkx  = x*(k(i+j)-k(i))
            co   = cos(dkx)
            si   = sin(dkx)
            P    = an*anpj + bn*bnpj
            Q    = an*bnpj - anpj*bn
            A    = P*co - Q*si
            B    = P*si + Q*co
            SA   = SA+A
            SB   = SB+B
            if (j.gt.0) then
            Cg   = j*dom/(k(i+j)-k(i))
            C    = (2*i+j)*dom/(k(i)+k(i+j))
            R    = (g2*Cg-gp5*C)/(C*Cg*Cg-gh*C)
            alfa = A*R
            beta = B*R
            Salfa= Salfa+alfa
            Sbeta= Sbeta+beta 
            endif
210      continue
         E(2*j+1) = SA/N
         E(2*j+2) = SB/N
         etabi(2*j+1) = Salfa/N
         etabi(2*j+2) = Sbeta/N 
         etalo(2*j+1) = etahi(2*j+1)
         etalo(2*j+2) = etahi(2*j+2)
         etahi(2*j+1) = 0.
         etahi(2*j+2) = 0.
200      continue
      write(*,'(6e12.5)')(E(ii),ii=1,12) 
      etabi(1) = 0.
      etabi(2) = 0.

      endif

      call realft (E    ,N,-1)
      call realft (etahi,N,-1)
      call realft (etalo,N,-1)
      call realft (etabi,N,-1)

      do 120 i=1,Nt
         E(i)     = E(i)*rhog/N
         etahi(i) = etahi(i)/N
         etalo(i) = etalo(i)/N
         etabi(i) = etabi(i)/N
120      continue
      do 121 i=Nt+1,2*N
         E(i)     = 0.
         etahi(i) = 0.
         etalo(i) = 0.
         etabi(i) = 0.
121      continue

      return
      end

      subroutine wavenr(h,tp,k)
      real k,k0,k1

      pi=3.141592654
      g=9.81

      omega=2.*pi/tp
      k0=omega*omega/g
      if (h.gt.(3.*tp*tp)) then
         k=k0
      elseif ((k0*h/2./pi).lt.0.0001) then
         k=2.*pi/(sqrt(g*h)*tp)
      else
         k1=k0/sqrt(tanh(k0*h))
         s=sinh(k1*h)
         c=cosh(k1*h)
         a1=h/(c*c)-k1*h*h/(c*c*c)*s
         a2=s/c+k1*h/(c*c)
         a3=k1*s/c-k0
         k=k1+(-a2+sqrt(a2*a2-4.*a1*a3))/(2.*a1)
      endif

      return
      end

c
c     Demean function y over N points      
c-----------------------------------------------------------------------
      subroutine demean (y,avrg,N)
      real y(N)

      avrg=0. 
      do 10 i=1,N  
         avrg=avrg+y(i)
10    continue
      avrg=avrg/N
      do 20 i=1,N
         y(i)=y(i)-avrg
20    continue

      end
c
c     Read time series from file 
c-----------------------------------------------------------------------
      subroutine geteta (filnam,y,N,Nsamp,Nstart)
      real y(N+Nstart),dum(10)
      character*78 filnam

      open(1,file=filnam,access='sequential',form='binary') 
      read(1)(y(i),(dum(is),is=1,Nsamp-1),i=Nstart+1,Nstart+N)

      close(1)
      end

c
c     Display a message on screen
c-----------------------------------------------------------------------
      subroutine mess (messg)
      character*(*) messg

      write(*,'(1x,78a1)')('-',i=1,78)
      write(*,'(1x,a)')messg 

      end 

c
c     Compute derivative of time-function
c-----------------------------------------------------------------------
      subroutine deriv (y,dt,N,yacc)
      real y(N),yacc(N)

      dti=1./dt
      dt2i=.5*dti 

      yacc(1)=dti*(y(2)-y(1))
      do 10 i=2,N-1
         yacc(i)=dt2i*(y(i+1)-y(i-1))
10    continue
      yacc(N)=dti*(y(N)-y(N-1))

      end

c
c     Average function y over N points      
c-----------------------------------------------------------------------
      function avrg (y,N)
      real y(N)

      avrg=0. 
      do 10 i=1,N  
         avrg=avrg+y(i)
10    continue
      avrg=avrg/N

      end

      function cgfun(f,h)
      real k,kh

            pi=3.141592654
            g=9.81
            t=1./f
            omega=2.*pi*f
            call wavenr(h,t,k)  
            kh=k*h
            tkh=tanh(kh)
            cgfun=.5*g/omega*(tkh+kh*(1-tkh*tkh))
 
      end

      SUBROUTINE CORREL(DATA1,DATA2,N,ANS)
      PARAMETER(NMAX=8192)
      DIMENSION DATA1(N),DATA2(N)
      COMPLEX FFT(NMAX),ANS(N)
      CALL TWOFFT(DATA1,DATA2,FFT,ANS,N)
      NO2=FLOAT(N)/2.0
      DO 11 I=1,N/2+1
        ANS(I)=FFT(I)*CONJG(ANS(I))/NO2
11    CONTINUE
      ANS(1)=CMPLX(REAL(ANS(1)),REAL(ANS(N/2+1)))
      CALL REALFT(ANS,N/2,-1)
      RETURN
      END
      SUBROUTINE REALFT(DATA,N,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      THETA=6.28318530717959D0/2.0D0/DBLE(N)
      C1=0.5
      IF (ISIGN.EQ.1) THEN
        C2=-0.5
        CALL FOUR1(DATA,N,+1)
      ELSE
        C2=0.5
        THETA=-THETA
      ENDIF
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.0D0+WPR
      WI=WPI
      N2P3=2*N+3
      DO 11 I=2,N/2+1
        I1=2*I-1
        I2=I1+1
        I3=N2P3-I2
        I4=I3+1
        WRS=SNGL(WR)
        WIS=SNGL(WI)
        H1R=C1*(DATA(I1)+DATA(I3))
        H1I=C1*(DATA(I2)-DATA(I4))
        H2R=-C2*(DATA(I2)+DATA(I4))
        H2I=C2*(DATA(I1)-DATA(I3))
        DATA(I1)=H1R+WRS*H2R-WIS*H2I
        DATA(I2)=H1I+WRS*H2I+WIS*H2R
        DATA(I3)=H1R-WRS*H2R+WIS*H2I
        DATA(I4)=-H1I+WRS*H2I+WIS*H2R
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
11    CONTINUE
      IF (ISIGN.EQ.1) THEN
        H1R=DATA(1)
        DATA(1)=H1R+DATA(2)
        DATA(2)=H1R-DATA(2)
      ELSE
        H1R=DATA(1)
        DATA(1)=C1*(H1R+DATA(2))
        DATA(2)=C1*(H1R-DATA(2))
        CALL FOUR1(DATA,N,-1)
      ENDIF
      RETURN
      END
      SUBROUTINE FOUR1(DATA,NN,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END

      SUBROUTINE TWOFFT(DATA1,DATA2,FFT1,FFT2,N)
      DIMENSION DATA1(N),DATA2(N)
      COMPLEX FFT1(N),FFT2(N),H1,H2,C1,C2
      C1=CMPLX(0.5,0.0)
      C2=CMPLX(0.0,-0.5)
      DO 11 J=1,N
        FFT1(J)=CMPLX(DATA1(J),DATA2(J))
11    CONTINUE
      CALL FOUR1(FFT1,N,1)
      FFT2(1)=CMPLX(AIMAG(FFT1(1)),0.0)
      FFT1(1)=CMPLX(REAL(FFT1(1)),0.0)
      N2=N+2
      DO 12 J=2,N/2+1
        H1=C1*(FFT1(J)+CONJG(FFT1(N2-J)))
        H2=C2*(FFT1(J)-CONJG(FFT1(N2-J)))
        FFT1(J)=H1
        FFT1(N2-J)=CONJG(H1)
        FFT2(J)=H2
        FFT2(N2-J)=CONJG(H2)
12    CONTINUE
      RETURN
      END
