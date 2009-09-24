      program genser
      parameter(Nmax=65536)
      real a(Nmax+2),S(Nmax+2)

      open(1,file='genser.inp')
      read(1,*) idum 
      read(1,*) Tmax
      read(1,*) dt
      read(1,*) fp
      read(1,*) gamma
      read(1,*) Hrms
      read(1,*) nstep

      twopi=8.*atan(1.)
      Nt=int(Tmax/dt)

      N=1
1        if (N.lt.Nt) then
         N=N+N
         goto 1
      endif
      N=N/2
      if (2*N.gt.Nmax) stop ' 2*N>Nmax!'
      df=1./(2*N*dt) 
      
      fnyq=1./2./dt
      ifnyq=fnyq/df 
      fpi=1./fp
      Sint=0.
      do 10 if=1,ifnyq
         f=if*df
         fpf=fp/f
         fpf4=fpf*fpf*fpf*fpf
         fpf5=fpf4*fpf
         ffp=1./fpf
         if (f.lt.fp) then
            sigma=0.07
         else
            sigma=0.09
         endif
         fac=(ffp-1.)/sigma
         S(if)=fpi*fpf5*exp(-1.25*fpf4)*gamma**exp(-.5*fac*fac)
         Sint=Sint+S(if)
10    continue
      Sint=Sint*df

      beta=.125/Sint
      betHr2=beta*Hrms*Hrms
      Sint=0.
      Sf=0.
      Sff=0.
      do 20 if=1,ifnyq
         f=if*df
         S(if)=S(if)*betHr2
         amp=sqrt(2.*S(if)*df)
         phi=ran2(idum)*twopi
         a(2*if+1)=  amp*cos(phi) 
         a(2*if+2)= -amp*sin(phi)
         Sint=Sint+S(if)
         Sf=Sf+S(if)*f
         Sff=Sff+S(if)*f*f

20    continue                  
      Sint=Sint*df
      Sf=Sf*df
      Sff=Sff*df
      fm01=Sf/Sint
      fm02=sqrt(Sff/Sint)

      write(*,*)' fm01 = ',fm01
      write(*,*)' fm02 = ',fm02
      write(*,*)' Hrms spec    = ',sqrt(8.*Sint)
      
      call realft (a,N,-1)
      open(2,file='gen.dat',access='sequential',form='binary')
      write(2)(a(it),it=1,Nt)
      call stdev(a,Nt,sigma)
      Hrms=sqrt(8.)*sigma  
      write(*,*)' Hrms st dev. = ',sqrt(8.*Sint)
       
      open(3,file='spec.tek')
      write(3,'(a4)')'BL01'
      write(3,*)Nf/nstep,2
      write(3,'(f10.4,f10.5)')((if-1)*df,S(if),if=1,ifNyq,nstep) 

      end 

      FUNCTION RAN2(IDUM)
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112E-6)
      DIMENSION IR(97)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IDUM=MOD(IC-IDUM,M)
        DO 11 J=1,97
          IDUM=MOD(IA*IDUM+IC,M)
          IR(J)=IDUM
11      CONTINUE
        IDUM=MOD(IA*IDUM+IC,M)
        IY=IDUM
      ENDIF
      J=1+(97*IY)/M
      IF(J.GT.97.OR.J.LT.1)PAUSE
      IY=IR(J)
      RAN2=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      RETURN
      END


      subroutine stdev (X,N,sigma)
      real X(N)
      SX=0.
      SX2=0.
      do 10 i=1,N
         SX=SX+X(i)
10    continue
      Xmean=SX/N
      do 30 i=1,N 
         SX2=SX2+(X(i)-Xmean)**2
30    continue
      sigma=sqrt(SX2/N)
      end 

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
