
C
C     SIMULATION OF N-ROOM QUANTUM MECHANICS
C
C (C) TOSHIAKI IITAKA 1994
C REFERENCE: T.Iitaka, Phys. Rev. E49 (1994) 4684.
C
      PROGRAM NROOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI=3.14159265358979323)
      PARAMETER(NMAX=100)
      COMPLEX*16 P(NMAX),Q(NMAX),R(NMAX),S(NMAX),T(NMAX),U(NMAX),V(NMAX)
      COMPLEX*16 W(NMAX),PH0(NMAX),IU,OVR
	      DIMENSION HA(NMAX),HB(NMAX-1),V0(NMAX)
C
C     INITIAL DATA
C
      OPEN(20,FILE='tst6.dat')
      WRITE(*,*) 'M, ALPHA'
      READ(*,*) M,ALPHA
      IU=(0.0,1.0)
      DX=1
      DX2=DX*DX
      EMAX=1/DX2
      DT=ALPHA/(EMAX*COS(M*PI/(NMAX+1)))
      MTIME=(2*PI/ALPHA) * 10
      MTIME=MIN(MTIME,10000)
C
C     CALCULATION OF HAMILTONIAN
C
      DO 50 N=1,NMAX
C        V0(N)=0
C        HA(N) = V0(N) + 1/DX2
      HA(N)=0
   50 CONTINUE
      DO 51 N=1,NMAX-1
          HB(N) = -0.5/DX2
   51 CONTINUE
C
C     SETTING INITIAL WAVE FUNCTION
C
      CALL EIGEN(M,DX,PH0,EM)
      DO 90 N=1,NMAX
      P(N)=PH0(N)
      Q(N)=PH0(N)*EXP(-IU*EM*DT)
      R(N)=PH0(N)*EXP(-IU*EM*2*DT)
      S(N)=PH0(N)*EXP(-IU*EM*3*DT)
      T(N)=PH0(N)*EXP(-IU*EM*4*DT)
      U(N)=PH0(N)*EXP(-IU*EM*5*DT)
   90 CONTINUE
C
C     TIME EVOLUTION LOOP 1000
C
      DO 1000 ITIME=5,MTIME,7
      TT = DT*ITIME
C
C      OUTPUT WAVE FUNCTION
C
      IF(MOD(ITIME,14).EQ. 5) THEN
      ANORM=0
      OVR=0.0
      DO 150 N=1,NMAX
        ANORM=ANORM+ABS(U(N))**2
        OVR=OVR+DCONJG(U(N))*PH0(N)*EXP(-IU*EM*TT)
  150   CONTINUE
      WRITE(*,100) ITIME,EMAX*TT,ANORM-1, 
     &1.- ABS(OVR),DATAN2(DIMAG(OVR),DBLE(OVR))
     &,ABS((EM*DT)**7*(41./840.)) * ITIME
      WRITE(20,100) ITIME,EMAX*TT,ABS(ANORM-1),
     &ABS(1- ABS(OVR)),ABS(DATAN2(DIMAG(OVR),DBLE(OVR)))
     &,ABS((EM*DT)**7*(41./840.)) * ITIME
  100 FORMAT(1H ,I4,10E15.7)
      ENDIF
C
C     TIME EVOLUTION 
C
      DO 80 N=1,NMAX
      W(N) = (11*Q(N)-14*R(N)+26*S(N)-14*T(N)+11*U(N))/20.0
   80 CONTINUE
      V(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + P(1)
      DO 70 N=2,NMAX-1
      V(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + P(N)
   70 CONTINUE
      V(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + P(NMAX)
C
      DO 81 N=1,NMAX
      W(N) = (11*R(N)-14*S(N)+26*T(N)-14*U(N)+11*V(N))/20.0
   81 CONTINUE
      P(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + Q(1)
      DO 71 N=2,NMAX-1
      P(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + Q(N)
   71 CONTINUE
      P(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + Q(NMAX)
C
      DO 82 N=1,NMAX
      W(N) = (11*S(N)-14*T(N)+26*U(N)-14*V(N)+11*P(N))/20.0
   82 CONTINUE
      Q(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + R(1)
      DO 72 N=2,NMAX-1
      Q(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + R(N)
   72 CONTINUE
      Q(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + R(NMAX)
C
      DO 83 N=1,NMAX
      W(N) = (11*T(N)-14*U(N)+26*V(N)-14*P(N)+11*Q(N))/20.0
   83 CONTINUE
      R(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + S(1)
      DO 73 N=2,NMAX-1
      R(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + S(N)
   73 CONTINUE
      R(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + S(NMAX)
C
      DO 84 N=1,NMAX
      W(N) = (11*U(N)-14*V(N)+26*P(N)-14*Q(N)+11*R(N))/20.0
   84 CONTINUE
      S(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + T(1)
      DO 74 N=2,NMAX-1
      S(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + T(N)
   74 CONTINUE
      S(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + T(NMAX)
C
      DO 85 N=1,NMAX
      W(N) = (11*V(N)-14*P(N)+26*Q(N)-14*R(N)+11*S(N))/20.0
   85 CONTINUE
      T(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + U(1)
      DO 75 N=2,NMAX-1
      T(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + U(N)
   75 CONTINUE
      T(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + U(NMAX)
C
      DO 86 N=1,NMAX
      W(N) = (11*P(N)-14*Q(N)+26*R(N)-14*S(N)+11*T(N))/20.0
   86 CONTINUE
      U(1)= -6*IU*DT*(HA(1)*W(1)+HB(1)*W(2)) + V(1)
      DO 76 N=2,NMAX-1
      U(N)= -6*IU*DT*(HB(N-1)*W(N-1)+HA(N)*W(N)+HB(N)*W(N+1)) + V(N)
   76 CONTINUE
      U(NMAX)=-6*IU*DT*(HB(NMAX-1)*W(NMAX-1)+HA(NMAX)*W(NMAX)) + V(NMAX)
 1000 CONTINUE
      STOP
      END
C
C     SETTING INITIAL WAVE FUNCTION
C
      SUBROUTINE EIGEN(M,DX,PH0,EM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI=3.14159265358979323)
      PARAMETER(NMAX=100,IMAX=1000)
      PARAMETER(EPS=1E-14)
      COMPLEX*16 PH0(NMAX),IU
      DIMENSION HA(NMAX),HB(NMAX-1),V(NMAX)
      IU=(0.0,1.0)
      DO 60 N=1,NMAX
        PH0(N)= SIN((PI*M*N)/(NMAX+1))
   60 CONTINUE
C     NORMALIZATION OF WAVE FUNCTION
      PA = 0.
      DO 100 N=1,NMAX
        PA = PA + ABS(PH0(N))**2
  100 CONTINUE
      PA = SQRT(PA)
      DO 110 N=1,NMAX
        PH0(N)=PH0(N)/PA
  110 CONTINUE
c      EM=-4*(-0.5/DX/DX) * SIN((PI*M)/(NMAX+1)/2.0)**2-1/dx/dx
      EM=2*(-0.5/DX/DX) * COS((PI*M)/(NMAX+1))
      RETURN
      END
