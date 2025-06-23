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
      COMPLEX*16 P(NMAX),Q(NMAX),R(NMAX),S(NMAX),T(NMAX),U(NMAX)
      COMPLEX*16 PH0(NMAX),IU,OVR
      DIMENSION HA(NMAX),HB(NMAX-1),V(NMAX)
      DIMENSION HAN(NMAX),HBN(NMAX-1)
C
C     INITIAL DATA
C
      OPEN(20,FILE='tst4.dat')
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
C        V(N)=0
C        HA(N) = V(N) + 1/DX2
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
   90 CONTINUE
C
C     TIME EVOLUTION LOOP 1000
C
      DO 1000 ITIME=3,MTIME,5
      TT = DT*ITIME
C
C      OUTPUT WAVE FUNCTION
C
      IF(MOD(ITIME,10).EQ. 3) THEN
      ANORM=0
      OVR=0.0
      DO 150 N=1,NMAX
        ANORM=ANORM+ABS(S(N))**2
        OVR=OVR+DCONJG(S(N))*PH0(N)*EXP(-IU*EM*TT)
  150   CONTINUE
      WRITE(*,100) ITIME,EMAX*TT,ANORM-1, 
     &1.- ABS(OVR),DATAN2(DIMAG(OVR),DBLE(OVR))
     &,ABS((EM*DT)**5*(7./90.)) * ITIME
      WRITE(20,100) ITIME,EMAX*TT,ABS(ANORM-1),
     &ABS(1- ABS(OVR)),ABS(DATAN2(DIMAG(OVR),DBLE(OVR)))
     &,ABS((EM*DT)**5*(7./90.)) * ITIME
  100 FORMAT(1H ,I4,10E15.7)
      ENDIF
C
C     TIME EVOLUTION 
C
      DO 80 N=1,NMAX
      U(N) = (2*Q(N)-R(N)+2*S(N))/3.0
   80 CONTINUE
      T(1)= -4*IU*DT*(HA(1)*U(1)+HB(1)*U(2)) + P(1)
      DO 70 N=2,NMAX-1
      T(N)= -4*IU*DT*(HB(N-1)*U(N-1)+HA(N)*U(N)+HB(N)*U(N+1)) + P(N)
   70 CONTINUE
      T(NMAX)=-4*IU*DT*(HB(NMAX-1)*U(NMAX-1)+HA(NMAX)*U(NMAX)) + P(NMAX)
C
      DO 81 N=1,NMAX
      U(N) = (2*R(N)-S(N)+2*T(N))/3.0
   81 CONTINUE
      P(1) = -4*IU*DT*(HA(1)*U(1)+HB(1)*U(2)) +Q(1)
      DO 71 N=2,NMAX-1
      P(N)= -4*IU*DT*(HB(N-1)*U(N-1)+HA(N)*U(N)+HB(N)*U(N+1))+Q(N)
   71 CONTINUE
      P(NMAX)=-4*IU*DT*(HB(NMAX-1)*U(NMAX-1)+HA(NMAX)*U(NMAX)) + Q(NMAX) 
C
C
      DO 82 N=1,NMAX
      U(N) = (2*S(N)-T(N)+2*P(N))/3.0
   82 CONTINUE
      Q(1) = -4*IU*DT*(HA(1)*U(1)+HB(1)*U(2)) + R(1)
      DO 72 N=2,NMAX-1
      Q(N)= -4*IU*DT*(HB(N-1)*U(N-1)+HA(N)*U(N)+HB(N)*U(N+1))+R(N)
   72 CONTINUE
      Q(NMAX)=-4*IU*DT*(HB(NMAX-1)*U(NMAX-1)+HA(NMAX)*U(NMAX))+R(NMAX)
C
C
      DO 83 N=1,NMAX
      U(N) = (2*T(N)-P(N)+2*Q(N))/3.0
   83 CONTINUE
      R(1) = -4*IU*DT*(HA(1)*U(1)+HB(1)*U(2))+S(1)
      DO 73 N=2,NMAX-1
      R(N)= -4*IU*DT*(HB(N-1)*U(N-1)+HA(N)*U(N)+HB(N)*U(N+1))+S(N) 
   73 CONTINUE
      R(NMAX)=-4*IU*DT*(HB(NMAX-1)*U(NMAX-1)+HA(NMAX)*U(NMAX))+S(NMAX)
C
C
      DO 84 N=1,NMAX
      U(N) = (2*P(N)-Q(N)+2*R(N))/3.0
   84 CONTINUE
      S(1) = -4*IU*DT*(HA(1)*U(1)+HB(1)*U(2))+T(1)
      DO 74 N=2,NMAX-1
      S(N)= -4*IU*DT*(HB(N-1)*U(N-1)+HA(N)*U(N)+HB(N)*U(N+1))+T(N)
   74 CONTINUE
      S(NMAX)=-4*IU*DT*(HB(NMAX-1)*U(NMAX-1)+HA(NMAX)*U(NMAX))+T(NMAX) 
C
C
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
C
C     NORMALIZATION OF WAVE FUNCTION
C
      PA = 0.
      DO 100 N=1,NMAX
        PA = PA + ABS(PH0(N))**2
  100 CONTINUE
      PA = SQRT(PA)
      DO 110 N=1,NMAX
        PH0(N)=PH0(N)/PA
  110 CONTINUE
C      EM=-4*(-0.5/DX/DX) * SIN((PI*M)/(NMAX+1)/2.0)**2
      EM=2*(-0.5/DX/DX) * COS((PI*M)/(NMAX+1))
      RETURN
      END
