C
C     SIMULATION OF N-ROOM QUANTUM MECHANICS
C
C (C) TOSHIAKI IITAKA 1994
C REFERENCEÅF"Introduction to Quatum Dynamics" by Toshiaki Iitaka.
C (Maruzen Publish. Co., 1994, Tokyo; Parity Physics Course, Close Up)
C 
C
      PROGRAM NROOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI=3.141592)
      PARAMETER(NMAX=100,MTIME=1500)
      COMPLEX P,Q,R,IU
      DIMENSION P(NMAX),Q(NMAX),R(NMAX),HA(NMAX),HB(NMAX-1)
C
C     I/O FILES
C
      OPEN(20,FILE='p202.dat')
C
C     INITIAL DATA
C
      IU=(0.0,1.0)
      DX=1
C
C     CALCULATION OF HAMILTONIAN
C
      A=1
      DT=1/A * 1E-1
      DO 50 N=1,NMAX
        HA(N) = 0
   50 CONTINUE
      DO 51 N=1,NMAX-1
          HB(N) = -A
   51 CONTINUE
C
C     SET INITIAL WAVE FUNCTION
C
C      WRITE(*,*)'INPUT THE CENTER OF WAVE PACKET 
C     &           AND ITS WIDTH AND MOMENTUM'
C      READ(*,*) X0,SG,PA
      X0=0.3*DX*NMAX
      SG=0.05*DX*NMAX
      PA=0.1*2*PI/DX
C
      DO 60 N=1,NMAX
        X=N*DX
        P(N)= EXP(-0.5*(X-X0)**2/SG**2) * EXP(IU*PA*(X-X0))
   60 CONTINUE
C
C     NORMALIZATION OF WAVE FUNCTION
C
      PA = 0.
      DO 100 N=1,NMAX
        PA = PA + CABS(P(N))**2
  100 CONTINUE
      PA = SQRT(PA)
      DO 110 N=1,NMAX
        P(N)=P(N)/PA
  110 CONTINUE
C
C  SECOND WAVE FUNCTION
C
      R(1)= -IU*DT*(HA(1)*P(1)+HB(1)*P(2))
      DO 61 N=2,NMAX-1
        R(N)= -IU*DT*( HB(N-1)*P(N-1) + HA(N)*P(N) + HB(N)*P(N+1) )
   61 CONTINUE
      R(NMAX)=-IU*DT*( HB(NMAX-1)*P(NMAX-1) + HA(NMAX)*P(NMAX) )
C
      DO 62 N=1,NMAX
       Q(N)=P(N)+R(N)
   62 CONTINUE
      Q(1)= Q(1)-0.5*IU*DT*(HA(1)*R(1)+HB(1)*R(2))
      DO 63 N=2,NMAX-1
        Q(N)= Q(N)-0.5*IU*DT*(HB(N-1)*R(N-1)+HA(N)*R(N)+HB(N)*R(N+1))
   63 CONTINUE
      Q(NMAX)=Q(NMAX)-0.5*IU*DT*(HB(NMAX-1)*R(NMAX-1)+HA(NMAX)*R(NMAX))
C
C     TIME EVOLUTION LOOP 1000
C
      DO 1000 ITIME=1,MTIME,3
      T= DT*ITIME
C
C      OUTPUT WAVE FUNCTION
C
      IF(MOD(ITIME,30).EQ. 1) THEN
      DO 150 N=1,NMAX
        WRITE(*,*)  N,CABS(Q(N))**2
        WRITE(20,*) CABS(Q(N))**2
  150   CONTINUE
        WRITE(20,*) ' '
      ENDIF
C
C     TIME EVOLUTION 
C
      R(1) = -2*IU*DT*(HA(1)*Q(1)+HB(1)*Q(2))+P(1)
      DO 70 N=2,NMAX-1
      R(N)=-2*IU*DT*(HB(N-1)*Q(N-1)+HA(N)*Q(N)+HB(N)*Q(N+1))+P(N)
   70 CONTINUE
      R(NMAX)=-2*IU*DT*(HB(NMAX-1)*Q(NMAX-1)+HA(NMAX)*Q(NMAX))+P(NMAX)
C
      P(1) = -2*IU*DT*(HA(1)*R(1)+HB(1)*R(2))+Q(1)
      DO 71 N=2,NMAX-1
      P(N)=-2*IU*DT*(HB(N-1)*R(N-1)+HA(N)*R(N)+HB(N)*R(N+1))+Q(N)
   71 CONTINUE
      P(NMAX)=-2*IU*DT*(HB(NMAX-1)*R(NMAX-1)+HA(NMAX)*R(NMAX))+Q(NMAX)
C
      Q(1) = -2*IU*DT*(HA(1)*P(1)+HB(1)*P(2))+R(1)
      DO 72 N=2,NMAX-1
      Q(N)=-2*IU*DT*(HB(N-1)*P(N-1)+HA(N)*P(N)+HB(N)*P(N+1))+R(N)
   72 CONTINUE
      Q(NMAX)=-2*IU*DT*(HB(NMAX-1)*P(NMAX-1)+HA(NMAX)*P(NMAX))+R(NMAX)
 1000 CONTINUE
      CLOSE(20)
      STOP
      END
