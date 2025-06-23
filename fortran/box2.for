C
C     SIMULATION OF 2-ROOM QUANTUM MECHANICS
C
C (C) TOSHIAKI IITAKA 1994
C REFERENCEÅF"Introduction to Quatum Dynamics" by Toshiaki Iitaka.
C (Maruzen Publish. Co., 1994,Tokyo; Parity Physics Course, Close Up)
C 
C
      PROGRAM TWOROOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI=3.14159265358979)
      PARAMETER(MTIME=10000)
      COMPLEX P,Q,R,IU
      DIMENSION P(2),Q(2),R(2),V(2),HA(2),HB(1)
C
C     I/O FILE
C
      OPEN(20,FILE='p1061.dat')
C
C     INITIAL DATA
C
      IU=(0.0,1.0)
      WRITE(*,*) 'INPUT PROBABILITY OF RIGHT AND LEFT ROOM'
      READ(*,*) PR,PL
      P(1)=DSQRT(PR)
      P(2)=DSQRT(PL)
C
C     SET HAMILTONIAN
C
      WRITE(*,*) 'INPUT ENERGY OF RIGHT AND LEFT ROOM'
      READ(*,*) HA(1),HA(2)
      WRITE(*,*) 'INPUT TRANSFER ENERGY A'
      READ(*,*) HB(1)
      DT=1E-3
C
C     NORMALIZATION OF WAVE FUNCTION
C
      PA = SQRT( CABS(P(1))**2 + CABS(P(2))**2)
      P(1)=P(1)/PA
      P(2)=P(2)/PA
C
C  SECOND WAVE FUNCTION
C
      R(1)= -IU*DT*( HA(1)*P(1) + HB(1)*P(2))
      R(2)= -IU*DT*( HB(1)*P(1) + HA(2)*P(2) )
      Q(1)= P(1)+R(1)
      Q(2)= P(2)+R(2)
      Q(1)= Q(1)-0.5*IU*DT*(HA(1)*R(1)+HB(1)*R(2))
      Q(2)= Q(2)-0.5*IU*DT*(HB(1)*R(1)+HA(2)*R(2))
C
C     TIME EVOLUTION LOOP 1000
C
      DO 1000 ITIME=1,MTIME,3
      T= DT*ITIME
C
C      OUTPUT
C
        PR=CABS(Q(1))**2
        PL=CABS(Q(2))**2
        WRITE(*,100)   T, Q(1), PR, -PR
        WRITE(20,100)  T, Q(1), PR, -PR
  100 FORMAT(1H ,10E15.7 )
C
C     TIME EVOLUTION * 3
C
      R(1) = -2*IU*DT*(HA(1)*Q(1)+HB(1)*Q(2))+P(1)
      R(2) = -2*IU*DT*(HB(1)*Q(1)+HA(2)*Q(2))+P(2)
      P(1) = -2*IU*DT*(HA(1)*R(1)+HB(1)*R(2))+Q(1)
      P(2) = -2*IU*DT*(HB(1)*R(1)+HA(2)*R(2))+Q(2)
      Q(1) = -2*IU*DT*(HA(1)*P(1)+HB(1)*P(2))+R(1)
      Q(2) = -2*IU*DT*(HB(1)*P(1)+HA(2)*P(2))+R(2)
 1000 CONTINUE
      END
