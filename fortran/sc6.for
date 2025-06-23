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
      COMPLEX*16 W(NMAX),PH0(NMAX),PH(NMAX),IU,OVR
      DIMENSION HA(NMAX),HB(NMAX-1),V0(NMAX),HAN(NMAX),HBN(NMAX)
C
C     INITIAL DATA
C
      OPEN(20,FILE='sc6.dat')
      WRITE(*,*) 'ALPHA'
      READ(*,*) ALPHA
      IU=(0.0,1.0)
      DX=1
      DX2=DX*DX
      EMAX=1/DX2+1
      DT=ALPHA/EMAX
      MTIME=(2*PI/ALPHA) * 40
      MTIME=MIN(MTIME,100000)
C
C     CALCULATION OF HAMILTONIAN
C
      DO 50 N=1,NMAX
      IF(45.LE.N .AND. N.LE.55) THEN
        V0(N)=0.05*COS((I-50)/5.0 * PI/2)**2
      ELSE
        V0(N)=0
      ENDIF
       HA(N) = V0(N)
c       HA(N) = V0(N) + 1/DX2
C      HA(N)=0
   50 CONTINUE
      DO 51 N=1,NMAX-1
          HB(N) = -0.5/DX2
   51 CONTINUE
C
C NORMALIZATION OF HAMILTONIAN
C
      VMAX=0
      VMIN=0
      DO 10 N=1,NMAX
      IF(V0(N).GT.VMAX) VMAX=V0(N)
      IF(V0(N).LT.VMIN) VMIN=V0(N)
   10 CONTINUE
C      AKMAX= 0.5 * (PI/DX)**2
C      EMIN=VMIN
C      EMAX=VMAX+AKMAX
      VMAX=1.0
      WRITE(*,*) 'VMAX=',VMAX
      EMAX=1/DX2 +VMAX
      EMIN=-1/DX2
      EGRID=EMAX-EMIN
      DO 20 N=1,NMAX
      HAN(N)=(HA(N)-EGRID/2-EMIN)*2/EGRID
      IF(N.LT.NMAX) HBN(N)=HB(N)*2/EGRID
   20 CONTINUE
      WRITE(*,*) 'EGRID=',EGRID
      WRITE(*,*) 'TIME,M'
C
C SET INITAIL FUNCTION
C
C
C     SETTING INITIAL WAVE FUNCTION
C
C      WRITE(*,*)'INPUT THE CENTER OF WAVE PACKET 
C     &           AND ITS WIDTH AND MOMENTUM'
C      READ(*,*) X0,SG,PA
      X0=0.25*DX*NMAX
      SG=0.1*DX*NMAX
      PA=0.1*2*PI/DX
      DO 60 N=1,NMAX
        X=N*DX
        PH0(N)= EXP(-0.5*(X-X0)**2/SG**2) * EXP(IU*PA*(X-X0))
   60 CONTINUE
C
C     NORMALIZATION OF WAVE FUNCTION
C
      PA = 0.
      DO 101 N=1,NMAX
        PA = PA + ABS(PH0(N))**2
  101 CONTINUE
      PA = SQRT(PA)
      DO 110 N=1,NMAX
        PH0(N)=PH0(N)/PA
        P(N)=PH0(N)
  110 CONTINUE
      CALL CH(HAN,HBN,EGRID,EMIN,PH0,DT,Q)
      CALL CH(HAN,HBN,EGRID,EMIN,PH0,2*DT,R)
      CALL CH(HAN,HBN,EGRID,EMIN,PH0,3*DT,S)
      CALL CH(HAN,HBN,EGRID,EMIN,PH0,4*DT,T)
      CALL CH(HAN,HBN,EGRID,EMIN,PH0,5*DT,U)
C
C     TIME EVOLUTION LOOP 1000
C
      DO 1000 ITIME=5,MTIME,7
      TT = DT*ITIME
C
C      OUTPUT WAVE FUNCTION
C
      IF(MOD(ITIME,28).EQ. 5) THEN
      CALL CH(HAN,HBN,EGRID,EMIN,PH0,TT,PH)
      ANORM=0
      OVR=0.0
      DO 150 N=1,NMAX
        ANORM=ANORM+ABS(U(N))**2
        OVR=OVR+DCONJG(U(N))*PH(N)
  150   CONTINUE
      WRITE(*,100) ITIME,EMAX*TT,ANORM-1, 
     &1.- ABS(OVR),DATAN2(DIMAG(OVR),DBLE(OVR))
      WRITE(20,100) ITIME,EMAX*TT,ABS(ANORM-1),
     &ABS(1- ABS(OVR)),ABS(DATAN2(DIMAG(OVR),DBLE(OVR)))
  100 FORMAT(1H ,I7,30E15.7)
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
C SUBROUTINE CH
C
      SUBROUTINE CH(HA,HB,EGRID,EMIN,PH0,T,PH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI=3.14159265358979323)
      PARAMETER(NMAX=100,IMAX=1000)
      PARAMETER(EPS=1E-16)
      COMPLEX*16 P(NMAX),Q(NMAX),R(NMAX),PH(NMAX),PH0(NMAX),IU
      DIMENSION HA(NMAX),HB(NMAX-1),V(NMAX)
      IU=(0.0,1.0)
C
      DO 20 N=1,NMAX
      P(N)=PH0(N)
   20 CONTINUE
      Q(1)   = -IU*(HA(1)*P(1)+HB(1)*P(2))
      DO 30 N=2,NMAX-1
      Q(N)   = -IU*(HB(N-1)*P(N-1)+HA(N)*P(N)+HB(N)*P(N+1))
   30 CONTINUE
      Q(NMAX)= -IU*(HB(NMAX-1)*P(NMAX-1)+HA(NMAX)*P(NMAX))
C
      AN=DJBES(0,EGRID*T/2)
      DO 40 N=1,NMAX
      PH(N)=AN*P(N)
   40 CONTINUE
      AN=2*DJBES(1,EGRID*T/2)
      DO 45 N=1,NMAX
      PH(N)=PH(N)+AN*Q(N)
   45 CONTINUE
      DO 1000 I=2, IMAX, 3
C
      AN=2*DJBES(I,EGRID*T/2)
      R(1) = -2*IU*(HA(1)*Q(1)+HB(1)*Q(2))+P(1)
      PH(1)=PH(1) +AN*R(1)
      DO 70 N=2,NMAX-1
      R(N)=-2*IU*(HB(N-1)*Q(N-1)+HA(N)*Q(N)+HB(N)*Q(N+1))+P(N)
      PH(N)=PH(N)+AN*R(N)
   70 CONTINUE
      R(NMAX)=-2*IU*(HB(NMAX-1)*Q(NMAX-1)+HA(NMAX)*Q(NMAX))+P(NMAX)
      PH(NMAX)=PH(NMAX)+AN*R(NMAX)
C
      AN=2*DJBES(I+1,EGRID*T/2)
      P(1) = -2*IU*(HA(1)*R(1)+HB(1)*R(2))+Q(1)
      PH(1)=PH(1) +AN*P(1)
      DO 71 N=2,NMAX-1
      P(N)=-2*IU*(HB(N-1)*R(N-1)+HA(N)*R(N)+HB(N)*R(N+1))+Q(N)
      PH(N)=PH(N)+AN*P(N)
   71 CONTINUE
      P(NMAX)=-2*IU*(HB(NMAX-1)*R(NMAX-1)+HA(NMAX)*R(NMAX))+Q(NMAX)
      PH(NMAX)=PH(NMAX)+AN*P(NMAX)
C
      AN=2*DJBES(I+2,EGRID*T/2)
      Q(1) = -2*IU*(HA(1)*P(1)+HB(1)*P(2))+R(1)
      PH(1)=PH(1) +AN*Q(1)
      DO 72 N=2,NMAX-1
      Q(N)=-2*IU*(HB(N-1)*P(N-1)+HA(N)*P(N)+HB(N)*P(N+1))+R(N)
      PH(N)=PH(N)+AN*Q(N)
   72 CONTINUE
      Q(NMAX)=-2*IU*(HB(NMAX-1)*P(NMAX-1)+HA(NMAX)*P(NMAX))+R(NMAX)
      PH(NMAX)=PH(NMAX)+AN*Q(NMAX)
C      WRITE(*,*) 'I=',I, AN
      IF(ABS(AN).LT.EPS) GOTO 2000
 1000 CONTINUE
      WRITE(*,*) 'CHEBYSHEV NOT CONVERGE'
 2000 CONTINUE
      DO 80 N=1,NMAX
      PH(N)=PH(N)*EXP(-IU*(EGRID/2+EMIN)*T)
   80 CONTINUE
      RETURN
      END
C
C  BESSEL FUNCTION DJBES(N,X)
C
      REAL FUNCTION DJBES * 8 (N,X)                                     
      REAL *8 X,Z,BJ,QJ,S,T1,T2,T3,ONE,D55                              
C      CALL OVERFL(KJ)                                                   
      ONE=1.0D0                                                         
      D55=1.0D-55                                                       
      NN=IABS(N)                                                        
      XA=DABS(X)                                                        
      IF(NN-30000)10,900,900                                            
   10 IF(XA-10.0)11,13,13                                               
   11 IF(XA-1.0)12,12,15                                                
   12 IF(XA-0.00002)1,1,16                                              
   13 IF(XA-100.0)17,14,14                                              
   14 IF(XA-30000.0)18,900,900                                          
    1 T1=0.5D0*X                                                        
      IF(NN) 3,2,3                                                      
    2 DJBES = ONE-T1*T1                                                 
      GO TO 1000                                                        
    3 IF(XA .LE. 1.0D-77) GO TO 500                                     
      T2=ONE                                                            
      T3=ONE                                                            
      DO 5 I=1,NN                                                       
      IF(DABS(T3) .LE. 1.0D-77*DABS(T2/T1)) GO TO 500                   
      T3=T3*T1/T2                                                       
      T2=T2+ONE                                                         
    5 CONTINUE                                                          
      BJ=T3*(ONE-T1*T1/T2)                                              
      GO TO 300                                                         
   15 L=1.4*XA+14.0                                                     
      GO TO 20                                                          
   16 L=14                                                              
      GO TO 20                                                          
   17 L=0.27*XA+27.0                                                    
      GO TO 20                                                          
   18 L=0.073*XA+47.0                                                   
   20 Z=2.0D0/X                                                         
      NM=MAX0(NN,IFIX(XA))+L                                            
      T3=0.0D0                                                          
      T2=1.0D-75                                                        
      S=0.0D0                                                           
      IF( MOD(NM,2)) 22,21,22                                           
   21 NM=NM+1                                                           
   22 DO 100 I=1,NM,2                                                   
      K=NM-I+1                                                          
      T1=DFLOAT(K+1)*T2*Z-T3                                            
      IF(NN-K) 40,30,40                                                 
   30 QJ=T1                                                             
   40 K=K-1                                                             
      T3=T2                                                             
      T2=T1                                                             
      T1=DFLOAT(K+1)*T2*Z-T3                                            
      IF(NN-K) 50,45,50                                                 
   45 QJ=T1                                                             
   50 S=S+T1                                                            
      IF( DABS(S)-1.0D55 ) 80,60,60                                     
   60 T1=T1*D55                                                         
      T2=T2*D55                                                         
      S =S*D55                                                          
      IF (NN-K) 80,70,70                                                
   70 QJ=QJ*D55                                                         
   80 T3=T2                                                             
      T2=T1                                                             
  100 CONTINUE                                                          
      S=S+S-T1                                                          
      BJ=QJ/S                                                           
  300 IF(N) 700,600,600                                                 
  500 DJBES = 0.0D0                                                     
      GO TO 1000                                                        
  600 DJBES = BJ                                                        
      GO TO 1000                                                        
  700 IF( MOD(NN,2) ) 800,600,800                                       
  800 DJBES = -BJ                                                       
      GO TO 1000                                                        
  900 WRITE(6,1001) N,X                                                 
      GO TO 500                                                         
 1000 RETURN                                                            
 1001 FORMAT( 1H ,5X,'THE VALUE OF DJBES IS NOT ACCURATE.  N=',I7,'  , X
     *=',D23.16 )                                                       
      END                                                               
