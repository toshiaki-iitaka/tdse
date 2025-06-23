C
C     SIMULATION OF 1D QUANTUM MECHANICS
C
C (C) TOSHIAKI IITAKA 1994
C REFERENCE: T.Iitaka, Phys. Rev. E49 (1994) 4684.
C
C   CHEBYSHEV SCHEME
C   REF: C.LEFORESTIER ET AL., J. COMP. PHYS. VOL.94 P.59 (1991)
C
      PROGRAM TEST
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI=3.14159265358979323)
      PARAMETER(NMAX=100)
      COMPLEX*16 PH0(NMAX),PH(NMAX),IU,OVR
      DIMENSION HA(NMAX),HB(NMAX-1),V(NMAX)
      DIMENSION HAN(NMAX),HBN(NMAX-1)
C
C     INITIAL DATA
C
      OPEN(20,FILE='tstch.dat')
      WRITE(*,*) 'M'
      READ(*,*) M
      IU=(0.0,1.0)
      DX=1
      DX2=DX*DX
      EMAX=1/DX2
C      DT=ALPHA/EMAX
C      MTIME=(2*PI/ALPHA) * 10
C      MTIME=MIN(MTIME,10000)
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
C NORMALIZATION OF HAMILTONIAN
C
      VMAX=0
      VMIN=0
      DO 10 N=1,NMAX
      IF(V(N).GT.VMAX) VMAX=V(N)
      IF(V(N).LT.VMIN) VMIN=V(N)
   10 CONTINUE
C      AKMAX= 0.5 * (PI/DX)**2
C      EMIN=VMIN
C      EMAX=VMAX+AKMAX
      EMAX=1/DX2
      EMIN=-1/DX2
      EGRID=EMAX-EMIN
      DO 20 N=1,NMAX
      HAN(N)=(HA(N)-EGRID/2-EMIN)*2/EGRID
      IF(N.LT.NMAX) HBN(N)=HB(N)*2/EGRID
   20 CONTINUE
      WRITE(*,*) 'EGRID=',EGRID
      WRITE(*,*) 'TIME,M'
      DO 1000 ITIME=1,40
      T=10**(ITIME/10.0)
      WRITE(*,*) 'T=',T
      WRITE(*,*) 'ESTIMATED STEP=',EGRID*T/2
      CALL EIGEN(M,DX,PH0,EM)
      WRITE(*,*) 'PERIOD=', 2*PI/EM
      CALL CH(HAN,HBN,EGRID,EMIN,PH0,T,PH)
      ANORM=0
      OVR=0.0
      DO 150 N=1,NMAX
        ANORM=ANORM+ABS(PH(N))**2
        OVR=OVR+DCONJG(PH(N))*PH0(N)*EXP(-IU*EM*T)
  150   CONTINUE
      WRITE(*,100) ITIME,EMAX*T,ANORM-1, 
     &1.- ABS(OVR),DATAN2(DIMAG(OVR),DBLE(OVR))
      WRITE(20,100) ITIME,EMAX*T,ABS(ANORM-1),
     &ABS(1- ABS(OVR)),ABS(DATAN2(DIMAG(OVR),DBLE(OVR)))
  100 FORMAT(1H ,I4,30E15.7)
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
      WRITE(*,*) 'I=',I, AN
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
