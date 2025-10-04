C   IMSL ROUTINE NAME   - FFTRC
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - HP9000/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE THE FAST FOURIER TRANSFORM OF A
C                           REAL VALUED SEQUENCE
C
C   USAGE               - CALL FFTRC (A,N,X,IWK,WK)
C
C   ARGUMENTS    A      - INPUT REAL VECTOR OF LENGTH N WHICH
C                           CONTAINS THE DATA TO BE TRANSFORMED.
C                N      - INPUT NUMBER OF DATA POINTS TO BE TRANSFORMED.
C                            N MUST BE A POSITIVE EVEN INTEGER.
C                X      - OUTPUT COMPLEX VECTOR OF LENGTH N/2+1
C                           CONTAINING THE FIRST N/2+1 COEFFICIENTS OF
C                           THE FOURIER TRANSFORM. THE REMAINING
C                           COEFFICIENTS MAY BE DETERMINED BY
C                             X(N+2-I) = CONJG(X(I)), FOR I=2,...,N/2.
C                IWK    - INTEGER WORK VECTOR.
C                           IF N IS A POWER OF 2, THEN IWK SHOULD BE OF
C                           LENGTH M WHERE N=2**M.
C                           OTHERWISE, IWK SHOULD BE OF LENGTH
C                           6*(N/2)+150.
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS)
C                WK     - REAL WORK VECTOR OF LENGTH 6*(N/2)+150.
C                           WK IS NOT USED IF N IS A POWER OF 2.
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS)
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - FFTCC,FFT2C
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  FFTRC COMPUTES THE FOURIER TRANSFORM, X, ACCORDING
C                TO THE FOLLOWING FORMULA;
C
C                  X(K+1) = SUM FROM J = 0 TO N-1 OF
C                           A(J+1)*CEXP((0.0,(2.0*PI*J*K)/N))
C                  FOR K=0,1,...,N/2 AND PI=3.1415...
C
C                THE USER CAN COMPUTE THE REMAINING X VALUES BY
C                PERFORMING THE FOLLOWING STEPS;
C
C                     ND2 = N/2
C                     DO 10 I=2,ND2
C                        X(N+2-I) = CONJG(X(I))
C                  10 CONTINUE
C            2.  FFTRC CAN BE USED TO COMPUTE
C
C                  X(K+1) = (1/N)*SUM FROM J = 0 TO N-1 OF
C                           A(J+1)*CEXP((0.0,(-2.0*PI*J*K)/N))
C                  FOR K=0,1,...,N/2 AND PI=3.1415...
C
C                BY PERFORMING THE FOLLOWING STEPS;
C
C                     CALL FFTRC (A,N,X,IWK,WK)
C                     ND2P1 = N/2+1
C                     DO 10 I=1,ND2P1
C                        X(I) = CONJG(X(I))/N
C                  10 CONTINUE
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FFTRC  (A,N,X,IWK,WK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IWK(1)
      DOUBLE PRECISION   A(N),WK(1)
      COMPLEX*16         X(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ND2P1,ND2,I,MTWO,M,IMAX,ND4,NP2,K,NMK,J
      DOUBLE PRECISION   RPI,ZERO,ONE,HALF,THETA,TP,G(2),B(2),Z(2),AI,
     1                   AR
      COMPLEX*16         XIMAG,ALPH,BETA,GAM,S1,ZD
      EQUIVALENCE        (GAM,G(1)),(ALPH,B(1)),(Z(1),AR),(Z(2),AI),
     1                   (ZD,Z(1))
      DATA               ZERO/0.0D0/,HALF/0.5D0/,ONE/1.0D0/,IMAX/24/
      DATA               RPI/3.141592653589793D0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (N .NE. 2) GO TO 5
C                                  N EQUAL TO 2
      ZD = DCMPLX(A(1),A(2))
      THETA = AR
      TP = AI
      X(2) = DCMPLX(THETA-TP,ZERO)
      X(1) = DCMPLX(THETA+TP,ZERO)
      GO TO 9005
    5 CONTINUE
C                                  N GREATER THAN 2
      ND2 = N/2
      ND2P1 = ND2+1
C                                  MOVE A TO X
      J = 1
      DO 6 I=1,ND2
         X(I) = DCMPLX(A(J),A(J+1))
         J = J+2
    6 CONTINUE
C                                  COMPUTE THE CENTER COEFFICIENT
      GAM = DCMPLX(ZERO,ZERO)
      DO 10 I=1,ND2
         GAM = GAM + X(I)
   10 CONTINUE
      TP = G(1)-G(2)
      GAM = DCMPLX(TP,ZERO)
C                                  DETERMINE THE SMALLEST M SUCH THAT
C                                  N IS LESS THAN OR EQUAL TO 2**M
      MTWO = 2
      M = 1
      DO 15 I=1,IMAX
         IF (ND2 .LE. MTWO) GO TO 20
         MTWO = MTWO+MTWO
         M = M+1
   15 CONTINUE
   20 IF (ND2 .EQ. MTWO) GO TO 25
C                                  N IS NOT A POWER OF TWO, CALL FFTCC
      CALL FFTCC (X,ND2,IWK,WK)
      GO TO 30
C                                  N IS A POWER OF TWO, CALL FFT2C
   25 CALL FFT2C (X,M,IWK)
   30 ALPH = X(1)
      X(1) = B(1) + B(2)
      ND4 = (ND2+1)/2
      IF (ND4 .LT. 2) GO TO 40
      NP2 = ND2 + 2
      THETA = RPI/ND2
      TP = THETA
      XIMAG = DCMPLX(ZERO,ONE)
C                                  DECOMPOSE THE COMPLEX VECTOR X
C                                  INTO THE COMPONENTS OF THE TRANSFORM
C                                  OF THE INPUT DATA.
      DO 35 K = 2,ND4
         NMK = NP2 - K
         S1 = DCONJG(X(NMK))
         ALPH = X(K) + S1
         BETA = XIMAG*(S1-X(K))
         S1 = DCMPLX(DCOS(THETA),DSIN(THETA))
         X(K) = (ALPH+BETA*S1)*HALF
         X(NMK) = DCONJG(ALPH-BETA*S1)*HALF
         THETA = THETA + TP
   35 CONTINUE
   40 CONTINUE
      X(ND2P1) = GAM
 9005 RETURN
      END
