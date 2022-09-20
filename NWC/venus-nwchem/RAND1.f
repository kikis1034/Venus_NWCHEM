      FUNCTION RAND1(ISEED)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         FUNCTION RAND1(ISEED)
C         RANDOM NUMBER GENERATION USING THE MULTIPLICATIVE
C         CONGRUENTIAL METHOD.
C
C         WRITTEN BY DAVID SCHWENKE (NASA AMES)
C         FOR MACHINES WITH A MAXIMUM INTEGER OF (2**32 - 1) OR LARGER
C         AND HAS PRECISION OF 1 IN 2**64
C
C         AMENDED BY KIERAN F. LIM (STANFORD uNIVERSITY)
C         FOR MACHINES WITH A MAXIMUM INTEGER OF (2**16 - 1) OR LARGER
C         NOTE:  DEC MICROVAX HAS A MAXIMUM INTEGER OF (2**31 - 1)
C
C         THIS ROUTINE HITS INTEGER OVERFLOW ONLY IF THE MAXIMUM
C         INTEGER IS (2**16 - 1) = 65535 OR SMALLER
C
C         IF RANLST(N) IS THE CURRENT SEED, THEN THE NEXT ONE RANLST(N+1)
C         IS GIVEN BY
C
C         RANLST(N+1)=(A*RANLST(N)+C)MOD(M)
C
C         AND THE RANDOM NUMBER IS RANLST(N+1)/M.
C         THE ARITHMETIC IS CARRIED OUT USING 8 "DIGIT" NUMBERS
C         OF BASE 2**8 (=256), THUS ISEED IS 8 ELEMENT ARRAY WITH
C         EACH ELEMENT A "DIGIT". THE UNITS "DIGIT" IS THE FIRST
C         ELEMENT OF ISEED.
C         A NINTH "DIGIT" IS PROVIDED FOR OVERFLOWS IN THE INTERMEDIATE
C         STEPS FOR ARITHMETIC OPERATIONS.
C
C         PARAMETERS:
C
C          M = 2**64 = (1) (0) (0) (0) (0) (0) (0) (0) (0) BASE 2**8
C          A = 6,364,136,223,846,793,005 BASE 10
C            = (22609) (62509) (19605) (322557) BASE 2**16
C            = (88) (81) (244) (45) (76) (149) (127) (45) BASE 256
C          C = 0
C
C         THESE PARAMETERS ARE FROM D.E. KNUTH, "THE ART OF COMPUTER
C         PROGRAMMING", VOL. 2 OF SEMINUMERICAL ALGORITHMS" (ADDISON-WESLEY
C         READING,MASS. 1981) 2ND EDITION,
C         SECTION 3.3.4, TABLE 1 AT PAGES 102-104
C         AND ARE ATTRIBUTED TO C.E. HAYNES
C
      COMMON/RANCOM/RANLST(100),ISEED3(8),IBFCTR
      DIMENSION ISEED(8),IA(8),IC(8),ID(16)
C
C         BASE 2**16 INFORMATION
C
C         DATA IA/32557,19605,62509,22609/
C         DATA IC/0,0,0,0/ 
C         DATA BI/1.5258789062D-5/
C
C         BASE 2**8 INFORMATION
C
      DATA IA/45,127,149,76,45,244,81,88/
      DATA IC/0,0,0,0,0,0,0,0/
      DATA BI/3.90625D-3/
C
C         ID WILL EQUAL ISEED*IA+IC.
C         SET IT EQUAL TO IC
C
      DO I=1,8
	   ID(I)=IC(I)
         ID(I+8)=0
      ENDDO
C
C         FORM IA*ISEED+IC
C
      DO J=1,8
         DO I=1,9-J
            K=J+I-1
C
C         IP IS UNNORMALIZED PRODUCT OF K DIGIT OF RESULT SO FAR.
C         IT SHOULD NEVER EXCEED IBFCTR*(IBFCTR-1).
C
            IP=IA(J)*ISEED(I)
C
C         NOW NORMALIZE THE RESULT (CARRY FORWARD TO NEXT DIGITS IF
C         NECESSARY)
C
  120       CONTINUE
            IP=IP+ID(K)
            ID(K)=MOD(IP,IBFCTR)
            IP=IP/IBFCTR
            IF (IP.EQ.0.OR.K.EQ.8) GOTO 130
            K=K+1
            GOTO 120
  130       CONTINUE
         ENDDO
      ENDDO
C
      DO I=1,8
         ISEED(I)=ID(I)
      ENDDO
C
C         NOW DETERMINE FLOATING POINT RANDOM NUMBER
C
      RAND1=DBLE(ISEED(1))
      DO I=2,8
         RAND1=DBLE(ISEED(I))+RAND1*BI
      ENDDO
      RAND1=RAND1*BI
      RETURN
      END     
