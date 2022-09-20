C
C         STATIONARY POINT SEARCH
C
C         DRIVER FOR BROYDEN'S METHOD.
C         IT IS A GLOBALLY CONVERGENT QUASI-NEWTON METHOD FOR FINDING
C         THE ROOTS OF NONLINEAR SYSTEMS OF EQUATIONS
C
C         ADAPTED FROM "NUMERICAL RECIPES, THE ART OF SCIENTIFIC 
C         COMPUTING, SECOND EDITION", by W.H. PRESS, S.A. TEUKOLSKY,
C         W.T. VETTERLING & B.P. FLANNERY, CAMBRIDGE UNIVERSITY PRESS.
C
      SUBROUTINE STATPT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      DIMENSION QOLD(NDA3),C(NDA3),D(NDA3),PDOLD(NDA3),TT(NDA3),
     *P(NDA3),QT(NDA3yf,NDA3yf),R(NDA3yf,NDA3yf),S(NDA3),G(NDA3),W(NDA3)
      LOGICAL RESTRT,SING,SKIP,CHECK,FIRST
      DATA EPS,TOLF,TOLMIN,TOLX,ALF,STPMX/1.0D-7,1.0D-3,1.0D-5,1.0D-7,
     *1.0D-4,100.0D0/FIRST/.TRUE./
      save ns,nip
C
    5 FORMAT(21X,'STATIONARY POINT SEARCH'///3X,'ENERGY ',
     *'GRADIENT ROOTS ARE SEARCHED BY THE MULTIDIMENSIONAL SECANT'/3X,
     *58HBROYDEN'S METHOD (IT'S A GLOBALLY CONVERGENT QUASI-NEWTON ,
     *'METHOD)'//5X,'MAXIMUM NUMBER OF ITERATIONS/STEPS : ',12X,
     *I8/5X,'NUMBER OF STEPS BETWEEN INTERMEDIATE PRINTOUT :  ',I8/
     *5X,'SELECTED CONVERGENCE CRITERION ON THE GRADIENT : ',G8.1//)
   15 FORMAT(/5X,'GRADIENT NORM : ',G15.8,/)
   25 FORMAT(5X,'MAXIMUM NUMBER OF ITERATIONS EXCEEDED WITHOUT'/
     *5X,'REACHING FULL CONVERGENCE TO A ZERO GRADIENT,'/5X,
     *'GIVEN THE SELECTED CONVERGENCE CRITERIA.')
   35 FORMAT(5X,'PARTIAL CONVERGENCE TO A ZERO GRADIENT WAS ACHIEVED',
     */,15X,'TREAT THE RESULTS WITH CARE')
   45 FORMAT(5X,'FULL CONVERGENCE TO A ZERO GRADIENT WAS ACHIEVED')
C
C         READ OPTIMIZATION PARAMETERS
C         AND CHECK THE COMPUTER LIMITING PRECISION 
C
      IF (FIRST) THEN
         NS=NC
         NIP=NX
         SCALE=TIME
         CALL PRECIS(DUM)
         IF (EPS*SCALE.LT.DUM) SCALE=DUM/EPS
         EPS=EPS*SCALE
         TOLF=TOLF*SCALE
         TOLMIN=TOLMIN*SCALE
         TOLX=TOLX*SCALE
         ALF=ALF*SCALE
         WRITE(6,5)NS,NIP,TOLF
         FIRST=.FALSE.
         RETURN
      ENDIF
C
C         INITIALIZE AND TEST THE INITIAL GUESS
C
      CHECK=.FALSE.
      CALL DVDQ
      SUM=0.0D0
      DO I=1,I3N
         SUM=SUM+PDOT(I)*PDOT(I)
      ENDDO
      F=0.5D0*SUM
      TEST=0.0D0
      DO I=1,I3N
         IF (ABS(PDOT(I)).GT.TEST) TEST=ABS(PDOT(I))
      ENDDO
      IF (TEST.LT.0.01D0*TOLF) GOTO 360
      SUM=0.0D0
      DO I=1,I3N
         SUM=SUM+Q(I)*Q(I)
      ENDDO
      STPMAX=STPMX*MAX(SQRT(SUM),DBLE(I3N))
      IF (NX.EQ.0) NX=NIP
      RESTRT=.TRUE.
C
C         START OF ITERATION LOOP
C
  350 NC=NC+1
C
C         INITIALIZE OR RE-INITIALIZE THE ESTIMATED JACOBIAN
C
      IF (RESTRT) THEN
         CALL FDJAC(R)
         CALL QRDCMP(R,I3N,NDA3,C,D,SING)
         IF (SING) PAUSE 'SINGULAR JACOBIAN IN BROYDN'
         DO I=1,I3N
            DO J=1,I3N
               QT(I,J)=0.0D0
            ENDDO
            QT(I,I)=1.0D0
         ENDDO
         DO K=1,I3N-1
            IF (C(K).NE.0.0D0) THEN
               DO J=1,I3N
                  SUM=0.0D0
                  DO I=K,I3N
                     SUM=SUM+R(I,K)*QT(I,J)
                  ENDDO
                  SUM=SUM/C(K)
                  DO I=K,I3N
                     QT(I,J)=QT(I,J)-SUM*R(I,K)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
         DO I=1,I3N
            R(I,I)=D(I)
            DO J=1,I-1
               R(I,J)=0.0D0
            ENDDO
         ENDDO
      ELSE
C
C         CARRY OUT BROYDEN UPDATE
C
         DO I=1,I3N
            S(I)=Q(I)-QOLD(I)
         ENDDO
         DO I=1,I3N
            SUM=0.0D0
            DO J=I,I3N
               SUM=SUM+R(I,J)*S(J)
            ENDDO
            TT(I)=SUM
         ENDDO
         SKIP=.TRUE.
         DO I=1,I3N
            SUM=0.0D0
            DO J=1,I3N
               SUM=SUM+QT(J,I)*TT(J)
            ENDDO
            W(I)=PDOT(I)-PDOLD(I)-SUM
            IF (ABS(W(I)).GE.EPS*(ABS(PDOT(I))+ABS(PDOLD(I)))) THEN
               SKIP=.FALSE.
            ELSE
               W(I)=0.0D0
            ENDIF
         ENDDO
         IF (.NOT.SKIP) THEN
            DO I=1,I3N
               SUM=0.0D0
               DO J=1,I3N
                  SUM=SUM+QT(I,J)*W(J)
               ENDDO
               TT(I)=SUM
            ENDDO
            DEN=0.0D0
            DO I=1,I3N
               DEN=DEN+S(I)*S(I)
            ENDDO
            DO I=1,I3N
               S(I)=S(I)/DEN
            ENDDO
            CALL QRUPDT(R,QT,I3N,NDA3,TT,S)
            DO I=1,I3N
               IF (R(I,I).EQ.0.0D0) PAUSE 'R SINGULAR IN BROYDN'
               D(I)=R(I,I)
            ENDDO
         ENDIF
      ENDIF
C
C         SETUP FOR THE LINE SEARCH
C
      DO I=1,I3N
         SUM=0.0D0
         DO J=1,I3N
            SUM=SUM+QT(I,J)*PDOT(J)
         ENDDO
         G(I)=SUM
      ENDDO
      DO I=I3N,1,-1
         SUM=0.0D0
         DO J=1,I
            SUM=SUM+R(J,I)*G(J)
         ENDDO
         G(I)=SUM
      ENDDO
      DO I=1,I3N
         QOLD(I)=Q(I)
         PDOLD(I)=PDOT(I)
      ENDDO
      FOLD=F
      DO I=1,I3N
         SUM=0.0D0
         DO J=1,I3N
            SUM=SUM+QT(I,J)*PDOT(J)
         ENDDO
         P(I)=-SUM
      ENDDO
C
C         SOLVE LINEAR EQUATIONS
C
      P(I3N)=P(I3N)/D(I3N)
      DO I=I3N-1,1,-1
         SUM=0.0D0
         DO J=I+1,I3N
            SUM=SUM+R(I,J)*P(J)
         ENDDO
         P(I)=(P(I)-SUM)/D(I)
      ENDDO
C
      CALL LNSRCH(QOLD,FOLD,G,P,F,STPMAX,ALF,TOLX,CHECK)
      TEST=0.0D0
      GNORM=0.0D0
      DO I=1,I3N
         GNORM=GNORM+PDOT(I)*PDOT(I)
         IF (ABS(PDOT(I)).GT.TEST) TEST=ABS(PDOT(I))
      ENDDO
      GNORM=SQRT(GNORM)
C
C         INTERMEDIATE PRINTOUT
C
      IF (NC.GE.NX) THEN
         WRITE(6,15)GNORM
         CALL ENERGY
         CALL GWRITE
         NX=NX+NIP
      ENDIF
C
C         CHECK FOR SUCCESS/FAILURE OF LINE SEARCH
C
      IF (TEST.LT.TOLF) THEN
         CHECK=.FALSE.
         GOTO 360
      ENDIF
      IF (CHECK) THEN
         IF (RESTRT) THEN
            GOTO 360
         ELSE
            TEST=0.0D0
            DEN=MAX(F,0.5D0*I3N)
            DO I=1,I3N
               TEMP=ABS(G(I))*MAX(ABS(Q(I)),1.0D0)/DEN
               IF (TEMP.GT.TEST) TEST=TEMP
            ENDDO
            IF (TEST.LT.TOLMIN) THEN
               GOTO 360
            ELSE
               RESTRT=.TRUE.
            ENDIF
         ENDIF
      ELSE
C
C         SUCCESSFUL STEP
C         CHECK FOR CONVERGENCE
C
         RESTRT=.FALSE.
         TEST=0.0D0
         DO I=1,I3N
            TEMP=(ABS(Q(I)-QOLD(I)))/MAX(ABS(Q(I)),1.0D0)
            IF (TEMP.GT.TEST) TEST=TEMP
          ENDDO
         IF (TEST.LT.TOLX) GOTO 360
      ENDIF
C
      IF (NC.LT.NS) GOTO 350
      WRITE(6,25)
      RETURN
  360 CONTINUE
C
      IF (CHECK) THEN
         WRITE(6,35)
      ELSE
         WRITE(6,45)
      ENDIF
C
      WRITE(6,15)GNORM
      CALL ENERGY
      CALL GWRITE
      CALL NMODE(NATOMS,0)
      RETURN
      END
C
C         COMPUTE THE FORWARD DIFFERENCE APPROXIMATION TO THE JACOBIAN
C
      SUBROUTINE FDJAC(DF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      PARAMETER(EPS=1.0D-7)
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      DIMENSION DF(NDA3yf,NDA3yf),F(NDA3)
      DO I=1,I3N
         F(I)=PDOT(I)
      ENDDO
      DO J=1,I3N
         TEMP=Q(J)
         H=EPS*ABS(TEMP)
         IF (H.EQ.0.0D0) H=EPS
         Q(J)=TEMP+H
         H=Q(J)-TEMP
         CALL DVDQ
         Q(J)=TEMP
         DO I=1,I3N
            DF(I,J)=(PDOT(I)-F(I))/H
         ENDDO
      ENDDO
      RETURN
      END
C
C         PERFORM LINE SEARCH
C
      SUBROUTINE LNSRCH(QOLD,FOLD,G,P,F,STPMAX,ALF,TOLX,CHECK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO
      LOGICAL CHECK
      DIMENSION G(NDA3),P(NDA3),QOLD(NDA3)
      CHECK=.FALSE.
C
C         DETERMINE THE SIZE OF THE NEWTON STEP
C
      SUM=0.0D0
      DO I=1,I3N
         SUM=SUM+P(I)*P(I)
      ENDDO
      SUM=SQRT(SUM)
      IF (SUM.GT.STPMAX) THEN
         DO I=1,I3N
            P(I)=P(I)*STPMAX/SUM
         ENDDO
      ENDIF
      SLOPE=0.0D0
      DO I=1,I3N
         SLOPE=SLOPE+G(I)*P(I)
      ENDDO
      TEST=0.0D0
      DO I=1,I3N
         TEMP=ABS(P(I))/MAX(ABS(QOLD(I)),1.0D0)
         IF (TEMP.GT.TEST) TEST=TEMP
      ENDDO
      ALAMIN=TOLX/TEST
      ALAM=1.0D0
C
C         START OF ITERATION LOOP
C
   50 CONTINUE
      DO I=1,I3N
         Q(I)=QOLD(I)+ALAM*P(I)
      ENDDO
      CALL DVDQ
      SUM=0.0D0
      DO I=1,I3N
         SUM=SUM+PDOT(I)*PDOT(I)
      ENDDO
      F=0.5D0*SUM
C
C         REFINE THE SIZE OF THE NEWTON STEP
C
      IF (ALAM.LT.ALAMIN) THEN
         DO I=1,I3N
            Q(I)=QOLD(I)
         ENDDO
         CHECK=.TRUE.
         RETURN
      ELSE IF(F.LE.FOLD+ALF*ALAM*SLOPE)THEN
         RETURN
      ELSE
         IF (ALAM.EQ.1.0D0) THEN
            TMPLAM=-SLOPE/(2.0D0*(F-FOLD-SLOPE))
         ELSE
            RHS1=F-FOLD-ALAM*SLOPE
            RHS2=F2-FOLD2-ALAM2*SLOPE
            A=(RHS1/ALAM**2.0D0-RHS2/ALAM2**2.0D0)/(ALAM-ALAM2)
            B=(-ALAM2*RHS1/ALAM**2.0D0+ALAM*RHS2/ALAM2**2.0D0)
     *         /(ALAM-ALAM2)
            IF (A.EQ.0.0D0) THEN
               TMPLAM=-SLOPE/(2.0D0*B)
            ELSE
               DISC=B*B-3.0D0*A*SLOPE
               TMPLAM=(-B+SQRT(DISC))/(3.0D0*A)
            ENDIF
            IF (TMPLAM.GT.0.5D0*ALAM) TMPLAM=0.5D0*ALAM
         ENDIF
      ENDIF
      ALAM2=ALAM
      F2=F
      FOLD2=FOLD
      ALAM=MAX(TMPLAM,0.1D0*ALAM)
      GOTO 50
      END
C
C         MATRIX UTILITY ROUTINE
C         CONSTRUCT THE QR DECOMPOSITION OF A MATRIX
C
      SUBROUTINE QRDCMP(A,N,NP,C,D,SING)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NP,NP),C(NP),D(NP)
      LOGICAL SING
      SING=.FALSE.
      SCALE=0.0D0
      DO K=1,N-1
         DO I=K,N
            SCALE=MAX(SCALE,ABS(A(I,K)))
         ENDDO
         IF (SCALE.EQ.0.0D0) THEN
            SING=.TRUE.
            C(K)=0.0D0
            D(K)=0.0D0
         ELSE
            DO I=K,N
               A(I,K)=A(I,K)/SCALE
            ENDDO
            SUM=0.0D0
            DO I=K,N
               SUM=SUM+A(I,K)*A(I,K)
            ENDDO
            SIGMA=SIGN(SQRT(SUM),A(K,K))
            A(K,K)=A(K,K)+SIGMA
            C(K)=SIGMA*A(K,K)
            D(K)=-SCALE*SIGMA
            DO J=K+1,N
               SUM=0.0D0
               DO I=K,N
                  SUM=SUM+A(I,K)*A(I,J)
               ENDDO
               TAU=SUM/C(K)
               DO I=K,N
                  A(I,J)=A(I,J)-TAU*A(I,K)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      D(N)=A(N,N)
      IF(D(N).EQ.0.0D0)SING=.TRUE.
      RETURN
      END
C
C         MATRIX UTILITY ROUTINE
C         UPDATE THE QR DECOMPOSITION OF A MATRIX
C
      SUBROUTINE QRUPDT(R,QT,N,NP,U,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NP,NP),QT(NP,NP),U(NP),V(NP)
      DO K=N,1,-1
         IF (U(K).NE.0.0D0) GOTO 20
      ENDDO
      K=1
   20 DO I=K-1,1,-1
         CALL JACROT(R,QT,N,NP,I,U(I),-U(I+1))
         IF (U(I).EQ.0.0D0) THEN
            U(I)=ABS(U(I+1))
         ELSEIF (ABS(U(I)).GT.ABS(U(I+1))) THEN
            U(I)=ABS(U(I))*SQRT(1.0D0+(U(I+1)/U(I))**2.0D0)
         ELSE
            U(I)=ABS(U(I+1))*SQRT(1.0D0+(U(I)/U(I+1))**2.0D0)
         ENDIF
      ENDDO
      DO J=1,N
         R(1,J)=R(1,J)+U(1)*V(J)
      ENDDO
      DO I=1,K-1
         CALL JACROT(R,QT,N,NP,I,R(I,I),-R(I+1,I))
      ENDDO
      RETURN
      END
C
C         MATRIX UTILITY ROUTINE
C         CARRY OUT A JACOBI ROTATION ON MATRIX ELEMENTS
C
      SUBROUTINE JACROT(R,QT,N,NP,I,A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NP,NP),QT(NP,NP)
      IF (A.EQ.0.0D0) THEN
         C=0.0D0
         S=SIGN(1.0D0,B)
      ELSEIF (ABS(A).GT.ABS(B)) THEN
         FACT=B/A
         C=SIGN(1.0D0/SQRT(1.0D0+FACT*FACT),A)
         S=FACT*C
      ELSE
         FACT=A/B
         S=SIGN(1.0D0/SQRT(1.0D0+FACT*FACT),B)
         C=FACT*S
      ENDIF
      DO J=I,N
         Y=R(I,J)
         W=R(I+1,J)
         R(I,J)=C*Y-S*W
         R(I+1,J)=S*Y+C*W
      ENDDO
      DO J=1,N
         Y=QT(I,J)
         W=QT(I+1,J)
         QT(I,J)=C*Y-S*W
         QT(I+1,J)=S*Y+C*W
      ENDDO
      RETURN
      END
C
C         DETERMINE THE COMPUTER LIMITING PRECISION 
C         THIS IS THE SMALLEST NUMBER EPS FOR WHICH 1+EPS.NE.1
C
      SUBROUTINE PRECIS(EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EPS=1.0D0
   10 IF (1.0D0+(0.5D0*EPS).EQ.1.0D0) RETURN
      EPS=0.5D0*EPS
      GOTO 10
      END
