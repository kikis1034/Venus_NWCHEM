      SUBROUTINE GLPAR(X1,X2,X,W,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'SIZES'
C
C         CALCULATE PARAMETERS FOR GAUSS-LEGENDRE QUADRATURE
C
C         THE BASIC OUTLINE OF THE ALGORITHM FOR CALCULATING THE
C         ROOTS OF LEGENDRE POLYNOMIALS (INITIAL GUESS AND NEWTON-
C         RAPHSON ROOT REFINEMENT) IS TAKEN FROM
C
C         "NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING",
C         W.H. PRESS, S.A. TEUKOLSKY, W.T. VETTERLING & B.P. FLANNERY,
C         CAMBRIDGE UNIVERSITY PRESS.
C
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      DIMENSION XX(NCF),DPL(NCF),X(NCF),W(NCF)
      LOGICAL FIRST
      DATA FIRST/.TRUE./,nold/0/
      SAVE FIRST,nold,xx,dpl
C
C         LEGENDRE POLYNOMIAL ROOTS
C
      M=(N+1)/2
      IF (FIRST.or.(nold.ne.n)) THEN
         DO I=1,M
            ICOUNT=0
            DUM=1.0D0
            XX(I)=COS(PI*(I-0.25D0)/(N+0.5D0))
C
C               NEWTON-RAPHSON FORMULA
C
            DO WHILE (ABS(DUM).GT.1.0D-13)
               DUM=PL(XX(I),N,DPL(I))
               XX(I)=XX(I)-DUM/DPL(I)
            ENDDO
         ENDDO
         FIRST=.FALSE.
         nold=n
      ENDIF
C
C         INTEGRATION POINTS AND WEIGHTS
C
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO I=1,M
         X(I)=XM-XL*XX(I)
         X(N+1-I)=XM+XL*XX(I)
         W(I)=XL/((1.0D0-XX(I)*XX(I))*DPL(I)*DPL(I))/0.5D0
         W(N+1-I)=W(I)
      ENDDO
C
      RETURN
      END
C
      FUNCTION PL(X,N,DPL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         LEGENDRE POLYNOMIAL OF ORDER N
C
      P2=X
      PL=1.5D0*X*X-0.5D0
      DO K=2,N-1
         P1=P2
         P2=PL
         PL=(DBLE(2*K+1)*P2*X-DBLE(K)*P1)/DBLE(K+1)
      ENDDO
      DPL=DBLE(N)*(X*PL-P2)/(X*X-1.0D0)
      RETURN
      END
