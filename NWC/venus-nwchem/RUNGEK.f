      SUBROUTINE RUNGEK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         PERFORM ONE CYCLE OF RUNGE-KUTTA-GILL INTEGRATION
C         OF THE EQUATIONS OF MOTION.
C
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/RKUTTA/RAA1,RA1,RA2,RA3,RB1,RB2,RB3,RC1,RC2
      COMMON/INTEGR/ATIME,NI,NID
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/TABLEB/TABLE(42*NDA)
      DIMENSION BQ(NDA3),BP(NDA3)
C
      DO I=1,NI
         RDUM=RAA1*TIME
         Q(I)=Q(I)+RDUM*QDOT(I)
         P(I)=P(I)+RDUM*PDOT(I)
         BQ(I)=QDOT(I)
         BP(I)=PDOT(I)
      ENDDO
      CALL PARTI
      DO I=1,NI
         RDUM=RA1*TIME
         Q(I)=Q(I)+(QDOT(I)-BQ(I))*RDUM
         P(I)=P(I)+(PDOT(I)-BP(I))*RDUM
         BQ(I)=(RA2*QDOT(I)-RA3*BQ(I))
         BP(I)=(RA2*PDOT(I)-RA3*BP(I))
      ENDDO
      CALL PARTI
      DO I=1,NI
         RDUM=RB1*TIME
         Q(I)=Q(I)+(QDOT(I)-BQ(I))*RDUM
         P(I)=P(I)+(PDOT(I)-BP(I))*RDUM
         BQ(I)=(RB2*QDOT(I)-RB3*BQ(I))
         BP(I)=(RB2*PDOT(I)-RB3*BP(I))
      ENDDO
      CALL PARTI
      DO I=1,NI
         Q(I)=Q(I)+(QDOT(I)*RC1-BQ(I)*RC2)*TIME
         P(I)=P(I)+(PDOT(I)*RC1-BP(I)*RC2)*TIME
      ENDDO
      CALL PARTI
      I=NC*NID
      DO K=1,NI
         M=I+K
         N=M+NI
         TABLE(M)=PDOT(K)
         TABLE(N)=QDOT(K)
      ENDDO
      RETURN
      END
