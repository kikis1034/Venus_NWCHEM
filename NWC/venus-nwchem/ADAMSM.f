      SUBROUTINE ADAMSM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         PERFORM ONE CYCLE OF SIXTH-ORDER ADAMS-MOULTON INTEGRATION
C         OF THE EQUATIONS OF MOTION.
C
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/INTEGR/ATIME,NI,NID
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
C     COMMON/PRLIST/T,V,H,TIME
      COMMON/TABLEB/TABLE(42*NDA)
C
C         CODE IS OPTIMIZED FOR VECTORIZATION
C
      DO I=1,NI
         J=I+NI
         TABLE(I)=P(I)
         TABLE(J)=Q(I)
      ENDDO
C
C         ADAMS-MOULTON PREDICTOR
C
      DO I=1,NI
         J=NID+I
         ASUM=-475.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM+2877.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM-7298.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM+9982.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM-7923.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM+4277.D0*TABLE(J)
         ASUM=ASUM*ATIME
         ASUM=ASUM+P(I)
         P(I)=ASUM
      ENDDO
      DO I=NI+1,NID
         J=NID+I
         ASUM=-475.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM+2877.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM-7298.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM+9982.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM-7923.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM+4277.D0*TABLE(J)
         ASUM=ASUM*ATIME
         L=I-NI
         ASUM=ASUM+Q(L)
         Q(L)=ASUM
      ENDDO
      CALL PARTI
C
C         ADAMS-MOULTON CORRECTOR
C
      DO I=1,NI
         J=2*NID+I
         ASUM=27.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM-173.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM+482.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM-798.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM+1427.D0*TABLE(J)
         ASUM=ASUM+475.D0*PDOT(I)
         ASUM=ASUM*ATIME+TABLE(I)
         P(I)=ASUM
      ENDDO
      DO I=NI+1,NID
         J=2*NID+I
         ASUM=27.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM-173.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM+482.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM-798.D0*TABLE(J)
         J=J+NID
         ASUM=ASUM+1427.D0*TABLE(J)
         L=I-NI
         ASUM=ASUM+475.0*QDOT(L)
         ASUM=ASUM*ATIME+TABLE(I)
         Q(L)=ASUM
      ENDDO
      CALL PARTI
      DO L=1,5
         I=NID*L
         DO K=1,NID
            M=I+K
            N=M+NID
            TABLE(M)=TABLE(N)
         ENDDO
      ENDDO
      L=NID*6
      DO I=1,NI
         J=L+I
         K=J+NI
         TABLE(J)=PDOT(I)
         TABLE(K)=QDOT(I)
      ENDDO
      RETURN
      END
