      SUBROUTINE WENMOD(T,E,ENSQ,NMA,N,NUM,NBOX,NI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         EVALUATE AND WRITE NORMAL MODE POPULATION ARRAYS
C
      DIMENSION NBOX(2000),PNT(2000),PLOG(2000),PT(2000),T(2000),
     *E(NMA,2000),ENSQ(NMA,2000),PE(NMA,2000),PENSQ(NMA,2000)
C
      ANUM=DBLE(NUM)
      BNUM=DBLE(NUM-1)
      DO I=1,NMA
         DO J=1,N
            IF (I.EQ.NI) THEN
               PNT(J)=DBLE(NBOX(J))/ANUM
               PLOG(J)=0.0D0
               IF (PNT(J).NE.0.0D0) PLOG(J)=LOG(PNT(J))
            ENDIF
            PT(J)=T(J)/100.0D0
            DUM1=(ENSQ(I,J)-E(I,J)*E(I,J)/ANUM)/BNUM
            IF (DUM1.LT.0.0D0) THEN
               WRITE(77,9997)DUM1
               DUM1=-DUM1
            ENDIF
            PENSQ(I,J)=SQRT(DUM1)
            PE(I,J)=E(I,J)/ANUM
         ENDDO
      ENDDO
      DO I=1,NMA
         WRITE(15,9990)I
         DO J=1,N
            IF (I.EQ.NI) THEN 
               WRITE(15,9996)PT(J),PE(I,J),PENSQ(I,J),PNT(J),PLOG(J)
            ELSE
               WRITE(15,9995)PT(J),PE(I,J),PENSQ(I,J)
            ENDIF
         ENDDO
      ENDDO
      WRITE(15,9991)NUM
      REWIND 15
 9991 FORMAT(I10)
 9990 FORMAT(5X,'NORMAL MODE COORDINATE ****',I4)
 9995 FORMAT((3(1PE15.4,',')))
 9996 FORMAT((5(1PE15.4,',')))
 9997 FORMAT(3X,'   DUM1=',1PE15.4)
      RETURN
      END
