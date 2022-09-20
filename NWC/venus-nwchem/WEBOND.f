      SUBROUTINE WEBOND(T,E,ESQ,NEBOX,N,NUM,NOUT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         EVALUATE AND WRITE EXCITED LOCAL MODE POPULATION ARRAYS
C
      DIMENSION T(2000),E(2000),ESQ(2000),NEBOX(2000),
     *PNT(2000),PT(2000),PESQ(2000),PE(2000),PLOG(2000)
C
      ANUM=DBLE(NUM)
      BNUM=DBLE(NUM-1)
      ANORM=DBLE(NUM)
      DO I=1,N
         PNT(I)=DBLE(NEBOX(I))/ANORM
         PT(I)=T(I)/100.0D0
         PESQ(I)=SQRT((ESQ(I)-E(I)*E(I)/ANUM)/BNUM)
         PE(I)=E(I)/ANUM
         PLOG(I)=0.0D0
         IF(PNT(I).NE.0)PLOG(I)=LOG(PNT(I))
      ENDDO
      IF (NOUT.EQ.0) IP=12
      IF (NOUT.EQ.1) IP=13
      IF (NOUT.EQ.-1) IP=14
      DO I=1,N
         WRITE(IP,9999)PT(I),PE(I),PESQ(I),PNT(I),PLOG(I)
      ENDDO
      WRITE(IP,9998) NUM
      REWIND IP
 9999 FORMAT(5(1PE15.4,','))
 9998 FORMAT(I10)
      RETURN
      END
