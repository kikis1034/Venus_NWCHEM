      SUBROUTINE EIGOUT(A,N,IP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         WRITE EIGENVALUES AND EIGENVECTORS IN UNIT IP
C
      COMMON/EIGVL/EIG(NDA3yf)
      DIMENSION A(NDA3yf,NDA3yf)
  101 FORMAT(//2X,6(9X,I3))
  102 FORMAT(5X,1P6E12.4)
  103 FORMAT(I3,1X,6F12.6)
  104 FORMAT(/)
      DO I=1,N
         DO J=1,I
            D=A(I,J)
            A(I,J)=A(J,I)
            A(J,I)=D
         ENDDO
      ENDDO
      K=0
    3 L=K+1
      K=K+6
      IF (N.LT.K) K=N
      WRITE(IP,101)(I,I=L,K)
c      WRITE(IP,104)
      WRITE(IP,102)(EIG(I),I=L,K)
      WRITE(IP,104)
      DO I=1,N
         WRITE(IP,103)I,(A(I,J),J=L,K)
      ENDDO
      IF(K.LT.N) GOTO 3
      DO I=1,N
         DO J=1,I
            D=A(I,J)
            A(I,J)=A(J,I)
            A(J,I)=D
         ENDDO
      ENDDO
      RETURN
      END
