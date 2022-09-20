      SUBROUTINE CENMAS(WT,QCM,VCM,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CALCULATE THE CENTER OF MASS MOMENTA AND COORDINATES
C
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,L(NDA),NAM
      DIMENSION QCM(3),VCM(3)
C
C         CENTER OF MASS COORDINATES AND MOMENTA ARE STORED IN
C         ARRAYS QQ AND PP
C
      DO I=1,3
         VCM(I)=0.0D0
         QCM(I)=0.0D0
      ENDDO
C
      DO I=1,N
         J=L(I)
         J3=3*J
         J2=J3-1
         J1=J2-1
         VCM(1)=VCM(1)+P(J1)
         VCM(2)=VCM(2)+P(J2)
         VCM(3)=VCM(3)+P(J3)
         QCM(1)=QCM(1)+W(J)*Q(J1)
         QCM(2)=QCM(2)+W(J)*Q(J2)
         QCM(3)=QCM(3)+W(J)*Q(J3)
      ENDDO
C
      DO I=1,3
         VCM(I)=VCM(I)/WT
         QCM(I)=QCM(I)/WT
      ENDDO
C
      DO I=1,N
         J=L(I)
         J3=3*J
         J2=J3-1
         J1=J2-1
         PP(J1)=P(J1)-W(J)*VCM(1)
         PP(J2)=P(J2)-W(J)*VCM(2)
         PP(J3)=P(J3)-W(J)*VCM(3)
         QQ(J1)=Q(J1)-QCM(1)
         QQ(J2)=Q(J2)-QCM(2)
         QQ(J3)=Q(J3)-QCM(3)
      ENDDO
C
      RETURN
      END
