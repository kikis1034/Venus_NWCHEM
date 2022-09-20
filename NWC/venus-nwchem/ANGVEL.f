      SUBROUTINE ANGVEL(N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         SUBTRACT OFF THE ANGULAR VELOCITY
C
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,L(NDA),NAM
C
      DO I=1,N
         J=L(I)
         J3=3*J
         J2=J3-1
         J1=J2-1
         P(J1)=P(J1)-(QQ(J3)*WY-QQ(J2)*WZ)*W(J)
         P(J2)=P(J2)-(QQ(J1)*WZ-QQ(J3)*WX)*W(J)
         P(J3)=P(J3)-(QQ(J2)*WX-QQ(J1)*WY)*W(J)
         PP(J1)=P(J1)
         PP(J2)=P(J2)
         PP(J3)=P(J3)
      ENDDO
      RETURN
      END
