      SUBROUTINE PARTI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         DRIVER ROUTINE FOR CALCULATING ENERGY DERIVATIVES
C
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      J=1
      DO I=1,NATOMS
         K=J+1
         L=J+2
         QDOT(J)=P(J)/W(I)
         QDOT(K)=P(K)/W(I)
         QDOT(L)=P(L)/W(I)
         J=J+3
      ENDDO
      CALL DVDQ
      RETURN
      END
