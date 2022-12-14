      SUBROUTINE RYDBG(NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CALCULATE RYDBERG POTENTIAL ENERGY AND DERIVATIVES
C
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/RYDBGB/RYDZ(100),DRYD(100),ARYD(100),VRYD(100),
     *N16J(100),N16K(100)
      COMMON/FORCES/NATOMS,I3N,NS,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
C
C         CALCULATE INDICES FOR COORDINATES
C
      J3=3*N16J(NL)
      J2=J3-1
      J1=J2-1
      K3=3*N16K(NL)
      K2=K3-1
      K1=K2-1
C
C        CALCULATE  INDEX FOR R
C
      JK=(N16J(NL)-1)*(2*NATOMS-N16J(NL))/2+N16K(NL)-N16J(NL)
C
C        CALCULATE RELATIVE COORDINATES AND POTENTIAL INCREMENT
C
      T1=Q(K1)-Q(J1)
      T2=Q(K2)-Q(J2)
      T3=Q(K3)-Q(J3)
      R(JK)=SQRT(T1*T1+T2*T2+T3*T3)
      ROJK=R(JK)/RYDZ(NL)-1.0D0
      DUM=-DRYD(NL)*EXP(-ARYD(NL)*ROJK)
      VRYD(NL)=DUM*(1.0D0+ARYD(NL)*ROJK)
      DUM=DUM*ARYD(NL)*ARYD(NL)*ROJK/R(JK)/RYDZ(NL)
C
C        CALCULATE (DV/DQ)'S
C
      TMP=DUM*T1
      PDOT(K1)=PDOT(K1)-TMP
      PDOT(J1)=PDOT(J1)+TMP
      TMP=DUM*T2
      PDOT(K2)=PDOT(K2)-TMP
      PDOT(J2)=PDOT(J2)+TMP
      TMP=DUM*T3
      PDOT(K3)=PDOT(K3)-TMP
      PDOT(J3)=PDOT(J3)+TMP
      RETURN
      END
