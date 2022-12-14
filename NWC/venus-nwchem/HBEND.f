      SUBROUTINE HBEND(NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CALCULATE HARMONIC BEND POTENTIAL ENERGY DERIVATIVES
C
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/BENDB/THETAZ(ND03),FBZ(ND03),CJ(ND03),CK(ND03),RJZ(ND03),
     *RKZ(ND03),FB(ND03),N3J(ND03),N3K(ND03),N3M(ND03)
      COMMON/FORCES/N,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
C
C         CALCULATE INDICES FOR COORDINATES
C
      J3=3*N3J(NL)
      J2=J3-1
      J1=J2-1
      M3=3*N3M(NL)
      M2=M3-1
      M1=M2-1
      K3=3*N3K(NL)
      K2=K3-1
      K1=K2-1
C
C         CALCULATE INDICES FOR R'S
C
      IF (N3M(NL).LE.N3K(NL)) THEN
         KM=(N3M(NL)-1)*(2*N-N3M(NL))/2+N3K(NL)-N3M(NL)
      ELSE
         KM=(N3K(NL)-1)*(2*N-N3K(NL))/2+N3M(NL)-N3K(NL)
      ENDIF
      IF (N3M(NL).LE.N3J(NL)) THEN
         JM=(N3M(NL)-1)*(2*N-N3M(NL))/2+N3J(NL)-N3M(NL)
      ELSE
         JM=(N3J(NL)-1)*(2*N-N3J(NL))/2+N3M(NL)-N3J(NL)
      ENDIF
C
C         CALCULATE RELATIVE COORDINATES
C
      T1=Q(J1)-Q(M1)
      T2=Q(J2)-Q(M2)
      T3=Q(J3)-Q(M3)
      T4=Q(K1)-Q(M1)
      T5=Q(K2)-Q(M2)
      T6=Q(K3)-Q(M3)
C
C         CALCULATE ANGLE
C
      CTHETA=(T1*T4+T2*T5+T3*T6)/R(JM)/R(KM)
      IF (CTHETA.GT.1.00D0) CTHETA=1.00D0
      IF (CTHETA.LT.-1.00D0) CTHETA=-1.00D0
      THETA(NL)=ACOS(CTHETA)
C
C         CALCULATE SWITCHING FUNCTIONS AND FORCE CONSTANTS
C
      SRJ=0.0D0
      DUM1=R(JM)-RJZ(NL)
      DUM2=CJ(NL)*DUM1**2
      IF (DUM2.LE.85.0D0) THEN
         SRJ=1.0D0
         IF (CJ(NL).GT.0.0D0.AND.DUM1.GT.0.0D0) SRJ=EXP(-DUM2)
      ENDIF
      SRK=0.0D0
      DUM1=R(KM)-RKZ(NL)
      DUM2=CK(NL)*DUM1**2
      IF (DUM2.LE.85.0D0) THEN
         SRK=1.0D0
         IF (CK(NL).GT.0.0D0.AND.DUM1.GT.0.0D0) SRK=EXP(-DUM2)
      ENDIF
      FB(NL)=SRJ*SRK*FBZ(NL)
C
C         CALCULATE (DV/DQ)'S
C
      DUM3=FB(NL)*(THETA(NL)-THETAZ(NL))
      IF (ABS(THETA(NL)-PI).GT.0.1D0) THEN
         DUM1=-1.0D0/SQRT(1.0D0-CTHETA**2)*DUM3
      ELSE
         DUM2=THETA(NL)-PI
         DUM1=FB(NL)/(1.D0-DUM2**2/6.D0+DUM2**4/120.D0-DUM2**6/5040.D0)
      ENDIF
      RJK=R(JM)*R(KM)
      RJ2=R(JM)**2
      RK2=R(KM)**2
      DUM2=CTHETA/RJ2
      PDOT(J1)=PDOT(J1)+DUM1*(T4/RJK-T1*DUM2)
      PDOT(J2)=PDOT(J2)+DUM1*(T5/RJK-T2*DUM2)
      PDOT(J3)=PDOT(J3)+DUM1*(T6/RJK-T3*DUM2)
      DUM2=CTHETA/RK2
      PDOT(K1)=PDOT(K1)+DUM1*(T1/RJK-T4*DUM2)
      PDOT(K2)=PDOT(K2)+DUM1*(T2/RJK-T5*DUM2)
      PDOT(K3)=PDOT(K3)+DUM1*(T3/RJK-T6*DUM2)
      PDOT(M1)=PDOT(M1)+DUM1*((-T1-T4)/RJK+CTHETA*(T1/RJ2+T4/RK2))
      PDOT(M2)=PDOT(M2)+DUM1*((-T2-T5)/RJK+CTHETA*(T2/RJ2+T5/RK2))
      PDOT(M3)=PDOT(M3)+DUM1*((-T3-T6)/RJK+CTHETA*(T3/RJ2+T6/RK2))
C
      DUM3=DUM3*(THETA(NL)-THETAZ(NL))
      DUM1=R(JM)-RJZ(NL)
      IF (CJ(NL).GE.0.0D0.AND.DUM1.GE.0.0D0) THEN
         DUM1=-DUM3*CJ(NL)*DUM1/R(JM)
         DUM2=T1*DUM1
         PDOT(J1)=PDOT(J1)+DUM2
         PDOT(M1)=PDOT(M1)-DUM2
         DUM2=T2*DUM1
         PDOT(J2)=PDOT(J2)+DUM2
         PDOT(M2)=PDOT(M2)-DUM2
         DUM2=T3*DUM1
         PDOT(J3)=PDOT(J3)+DUM2
         PDOT(M3)=PDOT(M3)-DUM2
      ENDIF
C
      DUM1=R(KM)-RKZ(NL)
      IF (CK(NL).GE.0.0D0.AND.DUM1.GE.0.0D0) THEN
         DUM1=-DUM3*CK(NL)*DUM1/R(KM)
         DUM2=T4*DUM1
         PDOT(K1)=PDOT(K1)+DUM2
         PDOT(M1)=PDOT(M1)-DUM2
         DUM2=T5*DUM1
         PDOT(K2)=PDOT(K2)+DUM2
         PDOT(M2)=PDOT(M2)-DUM2
         DUM2=T6*DUM1
         PDOT(K3)=PDOT(K3)+DUM2
         PDOT(M3)=PDOT(M3)-DUM2
      ENDIF
      RETURN
      END
