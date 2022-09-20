      SUBROUTINE LEPS1(NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CALCULATES GENERALIZED LEPS POTENTIAL ENERGY AND DERIVATIVES
C
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/LEPSA/RLZ1(100),RLZ2(100),RLZ3(100),BL1(100),BL2(100),
     *BL3(100),DL1(100),DL2(100),DL3(100),N18J1(100),N18K1(100),
     *N18J2(100),N18K2(100),N18J3(100),N18K3(100),DELTA1(100),
     *DELTA2(100),DELTA3(100),VLEPSA(100)
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
C
C         CALCULATE INDICES FOR COORDINATES
C
      J13=3*N18J1(NL)
      J12=J13-1
      J11=J12-1
      K13=3*N18K1(NL)
      K12=K13-1
      K11=K12-1
      J23=3*N18J2(NL)
      J22=J23-1
      J21=J22-1
      K23=3*N18K2(NL)
      K22=K23-1
      K21=K22-1
      J33=3*N18J3(NL)
      J32=J33-1
      J31=J32-1
      K33=3*N18K3(NL)
      K32=K33-1
      K31=K32-1
C
C         CALCULATE INDEX FOR R
C
      JK1=(N18J1(NL)-1)*(2*NATOMS-N18J1(NL))/2+N18K1(NL)-N18J1(NL)
      JK2=(N18J2(NL)-1)*(2*NATOMS-N18J2(NL))/2+N18K2(NL)-N18J2(NL)
      JK3=(N18J3(NL)-1)*(2*NATOMS-N18J3(NL))/2+N18K3(NL)-N18J3(NL)
C
C         CALCULATE RELATIVE COORDINATES AND R
C
      T11=Q(K11)-Q(J11)
      T12=Q(K12)-Q(J12)
      T13=Q(K13)-Q(J13)
      T21=Q(K21)-Q(J21)
      T22=Q(K22)-Q(J22)
      T23=Q(K23)-Q(J23)
      T31=Q(K31)-Q(J31)
      T32=Q(K32)-Q(J32)
      T33=Q(K33)-Q(J33)
      R(JK1)=SQRT(T11*T11+T12*T12+T13*T13)
      R(JK2)=SQRT(T21*T21+T22*T22+T23*T23)
      R(JK3)=SQRT(T31*T31+T32*T32+T33*T33)
C
C         CALCULATE THE POTENTIAL
C
      TERM1=EXP(-BL1(NL)*(R(JK1)-RLZ1(NL)))
      TERM2=EXP(-BL2(NL)*(R(JK2)-RLZ2(NL)))
      TERM3=EXP(-BL3(NL)*(R(JK3)-RLZ3(NL)))
      VM1=DL1(NL)*(TERM1*TERM1-2.0D0*TERM1)
      VM2=DL2(NL)*(TERM2*TERM2-2.0D0*TERM2)
      VM3=DL3(NL)*(TERM3*TERM3-2.0D0*TERM3)
      VAM1=0.5D0*DL1(NL)*(TERM1*TERM1+2.0D0*TERM1)
      VAM2=0.5D0*DL2(NL)*(TERM2*TERM2+2.0D0*TERM2)
      VAM3=0.5D0*DL3(NL)*(TERM3*TERM3+2.0D0*TERM3)
      CM1=1.0D0+DELTA1(NL)
      CM2=1.0D0+DELTA2(NL)
      CM3=1.0D0+DELTA3(NL)
      CAM1=1.0D0-DELTA1(NL)
      CAM2=1.0D0-DELTA2(NL)
      CAM3=1.0D0-DELTA3(NL)
      VQ1=(CM1*VM1+CAM1*VAM1)*0.5D0
      VQ2=(CM2*VM2+CAM2*VAM2)*0.5D0
      VQ3=(CM3*VM3+CAM3*VAM3)*0.5D0
      VJ1=(CM1*VM1-CAM1*VAM1)*0.5D0
      VJ2=(CM2*VM2-CAM2*VAM2)*0.5D0
      VJ3=(CM3*VM3-CAM3*VAM3)*0.5D0
      TERMJ=SQRT((VJ1/CM1)**2+(VJ2/CM2)**2+(VJ3/CM3)**2
     *-VJ1*VJ2/(CM1*CM2)-VJ2*VJ3/(CM2*CM3)-VJ3*VJ1/(CM3*CM1))
C
      VLEPSA(NL)=VQ1/CM1+VQ2/CM2+VQ3/CM3-TERMJ
C
C         CALCULATE (DV/DQ)'S
C
      DVMDR1=2.0D0*DL1(NL)*(1.0D0-TERM1)*BL1(NL)*TERM1/R(JK1)
      DVMDR2=2.0D0*DL2(NL)*(1.0D0-TERM2)*BL2(NL)*TERM2/R(JK2)
      DVMDR3=2.0D0*DL3(NL)*(1.0D0-TERM3)*BL3(NL)*TERM3/R(JK3)
      DVAMDR1=-DL1(NL)*(1.0D0+TERM1)*BL1(NL)*TERM1/R(JK1)
      DVAMDR2=-DL2(NL)*(1.0D0+TERM2)*BL2(NL)*TERM2/R(JK2)
      DVAMDR3=-DL3(NL)*(1.0D0+TERM3)*BL3(NL)*TERM3/R(JK3)
      DQDR1=(CM1*DVMDR1+CAM1*DVAMDR1)*0.5D0
      DQDR2=(CM2*DVMDR2+CAM2*DVAMDR2)*0.5D0
      DQDR3=(CM3*DVMDR3+CAM3*DVAMDR3)*0.5D0
      DJDR1=(CM1*DVMDR1-CAM1*DVAMDR1)*0.5D0
      DJDR2=(CM2*DVMDR2-CAM2*DVAMDR2)*0.5D0
      DJDR3=(CM3*DVMDR3-CAM3*DVAMDR3)*0.5D0
      DVDQ1=1.0D0/CM1
      DVDQ2=1.0D0/CM2
      DVDQ3=1.0D0/CM3
      DVDJ1=-0.5D0*(2.0D0*VJ1/CM1-VJ2/CM2-VJ3/CM3)/(TERMJ*CM1)
      DVDJ2=-0.5D0*(2.0D0*VJ2/CM2-VJ3/CM3-VJ1/CM1)/(TERMJ*CM2)
      DVDJ3=-0.5D0*(2.0D0*VJ3/CM3-VJ1/CM1-VJ2/CM2)/(TERMJ*CM3)
      DVDR1=DVDQ1*DQDR1+DVDJ1*DJDR1
      DVDR2=DVDQ2*DQDR2+DVDJ2*DJDR2
      DVDR3=DVDQ3*DQDR3+DVDJ3*DJDR3
C
      DUM2=T11*DVDR1
      PDOT(K11)=PDOT(K11)+DUM2
      PDOT(J11)=PDOT(J11)-DUM2
      DUM2=T12*DVDR1
      PDOT(K12)=PDOT(K12)+DUM2
      PDOT(J12)=PDOT(J12)-DUM2
      DUM2=T13*DVDR1
      PDOT(K13)=PDOT(K13)+DUM2
      PDOT(J13)=PDOT(J13)-DUM2
C
      DUM2=T21*DVDR2
      PDOT(K21)=PDOT(K21)+DUM2
      PDOT(J21)=PDOT(J21)-DUM2
      DUM2=T22*DVDR2
      PDOT(K22)=PDOT(K22)+DUM2
      PDOT(J22)=PDOT(J22)-DUM2
      DUM2=T23*DVDR2
      PDOT(K23)=PDOT(K23)+DUM2
      PDOT(J23)=PDOT(J23)-DUM2
C
      DUM2=T31*DVDR3
      PDOT(K31)=PDOT(K31)+DUM2
      PDOT(J31)=PDOT(J31)-DUM2
      DUM2=T32*DVDR3
      PDOT(K32)=PDOT(K32)+DUM2
      PDOT(J32)=PDOT(J32)-DUM2
      DUM2=T33*DVDR3
      PDOT(K33)=PDOT(K33)+DUM2
      PDOT(J33)=PDOT(J33)-DUM2
C
      RETURN
      END
