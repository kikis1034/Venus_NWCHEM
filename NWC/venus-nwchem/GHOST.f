      SUBROUTINE GHOST(NG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CALCULATE GHOST PAIR INTERACTION POTENTIAL ENERGY DERIVATIVES
C
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/GHOSTB/GC1(ND08),GEX1(ND08),GEX2(ND08),N8I(ND08),N8J(ND08),
     *N8K(ND08),N8L(ND08),N8M(ND08),N8N(ND08)
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
C
C         IJK ARE THE ATOMS OF THE FIRST WATER MOLECULE
C         I AND J REPRESENT HYDROGENS AND K IS OXYGEN
C         LMN REPRESENTS ATOMS FOR THE SECOND WATER MOLECULE
C         M IS OXYGEN AND L AND NATOMS ARE HYDROGENS
C
      I3=3*N8I(NG)
      I2=I3-1
      I1=I2-1
      J3=3*N8J(NG)
      J2=J3-1
      J1=J2-1
      K3=3*N8K(NG)
      K2=K3-1
      K1=K2-1
      L3=3*N8L(NG)
      L2=L3-1
      L1=L2-1
      M3=3*N8M(NG)
      M2=M3-1
      M1=M2-1
      N3=3*N8N(NG)
      N2=N3-1
      N1=N2-1
C
C         CALCULATE INDICES FOR COORDINATES FOR NORMAL ATOMS
C
      IF (N8K(NG).LE.N8J(NG)) THEN
         JK=(N8K(NG)-1)*(2*NATOMS-N8K(NG))/2+N8J(NG)-N8K(NG)
      ELSE
         JK=(N8J(NG)-1)*(2*NATOMS-N8J(NG))/2+N8K(NG)-N8J(NG)
      ENDIF
      IF (N8K(NG).LE.N8I(NG)) THEN
         IK=(N8K(NG)-1)*(2*NATOMS-N8K(NG))/2+N8I(NG)-N8K(NG)
      ELSE
         IK=(N8I(NG)-1)*(2*NATOMS-N8I(NG))/2+N8K(NG)-N8I(NG)
      ENDIF
C
C         CALCULATE RELATIVE COORDINATES
C
      T1=Q(I1)-Q(K1)
      T2=Q(I2)-Q(K2)
      T3=Q(I3)-Q(K3)
      T4=Q(J1)-Q(K1)
      T5=Q(J2)-Q(K2)
      T6=Q(J3)-Q(K3)
      TR1=T1/R(IK)+T4/R(JK)
      TR2=T2/R(IK)+T5/R(JK)
      TR3=T3/R(IK)+T6/R(JK)
C
C         CALCULATE VECTOR ALONG BISECTOR OF IJK
C
      R1=SQRT(TR1*TR1+TR2*TR2+TR3*TR3)
C
C         CONSIDER WATER WITH ATOMS LMN
C         CALCULATE INDICES FOR COORDINATES FOR NORMAL ATOMS
C
      IF (N8M(NG).LE.N8L(NG)) THEN
         LM=(N8M(NG)-1)*(2*NATOMS-N8M(NG))/2+N8L(NG)-N8M(NG)
      ELSE
         LM=(N8L(NG)-1)*(2*NATOMS-N8L(NG))/2+N8M(NG)-N8L(NG)
      ENDIF
      IF (N8M(NG).LE.N8N(NG)) THEN
         NM=(N8M(NG)-1)*(2*NATOMS-N8M(NG))/2+N8N(NG)-N8M(NG)
      ELSE
         NM=(N8N(NG)-1)*(2*NATOMS-N8N(NG))/2+N8M(NG)-N8N(NG)
      ENDIF
C
C         CALCULATE RELATIVE COORDINATES BETWEEN LMN
C
      T7=Q(L1)-Q(M1)
      T8=Q(L2)-Q(M2)
      T9=Q(L3)-Q(M3)
      T10=Q(N1)-Q(M1)
      T11=Q(N2)-Q(M2)
      T12=Q(N3)-Q(M3)
      TR4=T7/R(LM)+T10/R(NM)
      TR5=T8/R(LM)+T11/R(NM)
      TR6=T9/R(LM)+T12/R(NM)
C
C         CALCULATE VECTOR BISECTING LMN ATOMS
C
      R2=SQRT(TR4*TR4+TR5*TR5+TR6*TR6)
C
C         CALCULATE THE COORDINATES OF THE GHOST ATOM FOR IJK ATOMS
C
      X1=Q(K1)+GC1(NG)*TR1/R1
      Y1=Q(K2)+GC1(NG)*TR2/R1
      Z1=Q(K3)+GC1(NG)*TR3/R1
C
C         CALCULATE THE COORDINATES OF THE GHOST ATOM FOR LMN ATOMS
C
      X2=Q(M1)+GC1(NG)*TR4/R2
      Y2=Q(M2)+GC1(NG)*TR5/R2
      Z2=Q(M3)+GC1(NG)*TR6/R2
C
C         X1,Y1,Z1 ARE THE COORDINATES OF THE GHOST ATOM
C         LABELLED IG AND X2,Y2,Z2 ARE THE COORDINATES
C         OF THE GHOST ATOM LABELLED JG
C         GR(1): IS THE DISTANCE BETWEEN IG AND JG
C         GR(2): IS THE DISTANCE BETWEEN ATOM I AND JG
C         GR(3): IS THE DISTANCE BETWEEN ATOM J AND JG
C         GR(4): IS THE DISTANCE BETWEEN ATOM L AND IG
C         GR(5): ISS THE DISTANCE BETWEEN ATOM N AND IG
C         IG: LABEL OF THE GHOST ATOM CONNECTED WITH IJK MOLECULE
C         JG: LABEL OF THE GHOST ATOM CONNECTED WITH LMN MOLECULE
C
C         CALCULATE THE DISTANCES REQUIRED BETWEEN THE GHOST ATOMS
C         AND THE NORMAL ATOMS AND BETWEEN GHOST ATOMS
C
      GR1X=X1-X2
      GR1Y=Y1-Y2
      GR1Z=Z1-Z2
      GR(NG,1)=SQRT(GR1X*GR1X+GR1Y*GR1Y+GR1Z*GR1Z)
      GR2X=Q(I1)-X2
      GR2Y=Q(I2)-Y2
      GR2Z=Q(I3)-Z2
      GR(NG,2)=SQRT(GR2X*GR2X+GR2Y*GR2Y+GR2Z*GR2Z)
      GR3X=Q(J1)-X2
      GR3Y=Q(J2)-Y2
      GR3Z=Q(J3)-Z2
      GR(NG,3)=SQRT(GR3X*GR3X+GR3Y*GR3Y+GR3Z*GR3Z)
      GR4X=Q(L1)-X1
      GR4Y=Q(L2)-Y1
      GR4Z=Q(L3)-Z1
      GR(NG,4)=SQRT(GR4X*GR4X+GR4Y*GR4Y+GR4Z*GR4Z)
      GR5X=Q(N1)-X1
      GR5Y=Q(N2)-Y1
      GR5Z=Q(N3)-Z1
      GR(NG,5)=SQRT(GR5X*GR5X+GR5Y*GR5Y+GR5Z*GR5Z)
      CON1=GC1(NG)/R1
      CON2=CON1/R1
      CON3=GC1(NG)/R2
      CON4=CON3/R2
      SR1Q1=-T1/R(IK)
      SR1Q2=-T2/R(IK)
      SR1Q3=-T3/R(IK)
      SR1Q4=-SR1Q1
      SR1Q5=-SR1Q2
      SR1Q6=-SR1Q3
      SR2Q1=-T4/R(JK)
      SR2Q2=-T5/R(JK)
      SR2Q3=-T6/R(JK)
      SR2Q7=-SR2Q1
      SR2Q8=-SR2Q2
      SR2Q9=-SR2Q3
      SR3Q13=-T7/R(LM)
      SR3Q14=-T8/R(LM)
      SR3Q15=-T9/R(LM)
      SR3Q16=-SR3Q13
      SR3Q17=-SR3Q14
      SR3Q18=-SR3Q15
      SR4Q13=-T10/R(NM)
      SR4Q14=-T11/R(NM)
      SR4Q15=-T12/R(NM)
      SR4Q19=-SR4Q13
      SR4Q20=-SR4Q14
      SR4Q21=-SR4Q15
      CON5=-1.D0/R(IK)-1.D0/R(JK)
      CON6=-1.D0/R(LM)-1.D0/R(NM)
      R1Q011=CON5-SR1Q4*SR1Q1/R(IK)-SR2Q7*SR2Q1/R(JK)
      R1Q012=-SR1Q5*SR1Q1/R(IK)-SR2Q8*SR2Q1/R(JK)
      R1Q013=-SR1Q6*SR1Q1/R(IK)-SR2Q9*SR2Q1/R(JK)
      R1Q01=(R1Q011*TR1+R1Q012*TR2+R1Q013*TR3)/R1
      R1Q021=-SR1Q4*SR1Q2/R(IK)-SR2Q7*SR2Q2/R(JK)
      R1Q022=CON5-SR1Q5*SR1Q2/R(IK)-SR2Q8*SR2Q2/R(JK)
      R1Q023=-SR1Q6*SR1Q2/R(IK)-SR2Q9*SR2Q2/R(JK)
      R1Q02=(R1Q021*TR1+R1Q022*TR2+R1Q023*TR3)/R1
      R1Q031=-SR1Q4*SR1Q3/R(IK)-SR2Q7*SR2Q3/R(JK)
      R1Q032=-SR1Q5*SR1Q3/R(IK)-SR2Q8*SR2Q3/R(JK)
      R1Q033=CON5-SR1Q6*SR1Q3/R(IK)-SR2Q9*SR2Q3/R(JK)
      R1Q03=(R1Q031*TR1+R1Q032*TR2+R1Q033*TR3)/R1
      R1Q041=1.D0/R(IK)-SR1Q4*SR1Q4/R(IK)
      R1Q042=-SR1Q5*SR1Q4/R(IK)
      R1Q043=-SR1Q6*SR1Q4/R(IK)
      R1Q04=(R1Q041*TR1+R1Q042*TR2+R1Q043*TR3)/R1
      R1Q051=-SR1Q4*SR1Q5/R(IK)
      R1Q052=1.D0/R(IK)-SR1Q5*SR1Q5/R(IK)
      R1Q053=-SR1Q6*SR1Q5/R(IK)
      R1Q05=(R1Q051*TR1+R1Q052*TR2+R1Q053*TR3)/R1
      R1Q061=-SR1Q4*SR1Q6/R(IK)
      R1Q062=-SR1Q5*SR1Q6/R(IK)
      R1Q063=1.D0/R(IK)-SR1Q6*SR1Q6/R(IK)
      R1Q06=(R1Q061*TR1+R1Q062*TR2+R1Q063*TR3)/R1
      R1Q071=1.D0/R(JK)-SR2Q7*SR2Q7/R(JK)
      R1Q072=-SR2Q8*SR2Q7/R(JK)
      R1Q073=-SR2Q9*SR2Q7/R(JK)
      R1Q07=(R1Q071*TR1+R1Q072*TR2+R1Q073*TR3)/R1
      R1Q081=-SR2Q7*SR2Q8/R(JK)
      R1Q082=1.D0/R(JK)-SR2Q8*SR2Q8/R(JK)
      R1Q083=-SR2Q9*SR2Q8/R(JK)
      R1Q08=(R1Q081*TR1+R1Q082*TR2+R1Q083*TR3)/R1
      R1Q091=-SR2Q7*SR2Q9/R(JK)
      R1Q092=-SR2Q8*SR2Q9/R(JK)
      R1Q093=1.D0/R(JK)-SR2Q9*SR2Q9/R(JK)
      R1Q09=(R1Q091*TR1+R1Q092*TR2+R1Q093*TR3)/R1
      R2Q131=CON6-SR3Q16*SR3Q13/R(LM)-SR4Q19*SR4Q13/R(NM)
      R2Q132=-SR3Q17*SR3Q13/R(LM)-SR4Q20*SR4Q13/R(NM)
      R2Q133=-SR3Q18*SR3Q13/R(LM)-SR4Q21*SR4Q13/R(NM)
      R2Q13=(R2Q131*TR4+R2Q132*TR5+R2Q133*TR6)/R2
      R2Q141=-SR3Q16*SR3Q14/R(LM)-SR4Q19*SR4Q14/R(NM)
      R2Q142=CON6-SR3Q17*SR3Q14/R(LM)-SR4Q20*SR4Q14/R(NM)
      R2Q143=-SR3Q18*SR3Q14/R(LM)-SR4Q21*SR4Q14/R(NM)
      R2Q14=(R2Q141*TR4+R2Q142*TR5+R2Q143*TR6)/R2
      R2Q151=-SR3Q16*SR3Q15/R(LM)-SR4Q19*SR4Q15/R(NM)
      R2Q152=-SR3Q17*SR3Q15/R(LM)-SR4Q20*SR4Q15/R(NM)
      R2Q153=CON6-SR3Q18*SR3Q15/R(LM)-SR4Q21*SR4Q15/R(NM)
      R2Q15=(R2Q151*TR4+R2Q152*TR5+R2Q153*TR6)/R2
      R2Q161=1.D0/R(LM)-SR3Q16*SR3Q16/R(LM)
      R2Q162=-SR3Q17*SR3Q16/R(LM)
      R2Q163=-SR3Q18*SR3Q16/R(LM)
      R2Q16=(R2Q161*TR4+R2Q162*TR5+R2Q163*TR6)/R2
      R2Q171=-SR3Q16*SR3Q17/R(LM)
      R2Q172=1.D0/R(LM)-SR3Q17*SR3Q17/R(LM)
      R2Q173=-SR3Q18*SR3Q17/R(LM)
      R2Q17=(R2Q171*TR4+R2Q172*TR5+R2Q173*TR6)/R2
      R2Q181=-SR3Q16*SR3Q18/R(LM)
      R2Q182=-SR3Q17*SR3Q18/R(LM)
      R2Q183=1.D0/R(LM)-SR3Q18*SR3Q18/R(LM)
      R2Q18=(R2Q181*TR4+R2Q182*TR5+R2Q183*TR6)/R2
      R2Q191=1.D0/R(NM)-SR4Q19*SR4Q19/R(NM)
      R2Q192=-SR4Q20*SR4Q19/R(NM)
      R2Q193=-SR4Q21*SR4Q19/R(NM)
      R2Q19=(R2Q191*TR4+R2Q192*TR5+R2Q193*TR6)/R2
      R2Q201=-SR4Q19*SR4Q20/R(NM)
      R2Q202=1.D0/R(NM)-SR4Q20*SR4Q20/R(NM)
      R2Q203=-SR4Q21*SR4Q20/R(NM)
      R2Q20=(R2Q201*TR4+R2Q202*TR5+R2Q203*TR6)/R2
      R2Q211=-SR4Q19*SR4Q21/R(NM)
      R2Q212=-SR4Q20*SR4Q21/R(NM)
      R2Q213=-SR4Q21*SR4Q21/R(NM)+1.D0/R(NM)
      R2Q21=(R2Q211*TR4+R2Q212*TR5+R2Q213*TR6)/R2
      X8Q01=1.D0-CON2*R1Q01*TR1+R1Q011*CON1
      X8Q02=-CON2*R1Q02*TR1+CON1*R1Q021
      X8Q03=-CON2*R1Q03*TR1+CON1*R1Q031
      X8Q04=-CON2*R1Q04*TR1+CON1*R1Q041
      X8Q05=-CON2*R1Q05*TR1+CON1*R1Q051
      X8Q06=-CON2*R1Q06*TR1+CON1*R1Q061
      X8Q07=-CON2*R1Q07*TR1+CON1*R1Q071
      X8Q08=-CON2*R1Q08*TR1+CON1*R1Q081
      X8Q09=-CON2*R1Q09*TR1+CON1*R1Q091
      Y8Q01=-CON2*R1Q01*TR2+CON1*R1Q012
      Y8Q02=-CON2*R1Q02*TR2+CON1*R1Q022+1.D0
      Y8Q03=-CON2*R1Q03*TR2+CON1*R1Q032
      Y8Q04=-CON2*R1Q04*TR2+CON1*R1Q042
      Y8Q05=-CON2*R1Q05*TR2+CON1*R1Q052
      Y8Q06=-CON2*R1Q06*TR2+CON1*R1Q062
      Y8Q07=-CON2*R1Q07*TR2+CON1*R1Q072
      Y8Q08=-CON2*R1Q08*TR2+CON1*R1Q082
      Y8Q09=-CON2*R1Q09*TR2+CON1*R1Q092
      Z8Q01=-CON2*R1Q01*TR3+CON1*R1Q013
      Z8Q02=-CON2*R1Q02*TR3+CON1*R1Q023
      Z8Q03=-CON2*R1Q03*TR3+CON1*R1Q033+1.D0
      Z8Q04=-CON2*R1Q04*TR3+CON1*R1Q043
      Z8Q05=-CON2*R1Q05*TR3+CON1*R1Q053
      Z8Q06=-CON2*R1Q06*TR3+CON1*R1Q063
      Z8Q07=-CON2*R1Q07*TR3+CON1*R1Q073
      Z8Q08=-CON2*R1Q08*TR3+CON1*R1Q083
      Z8Q09=-CON2*R1Q09*TR3+CON1*R1Q093
      X9Q13=1.D0-CON4*R2Q13*TR4+CON3*R2Q131
      X9Q14=-CON4*R2Q14*TR4+CON3*R2Q141
      X9Q15=-CON4*R2Q15*TR4+CON3*R2Q151
      X9Q16=-CON4*R2Q16*TR4+CON3*R2Q161
      X9Q17=-CON4*R2Q17*TR4+CON3*R2Q171
      X9Q18=-CON4*R2Q18*TR4+CON3*R2Q181
      X9Q19=-CON4*R2Q19*TR4+CON3*R2Q191
      X9Q20=-CON4*R2Q20*TR4+CON3*R2Q201
      X9Q21=-CON4*R2Q21*TR4+CON3*R2Q211
      Y9Q13=-CON4*R2Q13*TR5+CON3*R2Q132
      Y9Q14=-CON4*R2Q14*TR5+CON3*R2Q142+1.D0
      Y9Q15=-CON4*R2Q15*TR5+CON3*R2Q152
      Y9Q16=-CON4*R2Q16*TR5+CON3*R2Q162
      Y9Q17=-CON4*R2Q17*TR5+CON3*R2Q172
      Y9Q18=-CON4*R2Q18*TR5+CON3*R2Q182
      Y9Q19=-CON4*R2Q19*TR5+CON3*R2Q192
      Y9Q20=-CON4*R2Q20*TR5+CON3*R2Q202
      Y9Q21=-CON4*R2Q21*TR5+CON3*R2Q212
      Z9Q13=-CON4*R2Q13*TR6+CON3*R2Q133
      Z9Q14=-CON4*R2Q14*TR6+CON3*R2Q143
      Z9Q15=-CON4*R2Q15*TR6+CON3*R2Q153+1.D0
      Z9Q16=-CON4*R2Q16*TR6+CON3*R2Q163
      Z9Q17=-CON4*R2Q17*TR6+CON3*R2Q173
      Z9Q18=-CON4*R2Q18*TR6+CON3*R2Q183
      Z9Q19=-CON4*R2Q19*TR6+CON3*R2Q193
      Z9Q20=-CON4*R2Q20*TR6+CON3*R2Q203
      Z9Q21=-CON4*R2Q21*TR6+CON3*R2Q213
      R89Q01=(GR1X*X8Q01+GR1Y*Y8Q01+GR1Z*Z8Q01)/GR(NG,1)
      R89Q02=(GR1X*X8Q02+GR1Y*Y8Q02+GR1Z*Z8Q02)/GR(NG,1)
      R89Q03=(GR1X*X8Q03+GR1Y*Y8Q03+GR1Z*Z8Q03)/GR(NG,1)
      R89Q04=(GR1X*X8Q04+GR1Y*Y8Q04+GR1Z*Z8Q04)/GR(NG,1)
      R89Q05=(GR1X*X8Q05+GR1Y*Y8Q05+GR1Z*Z8Q05)/GR(NG,1)
      R89Q06=(GR1X*X8Q06+GR1Y*Y8Q06+GR1Z*Z8Q06)/GR(NG,1)
      R89Q07=(GR1X*X8Q07+GR1Y*Y8Q07+GR1Z*Z8Q07)/GR(NG,1)
      R89Q08=(GR1X*X8Q08+GR1Y*Y8Q08+GR1Z*Z8Q08)/GR(NG,1)
      R89Q09=(GR1X*X8Q09+GR1Y*Y8Q09+GR1Z*Z8Q09)/GR(NG,1)
      R89Q13=(-GR1X*X9Q13-GR1Y*Y9Q13-GR1Z*Z9Q13)/GR(NG,1)
      R89Q14=(-GR1X*X9Q14-GR1Y*Y9Q14-GR1Z*Z9Q14)/GR(NG,1)
      R89Q15=(-GR1X*X9Q15-GR1Y*Y9Q15-GR1Z*Z9Q15)/GR(NG,1)
      R89Q16=(-GR1X*X9Q16-GR1Y*Y9Q16-GR1Z*Z9Q16)/GR(NG,1)
      R89Q17=(-GR1X*X9Q17-GR1Y*Y9Q17-GR1Z*Z9Q17)/GR(NG,1)
      R89Q18=(-GR1X*X9Q18-GR1Y*Y9Q18-GR1Z*Z9Q18)/GR(NG,1)
      R89Q19=(-GR1X*X9Q19-GR1Y*Y9Q19-GR1Z*Z9Q19)/GR(NG,1)
      R89Q20=(-GR1X*X9Q20-GR1Y*Y9Q20-GR1Z*Z9Q20)/GR(NG,1)
      R89Q21=(-GR1X*X9Q21-GR1Y*Y9Q21-GR1Z*Z9Q21)/GR(NG,1)
      R29Q04=GR2X/GR(NG,2)
      R29Q05=GR2Y/GR(NG,2)
      R29Q06=GR2Z/GR(NG,2)
      R29Q13=(-GR2X*X9Q13-GR2Y*Y9Q13-GR2Z*Z9Q13)/GR(NG,2)
      R29Q14=(-GR2X*X9Q14-GR2Y*Y9Q14-GR2Z*Z9Q14)/GR(NG,2)
      R29Q15=(-GR2X*X9Q15-GR2Y*Y9Q15-GR2Z*Z9Q15)/GR(NG,2)
      R29Q16=(-GR2X*X9Q16-GR2Y*Y9Q16-GR2Z*Z9Q16)/GR(NG,2)
      R29Q17=(-GR2X*X9Q17-GR2Y*Y9Q17-GR2Z*Z9Q17)/GR(NG,2)
      R29Q18=(-GR2X*X9Q18-GR2Y*Y9Q18-GR2Z*Z9Q18)/GR(NG,2)
      R29Q19=(-GR2X*X9Q19-GR2Y*Y9Q19-GR2Z*Z9Q19)/GR(NG,2)
      R29Q20=(-GR2X*X9Q20-GR2Y*Y9Q20-GR2Z*Z9Q20)/GR(NG,2)
      R29Q21=(-GR2X*X9Q21-GR2Y*Y9Q21-GR2Z*Z9Q21)/GR(NG,2)
      R39Q07=GR3X/GR(NG,3)
      R39Q08=GR3Y/GR(NG,3)
      R39Q09=GR3Z/GR(NG,3)
      R39Q13=(-GR3X*X9Q13-GR3Y*Y9Q13-GR3Z*Z9Q13)/GR(NG,3)
      R39Q14=(-GR3X*X9Q14-GR3Y*Y9Q14-GR3Z*Z9Q14)/GR(NG,3)
      R39Q15=(-GR3X*X9Q15-GR3Y*Y9Q15-GR3Z*Z9Q15)/GR(NG,3)
      R39Q16=(-GR3X*X9Q16-GR3Y*Y9Q16-GR3Z*Z9Q16)/GR(NG,3)
      R39Q17=(-GR3X*X9Q17-GR3Y*Y9Q17-GR3Z*Z9Q17)/GR(NG,3)
      R39Q18=(-GR3X*X9Q18-GR3Y*Y9Q18-GR3Z*Z9Q18)/GR(NG,3)
      R39Q19=(-GR3X*X9Q19-GR3Y*Y9Q19-GR3Z*Z9Q19)/GR(NG,3)
      R39Q20=(-GR3X*X9Q20-GR3Y*Y9Q20-GR3Z*Z9Q20)/GR(NG,3)
      R39Q21=(-GR3X*X9Q21-GR3Y*Y9Q21-GR3Z*Z9Q21)/GR(NG,3)
      R68Q01=(-GR4X*X8Q01-GR4Y*Y8Q01-GR4Z*Z8Q01)/GR(NG,4)
      R68Q02=(-GR4X*X8Q02-GR4Y*Y8Q02-GR4Z*Z8Q02)/GR(NG,4)
      R68Q03=(-GR4X*X8Q03-GR4Y*Y8Q03-GR4Z*Z8Q03)/GR(NG,4)
      R68Q04=(-GR4X*X8Q04-GR4Y*Y8Q04-GR4Z*Z8Q04)/GR(NG,4)
      R68Q05=(-GR4X*X8Q05-GR4Y*Y8Q05-GR4Z*Z8Q05)/GR(NG,4)
      R68Q06=(-GR4X*X8Q06-GR4Y*Y8Q06-GR4Z*Z8Q06)/GR(NG,4)
      R68Q07=(-GR4X*X8Q07-GR4Y*Y8Q07-GR4Z*Z8Q07)/GR(NG,4)
      R68Q08=(-GR4X*X8Q08-GR4Y*Y8Q08-GR4Z*Z8Q08)/GR(NG,4)
      R68Q09=(-GR4X*X8Q09-GR4Y*Y8Q09-GR4Z*Z8Q09)/GR(NG,4)
      R68Q16=GR4X/GR(NG,4)
      R68Q17=GR4Y/GR(NG,4)
      R68Q18=GR4Z/GR(NG,4)
      R78Q01=(-GR5X*X8Q01-GR5Y*Y8Q01-GR5Z*Z8Q01)/GR(NG,5)
      R78Q02=(-GR5X*X8Q02-GR5Y*Y8Q02-GR5Z*Z8Q02)/GR(NG,5)
      R78Q03=(-GR5X*X8Q03-GR5Y*Y8Q03-GR5Z*Z8Q03)/GR(NG,5)
      R78Q04=(-GR5X*X8Q04-GR5Y*Y8Q04-GR5Z*Z8Q04)/GR(NG,5)
      R78Q05=(-GR5X*X8Q05-GR5Y*Y8Q05-GR5Z*Z8Q05)/GR(NG,5)
      R78Q06=(-GR5X*X8Q06-GR5Y*Y8Q06-GR5Z*Z8Q06)/GR(NG,5)
      R78Q07=(-GR5X*X8Q07-GR5Y*Y8Q07-GR5Z*Z8Q07)/GR(NG,5)
      R78Q08=(-GR5X*X8Q08-GR5Y*Y8Q08-GR5Z*Z8Q08)/GR(NG,5)
      R78Q09=(-GR5X*X8Q09-GR5Y*Y8Q09-GR5Z*Z8Q09)/GR(NG,5)
      R78Q19=GR5X/GR(NG,5)
      R78Q20=GR5Y/GR(NG,5)
      R78Q21=GR5Z/GR(NG,5)
      DUM1=-GEX1(NG)/GR(NG,1)/GR(NG,1)
      DUM2=-GEX2(NG)/GR(NG,2)/GR(NG,2)
      DUM3=-GEX2(NG)/GR(NG,3)/GR(NG,3)
      DUM4=-GEX2(NG)/GR(NG,4)/GR(NG,4)
      DUM5=-GEX2(NG)/GR(NG,5)/GR(NG,5)
      PDOT(K1)=PDOT(K1)+DUM1*R89Q01+DUM4*R68Q01+DUM5*R78Q01
      PDOT(K2)=PDOT(K2)+DUM1*R89Q02+DUM4*R68Q02+DUM5*R78Q02
      PDOT(K3)=PDOT(K3)+DUM1*R89Q03+DUM4*R68Q03+DUM5*R78Q03
      PDOT(I1)=PDOT(I1)+DUM1*R89Q04+DUM4*R68Q04+DUM5*R78Q04
     *+DUM2*R29Q04
      PDOT(I2)=PDOT(I2)+DUM1*R89Q05+DUM4*R68Q05+DUM5*R78Q05
     *+DUM2*R29Q05
      PDOT(I3)=PDOT(I3)+DUM1*R89Q06+DUM4*R68Q06+DUM5*R78Q06
     *+DUM2*R29Q06
      PDOT(J1)=PDOT(J1)+DUM1*R89Q07+DUM4*R68Q07+DUM5*R78Q07
     *+DUM3*R39Q07
      PDOT(J2)=PDOT(J2)+DUM1*R89Q08+DUM4*R68Q08+DUM5*R78Q08
     *+DUM3*R39Q08
      PDOT(J3)=PDOT(J3)+DUM1*R89Q09+DUM4*R68Q09+DUM5*R78Q09
     *+DUM3*R39Q09
      PDOT(M1)=PDOT(M1)+DUM1*R89Q13+DUM2*R29Q13+DUM3*R39Q13
      PDOT(M2)=PDOT(M2)+DUM1*R89Q14+DUM2*R29Q14+DUM3*R39Q14
      PDOT(M3)=PDOT(M3)+DUM1*R89Q15+DUM2*R29Q15+DUM3*R39Q15
      PDOT(L1)=PDOT(L1)+DUM1*R89Q16+DUM2*R29Q16+DUM3*R39Q16
     *+DUM4*R68Q16
      PDOT(L2)=PDOT(L2)+DUM1*R89Q17+DUM2*R29Q17+DUM3*R39Q17
     *+DUM4*R68Q17
      PDOT(L3)=PDOT(L3)+DUM1*R89Q18+DUM2*R29Q18+DUM3*R39Q18
     *+DUM4*R68Q18
      PDOT(N1)=PDOT(N1)+DUM1*R89Q19+DUM2*R29Q19+DUM3*R39Q19
     *+DUM5*R78Q19
      PDOT(N2)=PDOT(N2)+DUM1*R89Q20+DUM2*R29Q20+DUM3*R39Q20
     *+DUM5*R78Q20
      PDOT(N3)=PDOT(N3)+DUM1*R89Q21+DUM2*R29Q21+DUM3*R39Q21
     *+DUM5*R78Q21
      RETURN
      END