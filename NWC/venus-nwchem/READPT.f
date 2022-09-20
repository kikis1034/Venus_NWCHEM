      SUBROUTINE READPT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C        READ IN PARAMETERS FOR ALL KINDS OF ANALYTICAL POTENTIALS
C
C      PARAMETER(ND1=2000)
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/STRETB/RSZ(ND01),FS(ND01),N1J(ND01),N1K(ND01)
      COMMON/MORSEB/RMZ(ND02),B(ND02),D(ND02),N2J(ND02),N2K(ND02),
     *CM1(ND02),CM2(ND02),CM3(ND02),CM4(ND02),CM5(ND02)
      COMMON/CRCO6B/a1,a2,a3,a4,a5,a6,a7,a8,dcrco6,rtmp(3*NDA),
     *dfixed(100),deltad
      COMMON/BENDB/THETAZ(ND03),FBZ(ND03),CJ(ND03),CK(ND03),RJZ(ND03),
     *RKZ(ND03),FB(ND03),N3J(ND03),N3K(ND03),N3M(ND03)
      COMMON/ALPHAB/FA(ND04),N4J(ND04),N4K(ND04),N4M(ND04),N4N(ND04)
      COMMON/LENJB/ALJ(ND05),BLJ(ND05),CLJ(ND05),N5J(ND05),N5K(ND05),
     *NREP(ND05),MREP(ND05),LREP(ND05)
      COMMON/TAUB/VZTAU(ND06),N6I(ND06),N6J(ND06),N6K(ND06),N6L(ND06),
     *N6M(ND06),N6N(ND06)
      COMMON/EXPB/AEX(ND07),BEX(ND07),CEX(ND07),N7J(ND07),N7K(ND07),
     *RNPOW(ND07)
      COMMON/GHOSTB/GC1(ND08),GEX1(ND08),GEX2(ND08),N8I(ND08),N8J(ND08),
     *N8K(ND08),N8L(ND08),N8M(ND08),N8N(ND08)
      COMMON/TETRAB/N9I(20),N9J(20),N9K(20),N9L(20),N9M(20),
     *               FT0(20,6),FT2(20,6),GT0(20,6),GT2(20,6),
     *               HT0(20,6),HT2(20,6),THT(20,6),R0(20,4),
     *               THT1(20,6),THT2(20,6),FD1(20,4),
     *               HD1(20,4),GN0(20,5),FT(20,6),GT(20,6),
     *               HT(20,6),FD(20,4),HD(20,4),DLTA(20,48),
     *               TETTST,SGN1,SGN2,SGN3,SGN4
      COMMON/VRRB/FKRRZ(ND10),FKRR(ND10),CIJ(ND10),CKL(ND10),
     *            RIJ0(ND10),RKL0(ND10),N10I(ND10),N10J(ND10),
     *            N10K(ND10),N10L(ND10)
      COMMON/VRTB/FKRTZ(ND11),FKRT(ND11),CRT(ND11),R110(ND11),
     *            N11I(ND11),N11J(ND11),N11B(ND11),NRT(ND11)
      COMMON/VTTB/FKTTZ(ND12),FKTT(ND12),N12B(ND12),N12BB(ND12),
     *            NTT(ND12)
      COMMON/ANGLEB/FDH(ND13I,ND13J),GDH(ND13I,ND13J),NDH(ND13I),
     *N13I(ND13I),N13J(ND13I),N13K(ND13I),N13L(ND13I)
      COMMON/AXTB/ZAXT(300),VAXT(300),N14I(300),N14J(300),N14K(300)
      COMMON/RYDBGB/RYDZ(100),DRYD(100),ARYD(100),VRYD(100),
     *N16J(100),N16K(100)
      COMMON/HFDB/AHFD(100),BHFD(100),RHFD(100),VHFD(100),C6HFD(100),
     *C8HFD(100),C10HFD(100),N17J(100),N17K(100)
      COMMON/LEPSA/RLZ1(100),RLZ2(100),RLZ3(100),BL1(100),BL2(100),
     *BL3(100),DL1(100),DL2(100),DL3(100),N18J1(100),N18K1(100),
     *N18J2(100),N18K2(100),N18J3(100),N18K3(100),DELTA1(100),
     *DELTA2(100),DELTA3(100),VLEPSA(100)
      COMMON/LEPSB/RLZS1(100),RLZS2(100),RLZS3(100),RLZT1(100),
     *RLZT2(100),RLZT3(100),BLS1(100),BLS2(100),BLS3(100),BLT1(100),
     *BLT2(100),BLT3(100),DLS1(100),DLS2(100),DLS3(100),DLT1(100),
     *DLT2(100),DLT3(100),VLEPSB(100),N19J1(100),N19K1(100),
     *N19J2(100),N19K2(100),N19J3(100),N19K3(100)
      COMMON/RELAXA/N21I(ND21),N21J(ND21),N21K(ND21),N21L(ND21),
     *N21M(ND21),N21N(ND21),N21O(ND21),N21P(ND21),N21Q(ND21),
     *N21R(ND21),N21S(ND21),N21T(ND21),N21U(ND21),N21V(ND21),
     *PH0D(3,ND21),PH0R(3,ND21),F0D(3,ND21),TH0D(3,ND21),TH0R(3,ND21),
     *FCCC(3,ND21),CC0D(3,ND21),CC0R(3,ND21),BET(3,ND21),DIS(3,ND21),
     *GA0D(6,ND21),GA0R(6,ND21),FCC1(6,ND21),FCC2(3,ND21)
      COMMON/NNB/N22J(ND22),N22K(ND22),ASR(ND22),BSR(ND22),CSR(ND22),
     *NDSR(ND22),ALR(ND22),BLR(ND22),CLR(ND22),NDLR(ND22),AAS(ND22),
     *AA0(ND22),NAA(ND22)
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
C
  801 FORMAT(5X,'  RZ=',F9.3,'  A=',G9.3,'  ALPHA=',F9.3/
     *       5X,'  C6=',G10.4,'  C8=',G10.4,'  C10=',G11.4)
  803 FORMAT('   NJ=',I4,'  NK=',I4,'  RZ=',F6.3,'  F=',F9.3)
  804 FORMAT('   NJ=',I4,'  NK=',I4,'  RZ=',F6.3,'  B=',F6.3,'  D=',F8.3
     *)
  805 FORMAT('   NJ=',I4,'  NK=',I4,'  NM=',I4,'  THETAZ=',F8.3,'  FZ=',
     *F9.3,'  RJZ=',F6.3,'  CJ=',F6.3,'  RKZ=',F6.3,'  CK=',F6.3)
  806 FORMAT(7X,'HARTREE-FOCK DISPERSION INTERACTIONS'/
     *7X,'BETWEEN ALL ATOMS FROM ',I4,' TO ',I4/
     *7X,'HAVE THE FOLLOWING PARAMETERS :'/)
  807 FORMAT(7X,'HARTREE-FOCK DISPERSION INTERACTIONS'/
     *7X,'BETWEEN ATOM ',I4,' AND ATOMS ',I4,' TO ',I4/
     *7X,'HAVE THE FOLLOWING PARAMETERS :'/)
  808 FORMAT('   NJ=',I4,'  NK=',I4,'  RZ=',F9.3,'  A=',G9.3,'  ALPHA=',
     *F9.3/19X,'  C6=',G10.4,'  C8=',G10.4,'  C10=',G11.4)
  815 FORMAT(/)
  830 FORMAT('   NJ=',I4,'  NK=',I4,'  NM=',I4,'  NN=',I4,'  F=',F7.3)
  831 FORMAT('   NJ=',I4,'  NK=',I4,'  ALJ=',1PE12.5,'  BLJ=',E12.5,
     *'  CLJ=',E12.5,'  NREP=',I4,'  MREP=',I4,'  LREP=',I4)
  834 FORMAT('   NI=',I4,'  NJ=',I4,'  NK=',I4,'  NL=',I4,'  NM=',
     *I4,'  NN=',I4,'  VZ=',F7.3)
  837 FORMAT('   NJ=',I4,'  NK=',I4,'  RNPOW=',f4.1,'  AEX=',1PE12.5,
     *'  BEX=',E12.5,'  CEX=',E12.5)
  849 FORMAT('   NJ=',I4,'  NK=',I4,'  RZ=',F6.3,'  B=',F6.3,'  D=',F8.3
     *,'   DELTA1=',F7.3)
  850 FORMAT('   NJ=',I4,'  NK=',I4,'  RZ(S)=',F6.3,'  B(S)=',F6.3,
     *'  D(S)=',F8.3/19X,'  RZ(T)=',F6.3,'  B(T)=',F6.3,'  D(T)=',
     *F8.3)
  854 FORMAT(' ATOMS FORMING THE BOND:'/
     *5X,' N21I=',I4,'   N21V=',I4/
     *5X,' N21J=',I4,'   N21M=',I4,'   N21N=',I4,'   N21O=',I4/
     *5X,' N21K=',I4,'   N21P=',I4,'   N21Q=',I4,'   N21R=',I4/
     *5X,' N21L=',I4,'   N21S=',I4,'   N21T=',I4,'   N21U=',I4)
  856 FORMAT(5X,' PH0D=',F10.5,5X,F10.5,5X,F10.5/
     *5X,' PH0R=',F10.5,5X,F10.5,5X,F10.5/
     *5X,'  F0D=',F10.5,5X,F10.5,5X,F10.5/
     *5X,' TH0D=',F10.5,5X,F10.5,5X,F10.5/
     *5X,' TH0R=',F10.5,5X,F10.5,5X,F10.5/
     *5X,' FCCC=',F10.5,5X,F10.5,5X,F10.5/
     *5X,' CC0D=',F10.5,5X,F10.5,5X,F10.5/
     *5X,' CC0R=',F10.5,5X,F10.5,5X,F10.5/
     *5X,'  BET=',F10.5,5X,F10.5,5X,F10.5/
     *5X,'  DIS=',F10.5,5X,F10.5,5X,F10.5/
     *5X,' GA0D=',F10.5,5X,F10.5,5X,F10.5/
     *        11X,F10.5,5X,F10.5,5X,F10.5/
     *5X,' GA0R=',F10.5,5X,F10.5,5X,F10.5/
     *        11X,F10.5,5X,F10.5,5X,F10.5/
     *5X,' FCC1=',F10.5,5X,F10.5,5X,F10.5/
     *        11X,F10.5,5X,F10.5,5X,F10.5/
     *5X,' FCC2=',F10.5,5X,F10.5,5X,F10.5/)
  858 FORMAT('   N20J=',I4,'  N20K=',I4,'  ASR=',1PE12.5,'  BSR=',E12.5,
     *'  CSR=',E12.5,'  NDSR=',I4,'  ALR=',E12.5,'  BLR=',E12.5,
     *'  CLR=',E12.5,'  NDLR=',I4,'  AAS=',E12.5,'  AA0=',E12.5,
     *'  NAA=',I4)
  898 FORMAT('     NB=',I4,' NBB=',I4,' NTT=',I4,'  FTT=',F9.5)
  899 FORMAT('   NJ=',I4,'  NK=',I4,'  NM=',I4,'  ZAXT=',1PE12.5)
  905 FORMAT('   CUBIC BETA:   CM1=',F10.6,'   CM2=',F10.6,'   CM3=',
     *F10.6,'   CM4=',F10.6,'   CM5=',F10.6)
  923 FORMAT('   NI=',I4,'  NJ=',I4,'  NK=',I4,'  NL=',I4,'  NM=',
     *I4,'  NN=',I4,'  GC1=',F10.5,'  GQ1=',1PE15.7,'  GQ2=',E15.7)
  926 FORMAT('   NI=',I4,'  NJ=',I4,'  NK=',I4,'  NL=',I4,
     *'  NM=',I4)
  927 FORMAT('   QUADRATIC FORCE CONSTANTS -- EQUILIBRIUM VALUES')
  928 FORMAT(5X,'FIJ=',F7.4,' FIK=',F7.4,' FIL=',F7.4,' FJK=',F7.4,
     *' FJL=',F7.4,' FKL=',F7.4)
  929 FORMAT('   QUADRATIC FORCE CONSTANTS -- ASYMPTOTIC VALUES')
  930 FORMAT('   CUBIC FORCE CONSTANTS -- EQUILIBRIUM VALUES')
  931 FORMAT(5X,'GIJ=',F7.4,' GIK=',F7.4,' GIL=',F7.4,' GJK=',F7.4,
     *' GJL=',F7.4,' GKL=',F7.4)
  932 FORMAT('   CUBIC FORCE CONSTANTS -- ASYMPTOTIC VALUES')
  933 FORMAT('   QUARTIC FORCE CONSTANTS -- EQUILIBRIUM VALUES')
  934 FORMAT(5X,'HIJ=',F7.4,' HIK=',F7.4,' HIL=',F7.4,' HJK=',F7.4,
     *' HJL=',F7.4,' HKL=',F7.4)
  935 FORMAT('   QUARTIC FORCE CONSTANTS -- ASYMPTOTIC VALUES')
  936 FORMAT('   EQUILIBRIUM ANGLE VALUES')
  937 FORMAT(5X,'THIJ=',F8.4,' THIK=',F8.4,' THIL=',F8.4,
     *' THJK=',F8.4,' THJL=',F8.4,' THKL=',F8.4)
  938 FORMAT('   ASYMPTOTIC ANGLE VALUES -- SET ONE')
  939 FORMAT('   ASYMPTOTIC ANGLE VALUES -- SET TWO')
  940 FORMAT('   OUT-OF-PLANE QUADRATIC AND QUARTIC ',
     *'FORCE CONSTANTS')
  941 FORMAT(5X,'F1=',F7.4,'  F2=',F7.4,'  F3=',F7.4,'  F4=',F7.4)
  942 FORMAT(5X,'H1=',F7.4,'  H2=',F7.4,'  H3=',F7.4,'  H4=',F7.4)
  943 FORMAT('   NON-DIAGONAL CUBIC FORCE CONSTANTS')
  944 FORMAT(5X,'GN1=',F7.4,' GN2=',F7.4,' GN3=',F7.4,
     *' GN4=',F7.4,' GN5=',F7.4)
  945 FORMAT('   EQUILIBRIUM BOND LENGTHS')
  946 FORMAT(5X,'R1=',F12.8,'  R2=',F12.8,'  R3=',F12.8,'  R4=',F12.8)
  949 FORMAT('   NI=',I4,'  NJ=',I4,'  NK=',I4,'  NL=',I4,' FRRZ=',F9.5,
     *' RIJ0=',F8.4,' CIJ=',F6.3,' RKL0=',F8.4,' CKL=',F6.3)
  950 FORMAT('   NI=',I4,'  NJ=',I4,'  NB=',I4,' NRT=',I4,'  FRT=',F9.5,
     *'  R110=',F8.4,' CIJ=',F6.3)
  951 FORMAT('   DIHEDRAL ANGLE',I3,',  NI=',I4,
     *'  NJ=',I4,'  NK=',I4,'  NL=',I4)
  952 FORMAT(10X,'N=  ',I2,'  FDH=',F7.3,'  GDH=',F8.3)
  960 FORMAT('   **  VZ LESS THAN ZERO INDICATES 3-FOLD TORSION  **')
C
C         READ POTENTIAL PARAMETERS
C
      IF (NST.NE.0) THEN
         DO I=1,NST
            READ(5,*)N1J(I),N1K(I),RSZ(I),FS(I)
            IF (RSZ(I).EQ.0.0D0) THEN
               RSZ(I)=RSZ(I-1)
               FS(I)=FS(I-1)
            ENDIF
            WRITE(6,803)N1J(I),N1K(I),RSZ(I),FS(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NM.NE.0) THEN
         DO I=1,NM
            READ(5,*)N2J(I),N2K(I),RMZ(I),B(I),D(I)
            IF (RMZ(I).EQ.0.0D0) THEN
               RMZ(I)=RMZ(I-1)
               B(I)=B(I-1)
               D(I)=D(I-1)
            ENDIF
            WRITE(6,804)N2J(I),N2K(I),RMZ(I),B(I),D(I)
            IF (B(I).LE.0.0D0) THEN
               READ(5,*)CM1(I),CM2(I),CM3(I),CM4(I),CM5(I)
               WRITE(6,905)CM1(I),CM2(I),CM3(I),CM4(I),CM5(I)
            ENDIF
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NB.NE.0) THEN
         DO I=1,NB
            READ(5,*)N3J(I),N3K(I),N3M(I),THETAZ(I),FBZ(I),RJZ(I),
     *               CJ(I),RKZ(I),CK(I)
            IF (THETAZ(I).EQ.0.0D0) THEN
               THETAZ(I)=THETAZ(I-1)
               FBZ(I)=FBZ(I-1)
               RJZ(I)=RJZ(I-1)
               CJ(I)=CJ(I-1)
               RKZ(I)=RKZ(I-1)
               CK(I)=CK(I-1)
            ENDIF
            WRITE(6,805)N3J(I),N3K(I),N3M(I),THETAZ(I),FBZ(I),
     *                  RJZ(I),CJ(I),RKZ(I),CK(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NA.NE.0) THEN 
         DO I=1,NA
            READ(5,*)N4J(I),N4K(I),N4M(I),N4N(I),FA(I)
            WRITE(6,830)N4J(I),N4K(I),N4M(I),N4N(I),FA(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NLJ.NE.0) THEN
         DO I=1,NLJ
            READ(5,*)N5J(I),N5K(I),NREP(I),MREP(I),LREP(I),
     *               ALJ(I),BLJ(I),CLJ(I)
            IF (NREP(I).EQ.0.AND.MREP(I).EQ.0) THEN
               NREP(I)=NREP(I-1)
               MREP(I)=MREP(I-1)
               LREP(I)=LREP(I-1)
               ALJ(I)=ALJ(I-1)
               BLJ(I)=BLJ(I-1)
               CLJ(I)=CLJ(I-1)
            ENDIF
            WRITE(6,831)N5J(I),N5K(I),ALJ(I),BLJ(I),CLJ(I),
     *                  NREP(I),MREP(I),LREP(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NTAU.NE.0) THEN
         WRITE(6,960)
         DO I=1,NTAU
            READ(5,*)N6I(I),N6J(I),N6K(I),N6L(I),N6M(I),N6N(I),VZTAU(I)
            WRITE(6,834)N6I(I),N6J(I),N6K(I),N6L(I),N6M(I),N6N(I),
     *                  VZTAU(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NEXP.NE.0) THEN
         DO I=1,NEXP
            READ(5,*)N7J(I),N7K(I),RNPOW(I),AEX(I),BEX(I),CEX(I)
            IF (RNPOW(I).EQ.0.AND.I.NE.1) THEN
               RNPOW(I)=RNPOW(I-1)
               AEX(I)=AEX(I-1)
               BEX(I)=BEX(I-1)
               CEX(I)=CEX(I-1)
             ENDIF
             WRITE(6,837)N7J(I),N7K(I),RNPOW(I),AEX(I),BEX(I),CEX(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NGHOST.NE.0) THEN
         DO I=1,NGHOST
            READ(5,*)N8I(I),N8J(I),N8K(I),N8L(I),N8M(I),N8N(I),
     *               GC1(I),GEX1(I),GEX2(I)
            WRITE(6,923)N8I(I),N8J(I),N8K(I),N8L(I),N8M(I),N8N(I),
     *                  GC1(I),GEX1(I),GEX2(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NTET.NE.0) THEN
         DO I=1,NTET
            READ(5,*)N9I(I),N9J(I),N9K(I),N9L(I),N9M(I)
            READ(5,*)(FT0(I,J),J=1,6)
            READ(5,*)(FT2(I,J),J=1,6)
            READ(5,*)(GT0(I,J),J=1,6)
            READ(5,*)(GT2(I,J),J=1,6)
            READ(5,*)(HT0(I,J),J=1,6)
            READ(5,*)(HT2(I,J),J=1,6)
            READ(5,*)(THT(I,J),J=1,6)
            READ(5,*)(THT1(I,J),J=1,6)
            READ(5,*)(THT2(I,J),J=1,6)
            READ(5,*)(FD1(I,J),J=1,4)
            READ(5,*)(HD1(I,J),J=1,4)
            READ(5,*)(GN0(I,J),J=1,5)
            READ(5,*)(R0(I,J),J=1,4)
            WRITE(6,926)N9I(I),N9J(I),N9K(I),N9L(I),N9M(I)
            WRITE(6,927)
            WRITE(6,928)(FT0(I,J),J=1,6)
            WRITE(6,929)
            WRITE(6,928)(FT2(I,J),J=1,6)
            WRITE(6,930)
            WRITE(6,931)(GT0(I,J),J=1,6)
            WRITE(6,932)
            WRITE(6,931)(GT2(I,J),J=1,6)
            WRITE(6,933)
            WRITE(6,934)(HT0(I,J),J=1,6)
            WRITE(6,935)
            WRITE(6,934)(HT2(I,J),J=1,6)
            WRITE(6,936)
            WRITE(6,937)(THT(I,J),J=1,6)
            WRITE(6,938)
            WRITE(6,937)(THT1(I,J),J=1,6)
            WRITE(6,939)
            WRITE(6,937)(THT2(I,J),J=1,6)
            WRITE(6,940)
            WRITE(6,941)(FD1(I,J),J=1,4)
            WRITE(6,942)(HD1(I,J),J=1,4)
            WRITE(6,943)
            WRITE(6,944)(GN0(I,J),J=1,5)
            WRITE(6,945)
            WRITE(6,946)(R0(I,J),J=1,4)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NVRR.NE.0) THEN
         DO I=1,NVRR
            READ(5,*)N10I(I),N10J(I),N10K(I),N10L(I),
     *           FKRRZ(I),RIJ0(I),CIJ(I),RKL0(I),CKL(I)
            WRITE(6,949)N10I(I),N10J(I),N10K(I),N10L(I),
     *           FKRRZ(I),RIJ0(I),CIJ(I),RKL0(I),CKL(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NVRT.NE.0) THEN
         DO I=1,NVRT
            READ(5,*)N11I(I),N11J(I),N11B(I),NRT(I),
     *               FKRTZ(I),R110(I),CRT(I)
            WRITE(6,950)N11I(I),N11J(I),N11B(I),NRT(I),
     *               FKRTZ(I),R110(I),CRT(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NVTT.NE.0) THEN
         DO I=1,NVTT
            READ(5,*)N12B(I),N12BB(I),NTT(I),FKTTZ(I)
            WRITE(6,898)N12B(I),N12BB(I),NTT(I),FKTTZ(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NANG.NE.0) THEN
         DO I=1,NANG
            READ(5,*)N13I(I),N13J(I),N13K(I),N13L(I)
            READ(5,*)NDH(I)
            K=NDH(I)
            WRITE(6,951)I,N13I(I),N13J(I),N13K(I),N13L(I)
            DO J=1,K
               READ(5,*) FDH(I,J),GDH(I,J)
               WRITE(6,952) J,FDH(I,J),GDH(I,J)
            ENDDO
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NAXT.NE.0) THEN
         DO I=1,NAXT
            READ(5,*)N14I(I),N14J(I),N14K(I),ZAXT(I)
            WRITE(6,899)N14I(I),N14J(I),N14K(I),ZAXT(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NRYD.NE.0) THEN
         DO I=1,NRYD
            READ(5,*)N16J(I),N16K(I),RYDZ(I),DRYD(I),ARYD(I)
            WRITE(6,804)N16J(I),N16K(I),RYDZ(I),ARYD(I),DRYD(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NHFD.NE.0) THEN
         IF (NHFD.LT.0) THEN
            I=0
            IF (NHFD.EQ.-1.OR.NHFD.EQ.-3) THEN
               READ(5,*)N17JA,N17KA,RHFDA,AHFDA,BHFDA,
     *         C6HFDA,C8HFDA,C10HFDA
               WRITE(6,806)N17JA,N17KA
               WRITE(6,801)RHFDA,AHFDA,BHFDA,C6HFDA,C8HFDA,C10HFDA
               DO J=N17JA,N17KA
                  DO K=J+1,N17KA
                     I=I+1
                     N17J(I)=J
                     N17K(I)=K
                     RHFD(I)=RHFDA
                     RHFD(I)=RHFDA
                     AHFD(I)=AHFDA
                     BHFD(I)=BHFDA
                     C6HFD(I)=C6HFDA
                     C8HFD(I)=C8HFDA
                     C10HFD(I)=C10HFDA
                  ENDDO
               ENDDO
            ENDIF
            IF (NHFD.EQ.-2.OR.NHFD.EQ.-3) THEN
               READ(5,*)N17JA,N17KA,RHFDA,AHFDA,BHFDA,
     *         C6HFDA,C8HFDA,C10HFDA,K
               IF (I.NE.0) WRITE(6,815)
               WRITE(6,807)K,N17JA,N17KA
               WRITE(6,801)RHFDA,AHFDA,BHFDA,C6HFDA,C8HFDA,C10HFDA
               DO J=N17JA,N17KA
                  I=I+1
                  N17J(I)=MIN0(J,K)
                  N17K(I)=MAX0(J,K)
                  RHFD(I)=RHFDA
                  RHFD(I)=RHFDA
                  AHFD(I)=AHFDA
                  BHFD(I)=BHFDA
                  C6HFD(I)=C6HFDA
                  C8HFD(I)=C8HFDA
                  C10HFD(I)=C10HFDA
               ENDDO
            ENDIF
            NHFD=I
         ELSE
            DO I=1,NHFD
               READ(5,*)N17J(I),N17K(I),RHFD(I),AHFD(I),BHFD(I),
     *         C6HFD(I),C8HFD(I),C10HFD(I)
               WRITE(6,808)N17J(I),N17K(I),RHFD(I),AHFD(I),BHFD(I),
     *         C6HFD(I),C8HFD(I),C10HFD(I)
            ENDDO
         ENDIF
         WRITE(6,815)
      ENDIF
C
      IF (NLEPSA.NE.0) THEN
         DO I=1,NLEPSA
            READ(5,*)N18J1(I),N18K1(I),RLZ1(I),BL1(I),DL1(I),DELTA1(I)
            WRITE(6,849)N18J1(I),N18K1(I),RLZ1(I),BL1(I),DL1(I),
     *                  DELTA1(I)
            READ(5,*)N18J2(I),N18K2(I),RLZ2(I),BL2(I),DL2(I),DELTA2(I)
            WRITE(6,849)N18J2(I),N18K2(I),RLZ2(I),BL2(I),DL2(I),
     *                  DELTA2(I)
            READ(5,*)N18J3(I),N18K3(I),RLZ3(I),BL3(I),DL3(I),DELTA3(I)
            WRITE(6,849)N18J3(I),N18K3(I),RLZ3(I),BL3(I),DL3(I),
     *                  DELTA3(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NLEPSB.NE.0) THEN
         DO I=1,NLEPSB
            READ(5,*)N19J1(I),N19K1(I),RLZS1(I),BLS1(I),DLS1(I),
     *               RLZT1(I),BLT1(I),DLT1(I)
            WRITE(6,850)N19J1(I),N19K1(I),RLZS1(I),BLS1(I),DLS1(I),
     *                  RLZT1(I),BLT1(I),DLT1(I)
            READ(5,*)N19J2(I),N19K2(I),RLZS2(I),BLS2(I),DLS2(I),
     *               RLZT2(I),BLT2(I),DLT2(I)
            WRITE(6,850)N19J2(I),N19K2(I),RLZS2(I),BLS2(I),DLS2(I),
     *                  RLZT2(I),BLT2(I),DLT2(I)
            READ(5,*)N19J3(I),N19K3(I),RLZS3(I),BLS3(I),DLS3(I),
     *               RLZT3(I),BLT3(I),DLT3(I)
            WRITE(6,850)N19J3(I),N19K3(I),RLZS3(I),BLS3(I),DLS3(I),
     *                  RLZT3(I),BLT3(I),DLT3(I)
         ENDDO
         WRITE(6,815)
      ENDIF
C
      IF (NRAX.NE.0) THEN
         DO J=1,NRAX
            READ(5,*)N21I(J),N21V(J)
            READ(5,*)N21J(J),N21M(J),N21N(J),N21O(J)
            READ(5,*)N21K(J),N21P(J),N21Q(J),N21R(J)
            READ(5,*)N21L(J),N21S(J),N21T(J),N21U(J)
            WRITE(6,854)N21I(J),N21V(J),
     *                  N21J(J),N21M(J),N21N(J),N21O(J),
     *                  N21K(J),N21P(J),N21Q(J),N21R(J),
     *                  N21L(J),N21S(J),N21T(J),N21U(J)
            READ(5,*)(PH0D(I,J),I=1,3)
            READ(5,*)(PH0R(I,J),I=1,3)
            READ(5,*)(F0D(I,J),I=1,3)
            READ(5,*)(TH0D(I,J),I=1,3)
            READ(5,*)(TH0R(I,J),I=1,3)
            READ(5,*)(FCCC(I,J),I=1,3)
            READ(5,*)(CC0D(I,J),I=1,3)
            READ(5,*)(CC0R(I,J),I=1,3)
            READ(5,*)(BET(I,J),I=1,3)
            READ(5,*)(DIS(I,J),I=1,3)
            READ(5,*)(GA0D(I,J),I=1,6)
            READ(5,*)(GA0R(I,J),I=1,6)
            READ(5,*)(FCC1(I,J),I=1,6)
            READ(5,*)(FCC2(I,J),I=1,3)

            WRITE(6,856)(PH0D(I,J),I=1,3),
     *                  (PH0R(I,J),I=1,3),
     *                  (F0D(I,J),I=1,3),
     *                  (TH0D(I,J),I=1,3),
     *                  (TH0R(I,J),I=1,3),
     *                  (FCCC(I,J),I=1,3),
     *                  (CC0D(I,J),I=1,3),
     *                  (CC0R(I,J),I=1,3),
     *                  (BET(I,J),I=1,3),
     *                  (DIS(I,J),I=1,3),
     *                  (GA0D(I,J),I=1,6),
     *                  (GA0R(I,J),I=1,6),
     *                  (FCC1(I,J),I=1,6),
     *                  (FCC2(I,J),I=1,3)
         ENDDO
      ENDIF
C
      IF (NONB.NE.0) THEN
         DO I=1,NONB
            READ(5,*)N22J(I),N22K(I),ASR(I),BSR(I),CSR(I),NDSR(I),
     *               ALR(I),BLR(I),CLR(I),NDLR(I)
            READ(5,*)AAS(I),AA0(I),NAA(I)
            WRITE(6,858)N22J(I),N22K(I),ASR(I),BSR(I),CSR(I),NDSR(I),
     *             ALR(I),BLR(I),CLR(I),NDLR(I),AAS(I),AA0(I),NAA(I)
         ENDDO
      ENDIF
C
C     ADD THIS CODE BLOCK FOR CRCO6 ON 11/19/2003 BY MING
C
      IF (NCRCO6.NE.0) THEN
         a1=5.20835
         a2=3.05986
         a3=4.42527
         a4=5.20312
         a5=3.01143
         a6=2.04798
         a7=0.0297789
         a8=0.478606

         dfixed(1)=32.3
         dfixed(2)=15.2
         dfixed(3)=13.6
         dfixed(4)=13.6
         dfixed(5)=22.6
         dfixed(6)=21.9       
         
C
C     Convert a1-a8 to venus unit
C     
         a1=a1*dsqrt(C1)
         a4=a4*dsqrt(C1)
         a5=a5*dsqrt(C1)
         
         DO i=1,6
             
             dfixed(i)=dfixed(i)*C1
             write(6,'(1x,a,i1,a,f12.6)')'D(',i,')=',dfixed(i)
         ENDDO
         WRITE(6,815)
      ENDIF      
C
C         CONVERT POTENTIAL PARAMETERS TO INTEGRATION UNITS
C
      DO I=1,NST
         FS(I)=FS(I)*C2
      ENDDO
C
      DO I=1,NM
         D(I)=D(I)*C1
      ENDDO
C
      DO I=1,NB
         FBZ(I)=FBZ(I)*C3
         THETAZ(I)=THETAZ(I)*C4
      ENDDO
C
      DO I=1,NA
         FA(I)=FA(I)*C3
      ENDDO
C
      DO I=1,NLJ
         ALJ(I)=ALJ(I)*C1
         BLJ(I)=BLJ(I)*C1
         CLJ(I)=CLJ(I)*C1
      ENDDO
C
      DO I=1,NTAU
         VZTAU(I)=VZTAU(I)*C1
      ENDDO
C
      DO I=1,NEXP
         AEX(I)=AEX(I)*C1
         CEX(I)=CEX(I)*C1
      ENDDO
C
      DO I=1,NGHOST
         GEX1(I)=GEX1(I)*C1
         GEX2(I)=GEX2(I)*C1
      ENDDO
C
      DO I=1,NTET
         DO J=1,6
            FT0(I,J)=FT0(I,J)*C3
            FT2(I,J)=FT2(I,J)*C3
            GT0(I,J)=GT0(I,J)*C3
            GT2(I,J)=GT2(I,J)*C3
            HT0(I,J)=HT0(I,J)*C3
            HT2(I,J)=HT2(I,J)*C3
            THT(I,J)=THT(I,J)*C4
            THT1(I,J)=THT1(I,J)*C4
            THT2(I,J)=THT2(I,J)*C4
         ENDDO
         DO J=1,4
            FD1(I,J)=FD1(I,J)*C3
            HD1(I,J)=HD1(I,J)*C3
         ENDDO
         DO J=1,5
            GN0(I,J)=GN0(I,J)*C3
         ENDDO
      ENDDO
C
      DO I=1,NVRR
         FKRRZ(I)=FKRRZ(I)*C2
      ENDDO
C
      DO I=1,NVRT
         FKRTZ(I)=FKRTZ(I)*C2
      ENDDO
C
      DO I=1,NVTT
         FKTTZ(I)=FKTTZ(I)*C3
      ENDDO
C
      DO I=1,NANG
         K=NDH(I)
         DO J=1,K
            FDH(I,J)=FDH(I,J)*C1
            GDH(I,J)=GDH(I,J)*C4
         ENDDO
      ENDDO
C
      DO I=1,NAXT
         ZAXT(I)=ZAXT(I)*C1
      ENDDO
C
      DO I=1,NRYD
         DRYD(I)=DRYD(I)*C1
      ENDDO
C
      DO I=1,NHFD
         AHFD(I)=AHFD(I)*C1
         C6HFD(I)=C6HFD(I)*C1
         C8HFD(I)=C8HFD(I)*C1
         C10HFD(I)=C10HFD(I)*C1
      ENDDO
C
      DO I=1,NLEPSA
         DL1(I)=DL1(I)*C1
         DL2(I)=DL2(I)*C1
         DL3(I)=DL3(I)*C1
      ENDDO
C
      DO I=1,NLEPSB
         DLS1(I)=DLS1(I)*C1
         DLS2(I)=DLS2(I)*C1
         DLS3(I)=DLS3(I)*C1
         DLT1(I)=DLT1(I)*C1
         DLT2(I)=DLT2(I)*C1
         DLT3(I)=DLT3(I)*C1
      ENDDO
C
      DO J=1,NRAX
         DO I=1,3
            PH0D(I,J)=PH0D(I,J)*C4
            PH0R(I,J)=PH0R(I,J)*C4
            F0D(I,J)=F0D(I,J)*C3
            TH0D(I,J)=TH0D(I,J)*C4
            TH0R(I,J)=TH0R(I,J)*C4
            FCCC(I,J)=FCCC(I,J)*C3
            FCC2(I,J)=FCC2(I,J)*C3
            DIS(I,J)=DIS(I,J)*C1
         ENDDO
         DO I=1,6
            GA0D(I,J)=GA0D(I,J)*C4
            GA0R(I,J)=GA0R(I,J)*C4
            FCC1(I,J)=FCC1(I,J)*C3
         ENDDO
      ENDDO
C
      DO I=1,NONB
         ASR(I)=ASR(I)*C1
         CSR(I)=CSR(I)*C1
         ALR(I)=ALR(I)*C1
         CLR(I)=CLR(I)*C1
      ENDDO
C
      RETURN
      END
