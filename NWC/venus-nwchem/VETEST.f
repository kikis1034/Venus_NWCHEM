      SUBROUTINE VETEST
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CHECK FOR INTERMEDIATE AND FINAL EVENTS
C
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/PSN2/PESN2,GA,RA,RB
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,L(NDA),NAM
      COMMON/TESTIN/VRELO,INTST
      COMMON/TESTSN2/GAO,NSAD,NCBA,NCAB,IBAR
      COMMON/FINALB/EROTA,EROTB,EA(3),EB(3),AMA(4),AMB(4),AN,AJ,BN,BJ,
     *OAM(4),EREL,ERELSQ,BF,SDA,SDB,DELH(NDP),ANG(NDG),NFINAL
      COMMON/VMAXB/QVMAX(NDA3),PVMAX(NDA3),VMAX,NCVMAX      
c      COMMON/SELTB/QZ(nda3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      DIMENSION QCMA(3),VCMA(3),QCMB(3),VCMB(3),QR(3),VR(3)
      CHARACTER*10 TYPE
      CHARACTER*15 COMP
C
 900  FORMAT(4X,' TURNING POINT #  ','  CYCLE ','  RCM(A)',
     *       '    EA     ','    EB     ','    EROTA  ','    EROTB  ',
     *       '    JA     ','    JB     ','    L      ')
 903  FORMAT(8X,A10,I4,I8,F8.3,1P7D11.4)
 905  FORMAT(/5X,'$$$$BARRIER CROSSING NUMBER$$$$ ',I6,
     &       '  AT CYCLE',I8)
 906  FORMAT(5X,'$$$$BARRIER CROSSING FROM B TO A $$$$')
 907  FORMAT(5X,'$$$$BARRIER CROSSING FROM A TO B $$$$')
 935  FORMAT(7X,'RA= ',F7.3,3X,'RB= ',F7.3,3X,'GA= ',F7.3)
 910  FORMAT(4X,' TURNING POINT #  ',5X,' COMPLEX ',6X,
     *'  CYCLE ','  RCM(A)',
     *       '    EA     ','    EB     ','    EROTA  ','    EROTB  ',
     *       '    JA     ','    JB     ','    L      ')
 913  FORMAT(8X,A10,I4,2X,A17,I8,F8.3,1P7D11.4)
C
C         TEST FOR MAXIMUM IN POTENTIAL ENERGY
C
      CALL ENERGY      
      IF (V.GT.VMAX) THEN
	 VMAX=V
	 NCVMAX=NC
	 DO I=1,NDA3
	 QVMAX(I)=Q(I)
	 PVMAX(I)=P(I)
	 ENDDO
      ENDIF
C
C         COUNT INNER TURNING POINTS BETWEEN THE TWO REACTANT'S
C         CENTERS OF MASS.
C
      IF (NSN2.EQ.0) THEN
         IF (NABJ(1).EQ.NABK(1)) GOTO 1010
C
C         CALCULATE CENTER OF MASS COORDINATES AND VELOCITIES FOR A
C
         WT=WTA(1)
         NN=NATOMA(1)
         DO I=1,NN
            L(I)=LA(1,I)
         ENDDO
         CALL CENMAS(WT,QCMA,VCMA,NN)
C
C         CALCULATE CENTER OF MASS COORDINATES AND VELOCITIES FOR B
C
         WT=WTB(1)
         NN=NATOMB(1)
         DO I=1,NN
            L(I)=LB(1,I)
         ENDDO
         CALL CENMAS(WT,QCMB,VCMB,NN)
C
C         CALCULATE INTERNUCLEAR SEPARATION BETWEEN THE CENTERS-OF-
C         MASS, AND THE RELATIVE VELOCITY ALONG THE INTERNUCLEAR AXIS.
C
         DO I=1,3
            QR(I)=QCMA(I)-QCMB(I)
            VR(I)=VCMA(I)-VCMB(I)
         ENDDO
         RCM=0.0D0
         VREL=0.0D0
         DO I=1,3
            RCM=RCM+QR(I)**2
            VREL=VREL+QR(I)*VR(I)
         ENDDO
         RCM=SQRT(RCM)
         VREL=VREL/RCM
C
         SGN=VRELO*VREL
         IF (SGN.LT.0.0D0) THEN
            INTST=INTST+1
            IF (INTST.EQ.1) WRITE(6,900)
C
C         WRITE ENERGIES AND ANGULAR MOMENTUM FOR PATH=1 AT THE
C         TURNING POINT
C
            NPATH=1
            CALL FINAL
            IF (VRELO.LT.0.0D0) TYPE=' **INNER** '
            IF (VRELO.GT.0.0D0) TYPE=' **OUTER** '
            WRITE(6,903)TYPE,INTST,NC,RCM,EA(3),EB(3),EROTA,
     *                  EROTB,AMA(4),AMB(4),OAM(4)
            NFINAL=0
         ENDIF
         VRELO=VREL
C 
C         CODE FOR SN2 DYNAMICS
C
      ELSE
         IF (NSN2.LE.6) GATS=0.0D0
         IF (NSN2.EQ.7.OR.NSN2.EQ.8) GATS=0.008D0
         CBAR=GAO*(GA-GATS)
         IF (CBAR.LT.0.0D0) THEN
            IBAR=IBAR+1
            WRITE(6,905)IBAR,NC
            WRITE(6,935)RA,RB,GA
            IF (GAO.LT.0.0D0) THEN
               NCBA=NCBA+1
               WRITE(6,906)
            ENDIF
            IF (GAO.GT.0.0D0) THEN
               NCAB=NCAB+1
               WRITE(6,907)
            ENDIF
         ENDIF 
C
         GAO=GA-GATS
         IF (NABJ(1).EQ.NABK(1)) THEN
            IF (GA.GE.GATS) NCPL=2
            IF (GA.LT.GATS) NCPL=3
         ELSE
            IF (GA.GE.GATS) NCPL=1
            IF (GA.LT.GATS) NCPL=2
         ENDIF
         WT=WTA(NCPL)
         NN=NATOMA(NCPL)
         DO I=1,NN
            L(I)=LA(NCPL,I)
         ENDDO
         CALL CENMAS(WT,QCMA,VCMA,NN)
C
C         CALCULATE CENTER OF MASS COORDINATES AND VELOCITIES FOR B
C
         WT=WTB(NCPL)
         NN=NATOMB(NCPL)
         DO I=1,NN
            L(I)=LB(NCPL,I)
         ENDDO
         CALL CENMAS(WT,QCMB,VCMB,NN)
C
C         CALCULATE INTERNUCLEAR SEPARATION BETWEEN THE CENTERS OF 
C         MASS, AND THE RELATIVE VELOCITY ALONG THE INTERNUCLEAR AXIS.
C
         DO I=1,3
            QR(I)=QCMA(I)-QCMB(I)
            VR(I)=VCMA(I)-VCMB(I)
         ENDDO
         RCM=0.0D0
         VREL=0.0D0
         DO I=1,3
            RCM=RCM+QR(I)*QR(I)
            VREL=VREL+QR(I)*VR(I)
         ENDDO
         RCM=SQRT(RCM)
         VREL=VREL/RCM
C
         SGN=VRELO*VREL
         IF (SGN.LT.0.0D0) THEN
            INTST=INTST+1
            IF (INTST.EQ.1) WRITE(6,910)
C
C         WRITE ENERGIES AND ANGULAR MOMENTUM FOR PATH=NCPL AT THE
C         TURNING POINT
C
            NPATH=NCPL
            CALL FINAL
            IF (NABJ(1).EQ.NABK(1)) THEN
               IF (NCPL.EQ.2) COMP='## X---CH3Y ## '
               IF (NCPL.EQ.3) COMP='## XCH3---Y ## '
            ELSE
               IF (NCPL.EQ.1) COMP='## X---CH3Y ## '
               IF (NCPL.EQ.2) COMP='## XCH3---Y ## '
            ENDIF
            IF (VRELO.LT.0.0D0) TYPE=' **INNER** '
            IF (VRELO.GT.0.0D0) TYPE=' **OUTER** '
            WRITE(6,913)TYPE,INTST,COMP,NC,RCM,EA(3),EB(3),
     *                  EROTA,EROTB,AMA(4),AMB(4),OAM(4)
            NFINAL=0
         ENDIF
         VRELO=VREL
      ENDIF
C
C         END OF CODE WHICH TESTS FOR INNER TURNING POINTS
C
 1010 CONTINUE
C
C         DEALING WITH TURNING POINT OF HEIGHT ABOVE SURFACE 
C
c      IF (NSURF.NE.0) CALL HTURN(1,13.0)
c
c         TEST FOR REACHING RBAR(I) OR RMAX(I)
c
      NTEST=0
      M=NPATHS+1
      DO I=1,M
         NPATH=I
         J3=3*NABJ(I)
         J2=J3-1
         J1=J2-1
         K3=3*NABK(I)
         IF (J3.NE.K3) THEN
            K2=K3-1
            K1=K2-1
            JK=(NABJ(I)-1)*(2*NATOMS-NABJ(I))/2+NABK(I)-NABJ(I)
            T1=Q(K1)-Q(J1)
            T2=Q(K2)-Q(J2)
            T3=Q(K3)-Q(J3)
            R(JK)=SQRT(T1*T1+T2*T2+T3*T3)
            IF (NABL(I).GT.0)THEN
                J3=3*NABL(I)
                J2=J3-1
                J1=J2-1
                K3=3*NABM(I)
                K2=K3-1
                K1=K2-1
                LM=(NABL(I)-1)*(2*NATOMS-NABL(I))/2+NABM(I)-NABL(I)
                T1=Q(K1)-Q(J1)
                T2=Q(K2)-Q(J2)
                T3=Q(K3)-Q(J3)
                R(LM)=SQRT(T1*T1+T2*T2+T3*T3)
                IF((R(JK).GT.RBAR(I)).AND.(R(LM).GT.RBAR(I))) NTEST =1
                IF((R(JK).GT.RMAX(I)).AND.(R(LM).GT.RMAX(I))) NTEST =2 
                IF (NTEST.GT.0) GOTO 3
            ELSE
                IF (R(JK).GT.RBAR(I)) NTEST=1
                IF (R(JK).GT.RMAX(I)) NTEST=2
                IF (NTEST.GT.0) GOTO 3
            ENDIF
             
         ENDIF
      ENDDO
      NPATH=1
    3 RETURN
      END
