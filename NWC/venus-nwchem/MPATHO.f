C
C         PERFORM REACTION PATH FOLLOWING
C
      SUBROUTINE MPATHO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      character*80 title1,title2
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,LL(NDA),NAM
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/PRFLAG/NFQP,NCOOR,NFR,NUMR,NFB,NUMB,NFA,NUMA,NFTAU,NUMTAU,
     *NFTET,NUMTET,NFDH,NUMDH,NFHT,NUMHT
      COMMON/PARRAY/KR(300),JR(300),KB(300),MB(300),IB(300),IA(300),
     *ITAU(300),ITET(300),IDH(300),IHT(300)
      COMMON/INERT/UXX,UXY,UXZ,UYY,UYZ,UZZ,AIXX,AIXY,AIXZ,AIYY,AIYZ,
     *AIZZ
      COMMON/FR2/DG(NDA3,2),DIM(NDA3)
      COMMON/ARRAYS/FA(NDA3yf,NDA3yf),DA(NDA3),B(NDA3yf,NDA3yf),DB(NDA3)
      COMMON/RSTART/HINC,NPTS
      COMMON/EIGVL/EIG(NDA3yf)
      COMMON/RKUTTA/RAA1,RA1,RA2,RA3,RB1,RB2,RB3,RC1,RC2
      COMMON/GPATHB/WM(NDA3),TEMP(NDP),AI1D(5),AAI(2),BBI(2),SYMM(5),
     *SYMA,SYMB,GTEMP(NDP),NFLAG(NDP),N1DR,N2DR
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/SYBB/TITLE1,TITLE2,SYBTI
      DIMENSION BQ(nda3),PJ(nda3yf,nda3yf),FAP(nda3yf,nda3yf),
     * QM(nda3),GZ(nda3),QCM(3),VCM(3),AM(4)
      DATA RHO/1.0D-08/
      DATA TFACT/2.81837552D05/
C
  150 FORMAT(1X,'XXXXXXXXXXXXXXXXXXXXX REACTION PATH  NC = ',I7,
     *' XXXXXXXXXXXXXXXXXXXXXXX')
  151 FORMAT(1H ,//,10X,33H*****CALCULATE REACTION PATH*****,//,
     *4X,3HNS=,I10,6H  NIP=,I10,5H  DS=,F13.10,/)
  152 FORMAT(/'    NC=',I8,5X,'POTENTIAL ENERGY(KCAL/MOLE)=',F9.4)
  153 FORMAT(10X,'CARTESIAN COORDINATES')
  154 FORMAT(3X,1p3D18.10,3X,3D18.10)
  155 FORMAT(10X,'EIGENVALUES')
  156 FORMAT(3X,1p7D14.6)
  157 FORMAT(10X,'PRINCIPAL MOMENTS OF INERTIA: IX,IY,IZ')
  158 FORMAT(4X,' ATOMS  ',5X,'BOND LENGTH(A)')
  159 FORMAT(2X,2I4,10X,1pe15.8)
  160 FORMAT(1X,3F11.7)
  161 FORMAT(10X,'REACTANT VIBRATIONAL FREQUENCIES(CM-1)')
  162 FORMAT('    TEMPERATURES=',10F9.2)
  164 FORMAT('    NUMBER OF 1-DIMENSIONAL ROTORS =',I2)
  165 FORMAT(10X,'SYMA=',F4.1)
  166 FORMAT(10X,'SYMB=',F4.1)
  167 FORMAT(10X,'AI =',F8.3,'SYMM =',F5.1)
  168 FORMAT('    NUMBER OF 2-DIMENSIONAL ROTORS =',I2)
  171 FORMAT(10X,'MOMENTS OF INERTIA FOR FRAGMENT A: IX,IY,IZ')
  172 FORMAT(10X,'REDUCED MOMENTS OF INERTIA FOR 2-D ROTOR:',
     *'  IX=',1pe14.6,'  IZ=',e14.6)
  173 FORMAT(5X,'ATOMS',7X,'ANGLE (DEGREES)')
  174 FORMAT(2X,3I3,8X,F8.3)
  175 FORMAT(10X,'MOMENTS OF INERTIA FOR FRAGMENT B: IX,IY,IZ')
C
C         READ CYCLE COUNTS AND INTEGRATION PARAMETERS
C
      READ(5,*)NS,NIP,DS
      WRITE(6,151)NS,NIP,DS
C
C         READ REACTANT VIBRATIONAL FREQUENCIES(CM-1)
C
      J=I3N-6
      READ(5,*)(WM(I),I=1,J)
      WRITE(6,161)
      WRITE(6,156)(WM(I),I=1,J)
C
C         READ TEMPERATURES(K) FOR CALCULATING THE FREE ENERGY
C         ALONG THE REACTION PATH
C
      READ(5,*)NTEMP
      READ(5,*)(TEMP(I),I=1,NTEMP)
      WRITE(6,162)(TEMP(I),I=1,NTEMP)
C
C         READ PARAMETERS FOR CALCULATING 1- AND 2-DIMENSIONAL INTERNAL
C         ROTOR PARTITION FUNCTIONS.  N2DR MAY EQUAL 0, 1, OR 2.
C
      READ(5,*)N1DR,N2DR
      WRITE(6,164)N1DR
      IF (N1DR.NE.0) THEN
         READ(5,*)(AI1D(I),I=1,N1DR)
         READ(5,*)(SYMM(I),I=1,N1DR)
         WRITE(6,167)(AI1D(I),SYMM(I),I=1,N1DR)
      ENDIF
      WRITE(6,168)N2DR
      IF (N2DR.NE.0) THEN
         READ(5,*)SYMA
         WRITE(6,165)SYMA
         IF (N2DR.NE.1) THEN
            READ(5,*)SYMB
            WRITE(6,166)SYMB
         ENDIF
      ENDIF
C
C         CALCULATE DIM ARRAY USED FOR MASS-WEIGHTING AND TOTAL MASS
C
      K=0
      WTM=0.0D0
      DO I=1,NATOMS
         WTM=WTM+W(I)
         DO J=1,3
            K=K+1
            DIM(K)=1.D0/SQRT(W(I))
         ENDDO
      ENDDO
C
C         SET INDEX FOR LL ARRAY
C
      DO I=1,NATOMS
         LL(I)=I
      ENDDO
C
C         CALCULATE CENTER OF MASS COORDINATES
C
      CALL CENMAS(WTM,QCM,VCM,NATOMS)
      DO I=1,I3N
         Q(I)=QQ(I)
      ENDDO 
      IF (NCOOR.EQ.1) THEN
         SYBTI=-DS*DBLE(NIP)
         CALL SYBMOL
      ENDIF
C
C         CONVERT GRADIENT AND CENTER OF MASS COORDINATES TO
C         MASS WEIGHTED COORDINATES
C
      DUM=0.0D0
      DO I=1,I3N
         QM(I)=QQ(I)/DIM(I)
         GZ(I)=-PDOT(I)*DIM(I)
         DUM=DUM+GZ(I)*GZ(I)
      ENDDO
      DUM=SQRT(DUM)
      DO I=1,I3N
         GZ(I)=GZ(I)/DUM
      ENDDO
C
C         CALCULATE CARTESIAN FORCE CONSTANT MATRIX
C
  100 CONTINUE
      REWIND 7
      I=0
      CALL NMODE(NATOMS,I)
      IF (NC.NE.1) NX=NX+NIP
C
C         CALCULATE MASS-WEIGHTED CARTESIAN FORCE CONSTANT MATRIX
C
      DO I=1,I3N
         DO J=1,I
            FA(I,J)=DIM(I)*FA(I,J)*DIM(J)
            FA(J,I)=FA(I,J)
         ENDDO
      ENDDO
C
C         CALCULATE INERTIA MATRIX AND INVERSE OF INERTIA MATRIX
C         THIS REQUIRES QQ ARRAY
C
      CALL ROTN(AM,EROT,NATOMS)
C
C         CALCULATE PROJECTOR COMPONENTS
C
   21 CONTINUE
      DO I=1,NATOMS
         I3=3*I
         I2=I3-1
         I1=I3-2
         DUM1=1.0D0/WTM/DIM(I1)/DIM(I1)
         PJ(I1,I1)=1.D0-(GZ(I1)**2+DUM1+QM(I2)**2*UZZ
     *             -2.D0*QM(I2)*QM(I3)*UYZ+QM(I3)**2*UYY)
         PJ(I2,I2)=1.D0-(GZ(I2)**2+DUM1+QM(I1)**2*UZZ
     *             -2.D0*QM(I1)*QM(I3)*UXZ+QM(I3)**2*UXX)
         PJ(I3,I3)=1.D0-(GZ(I3)**2+DUM1+QM(I1)**2*UYY
     *             -2.D0*QM(I1)*QM(I2)*UXY+QM(I2)**2*UXX)
         PJ(I1,I2)=-(GZ(I1)*GZ(I2)-QM(I1)*QM(I2)*UZZ+QM(I2)*QM(I3)*UXZ
     *             +QM(I1)*QM(I3)*UYZ-QM(I3)*QM(I3)*UXY)
         PJ(I2,I1)=PJ(I1,I2)
         PJ(I1,I3)=-(GZ(I1)*GZ(I3)+QM(I1)*QM(I2)*UYZ-QM(I2)*QM(I2)*UXZ
     *             -QM(I1)*QM(I3)*UYY+QM(I2)*QM(I3)*UXY)
         PJ(I3,I1)=PJ(I1,I3)
         PJ(I2,I3)=-(GZ(I2)*GZ(I3)-QM(I1)*QM(I1)*UYZ+QM(I1)*QM(I2)*UXZ
     *             +QM(I1)*QM(I3)*UXY-QM(I2)*QM(I3)*UXX)
         PJ(I3,I2)=PJ(I2,I3)
C
         L=I+1
         DO J=L,NATOMS
            J3=3*J
            J2=J3-1
            J1=J3-2
            DUM1=1.0D0/WTM/DIM(I1)/DIM(J1)
            PJ(I1,J1)=-(GZ(I1)*GZ(J1)+DUM1+QM(I2)*QM(J2)*UZZ
     *                -QM(I2)*QM(J3)*UYZ-QM(I3)*QM(J2)*UYZ
     *                +QM(I3)*QM(J3)*UYY)
            PJ(J1,I1)=PJ(I1,J1)
            PJ(I2,J2)=-(GZ(I2)*GZ(J2)+DUM1+QM(I1)*QM(J1)*UZZ
     *                -QM(I1)*QM(J3)*UXZ-QM(I3)*QM(J1)*UXZ
     *                +QM(I3)*QM(J3)*UXX)
            PJ(J2,I2)=PJ(I2,J2)
            PJ(I3,J3)=-(GZ(I3)*GZ(J3)+DUM1+QM(I1)*QM(J1)*UYY
     *                -QM(I1)*QM(J2)*UXY-QM(I2)*QM(J1)*UXY
     *                +QM(I2)*QM(J2)*UXX)
            PJ(J3,I3)=PJ(I3,J3)
            PJ(I1,J2)=-(GZ(I1)*GZ(J2)-QM(I2)*QM(J1)*UZZ
     *                +QM(I2)*QM(J3)*UXZ+QM(I3)*QM(J1)*UYZ
     *                -QM(I3)*QM(J3)*UXY)
            PJ(J2,I1)=PJ(I1,J2)
            PJ(I1,J3)=-(GZ(I1)*GZ(J3)+QM(I2)*QM(J1)*UYZ
     *                -QM(I2)*QM(J2)*UXZ-QM(I3)*QM(J1)*UYY
     *                +QM(I3)*QM(J2)*UXY)
            PJ(J3,I1)=PJ(I1,J3)
            PJ(I2,J1)=-(GZ(I2)*GZ(J1)-QM(I1)*QM(J2)*UZZ
     *                +QM(I1)*QM(J3)*UYZ+QM(I3)*QM(J2)*UXZ
     *                -QM(I3)*QM(J3)*UXY)
            PJ(J1,I2)=PJ(I2,J1)
            PJ(I2,J3)=-(GZ(I2)*GZ(J3)-QM(I1)*QM(J1)*UYZ
     *                +QM(I1)*QM(J2)*UXZ+QM(I3)*QM(J1)*UXY
     *                -QM(I3)*QM(J2)*UXX)
            PJ(J3,I2)=PJ(I2,J3)
            PJ(I3,J1)=-(GZ(I3)*GZ(J1)+QM(I1)*QM(J2)*UYZ
     *                -QM(I1)*QM(J3)*UYY-QM(I2)*QM(J2)*UXZ
     *                +QM(I2)*QM(J3)*UXY)
            PJ(J1,I3)=PJ(I3,J1)
            PJ(I3,J2)=-(GZ(I3)*GZ(J2)-QM(I1)*QM(J1)*UYZ
     *                +QM(I1)*QM(J3)*UXY+QM(I2)*QM(J1)*UXZ
     *                -QM(I2)*QM(J3)*UXX)
            PJ(J2,I3)=PJ(I3,J2)
         ENDDO
      ENDDO
C
C         CALCULATE PROJECTED FORCE CONSTANT MATRIX: (1-P)FA(1-P)
C
      DO I=1,I3N
         DO J=1,I
            FAP(I,J)=0.0D0
            DO L=1,I3N
               DUM1=0.0D0
               DO K=1,I3N
                  DUM1=DUM1+FA(L,K)*PJ(K,J)
               ENDDO
               FAP(I,J)=FAP(I,J)+PJ(I,L)*DUM1
            ENDDO
            FAP(J,I)=FAP(I,J)
         ENDDO
      ENDDO
C
C         DIAGONALIZE PROJECTED FORCE CONSTANT MATRIX
C
      CALL EIGN(FAP,FA,I3N,RHO)
C
C         CALCULATE EIGENVALUES
C
      DO I=1,I3N
         EIG(I)=SQRT(ABS(TFACT*EIG(I)))
      ENDDO
C
      WRITE(6,152)NC,V
      WRITE(6,153)
      WRITE(6,154)(Q(I),I=1,I3N)
C
C         WRITE COORDINATES FOR GRAPHICS
C
      IF (NCOOR.EQ.1.AND.NC.NE.1) THEN
         WRITE(8,150)NC
         WRITE(8,160)(Q(I),I=1,I3N)
      ENDIF
C
C         CALCULATE AND WRITE POSSIBLE INTERATOMIC DISTANCES
C
      IF (NFR.NE.0) THEN
         WRITE(6,158)
         DO I=1,NUMR
            J3=3*JR(I)
            J2=J3-1
            J1=J2-1
            K3=3*KR(I)
            K2=K3-1
            K1=K2-1
            T1=Q(K1)-Q(J1)
            T2=Q(K2)-Q(J2)
            T3=Q(K3)-Q(J3)
            RR=SQRT(T1*T1+T2*T2+T3*T3)
            WRITE(6,159)JR(I),KR(I),RR
         ENDDO
      ENDIF
C
C         CALCULATE AND WRITE POSSIBLE ANGLES 
C
      IF (NFB.NE.0) THEN
         WRITE(6,173)
         DO I=1,NUMB
            K3=3*KB(I)
            K2=K3-1
            K1=K2-1
            M3=3*MB(I)
            M2=M3-1
            M1=M2-1
            I3=3*IB(I)
            I2=I3-1
            I1=I2-1
            T1=Q(I1)-Q(M1)
            T2=Q(I2)-Q(M2)
            T3=Q(I3)-Q(M3)
            T4=Q(K1)-Q(M1)
            T5=Q(K2)-Q(M2)
            T6=Q(K3)-Q(M3)
            R1=SQRT(T1*T1+T2*T2+T3*T3)
            R2=SQRT(T4*T4+T5*T5+T6*T6)
            CTHETA=(T1*T4+T2*T5+T3*T6)/R1/R2
            IF (CTHETA.GT. 1.00D0) CTHETA= 1.00D0
            IF (CTHETA.LT.-1.00D0) CTHETA=-1.00D0
            DUM=ACOS(CTHETA)/C4
            WRITE(6,174)KB(I),MB(I),IB(I),DUM
         ENDDO
      ENDIF
C
      WRITE(6,155)
      WRITE(6,156)(EIG(I),I=1,I3N)
C
C         CALCULATE PRINCIPAL MOMENTS OF INERTIA
C
      FAP(1,1)=AIXX
      FAP(2,1)=-AIXY
      FAP(2,2)=AIYY
      FAP(3,1)=-AIXZ
      FAP(3,2)=-AIYZ
      FAP(3,3)=AIZZ
      I=3
      CALL EIGN(FAP,FA,I,RHO)
      EIG(4)=EIG(1)
      EIG(5)=EIG(2)
      EIG(6)=EIG(3)
C
C             CALCULATE THE TWO MOMENTS OF INERTIA WHICH ARE NOT ABOUT
C             THE Y-COORDINATE, WHICH LIES ALONG THE REACTION PATH.
C
      DO I=1,3
         IF (FA(I,1)**2.GT.(FA(I,2)**2+FA(I,3)**2)) THEN
            DUM1=EIG(I)
         ELSEIF (FA(I,2)**2.GT.(FA(I,1)**2+FA(I,3)**2)) THEN
            DUM2=EIG(I)
         ELSE
            DUM3=EIG(I)
         ENDIF
      ENDDO
      WRITE(6,157)
      WRITE(6,156)DUM1,DUM2,DUM3
C
C         CALCULATE REDUCED MOMENTS OF INERTIA FOR 2-D INTERNAL ROTORS.
C
      IF (N2DR.NE.0) THEN
C
C             REDUCED MOMENTS OF INERTIA FOR REACTANT A.
C
         WT=WTA(1)
         N=NATOMA(1)
         CALL CENMAS(WT,QCM,VCM,N)
         CALL ROTN(AM,EROT,N)
         FAP(1,1)=AIXX
         FAP(2,1)=-AIXY
         FAP(2,2)=AIYY
         FAP(3,1)=-AIXZ
         FAP(3,2)=-AIYZ
         FAP(3,3)=AIZZ
         I=3
         CALL EIGN(FAP,FA,I,RHO)
C
         DO I=1,3
            IF (FA(I,1)**2.GT.(FA(I,2)**2+FA(I,3)**2)) THEN
               DUM4=EIG(I)
            ELSEIF (FA(I,2)**2.GT.(FA(I,1)**2+FA(I,3)**2)) THEN
               DUM5=EIG(I)
            ELSE
               DUM6=EIG(I)
            ENDIF
         ENDDO
         WRITE(6,171)
         WRITE(6,156)DUM4,DUM5,DUM6
         AAI(1)=DUM4*(1.0D0-DUM4/DUM1)
         AAI(2)=DUM6*(1.0D0-DUM6/DUM3)
         WRITE(6,172)AAI(1),AAI(2)
C
C             REDUCED MOMENTS OF INERTIA FOR REACTANT B.
C
         IF (N2DR.NE.1) THEN
            WT=WTB(1)
            N=NATOMB(1)
            DO I=1,N
               LL(I)=LB(1,I)
            ENDDO
            CALL CENMAS(WT,QCM,VCM,N)
            CALL ROTN(AM,EROT,N)
            FAP(1,1)=AIXX
            FAP(2,1)=-AIXY
            FAP(2,2)=AIYY
            FAP(3,1)=-AIXZ
            FAP(3,2)=-AIYZ
            FAP(3,3)=AIZZ
            I=3
            CALL EIGN(FAP,FA,I,RHO)
C
            DO I=1,3
               IF (FA(I,1)**2.GT.(FA(I,2)**2+FA(I,3)**2)) THEN
                  DUM4=EIG(I)
               ELSEIF (FA(I,2)**2.GT.(FA(I,1)**2+FA(I,3)**2)) THEN
                  DUM5=EIG(I)
               ELSE
                  DUM6=EIG(I)
               ENDIF
            ENDDO
            WRITE(6,175)
            WRITE(6,156)DUM4,DUM5,DUM6
            BBI(1)=DUM4*(1.0D0-DUM4/DUM1)
            BBI(2)=DUM6*(1.0D0-DUM6/DUM3)
            WRITE(6,172)BBI(1),BBI(2)
C
C             RESET LL ARRAY
C
            DO I=1,NATOMS
               LL(I)=I
            ENDDO
         ENDIF
      ENDIF
C
C         DETERMINE MAXIMUM AND MINIMUM IN FREE ENERGY ALONG
C         THE REACTION PATH(NONLINEAR GEOMETRY)
C
      VZ=V
      DO I=1,NTEMP
         CALL GPATH(G,VZ,I,NC)
      ENDDO
C
C         SOLVE FOR REACTION PATH
C
   47 DO I=1,I3N
         RDUM=RAA1*DS
         QM(I)=QM(I)+RDUM*GZ(I)
         BQ(I)=GZ(I)
      ENDDO
      call flush(6)
      CALL GRCONV(QM,GZ)
      DO I=1,I3N
         RDUM=RA1*DS
         QM(I)=QM(I)+(GZ(I)-BQ(I))*RDUM
         BQ(I)=(RA2*GZ(I)-RA3*BQ(I))
      ENDDO
      CALL GRCONV(QM,GZ)
      DO I=1,I3N
         RDUM=RB1*DS
         QM(I)=QM(I)+(GZ(I)-BQ(I))*RDUM
         BQ(I)=(RB2*GZ(I)-RB3*BQ(I))
      ENDDO
      CALL GRCONV(QM,GZ)
      DO I=1,I3N
         QM(I)=QM(I)+(GZ(I)*RC1-BQ(I)*RC2)*DS
      ENDDO
      CALL GRCONV(QM,GZ)
C
C         TEST FOR COMPLETION OF INTEGRATION AND NORMAL MODE ANALYSIS
C
      NC=NC+1
      IF (NC.GT.NS) STOP
      IF (NC.NE.1) THEN
         IF (NC.LT.NX) GOTO 47
      ENDIF
      CALL ENERGY
C
C         SAVE GZ, AND EQUATE QQ AND Q.
C
      DO I=1,I3N
         QQ(I)=Q(I)
      ENDDO
      GOTO 100
      END
