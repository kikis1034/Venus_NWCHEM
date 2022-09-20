      SUBROUTINE TETRA(NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CALCULATE TETRAHEDRAL CENTER POTENTIAL ENERGY DERIVATIVES
C
C************************************************************
C                                                           *
C      WRITTEN BY:  R. J. DUCHOVIC                          *
C      DATE:  JULY, 1983                                    *
C      MODIFIED:  MARCH,1986                                *
C                                                           *
C         MODIFICATIONS:  NEW AB INITIO DATA FROM DAVID     *
C                         HIRST (UNIVERSITY OF WARWICK,     *
C                         ENGLAND) HAVE RESULTED IN A NEW   *
C                         FIT FOR THE POTENTIAL PARAMETERS  *
C                         WHICH REQUIRES A CHANGE IN THE    *
C                         ANGULAR SWITCHING FUNCTIONS.  IN  *
C                         ADDITION TO NEW VALUES FOR THE    *
C                         SWITCHING FUNCTION PARAMETERS,    *
C                         A NEW SWITCHING FUNCTION, S4,     *
C                         HAS BEEN ADDED FOR THE CUBIC AND  *
C                         QUARTIC PORTIONS OF THE PHI-TYPE  *
C                         ANGULAR POTENTIAL.  FURTHER, THE  *
C                         EXPONENT IN S1 HAS BEEN MADE A    *
C                         VARIABLE PARAMETER AS IS THE      *
C                         EXPONENT IN THE FUNCTION S4.      *
C                                                           *
C      MODIFIED BY:  XICHE HU                               *
C      DATE:  MARCH 12, 1990                                *
C                                                           *
C         MODIFICATIONS:  THE SWITCHING FUNCTIONS ST AND    *
C                         SP FOR THE THETA AND PHI EQ.      *
C                         ANGLES HAVE BEEN CHANGED.  THE    *
C                         FUNCTIONAL FORMS FOR THE THETA(   *
C                         IJ) QUADRATIC, CUBIC, & QUARTIC   *
C                         FORCE CONSTANTS, AND THE NON-     *
C                         DIAGONAL CUBIC GN4 FORCE CONST.   *
C                         ARE CHANGED.  THESE CHANGES WILL  *
C                         BE PUBLISHED.                     *
C                                                           *
C   THIS SUBROUTINE IS DESIGNED TO COMPUTE DV/DQ FOR THE    *
C   ANGULAR MOTIONS OF A TETRAHEDRAL CENTER.  THE BENDS     *
C   ARE EXPANDED THROUGH QUARTIC TERMS FOR THE DIAGONAL     *
C   FORCES.  THE SUBROUTINE INCLUDES CODE FOR THE OUT-OF-   *
C   PLANE MOTION RESULTING FROM THE RUPTURE OF A SINGLE     *
C   BOND AT THE TETRAHEDRAL CENTER.                         *
C                                                           *
C   DEFINITIONS OF IMPORTANT VARIABLES                      *
C                                                           *
C   ARRAY IQ(15)                                            *
C    -- ARRAY OF INDICES FOR THE COORDINATES OF EACH ATOM   *
C                                                           *
C      IQ(1) => X COORDINATE ATOM I                         *
C      IQ(2) => Y COORDINATE ATOM I                         *
C      IQ(3) => Z COORDINATE ATOM I                         *
C      IQ(4) => X COORDINATE ATOM J                         *
C      IQ(5) => Y COORDINATE ATOM J                         *
C      IQ(6) => Z COORDINATE ATOM J                         *
C      IQ(7) => X COORDINATE ATOM K                         *
C      IQ(8) => Y COORDINATE ATOM K                         *
C      IQ(9) => Z COORDINATE ATOM K                         *
C      IQ(10) => X COORDINATE ATOM L                        *
C      IQ(11) => Y COORDINATE ATOM L                        *
C      IQ(12) => Z COORDINATE ATOM L                        *
C      IQ(13) => X COORDINATE ATOM M                        *
C      IQ(14) => Y COORDINATE ATOM M                        *
C      IQ(15) => Z COORDINATE ATOM M                        *
C      NOTE: ATOM M IS THE TETRAHEDRAL CENTER               *
C                                                           *
C   ARRAY RC(12)                                            *
C    -- ARRAY OF RELATIVE COORDINATES                       *
C      RC(1) => RELATIVE X COORDINATE  XI-XM                *
C      RC(2) => RELATIVE Y COORDINATE  YI-YM                *
C      RC(3) => RELATIVE Z COORDINATE  ZI-ZM                *
C      RC(4) => RELATIVE X COORDINATE  XJ-XM                *
C      RC(5) => RELATIVE Y COORDINATE  YJ-YM                *
C      RC(6) => RELATIVE Z COORDINATE  ZJ-ZM                *
C      RC(7) => RELATIVE X COORDINATE  XK-XM                *
C      RC(8) => RELATIVE Y COORDINATE  YK-YM                *
C      RC(9) => RELATIVE Z COORDINATE  ZK-ZM                *
C      RC(10) => RELATIVE X COORDINATE  XL-XM               *
C      RC(11) => RELATIVE Y COORDINATE  YL-YM               *
C      RC(12) => RELATIVE Z COORDINATE  ZL-ZM               *
C                                                           *
C   ARRAYS XDF(10), YDF(10), ZDF(10)                        *
C    -- ARRAYS OF COORDINATE DIFFERENCES                    *
C                                                           *
C   XDF(1) => XI-XM  YDF(1) => YI-YM  ZDF(1) => ZI-ZM       *
C   XDF(2) => XJ-XM  YDF(2) => YJ-YM  ZDF(2) => ZJ-ZM       *
C   XDF(3) => XJ-XI  YDF(3) => YJ-YI  ZDF(3) => ZJ-ZI       *
C   XDF(4) => XK-XM  YDF(4) => YK-YM  ZDF(4) => ZK-ZM       *
C   XDF(5) => XK-XI  YDF(5) => YK-YI  ZDF(5) => ZK-ZI       *
C   XDF(6) => XK-XJ  YDF(6) => YK-YJ  ZDF(6) => ZK-ZJ       *
C   XDF(7) => XL-XM  YDF(7) => YL-YM  ZDF(7) => ZL-ZM       *
C   XDF(8) => XL-XI  YDF(8) => YL-YI  ZDF(8) => ZL-ZI       *
C   XDF(9) => XL-XJ  YDF(9) => YL-YJ  ZDF(9) => ZL-ZI       *
C   XDF(10) => XL-XK  YDF(10) => YL-YK  ZDF(10) => ZL-ZK    *
C                                                           *
C   ARRAYS TT(NL,6) AND CTT(6)                              *
C   -- T ARRAY HOLDS THETA ANGLES                           *
C   -- CCT ARRAY HOLDS COSINES OF THETA ANGLES              *
C                                                           *
C    TT(NL,1) => THETA(I,J)  CTT(1) => COS(THETA(I,J))      *
C    TT(NL,2) => THETA(I,K)  CTT(2) => COS(THETA(I,K))      *
C    TT(NL,3) => THETA(I,L)  CTT(3) => COS(THETA(I,L))      *
C    TT(NL,4) => THETA(J,K)  CTT(4) => COS(THETA(J,K))      *
C    TT(NL,5) => THETA(J,L)  CTT(5) => COS(THETA(J,L))      *
C    TT(NL,6) => THETA(K,L)  CTT(6) => COS(THEAT(K,L))      *
C                                                           *
C   ARRAY IB(4)                                             *
C   -- ARRAY OF BOND INDICES                                *
C                                                           *
C      IB(1) => INDEX FOR I-M BOND                          *
C      IB(2) => INDEX FOR J-M BOND                          *
C      IB(3) => INDEX FOR K-M BOND                          *
C      IB(4) => INDEX FOR L-M BOND                          *
C                                                           *
C   ARRAY R0(NL,4)                                          *
C   -- ARRAY OF EQUILIBRIUM BOND LENGTHS                    *
C                                                           *
C      R0(NL,1) => BOND I                                   *
C      R0(NL,2) => BOND J                                   *
C      R0(NL,3) => BOND K                                   *
C      R0(NL,4) => BOND L                                   *
C                                                           *
C   ARRAY DLTA(NL,48)                                       *
C   -- ARRAY OF DELTA ANGLES AND HIGHER POWERS OF THESE     *
C      ANGLES                                               *
C                                                           *
C DLTA(NL,1) => DELTA(I,J)    DLTA(NL,25) => DELTA(K,I)     *
C DLTA(NL,2) => DELTA(I,J)**2 DLTA(NL,26) => DELTA(K,I)**2  *
C DLTA(NL,3) => DELTA(I,J)**3 DLTA(NL,27) => DELTA(K,I)**3  *
C DLTA(NL,4) => DELTA(I,J)**4 DLTA(NL,28) => DELTA(K,I)**4  *
C DLTA(NL,5) => DELTA(I,K)    DLTA(NL,29) => DELTA(K,J)     *
C DLTA(NL,6) => DELTA(I,K)**2 DLTA(NL,30) => DELTA(K,J)**2  *
C DLTA(NL,7) => DELTA(I,K)**3 DLTA(NL,31) => DELTA(K,J)**3  *
C DLTA(NL,8) => DELTA(I,K)**4 DLTA(NL,32) => DELTA(K,J)**4  *
C DLTA(NL,9) => DELTA(I,L)    DLTA(NL,33) => DELTA(K,L)     *
C DLTA(NL,10) => DELTA(I,L)**2 DLTA(NL,34) => DELTA(K,L)**2 *
C DLTA(NL,11) => DELTA(I,L)**3 DLTA(NL,35) => DELTA(K,L)**3 *
C DLTA(NL,12) => DELTA(I,L)**4 DLTA(NL,36) => DELTA(K,L)**4 *
C DLTA(NL,13) => DELTA(J,I)    DLTA(NL,37) => DELTA(L,I)    *
C DLTA(NL,14) => DELTA(J,I)**2 DLTA(NL,38) => DELTA(L,I)**2 *
C DLTA(NL,15) => DELTA(J,I)**3 DLTA(NL,39) => DELTA(L,I)**3 *
C DLTA(NL,16) => DELTA(J,I)**4 DLTA(NL,40) => DELTA(L,I)**4 *
C DLTA(NL,17) => DELTA(J,K)    DLTA(NL,41) => DELTA(L,J)    *
C DLTA(NL,18) => DELTA(J,K)**2 DLTA(NL,42) => DELTA(L,J)**2 *
C DLTA(NL,19) => DELTA(J,K)**3 DLTA(NL,43) => DELTA(L,J)**3 *
C DLTA(NL,20) => DELTA(J,K)**4 DLTA(NL,44) => DELTA(L,J)**4 *
C DLTA(NL,21) => DELTA(J,L)    DLTA(NL,45) => DELTA(L,K)    *
C DLTA(NL,22) => DELTA(J,L)**2 DLTA(NL,46) => DELTA(L,K)**2 *
C DLTA(NL,23) => DELTA(J,L)**3 DLTA(NL,47) => DELTA(L,K)**3 *
C DLTA(NL,24) => DELTA(J,L)**4 DLTA(NL,48) => DELTA(L,K)**4 *
C                                                           *
C   ARRAY CDLTA(12)                                         *
C   -- ARRAY OF COSINES FOR OUT-OF-PLANE ANGLES             *
C                                                           *
C   CDLTA(1) => COS(DELTA(I,J)+THETAZERO(I,J))              *
C   CDLTA(2) => COS(DELTA(I,K)+THETAZERO(I,K))              *
C   CDLTA(3) => COS(DELTA(I,L)+THETAZERO(I,L))              *
C   CDLTA(4) => COS(DELTA(J,I)+THETAZERO(I,J))              *
C   CDLTA(5) => COS(DELTA(J,K)+THETAZERO(J,K))              *
C   CDLTA(6) => COS(DELTA(J,L)+THETAZERO(J,L))              *
C   CDLTA(7) => COS(DELTA(K,I)+THETAZERO(I,K))              *
C   CDLTA(8) => COS(DELTA(K,J)+THETAZERO(J,K))              *
C   CDLTA(9) => COS(DELTA(K,L)+THETAZERO(K,L))              *
C   CDLTA(10) => COS(DELTA(L,I)+THETAZERO(I,L))             *
C   CDLTA(11) => COS(DELTA(L,J)+THETAZERO(J,L))             *
C   CDLTA(12) => COS(DELTA(L,K)+THETAZERO(K,L))             *
C                                                           *
C   SPECIAL INDEXING ARRAYS J1(6),J2(6),J3(6),J4(6)         *
C                                                           *
C   J1(1)=1   J2(1)=2   J3(1)=3   J4(1)=4                   *
C   J1(2)=1   J2(2)=3   J3(2)=2   J4(2)=4                   *
C   J1(3)=1   J2(3)=4   J3(3)=2   J4(3)=3                   *
C   J1(4)=2   J2(4)=3   J3(4)=1   J4(4)=4                   *
C   J1(5)=2   J2(5)=4   J3(5)=1   J4(5)=3                   *
C   J1(6)=3   J2(6)=4   J3(6)=1   J4(6)=2                   *
C                                                           *
C   ADDITIONAL SPECIAL INDEXING ARRAYS:  K1(4),K2(4),K3(4)  *
C      K4(4),K5(4),K6(4),K7(4),K8(4),K9(4)                  *
C                                                           *
C   K1(1)=2  K2(1)=3  K3(1)=4  K4(1)=1  K5(1)=2             *
C   K1(2)=1  K2(2)=3  K3(2)=4  K4(2)=1  K5(2)=4             *
C   K1(3)=1  K2(3)=2  K3(3)=4  K4(3)=2  K5(3)=4             *
C   K1(4)=1  K2(4)=2  K3(4)=3  K4(4)=3  K5(4)=5             *
C                                                           *
C   K6(1)=3  K7(1)=4  K8(1)=5  K9(1)=6                      *
C   K6(2)=5  K7(2)=2  K8(2)=3  K9(2)=6                      *
C   K6(3)=6  K7(3)=1  K8(3)=3  K9(3)=5                      *
C   K6(4)=6  K7(4)=1  K8(4)=2  K9(4)=4                      *
C                                                           *
C   ARRAY CBIC(15,6)                                        *
C   -- ARRAY USED TO SAVE THE DERIVATIVES OF THE THETA      *
C      ANGLES FOR USE IN THE NON-DIAGONAL CUBIC             *
C      SUBROUTINE CUBEND(NL,N1,N2,D).  EACH ROW OF THE      *
C      ARRAY STORES DATA FOR A CARTESIAN COORDINATE         *
C      (A TOTAL OF 15) WHILE THE COLUMNS ARE THE            *
C      DIFFERENCES (DTHETA/DXI - DTHETAZERO/DXI) FOR THE    *
C      SIX THETA ANGLES IN METHANE.                         *
C                                                           *
C   ARRAY OUTPL(15,6)                                       *
C   -- ARRAY USED TO SAVE THE DERIVATIVES OF THE            *
C      THETAZERO ANGLES FOR USE IN THE OUT-OF-PLANE BEND    *
C      CALCULATION.  EACH ROW OF THE ARRAY STORES DATA      *
C      FOR A CARTESIAN COORDINATE (A TOTAL OF 15) AND       *
C      THE QUANTITIES STORED ARE: DTHETAZERO/DXI            *
C                                                           *
C   ARRAY DERV(5,4)                                         *
C   -- ARRAY USED TO HOLD THE DERIVATITVES WITH RESEPECT    *
C      TO THE CARTESIAN COORDINATES OF THE OUT-OF-PLANE     *
C      SWITCHING FUNCTION                                   *
C                                                           *
C************************************************************
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/TETRAB/ N9I(20),N9J(20),N9K(20),N9L(20),N9M(20),
     *               FT0(20,6),FT2(20,6),GT0(20,6),GT2(20,6),
     *               HT0(20,6),HT2(20,6),THT(20,6),R0(20,4),
     *               THT1(20,6),THT2(20,6),FD1(20,4),
     *               HD1(20,4),GN0(20,5),FT(20,6),GT(20,6),
     *               HT(20,6),FD(20,4),HD(20,4),DLTA(20,48),
     *               TETTST,SGN1,SGN2,SGN3,SGN4
      COMMON/CUBEB/ S3(4),DS3(4),CBIC(15,6),ANG1(20,6,4),GN4(20)
      DIMENSION IB(4),IQ(15),J1(6),J2(6),J3(6),J4(6),K1(4),
     *          K2(4),K3(4),K4(4),K5(4),K6(4),K7(4),K8(4),
     *          K9(4),RC(12),CDLTA(12),CTT(6),DERV(5,4),
     *          TN(6),S1(4),S2(4),SP(4),ST(4),SWF(6,2),
     *          SWG(6,2),SWH(6,2),EP(4),ET(4),DS1(4),DS2(4),
     *          DSP(4),DST(4),DF1(4),DF2(4),DF3(4),DF4(6),
     *          DFN1(15),DFN2(15),DFN3(15),DFN4(15),
     *          TRMA(2),TRMB(2),TRMC(2),TRMD(2),
     *          DA1(15),DA2(15),DB1(15),DB2(15),
     *          DC1(15),DC2(15),DD1(15),DD2(15),OUTPL(15,6),
     *          DANG1(12),XDF(10),YDF(10),ZDF(10),BNG(6),
     *          CNG(6),FTD(6),GTD(6),HTD(6),DUMA(6),DUMT(6)
     *          ,S4(4),DS4(4)
      DATA J1/1,1,1,2,2,3/,J2/2,3,4,3,4,4/
      DATA J3/3,2,2,1,1,1/,J4/4,4,3,4,3,2/
      DATA K1/2,1,1,1/,K2/3,3,2,2/,K3/4,4,4,3/
      DATA K4/1,1,2,3/,K5/2,4,4,5/,K6/3,5,6,6/
      DATA K7/4,2,1,1/,K8/5,3,3,2/,K9/6,6,5,4/
C***********************************************************
C                                                          *
C   DEFINE CONSTANTS FOR SWITCHING FUNCTIONS               *
C                                                          *
C***********************************************************
      A1=0.5045375D0
      B1=0.4191092D0
      C1=0.6989260D0
      A2=1.0147402D-07
      B2=-1.2362798D01
      C2=6.0D0
      A3=1.4191474D-01
      B3=-3.0684503D-01
      C3=2.0D0
      A4=0.38387689D0
      B4=-0.169915274D0
      C4=0.971186697D0
      AP=0.224147141D0
      BP=-0.990736407D-05
      CP=0.238179398D+02
      AT=0.330879271D0
      BT=-0.12408758D-02
      CT=0.880601978D+01
C***********************************************************
C                                                          *
C   ZERO OUT CUBIC ARRAY FOR LATER USE                     *
C                                                          *
C***********************************************************
      DO I=1,15
        DO J=1,6
          CBIC(I,J)=0.0D0
          OUTPL(I,J)=0.0D0
        ENDDO
      ENDDO
C***********************************************************
C                                                          *
C   CALCULATE INDICES FOR COORDINATES                      *
C                                                          *
C***********************************************************
      IQ(3)=3*N9I(NL)
      IQ(2)=IQ(3)-1
      IQ(1)=IQ(2)-1
      IQ(6)=3*N9J(NL)
      IQ(5)=IQ(6)-1
      IQ(4)=IQ(5)-1
      IQ(9)=3*N9K(NL)
      IQ(8)=IQ(9)-1
      IQ(7)=IQ(8)-1
      IQ(12)=3*N9L(NL)
      IQ(11)=IQ(12)-1
      IQ(10)=IQ(11)-1
      IQ(15)=3*N9M(NL)
      IQ(14)=IQ(15)-1
      IQ(13)=IQ(14)-1
C***********************************************************
C                                                          *
C   CALCULATE INDICES FOR R'S                              *
C                                                          *
C***********************************************************
      NUM1=(N9M(NL)-1)*(2*NATOMS-N9M(NL))/2
      IB(1)=NUM1+N9I(NL)-N9M(NL)
      IF(N9M(NL).GT.N9I(NL))
     +IB(1)=(N9I(NL)-1)*(2*NATOMS-N9I(NL))/2+N9M(NL)-N9I(NL)
      IB(2)=NUM1+N9J(NL)-N9M(NL)
      IF(N9M(NL).GT.N9J(NL))
     +IB(2)=(N9J(NL)-1)*(2*NATOMS-N9J(NL))/2+N9M(NL)-N9J(NL)
      IB(3)=NUM1+N9K(NL)-N9M(NL)
      IF(N9M(NL).GT.N9K(NL))
     +IB(3)=(N9K(NL)-1)*(2*NATOMS-N9K(NL))/2+N9M(NL)-N9K(NL)
      IB(4)=NUM1+N9L(NL)-N9M(NL)
      IF(N9M(NL).GT.N9L(NL))
     +IB(4)=(N9L(NL)-1)*(2*NATOMS-N9L(NL))/2+N9M(NL)-N9L(NL)
C***********************************************************
C                                                          *
C   CALCULATE RELATIVE COORDINATES                         *
C                                                          *
C***********************************************************
      DO I=1,4
         J=3*I-2
         RC(J)=Q(IQ(J))-Q(IQ(13))
         RC(J+1)=Q(IQ(J+1))-Q(IQ(14))
         RC(J+2)=Q(IQ(J+2))-Q(IQ(15))
      ENDDO
C***********************************************************
C                                                          *
C   CALCULATE VALUES OF THE SWITCHING FUNCTIONS            *
C      CALCULATION OF SP AND ST MODIFIED                   *
C      DECEMBER 1985                                       *
C                                                          *
C***********************************************************
      DO I=1,4
         S1(I)=1.0D0-TANH(A1*(R(IB(I))-R0(NL,I))*(R(IB(I))-B1)**C1)
         S2(I)=1.0D0-TANH(A2*(R(IB(I))-R0(NL,I))*(R(IB(I))-B2)**C2)
         S3(I)=1.0D0-TANH(A3*(R(IB(I))-R0(NL,I))*(R(IB(I))-B3)**C3)
         S4(I)=1.0D0-TANH(A4*(R(IB(I))-R0(NL,I))*(R(IB(I))-B4)**C4)
         EP(I)=5.54D+34
         SP(I)=0.0D0
         DUM1=BP*(R(IB(I))-CP)**3
         IF (DUM1.LE.80.0D0) THEN
            EP(I)=EXP(DUM1)
            SP(I)=1.0D0-TANH(AP*(R(IB(I))-R0(NL,I))*(1+EP(I)))
         ENDIF
         ET(I)=5.54D+34
         ST(I)=0.0D0
         DUM1=BT*(R(IB(I))-CT)**3
         IF (DUM1.LE.80.0D0) THEN
            ET(I)=EXP(DUM1)
            ST(I)=1.0D0-TANH(AT*(R(IB(I))-R0(NL,I))*(1+ET(I)))
         ENDIF
      ENDDO
C***********************************************************
C                                                          *
C   CALCULATE THETA(I,J) AND COS(THETA(I,J)) FOR ALL I,J   *
C                                                          *
C***********************************************************
      DUM1=1.0D0/R(IB(1))
      DUM2=1.0D0/R(IB(2))
      DUM3=1.0D0/R(IB(3))
      DUM4=1.0D0/R(IB(4))
      CTT(1)=(RC(1)*RC(4)+RC(2)*RC(5)+RC(3)*RC(6))*DUM1*DUM2
      CTT(2)=(RC(1)*RC(7)+RC(2)*RC(8)+RC(3)*RC(9))*DUM1*DUM3
      CTT(3)=(RC(1)*RC(10)+RC(2)*RC(11)+RC(3)*RC(12))*DUM1*DUM4
      CTT(4)=(RC(4)*RC(7)+RC(5)*RC(8)+RC(6)*RC(9))*DUM2*DUM3
      CTT(5)=(RC(4)*RC(10)+RC(5)*RC(11)+RC(6)*RC(12))*DUM2*DUM4
      CTT(6)=(RC(7)*RC(10)+RC(8)*RC(11)+RC(9)*RC(12))*DUM3*DUM4
      DO I=1,6
         TT(NL,I)=ACOS(CTT(I))
      ENDDO
C***********************************************************
C                                                          *
C   CALCULATE THETAZERO(I,J) FOR ALL I,J                   *
C                                                          *
C***********************************************************
      DO I=1,6
         TN(I)=THT(NL,I)+(THT(NL,I)-THT1(NL,I))
     *         *(SP(J1(I))*SP(J2(I))-1.0D0)+(THT(NL,I)-THT2(NL,I))
     *         *(ST(J3(I))*ST(J4(I))-1.0D0)
      ENDDO
C***********************************************************
C                                                          *
C   NOW DEFINE SWITCHING FUNCTION ARRAYS NEEDED FOR THE    *
C   CALCULATION OF THE FORCE CONSTANTS                     *
C                                                          *
C***********************************************************
      DO I=1,6
         SWF(I,1)=S1(J1(I))*S1(J2(I))-1.0D0
         SWF(I,2)=S2(J3(I))*S2(J4(I))-1.0D0
         SWG(I,1)=S4(J1(I))*S4(J2(I))-1.0D0
         SWG(I,2)=S2(J3(I))*S2(J4(I))-1.0D0
         SWH(I,1)=SWG(I,1)
         SWH(I,2)=SWG(I,2)
      ENDDO
C***********************************************************
C                                                          *
C   CALCULATE THE DIAGONAL QUADRATIC, CUBIC, AND QUARTIC   *
C   FORCE CONSTANTS                                        *
C                                                          *
C***********************************************************
      DO I=1,6
         FT(NL,I)=FT0(NL,I)+FT0(NL,I)*SWF(I,1)
     *           +(FT0(NL,I)-FT2(NL,I))*SWF(I,2)*(SWF(I,1)+1.0D0)
         GT(NL,I)=GT0(NL,I)+GT0(NL,I)*SWG(I,1)
     *           +(GT0(NL,I)-GT2(NL,I))*SWG(I,2)*(SWG(I,1)+1.0D0)
         HT(NL,I)=HT0(NL,I)+HT0(NL,I)*SWH(I,1)
     *           +(HT0(NL,I)-HT2(NL,I))*SWH(I,2)*(SWH(I,1)+1.0D0)
      ENDDO
C***********************************************************
C                                                          *
C   WE ARE NOW READY TO CALCULATE DV/DQ FOR THE DIAGONAL   *
C   QUADRATIC, CUBIC, AND QUARTIC ANGULAR TERMS OF THE     *
C   POTENTIAL.                                             *
C                                                          *
C***********************************************************
      DO I=1,6
         ANG1(NL,I,1)=TT(NL,I)-TN(I)
         ANG1(NL,I,2)=ANG1(NL,I,1)*ANG1(NL,I,1)
         ANG1(NL,I,3)=ANG1(NL,I,2)*ANG1(NL,I,1)
         ANG1(NL,I,4)=ANG1(NL,I,3)*ANG1(NL,I,1)
         BNG(I)=THT(NL,I)-THT1(NL,I)
         CNG(I)=THT(NL,I)-THT2(NL,I)
         FTD(I)=FT0(NL,I)-FT2(NL,I)
         GTD(I)=GT0(NL,I)-GT2(NL,I)
         HTD(I)=HT0(NL,I)-HT2(NL,I)
      ENDDO
C***********************************************************
C                                                          *
C   BEGIN MAIN LOOP FOR DV/DQ CALCULATION FOR QUADRATIC,   *
C   CUBIC, AND QUARTIC ANGULAR TERMS                       *
C                                                          *
C***********************************************************
      DO I=1,3
         DO K=1,4
            J=I+3*(K-1)
            D1=1.0D0/R(IB(K))
            D2=R0(NL,K)*D1
            TERM1=(1.0D0-S1(K))**2-1.0D0
            TERM2=A1*RC(J)*D1*(R(IB(K))-B1)**C1
            TERM3=C1*A1*RC(J)*(R(IB(K))-B1)**(C1-1.0D0)*(1.0D0-D2)
            DS1(K)=TERM1*(TERM2+TERM3)
            TERM1=(1.0D0-S2(K))**2-1.0D0
            TERM2=A2*RC(J)*D1*(R(IB(K))-B2)**C2
            TERM3=C2*A2*RC(J)*(R(IB(K))-B2)**(C2-1.0D0)*(1.0D0-D2)
            DS2(K)=TERM1*(TERM2+TERM3)
            TERM1=(1.0D0-S4(K))**2-1.0D0
            TERM2=A4*RC(J)*D1*(R(IB(K))-B4)**C4
            TERM3=C4*A4*RC(J)*(R(IB(K))-B4)**(C4-1.0D0)*(1.0D0-D2)
            DS4(K)=TERM1*(TERM2+TERM3)
            TERM1=(1.0D0-SP(K))**2-1.0D0
            TERM2=AP*RC(J)*EP(K)
            TERM3=D1+3.0D0*BP*(R(IB(K))-CP)**2*(1.0D0-D2)
            DSP(K)=TERM1*(TERM2*TERM3+AP*RC(J)/R(IB(K)))
            TERM1=(1.0D0-ST(K))**2-1.0D0
            TERM2=AT*RC(J)*ET(K)
            TERM3=D1+3.0D0*BT*(R(IB(K))-CT)**2*(1.0D0-D2)
            DST(K)=TERM1*(TERM2*TERM3+AT*RC(J)/R(IB(K)))
         ENDDO
C***********************************************************
C                                                          *
C   WE HAVE CALCULATED THE DERIVATIVES OF THE SWITCHING    *
C   FUNCTIONS S1, S2, SP, ST.  WE NOW CALCULATE THE        *
C   FRACTIONAL TERMS IN EACH DERIVATIVE WITH RESPECT TO    *
C   ATOMS I, J, K, L.                                      *
C                                                          *
C***********************************************************
         R11=1.0D0/R(IB(1))
         R12=1.0D0/R(IB(2))
         R13=1.0D0/R(IB(3))
         R14=1.0D0/R(IB(4))
         S11=-1.0D0/SQRT(1.0D0-CTT(1)**2)
         S12=-1.0D0/SQRT(1.0D0-CTT(2)**2)
         S13=-1.0D0/SQRT(1.0D0-CTT(3)**2)
         S14=-1.0D0/SQRT(1.0D0-CTT(4)**2)
         S15=-1.0D0/SQRT(1.0D0-CTT(5)**2)
         S16=-1.0D0/SQRT(1.0D0-CTT(6)**2)
         R21=R11**2
         FRC1=RC(I+3)*R11*R12
         FRC2=RC(I)*CTT(1)*R21
         DF1(1)=S11*(FRC1-FRC2)
         FRC1=RC(I+6)*R11*R13
         FRC2=RC(I)*CTT(2)*R21
         DF2(1)=S12*(FRC1-FRC2)
         FRC1=RC(I+9)*R11*R14
         FRC2=RC(I)*CTT(3)*R21
         DF3(1)=S13*(FRC1-FRC2)
         R22=R12**2
         FRC1=RC(I)*R11*R12
         FRC2=RC(I+3)*CTT(1)*R22
         DF1(2)=S11*(FRC1-FRC2)
         FRC1=RC(I+6)*R12*R13
         FRC2=RC(I+3)*CTT(4)*R22
         DF2(2)=S14*(FRC1-FRC2)
         FRC1=RC(I+9)*R12*R14
         FRC2=RC(I+3)*CTT(5)*R22
         DF3(2)=S15*(FRC1-FRC2)
         R23=R13**2
         FRC1=RC(I)*R11*R13
         FRC2=RC(I+6)*CTT(2)*R23
         DF1(3)=S12*(FRC1-FRC2)
         FRC1=RC(I+3)*R12*R13
         FRC2=RC(I+6)*CTT(4)*R23
         DF2(3)=S14*(FRC1-FRC2)
         FRC1=RC(I+9)*R13*R14
         FRC2=RC(I+6)*CTT(6)*R23
         DF3(3)=S16*(FRC1-FRC2)
         R24=R14**2
         FRC1=RC(I)*R11*R14
         FRC2=RC(I+9)*CTT(3)*R24
         DF1(4)=S13*(FRC1-FRC2)
         FRC1=RC(I+3)*R12*R14
         FRC2=RC(I+9)*CTT(5)*R24
         DF2(4)=S15*(FRC1-FRC2)
         FRC1=RC(I+6)*R13*R14
         FRC2=RC(I+9)*CTT(6)*R24
         DF3(4)=S16*(FRC1-FRC2)
C***********************************************************
C                                                          *
C   AT THIS POINT WE ONLY NEED TO CALCULATE THE FRACTIONAL *
C   TERMS IN THE DERIVATIVE WITH RESPECT TO ATOM M         *
C   (THE TETRAHEDRAL CENTER)                               *
C                                                          *
C***********************************************************
         Q1=2.0D0*Q(IQ(I+12))
         FRC1=(Q1-Q(IQ(I))-Q(IQ(I+3)))*R11*R12
         FRC2=RC(I)*CTT(1)*R21
         FRC3=RC(I+3)*CTT(1)*R22
         DF4(1)=S11*(FRC1+FRC2+FRC3)
         FRC1=(Q1-Q(IQ(I))-Q(IQ(I+6)))*R11*R13
         FRC2=RC(I)*CTT(2)*R21
         FRC3=RC(I+6)*CTT(2)*R23
         DF4(2)=S12*(FRC1+FRC2+FRC3)
         FRC1=(Q1-Q(IQ(I))-Q(IQ(I+9)))*R11*R14
         FRC2=RC(I)*CTT(3)*R21
         FRC3=RC(I+9)*CTT(3)*R24
         DF4(3)=S13*(FRC1+FRC2+FRC3)
         FRC1=(Q1-Q(IQ(I+3))-Q(IQ(I+6)))*R12*R13
         FRC2=RC(I+3)*CTT(4)*R22
         FRC3=RC(I+6)*CTT(4)*R23
         DF4(4)=S14*(FRC1+FRC2+FRC3)
         FRC1=(Q1-Q(IQ(I+3))-Q(IQ(I+9)))*R12*R14
         FRC2=RC(I+3)*CTT(5)*R22
         FRC3=RC(I+9)*CTT(5)*R24
         DF4(5)=S15*(FRC1+FRC2+FRC3)
         FRC1=(Q1-Q(IQ(I+6))-Q(IQ(I+9)))*R13*R14
         FRC2=RC(I+6)*CTT(6)*R23
         FRC3=RC(I+9)*CTT(6)*R24
         DF4(6)=S16*(FRC1+FRC2+FRC3)
C***********************************************************
C                                                          *
C   WE ARE NOW READY TO EVALUATE DV/DQ FOR THE QUADRATIC,  *
C   CUBIC, AND QUARTIC TERMS OF THE POTENTIAL.             *
C                                                          *
C***********************************************************
         DO J=1,4
            L=I+3*(J-1)
C***********************************************************
C                                                          *
C   QUADRATIC TERMS  ATOMS I,J,K,L                         *
C   INDEX J REFERS TO THE CARTESIAN COORDINATES OF ATOMS   *
C   I,J,K,L (XI,XJ,XK,XL, ETC)                             *
C   INDEX L ADDRESSES THE APPROPRIATE LOCATION IN          *
C   THE STORAGE ARRAYS                                     *
C                                                          *
C***********************************************************
            TEMP1=FT0(NL,K4(J))+FTD(K4(J))*(S2(K3(J))*S2(K2(J))-1.0D0)
            TEMP2=FT0(NL,K5(J))+FTD(K5(J))*(S2(K3(J))*S2(K1(J))-1.0D0)
            TEMP3=FT0(NL,K6(J))+FTD(K6(J))*(S2(K1(J))*S2(K2(J))-1.0D0)
            DUM1=TEMP1*S1(K1(J))*ANG1(NL,K4(J),2)
     *           +TEMP2*S1(K2(J))*ANG1(NL,K5(J),2)
     *           +TEMP3*S1(K3(J))*ANG1(NL,K6(J),2)
            DUM2=FTD(K7(J))*S2(K3(J))*ANG1(NL,K7(J),2)
     *                     *S1(K1(J))*S1(K2(J))
     *          +FTD(K8(J))*S2(K2(J))*ANG1(NL,K8(J),2)
     *                     *S1(K1(J))*S1(K3(J))
     *          +FTD(K9(J))*S2(K1(J))*ANG1(NL,K9(J),2)
     *                     *S1(K2(J))*S1(K3(J))
            TEMP1=BNG(K4(J))
            TEMP2=BNG(K5(J))
            TEMP3=BNG(K6(J))
            DUM3=FT(NL,K4(J))*ANG1(NL,K4(J),1)*TEMP1*SP(K1(J))
     *           +FT(NL,K5(J))*ANG1(NL,K5(J),1)*TEMP2*SP(K2(J))
     *           +FT(NL,K6(J))*ANG1(NL,K6(J),1)*TEMP3*SP(K3(J))
            DUM4=FT(NL,K7(J))*ANG1(NL,K7(J),1)*CNG(K7(J))*ST(K3(J))
     *           +FT(NL,K8(J))*ANG1(NL,K8(J),1)
     *                        *CNG(K8(J))*ST(K2(J))
     *           +FT(NL,K9(J))*ANG1(NL,K9(J),1)
     *                        *CNG(K9(J))*ST(K1(J))
            PDOT(IQ(L))=PDOT(IQ(L))+0.5D0*DUM1*DS1(J)
     *                  +0.5D0*DUM2*DS2(J)
     *                  +FT(NL,K4(J))*ANG1(NL,K4(J),1)*DF1(J)
     *                  +FT(NL,K5(J))*ANG1(NL,K5(J),1)*DF2(J)
     *                  +FT(NL,K6(J))*ANG1(NL,K6(J),1)*DF3(J)
     *                  -DUM3*DSP(J)-DUM4*DST(J)
C***********************************************************
C                                                          *
C   CODE NEEDED TO SAVE DATA FOR OUT-OF-PLANE BEND         *
C   CALCULATION                                            *
C                                                          *
C***********************************************************
            OUTPL(L,K4(J))=TEMP1*SP(K1(J))*DSP(J)
            OUTPL(L,K5(J))=TEMP2*SP(K2(J))*DSP(J)
            OUTPL(L,K6(J))=TEMP3*SP(K3(J))*DSP(J)
            OUTPL(L,K7(J))=CNG(K7(J))*ST(K3(J))*DST(J)
            OUTPL(L,K8(J))=CNG(K8(J))*ST(K2(J))*DST(J)
            OUTPL(L,K9(J))=CNG(K9(J))*ST(K1(J))*DST(J)
C***********************************************************
C                                                          *
C   ADDITIONAL CODE REQUIRED TO SAVE ANGULAR DERIVATIVES   *
C   FOR THE NON-DIAGONAL CUBIC CALCULATION                 *
C   SAVE THESE QUANTITIES IN THE ARRAY CBIC                *
C                                                          *
C***********************************************************
C
            CBIC(L,K4(J))=DF1(J)-OUTPL(L,K4(J))
            CBIC(L,K5(J))=DF2(J)-OUTPL(L,K5(J))
            CBIC(L,K6(J))=DF3(J)-OUTPL(L,K6(J))
            CBIC(L,K7(J))=-OUTPL(L,K7(J))
            CBIC(L,K8(J))=-OUTPL(L,K8(J))
            CBIC(L,K9(J))=-OUTPL(L,K9(J))
C***********************************************************
C                                                          *
C   CUBIC TERMS  ATOMS I, J, K, L                          *
C                                                          *
C***********************************************************
            TEMP31=GT0(NL,K4(J))+GTD(K4(J))*(S2(K3(J))*S2(K2(J))-1.0D0)
            TEMP32=GT0(NL,K5(J))+GTD(K5(J))*(S2(K3(J))*S2(K1(J))-1.0D0)
            TEMP33=GT0(NL,K6(J))+GTD(K6(J))*(S2(K1(J))*S2(K2(J))-1.0D0)
            DUM1=TEMP31*S4(K1(J))*ANG1(NL,K4(J),3)
     *           +TEMP32*S4(K2(J))*ANG1(NL,K5(J),3)
     *           +TEMP33*S4(K3(J))*ANG1(NL,K6(J),3)
            DUM2=GTD(K7(J))*S2(K3(J))*ANG1(NL,K7(J),3)
     *                     *S4(K1(J))*S4(K2(J))
     *           +GTD(K8(J))*S2(K2(J))*ANG1(NL,K8(J),3)
     *                     *S4(K1(J))*S4(K3(J))
     *           +GTD(K9(J))*S2(K1(J))*ANG1(NL,K9(J),3)
     *                     *S4(K2(J))*S4(K3(J))
            DUM3=GT(NL,K4(J))*ANG1(NL,K4(J),2)*TEMP1*SP(K1(J))
     *           +GT(NL,K5(J))*ANG1(NL,K5(J),2)*TEMP2*SP(K2(J))
     *           +GT(NL,K6(J))*ANG1(NL,K6(J),2)*TEMP3*SP(K3(J))
            DUM4=GT(NL,K7(J))*ANG1(NL,K7(J),2)*CNG(K7(J))*ST(K3(J))
     *           +GT(NL,K8(J))*ANG1(NL,K8(J),2)*CNG(K8(J))*ST(K2(J))
     *           +GT(NL,K9(J))*ANG1(NL,K9(J),2)*CNG(K9(J))*ST(K1(J))
            PDOT(IQ(L))=PDOT(IQ(L))+DUM1*DS4(J)+DUM2*DS2(J)
     *                  +3.0D0*(GT(NL,K4(J))*ANG1(NL,K4(J),2)*DF1(J)
     *                  +GT(NL,K5(J))*ANG1(NL,K5(J),2)*DF2(J)
     *                  +GT(NL,K6(J))*ANG1(NL,K6(J),2)*DF3(J)
     *                  -DUM3*DSP(J)-DUM4*DST(J))
C***********************************************************
C                                                          *
C   QUARTIC TERMS  ATOMS I, J, K, L                        *
C                                                          *
C***********************************************************
            TEMP41=HT0(NL,K4(J))+HTD(K4(J))*(S2(K3(J))*S2(K2(J))-1.0D0)
            TEMP42=HT0(NL,K5(J))+HTD(K5(J))*(S2(K3(J))*S2(K1(J))-1.0D0)
            TEMP43=HT0(NL,K6(J))+HTD(K6(J))*(S2(K1(J))*S2(K2(J))-1.0D0)
            DUM1=TEMP41*S4(K1(J))*ANG1(NL,K4(J),4)
     *           +TEMP42*S4(K2(J))*ANG1(NL,K5(J),4)
     *           +TEMP43*S4(K3(J))*ANG1(NL,K6(J),4)
            DUM2=HTD(K7(J))*S2(K3(J))*ANG1(NL,K7(J),4)
     *                     *S4(K1(J))*S4(K2(J))
     *           +HTD(K8(J))*S2(K2(J))*ANG1(NL,K8(J),4)
     *                     *S4(K1(J))*S4(K3(J))
     *           +HTD(K9(J))*S2(K1(J))*ANG1(NL,K9(J),4)
     *                     *S4(K2(J))*S4(K3(J))
            DUM3=HT(NL,K4(J))*ANG1(NL,K4(J),3)*TEMP1*SP(K1(J))
     *           +HT(NL,K5(J))*ANG1(NL,K5(J),3)*TEMP2*SP(K2(J))
     *           +HT(NL,K6(J))*ANG1(NL,K6(J),3)*TEMP3*SP(K3(J))
            DUM4=HT(NL,K7(J))*ANG1(NL,K7(J),3)*CNG(K7(J))*ST(K3(J))
     *           +HT(NL,K8(J))*ANG1(NL,K8(J),3)*CNG(K8(J))*ST(K2(J))
     *           +HT(NL,K9(J))*ANG1(NL,K9(J),3)*CNG(K9(J))*ST(K1(J))
            PDOT(IQ(L))=PDOT(IQ(L))+DUM1*DS4(J)+DUM2*DS2(J)
     *                  +4.0D0*(HT(NL,K4(J))*ANG1(NL,K4(J),3)*DF1(J)
     *                  +HT(NL,K5(J))*ANG1(NL,K5(J),3)*DF2(J)
     *                  +HT(NL,K6(J))*ANG1(NL,K6(J),3)*DF3(J)
     *                  -DUM3*DSP(J)-DUM4*DST(J))
         ENDDO
C***********************************************************
C                                                          *
C   NOW COMPUTE THE DERIVATIVE OF THE QUADRATIC TERM       *
C   WITH RESPECT TO ATOM M (THE TETRAHEDRAL CENTER)        *
C                                                          *
C***********************************************************
         DO K=1,6
            DUMA(K)=-((DS1(J1(K))*S1(J2(K))+S1(J1(K))*DS1(J2(K)))
     *              *(FT0(NL,K)+FTD(K)*(S2(J3(K))*S2(J4(K))-1.0D0))
     *              +(DS2(J3(K))*S2(J4(K))+S2(J3(K))*DS2(J4(K)))
     *              *FTD(K)*S1(J1(K))*S1(J2(K)))
            DUMT(K)=-((DSP(J1(K))*SP(J2(K))+SP(J1(K))*DSP(J2(K)))
     *              *BNG(K)
     *              +(DST(J3(K))*ST(J4(K))+ST(J3(K))*DST(J4(K)))
     *              *CNG(K))
         ENDDO
         DO K=1,6
            CBIC(I+12,K)=DF4(K)-DUMT(K)
            OUTPL(I+12,K)=DUMT(K)
            PDOT(IQ(I+12))=PDOT(IQ(I+12))+0.5D0*DUMA(K)*ANG1(NL,K,2)
     *                    +FT(NL,K)*ANG1(NL,K,1)*CBIC(I+12,K)
         ENDDO
C***********************************************************
C                                                          *
C   COMPUTE THE DERIVATIVE OF THE CUBIC TERM WITH          *
C   RESPECT TO ATOM M                                      *
C                                                          *
C***********************************************************
         DO K=1,6
            DUMA(K)=-((DS4(J1(K))*S4(J2(K))+S4(J1(K))*DS4(J2(K)))
     *              *(GT0(NL,K)+GTD(K)*(S2(J3(K))*S2(J4(K))-1.0D0))
     *              +(DS2(J3(K))*S2(J4(K))+S2(J3(K))*DS2(J4(K)))
     *              *GTD(K)*S4(J1(K))*S4(J2(K)))
         ENDDO
         DO K=1,6
            PDOT(IQ(I+12))=PDOT(IQ(I+12))+DUMA(K)*ANG1(NL,K,3)
     *               +3.0D0*GT(NL,K)*ANG1(NL,K,2)*(DF4(K)-DUMT(K))
         ENDDO
C***********************************************************
C                                                          *
C   COMPUTE THE DERIVATIVE OF THE QUARTIC TERM WITH        *
C   RESPECT TO ATOM M                                      *
C                                                          *
C***********************************************************
         DO K=1,6
            DUMA(K)=-((DS4(J1(K))*S4(J2(K))+S4(J1(K))*DS4(J2(K)))
     *              *(HT0(NL,K)+HTD(K)*(S2(J3(K))*S2(J4(K))-1.0D0))
     *              +(DS2(J3(K))*S2(J4(K))+S2(J3(K))*DS2(J4(K)))
     *              *HTD(K)*S4(J1(K))*S4(J2(K)))
         ENDDO
         DO K=1,6
            PDOT(IQ(I+12))=PDOT(IQ(I+12))+DUMA(K)*ANG1(NL,K,4)
     *            +4.0D0*HT(NL,K)*ANG1(NL,K,3)*(DF4(K)-DUMT(K))
         ENDDO
      ENDDO
C***********************************************************
C                                                          *
C   WE HAVE NOW COMPLETED THE CALCULATION OF THE           *
C   DERIVATIVES FOR THE DIAGONAL QUADRATIC, CUBIC, AND     *
C   QUARTIC ANGULAR TERMS OF THE POTENTIAL                 *
C      I=1 => X COORDINATE                                 *
C      I=2 => Y COORDINATE                                 *
C      I=3 => Z COORDINATE                                 *
C                                                          *
C   NOW BEGIN THE OUT-OF-PLANE CALCULATION                 *
C   THE FIRST TASK IS TO COMPUTE DELTA(I,J) (A TOTAL OF    *
C   12 TERMS) AND COS(DELTA(I,J)) ALONG WITH THE HIGHER    *
C   ORDER POWERS OF DELTA(I,J)                             *
C   WE BEGIN BY CALCULATING COORDINATE DIFFERENCES AND     *
C   STORING THEM IN THE ARRAYS XDF, YDF, ZDF               *
C                                                          *
C***********************************************************
      XDF(1)=RC(1)
      XDF(2)=RC(4)
      XDF(3)=Q(IQ(4))-Q(IQ(1))
      XDF(4)=RC(7)
      XDF(5)=Q(IQ(7))-Q(IQ(1))
      XDF(6)=Q(IQ(7))-Q(IQ(4))
      XDF(7)=RC(10)
      XDF(8)=Q(IQ(10))-Q(IQ(1))
      XDF(9)=Q(IQ(10))-Q(IQ(4))
      XDF(10)=Q(IQ(10))-Q(IQ(7))
      YDF(1)=RC(2)
      YDF(2)=RC(5)
      YDF(3)=Q(IQ(5))-Q(IQ(2))
      YDF(4)=RC(8)
      YDF(5)=Q(IQ(8))-Q(IQ(2))
      YDF(6)=Q(IQ(8))-Q(IQ(5))
      YDF(7)=RC(11)
      YDF(8)=Q(IQ(11))-Q(IQ(2))
      YDF(9)=Q(IQ(11))-Q(IQ(5))
      YDF(10)=Q(IQ(11))-Q(IQ(8))
      ZDF(1)=RC(3)
      ZDF(2)=RC(6)
      ZDF(3)=Q(IQ(6))-Q(IQ(3))
      ZDF(4)=RC(9)
      ZDF(5)=Q(IQ(9))-Q(IQ(3))
      ZDF(6)=Q(IQ(9))-Q(IQ(6))
      ZDF(7)=RC(12)
      ZDF(8)=Q(IQ(12))-Q(IQ(3))
      ZDF(9)=Q(IQ(12))-Q(IQ(6))
      ZDF(10)=Q(IQ(12))-Q(IQ(9))
      DUM1=YDF(4)*Q(IQ(12))-YDF(7)*Q(IQ(9))+YDF(10)*Q(IQ(15))
      DUM2=-XDF(4)*Q(IQ(12))+XDF(7)*Q(IQ(9))-XDF(10)*Q(IQ(15))
      DUM3=XDF(4)*Q(IQ(11))-XDF(7)*Q(IQ(8))+XDF(10)*Q(IQ(14))
C*************************************************************
C                                                            *
C   NEW CODE INSERTED TO DETERMINE THE CORRECT NORMAL VECTOR *
C                                                            *
C*************************************************************
      IF (TETTST.NE.1) THEN
         TSTN1=-(DUM1*XDF(2)+DUM2*YDF(2)+DUM3*ZDF(2))
         DST1=SQRT(DUM1**2+DUM2**2+DUM3**2)
         COSN1=TSTN1/(DST1*R(IB(2)))
         TSTN2=-(DUM1*XDF(1)+DUM2*YDF(1)+DUM3*ZDF(1))
         COSN2=TSTN2/(DST1*R(IB(1)))
         SGN1=1.0D0
         IF ((COSN1.GE.0.0D0).AND.(COSN2.LE.0.0D0)) SGN1=-1.0D0
         IF ((COSN1.GE.0.0D0).AND.(COSN1.GE.COSN2)) SGN1=-1.0D0
         SGN2=SGN1
      ENDIF
      TRMA(1)=-(DUM1*XDF(2)+DUM2*YDF(2)+DUM3*ZDF(2))
      TRMB(1)=DUM1*XDF(1)+DUM2*YDF(1)+DUM3*ZDF(1)
      SQ1=XDF(9)**2
      SQ2=YDF(9)**2
      SQ3=ZDF(9)**2
      DUM1=XDF(6)**2*(SQ2+SQ3)
      DUM2=YDF(6)**2*(SQ1+SQ3)
      DUM3=ZDF(6)**2*(SQ1+SQ2)
      DUM4=XDF(6)*XDF(9)
      DUM5=YDF(6)*YDF(9)
      DUM6=ZDF(6)*ZDF(9)
      TRMA(2)=SQRT(DUM1+DUM2+DUM3
     +       -2.0D0*(DUM5*DUM6+DUM4*DUM6+DUM4*DUM5))
      SQ1=XDF(8)**2
      SQ2=YDF(8)**2
      SQ3=ZDF(8)**2
      DUM1=XDF(5)**2*(SQ2+SQ3)
      DUM2=YDF(5)**2*(SQ1+SQ3)
      DUM3=ZDF(5)**2*(SQ1+SQ2)
      DUM4=XDF(5)*XDF(8)
      DUM5=YDF(5)*YDF(8)
      DUM6=ZDF(5)*ZDF(8)
      TRMB(2)=SQRT(DUM1+DUM2+DUM3
     +       -2.0D0*(DUM5*DUM6+DUM4*DUM6+DUM4*DUM5))
      DUM1=YDF(1)*Q(IQ(6))-YDF(2)*Q(IQ(3))+YDF(3)*Q(IQ(15))
      DUM2=-XDF(1)*Q(IQ(6))+XDF(2)*Q(IQ(3))-XDF(3)*Q(IQ(15))
      DUM3=XDF(1)*Q(IQ(5))-XDF(2)*Q(IQ(2))+XDF(3)*Q(IQ(14))
C**********************************************************
C                                                         *
C   NEW CODE ADDED TO DETERMINE THE CORRECT NORMAL VECTOR *
C                                                         *
C**********************************************************
      IF (TETTST.NE.1) THEN
         TETTST=1
         TSTN3=-(DUM1*XDF(7)+DUM2*YDF(7)+DUM3*ZDF(7))
         DST2=SQRT(DUM1**2+DUM2**2+DUM3**2)
         COSN3=TSTN3/(DST2*R(IB(4)))
         TSTN4=-(DUM1*XDF(4)+DUM2*YDF(4)+DUM3*ZDF(4))
         COSN4=TSTN4/(DST2*R(IB(3)))
         SGN3=1.0D0
         IF ((COSN3.GE.0.0D0).AND.(COSN4.LE.0.0D0)) SGN3=-1.0D0
         IF ((COSN3.GE.0.0D0).AND.(COSN3.GE.COSN4)) SGN3=-1.0D0
         SGN4=SGN3
      ENDIF
      TRMC(1)=-(DUM1*XDF(7)+DUM2*YDF(7)+DUM3*ZDF(7))
      TRMD(1)=DUM1*XDF(4)+DUM2*YDF(4)+DUM3*ZDF(4)
      TRMA(1)=TRMA(1)*SGN1
      TRMB(1)=TRMB(1)*SGN2
      TRMC(1)=TRMC(1)*SGN3
      TRMD(1)=TRMD(1)*SGN4
C***********************************************************
C                                                          *
C   NOTE:  WE DO NOT NEED TO RECOMPUTE SQ1, SQ2, SQ3       *
C                                                          *
C***********************************************************
      DUM1=XDF(3)**2*(SQ2+SQ3)
      DUM2=YDF(3)**2*(SQ1+SQ3)
      DUM3=ZDF(3)**2*(SQ1+SQ2)
      DUM4=XDF(3)*XDF(8)
      DUM5=YDF(3)*YDF(8)
      DUM6=ZDF(3)*ZDF(8)
      TRMC(2)=SQRT(DUM1+DUM2+DUM3
     *       -2.0D0*(DUM5*DUM6+DUM4*DUM6+DUM4*DUM5))
      SQ1=XDF(5)**2
      SQ2=YDF(5)**2
      SQ3=ZDF(5)**2
      DUM1=XDF(3)**2*(SQ2+SQ3)
      DUM2=YDF(3)**2*(SQ1+SQ3)
      DUM3=ZDF(3)**2*(SQ1+SQ2)
      DUM4=XDF(3)*XDF(5)
      DUM5=YDF(3)*YDF(5)
      DUM6=ZDF(3)*ZDF(5)
      TRMD(2)=SQRT(DUM1+DUM2+DUM3
     *       -2.0D0*(DUM5*DUM6+DUM4*DUM6+DUM4*DUM5))
C***********************************************************
C                                                          *
C   CALCULATION OF DELTA(I,J) AND COS(DELTA(I,J))          *
C                                                          *
C***********************************************************
      DUMA(1)=1.0D0/R(IB(1))
      DUMA(2)=1.0D0/R(IB(2))
      DUMA(3)=1.0D0/R(IB(3))
      DUMA(4)=1.0D0/R(IB(4))
      DUMT(1)=TRMA(1)/TRMA(2)
      DUMT(2)=TRMB(1)/TRMB(2)
      DUMT(3)=TRMC(1)/TRMC(2)
      DUMT(4)=TRMD(1)/TRMD(2)
      DO I=1,4
         J=1+3*(I-1)
         CDLTA(J)=DUMT(I)*DUMA(K1(I))
         CDLTA(J+1)=DUMT(I)*DUMA(K2(I))
         CDLTA(J+2)=DUMT(I)*DUMA(K3(I))
      ENDDO
      DO I=1,4
         J=1+12*(I-1)
         K=1+3*(I-1)
         DLTA(NL,J)=ACOS(CDLTA(K))-TN(K4(I))
         DLTA(NL,J+4)=ACOS(CDLTA(K+1))-TN(K5(I))
         DLTA(NL,J+8)=ACOS(CDLTA(K+2))-TN(K6(I))
      ENDDO
      DO I=1,12
         J=4*I-3
         DLTA(NL,J+1)=DLTA(NL,J)*DLTA(NL,J)
         DLTA(NL,J+2)=DLTA(NL,J+1)*DLTA(NL,J)
         DLTA(NL,J+3)=DLTA(NL,J+2)*DLTA(NL,J)
      ENDDO
C***********************************************************
C                                                          *
C   CALCULATE THE QUADRATIC AND QUARTIC FORCE CONSTANTS    *
C   FOR THE OUT-OF-PLANE BEND                              *
C                                                          *
C***********************************************************
      FD(NL,1)=(1.0D0-S3(1))*S3(2)*S3(3)*S3(4)*FD1(NL,1)
      FD(NL,2)=(1.0D0-S3(2))*S3(1)*S3(3)*S3(4)*FD1(NL,2)
      FD(NL,3)=(1.0D0-S3(3))*S3(1)*S3(2)*S3(4)*FD1(NL,3)
      FD(NL,4)=(1.0D0-S3(4))*S3(1)*S3(2)*S3(3)*FD1(NL,4)
      HD(NL,1)=(1.0D0-S3(1))*S3(2)*S3(3)*S3(4)*HD1(NL,1)
      HD(NL,2)=(1.0D0-S3(2))*S3(1)*S3(3)*S3(4)*HD1(NL,2)
      HD(NL,3)=(1.0D0-S3(3))*S3(1)*S3(2)*S3(4)*HD1(NL,3)
      HD(NL,4)=(1.0D0-S3(4))*S3(1)*S3(2)*S3(3)*HD1(NL,4)
C***********************************************************
C                                                          *
C   WE CAN NOW BEGIN CALCULATING THE NECESSARY DERIVATIVES *
C   OF THE DELTA(I,J).  THESE ARE COMPUTED USING THE       *
C   QUANTITIES TRMA, TRMB, TRMC, TRMD, AND THE EXPRESSIONS *
C   FOR THE LENGTHS OF THE FOUR BONDS.  WE ALSO NEED THE   *
C   DERIVATIVES OF THE S3 SWITCHING FUNCTION.  BEGIN BY    *
C   CALCULATING THE DERIVATIVES OF F1, F2, F3, F4 (THE     *
C   EXPRESSIONS FOR THE BOND LENGTHS).                     *
C                                                          *
C***********************************************************
      DO I=1,15
         DFN1(I)=0.0D0
         DFN2(I)=0.0D0
         DFN3(I)=0.0D0
         DFN4(I)=0.0D0
      ENDDO
      D1=1.0D0/R(IB(1))
      DFN1(1)=XDF(1)*D1
      DFN1(2)=YDF(1)*D1
      DFN1(3)=ZDF(1)*D1
      DFN1(13)=-DFN1(1)
      DFN1(14)=-DFN1(2)
      DFN1(15)=-DFN1(3)
      D1=1.0D0/R(IB(2))
      DFN2(4)=XDF(2)*D1
      DFN2(5)=YDF(2)*D1
      DFN2(6)=ZDF(2)*D1
      DFN2(13)=-DFN2(4)
      DFN2(14)=-DFN2(5)
      DFN2(15)=-DFN2(6)
      D1=1.0D0/R(IB(3))
      DFN3(7)=XDF(4)*D1
      DFN3(8)=YDF(4)*D1
      DFN3(9)=ZDF(4)*D1
      DFN3(13)=-DFN3(7)
      DFN3(14)=-DFN3(8)
      DFN3(15)=-DFN3(9)
      D1=1.0D0/R(IB(4))
      DFN4(10)=XDF(7)*D1
      DFN4(11)=YDF(7)*D1
      DFN4(12)=ZDF(7)*D1
      DFN4(13)=-DFN4(10)
      DFN4(14)=-DFN4(11)
      DFN4(15)=-DFN4(12)
C***********************************************************
C                                                          *
C   CALCULATION OF DERIVATIVES OF TRMA(1)                  *
C   BEGIN BY FILLING THE DERIVATIVE STORAGE ARRAYS WITH    *
C   ZEROS.                                                 *
C                                                          *
C***********************************************************
      DO I=1,15
         DA1(I)=0.0D0
         DA2(I)=0.0D0
         DB1(I)=0.0D0
         DB2(I)=0.0D0
         DC1(I)=0.0D0
         DC2(I)=0.0D0
         DD1(I)=0.0D0
         DD2(I)=0.0D0
      ENDDO
      DA1(4)=-YDF(4)*Q(IQ(12))+YDF(7)*Q(IQ(9))-YDF(10)*Q(IQ(15))
      DA1(5)=XDF(4)*Q(IQ(12))-XDF(7)*Q(IQ(9))+XDF(10)*Q(IQ(15))
      DA1(6)=-XDF(4)*Q(IQ(11))+XDF(7)*Q(IQ(8))-XDF(10)*Q(IQ(14))
      DA1(7)=ZDF(7)*YDF(2)-YDF(7)*ZDF(2)
      DA1(8)=XDF(7)*ZDF(2)-ZDF(7)*XDF(2)
      DA1(9)=YDF(7)*XDF(2)-XDF(7)*YDF(2)
      DA1(10)=YDF(4)*ZDF(2)-ZDF(4)*YDF(2)
      DA1(11)=ZDF(4)*XDF(2)-XDF(4)*ZDF(2)
      DA1(12)=XDF(4)*YDF(2)-YDF(4)*XDF(2)
      DA1(13)=YDF(10)*ZDF(2)-ZDF(10)*YDF(2)-DA1(4)
      DA1(14)=ZDF(10)*XDF(2)-XDF(10)*ZDF(2)-DA1(5)
      DA1(15)=XDF(10)*YDF(2)-YDF(10)*XDF(2)-DA1(6)
      DO I=1,15
         DA1(I)=DA1(I)*SGN1
      ENDDO
C***********************************************************
C                                                          *
C   CALCULATE DERIVATIVES FOR TRMA(2)                      *
C                                                          *
C***********************************************************
      Q1=1.0D0/TRMA(2)
      DUM1=XDF(6)**2
      DUM2=YDF(6)**2
      DUM3=ZDF(6)**2
      DUM4=XDF(9)**2
      DUM5=YDF(9)**2
      DUM6=ZDF(9)**2
      TERM1=XDF(6)*XDF(9)
      TERM2=YDF(6)*YDF(9)
      TERM3=ZDF(6)*ZDF(9)
      DUM45=DUM4+DUM5
      DUM46=DUM4+DUM6
      DUM56=DUM5+DUM6
      DUM12=DUM1+DUM2
      DUM13=DUM1+DUM3
      DUM23=DUM2+DUM3
      TERM12=TERM1+TERM2
      TERM13=TERM1+TERM3
      TERM23=TERM2+TERM3
      DA2(4)=-Q1*(XDF(6)*DUM56+XDF(9)*DUM23
     +      +(2.0D0*Q(IQ(4))-Q(IQ(7))-Q(IQ(10)))*TERM23)
      DA2(5)=-Q1*(YDF(6)*DUM46+YDF(9)*DUM13
     +      +(2.0D0*Q(IQ(5))-Q(IQ(8))-Q(IQ(11)))*TERM13)
      DA2(6)=-Q1*(ZDF(6)*DUM45+ZDF(9)*DUM12
     +      +(2.0D0*Q(IQ(6))-Q(IQ(9))-Q(IQ(12)))*TERM12)
      DA2(7)=Q1*(XDF(6)*DUM56-XDF(9)*TERM23)
      DA2(8)=Q1*(YDF(6)*DUM46-YDF(9)*TERM13)
      DA2(9)=Q1*(ZDF(6)*DUM45-ZDF(9)*TERM12)
      DA2(10)=Q1*(XDF(9)*DUM23-XDF(6)*TERM23)
      DA2(11)=Q1*(YDF(9)*DUM13-YDF(6)*TERM13)
      DA2(12)=Q1*(ZDF(9)*DUM12-ZDF(6)*TERM12)
C***********************************************************
C                                                          *
C   CALCULATION OF DERIVATIVES OF TRMB(1)                  *
C                                                          *
C***********************************************************
      FACT=1.0D0/SGN1
      DB1(1)=-DA1(4)*FACT
      DB1(2)=-DA1(5)*FACT
      DB1(3)=-DA1(6)*FACT
      DB1(7)=YDF(7)*ZDF(1)-ZDF(7)*YDF(1)
      DB1(8)=ZDF(7)*XDF(1)-XDF(7)*ZDF(1)
      DB1(9)=XDF(7)*YDF(1)-YDF(7)*XDF(1)
      DB1(10)=ZDF(4)*YDF(1)-YDF(4)*ZDF(1)
      DB1(11)=XDF(4)*ZDF(1)-ZDF(4)*XDF(1)
      DB1(12)=YDF(4)*XDF(1)-XDF(4)*YDF(1)
      DB1(13)=ZDF(10)*YDF(1)-YDF(10)*ZDF(1)+DA1(4)*FACT
      DB1(14)=XDF(10)*ZDF(1)-ZDF(10)*XDF(1)+DA1(5)*FACT
      DB1(15)=YDF(10)*XDF(1)-XDF(10)*YDF(1)+DA1(6)*FACT
      DO I=1,15
         DB1(I)=DB1(I)*SGN2
      ENDDO
C***********************************************************
C                                                          *
C   CALCULATION OF DERIVATIVES FOR TRMB(2)                 *
C                                                          *
C***********************************************************
      Q1=1.0D0/TRMB(2)
      DUM1=XDF(5)**2
      DUM2=YDF(5)**2
      DUM3=ZDF(5)**2
      DUM4=XDF(8)**2
      DUM5=YDF(8)**2
      DUM6=ZDF(8)**2
      TERM1=XDF(5)*XDF(8)
      TERM2=YDF(5)*YDF(8)
      TERM3=ZDF(5)*ZDF(8)
      DUM45=DUM4+DUM5
      DUM46=DUM4+DUM6
      DUM56=DUM5+DUM6
      DUM12=DUM1+DUM2
      DUM13=DUM1+DUM3
      DUM23=DUM2+DUM3
      TERM12=TERM1+TERM2
      TERM13=TERM1+TERM3
      TERM23=TERM2+TERM3
      DB2(1)=-Q1*(XDF(5)*DUM56+XDF(8)*DUM23
     +      +(2.0D0*Q(IQ(1))-Q(IQ(7))-Q(IQ(10)))*TERM23)
      DB2(2)=-Q1*(YDF(5)*DUM46+YDF(8)*DUM13
     +      +(2.0D0*Q(IQ(2))-Q(IQ(8))-Q(IQ(11)))*TERM13)
      DB2(3)=-Q1*(ZDF(5)*DUM45+ZDF(8)*DUM12
     +      +(2.0D0*Q(IQ(3))-Q(IQ(9))-Q(IQ(12)))*TERM12)
      DB2(7)=Q1*(XDF(5)*DUM56-XDF(8)*TERM23)
      DB2(8)=Q1*(YDF(5)*DUM46-YDF(8)*TERM13)
      DB2(9)=Q1*(ZDF(5)*DUM45-ZDF(8)*TERM12)
      DB2(10)=Q1*(XDF(8)*DUM23-XDF(5)*TERM23)
      DB2(11)=Q1*(YDF(8)*DUM13-YDF(5)*TERM13)
      DB2(12)=Q1*(ZDF(8)*DUM12-ZDF(5)*TERM12)
C***********************************************************
C                                                          *
C   CALCULATION OF DERIVATIVES FOR TRMC(1)                 *
C                                                          *
C***********************************************************
      DC1(1)=ZDF(2)*YDF(7)-YDF(2)*ZDF(7)
      DC1(2)=XDF(2)*ZDF(7)-ZDF(2)*XDF(7)
      DC1(3)=YDF(2)*XDF(7)-XDF(2)*YDF(7)
      DC1(4)=YDF(1)*ZDF(7)-ZDF(1)*YDF(7)
      DC1(5)=ZDF(1)*XDF(7)-XDF(1)*ZDF(7)
      DC1(6)=XDF(1)*YDF(7)-YDF(1)*XDF(7)
      DC1(10)=-YDF(1)*Q(IQ(6))+YDF(2)*Q(IQ(3))-YDF(3)*Q(IQ(15))
      DC1(11)=XDF(1)*Q(IQ(6))-XDF(2)*Q(IQ(3))+XDF(3)*Q(IQ(15))
      DC1(12)=-XDF(1)*Q(IQ(5))+XDF(2)*Q(IQ(2))-XDF(3)*Q(IQ(14))
      DC1(13)=YDF(3)*ZDF(7)-ZDF(3)*YDF(7)-DC1(10)
      DC1(14)=ZDF(3)*XDF(7)-XDF(3)*ZDF(7)-DC1(11)
      DC1(15)=XDF(3)*YDF(7)-YDF(3)*XDF(7)-DC1(12)
      DO I=1,15
         DC1(I)=DC1(I)*SGN3
      ENDDO
C***********************************************************
C                                                          *
C   CALCULATION OF DERIVATIVES FOR TRMC(2)                 *
C   NOTE: WE DO NOT NEED TO RECOMPUTE DUM4, DUM5, DUM6,    *
C         DUM45, DUM46, DUM56                              *
C                                                          *
C***********************************************************
      Q1=1.0D0/TRMC(2)
      DUM1=XDF(3)**2
      DUM2=YDF(3)**2
      DUM3=ZDF(3)**2
      TERM1=XDF(3)*XDF(8)
      TERM2=YDF(3)*YDF(8)
      TERM3=ZDF(3)*ZDF(8)
      DUM12=DUM1+DUM2
      DUM13=DUM1+DUM3
      DUM23=DUM2+DUM3
      TERM12=TERM1+TERM2
      TERM13=TERM1+TERM3
      TERM23=TERM2+TERM3
      DC2(1)=-Q1*(XDF(3)*DUM56+XDF(8)*DUM23
     +      +(2.0D0*Q(IQ(1))-Q(IQ(4))-Q(IQ(10)))*TERM23)
      DC2(2)=-Q1*(YDF(3)*DUM46+YDF(8)*DUM13
     +      +(2.0D0*Q(IQ(2))-Q(IQ(5))-Q(IQ(11)))*TERM13)
      DC2(3)=-Q1*(ZDF(3)*DUM45+ZDF(8)*DUM12
     +      +(2.0D0*Q(IQ(3))-Q(IQ(6))-Q(IQ(12)))*TERM12)
      DC2(4)=Q1*(XDF(3)*DUM56-XDF(8)*TERM23)
      DC2(5)=Q1*(YDF(3)*DUM46-YDF(8)*TERM13)
      DC2(6)=Q1*(ZDF(3)*DUM45-ZDF(8)*TERM12)
      DC2(10)=Q1*(XDF(8)*DUM23-XDF(3)*TERM23)
      DC2(11)=Q1*(YDF(8)*DUM13-YDF(3)*TERM13)
      DC2(12)=Q1*(ZDF(8)*DUM12-ZDF(3)*TERM12)
C***********************************************************
C                                                          *
C   CALCULATION OF DERIVATIVES FOR TRMD(1)                 *
C                                                          *
C***********************************************************
      FACT=1.0D0/SGN3
      DD1(1)=YDF(2)*ZDF(4)-ZDF(2)*YDF(4)
      DD1(2)=ZDF(2)*XDF(4)-XDF(2)*ZDF(4)
      DD1(3)=XDF(2)*YDF(4)-YDF(2)*XDF(4)
      DD1(4)=ZDF(1)*YDF(4)-YDF(1)*ZDF(4)
      DD1(5)=XDF(1)*ZDF(4)-ZDF(1)*XDF(4)
      DD1(6)=YDF(1)*XDF(4)-XDF(1)*YDF(4)
      DD1(7)=-DC1(10)*FACT
      DD1(8)=-DC1(11)*FACT
      DD1(9)=-DC1(12)*FACT
      DD1(13)=ZDF(3)*YDF(4)-YDF(3)*ZDF(4)+DC1(10)*FACT
      DD1(14)=XDF(3)*ZDF(4)-ZDF(3)*XDF(4)+DC1(11)*FACT
      DD1(15)=YDF(3)*XDF(4)-XDF(3)*YDF(4)+DC1(12)*FACT
      DO I=1,15
         DD1(I)=DD1(I)*SGN4
      ENDDO
C***********************************************************
C                                                          *
C   CALCULATION OF DERIVATIVES FOR TRMD(2)                 *
C   NOTE:  WE NEED NOT RECOMPUTE DUM1, DUM2, DUM3,         *
C          DUM12, DUM13, DUM23                             *
C                                                          *
C***********************************************************
      Q1=1.0D0/TRMD(2)
      DUM4=XDF(5)**2
      DUM5=YDF(5)**2
      DUM6=ZDF(5)**2
      TERM1=XDF(3)*XDF(5)
      TERM2=YDF(3)*YDF(5)
      TERM3=ZDF(3)*ZDF(5)
      DUM45=DUM4+DUM5
      DUM46=DUM4+DUM6
      DUM56=DUM5+DUM6
      TERM12=TERM1+TERM2
      TERM13=TERM1+TERM3
      TERM23=TERM2+TERM3
      DD2(1)=-Q1*(XDF(3)*DUM56+XDF(5)*DUM23
     +      +(2.0D0*Q(IQ(1))-Q(IQ(4))-Q(IQ(7)))*TERM23)
      DD2(2)=-Q1*(YDF(3)*DUM46+YDF(5)*DUM13
     +      +(2.0D0*Q(IQ(2))-Q(IQ(5))-Q(IQ(8)))*TERM13)
      DD2(3)=-Q1*(ZDF(3)*DUM45+ZDF(5)*DUM12
     +      +(2.0D0*Q(IQ(3))-Q(IQ(6))-Q(IQ(9)))*TERM12)
      DD2(4)=Q1*(XDF(3)*DUM56-XDF(5)*TERM23)
      DD2(5)=Q1*(YDF(3)*DUM46-YDF(5)*TERM13)
      DD2(6)=Q1*(ZDF(3)*DUM45-ZDF(5)*TERM12)
      DD2(7)=Q1*(XDF(5)*DUM23-XDF(3)*TERM23)
      DD2(8)=Q1*(YDF(5)*DUM13-YDF(3)*TERM13)
      DD2(9)=Q1*(ZDF(5)*DUM12-ZDF(3)*TERM12)
C***********************************************************
C                                                          *
C   WE ARE NOW READY TO COMPUTE DV/DQ FOR THE OUT-OF-PLANE *
C   BEND.  TO DO THIS WE NEED D(DELTA(I,J))/DXI AND        *
C   DS3/SXI.                                               *
C                                                          *
C   BEGIN MAIN LOOP FOR OUT-OF-PLANE BEND CALCULATION      *
C                                                          *
C***********************************************************
      DO 700 I=1,3
C***********************************************************
C                                                          *
C   FIRST COMPUTE THE DERIVATIVE OF THE S3 SWITCHING       *
C   FUNCTION FOR BONDS 1, 2, 3, 4                          *
C                                                          *
C***********************************************************
         DO K=1,4
            L=I+3*(K-1)
            D1=1.0D0/R(IB(K))
            D2=R0(NL,K)*D1
            TERM1=(1.0D0-S3(K))**2-1.0D0
            TERM2=A3*RC(L)*D1*(R(IB(K))-B3)**C3
            TERM3=C3*A3*RC(L)*(R(IB(K))-B3)**(C3-1.0D0)*(1.0D0-D2)
            DS3(K)=TERM1*(TERM2+TERM3)
         ENDDO
C************************************************************
C                                                           *
C   NOW CALCULATE THE DERIVATIVES USING THE QUANTITIES      *
C   TRMA, TRMB, TRMC, TRMD, F1, F2, F3, F4 AND THEIR        *
C   DERIVATIVES.  WE USE THE ARRAYS DUMA AND DUMT TO        *
C   HOLD CALCULATED QUANTITIES WHICH WILL BE USED           *
C   REPETITIVELY.                                           *
C                                                           *
C************************************************************
         DO K=1,6
            DUMA(K)=-1.0D0/SQRT(1.0D0-CDLTA(K)**2)
            DUMT(K)=-1.0D0/SQRT(1.0D0-CDLTA(K+6)**2)
         ENDDO
         DERV(1,1)=-DS3(1)*S3(2)*S3(3)*S3(4)
         DERV(1,2)=(1.0D0-S3(2))*DS3(1)*S3(3)*S3(4)
         DERV(1,3)=(1.0D0-S3(3))*DS3(1)*S3(2)*S3(4)
         DERV(1,4)=(1.0D0-S3(4))*DS3(1)*S3(2)*S3(3)
         DERV(2,1)=(1.0D0-S3(1))*DS3(2)*S3(3)*S3(4)
         DERV(2,2)=-DS3(2)*S3(1)*S3(3)*S3(4)
         DERV(2,3)=(1.0D0-S3(3))*DS3(2)*S3(1)*S3(4)
         DERV(2,4)=(1.0D0-S3(4))*DS3(2)*S3(1)*S3(3)
         DERV(3,1)=(1.0D0-S3(1))*DS3(3)*S3(2)*S3(4)
         DERV(3,2)=(1.0D0-S3(2))*DS3(3)*S3(1)*S3(4)
         DERV(3,3)=-DS3(3)*S3(1)*S3(2)*S3(4)
         DERV(3,4)=(1.0D0-S3(4))*DS3(3)*S3(1)*S3(2)
         DERV(4,1)=(1.0D0-S3(1))*DS3(4)*S3(2)*S3(3)
         DERV(4,2)=(1.0D0-S3(2))*DS3(4)*S3(1)*S3(3)
         DERV(4,3)=(1.0D0-S3(3))*DS3(4)*S3(1)*S3(2)
         DERV(4,4)=-DS3(4)*S3(1)*S3(2)*S3(3)
         DERV(5,1)=DS3(1)*S3(2)*S3(3)*S3(4)
     *             -(1.0D0-S3(1))*DS3(2)*S3(3)*S3(4)
     *             -(1.0D0-S3(1))*S3(2)*DS3(3)*S3(4)
     *             -(1.0D0-S3(1))*S3(2)*S3(3)*DS3(4)
         DERV(5,2)=DS3(2)*S3(1)*S3(3)*S3(4)
     *             -(1.0D0-S3(2))*DS3(1)*S3(3)*S3(4)
     *             -(1.0D0-S3(2))*S3(1)*DS3(3)*S3(4)
     *             -(1.0D0-S3(2))*S3(1)*S3(3)*DS3(4)
         DERV(5,3)=DS3(3)*S3(1)*S3(2)*S3(4)
     *             -(1.0D0-S3(3))*DS3(1)*S3(2)*S3(4)
     *             -(1.0D0-S3(3))*S3(1)*DS3(2)*S3(4)
     *             -(1.0D0-S3(3))*S3(1)*S3(2)*DS3(4)
         DERV(5,4)=DS3(4)*S3(1)*S3(2)*S3(3)
     *             -(1.0D0-S3(4))*DS3(1)*S3(2)*S3(3)
     *             -(1.0D0-S3(4))*S3(1)*DS3(2)*S3(3)
     *             -(1.0D0-S3(4))*S3(1)*S3(2)*DS3(3)
         DO 650 J=1,5
            L=I+3*(J-1)
            DUM1=TRMA(2)*R(IB(2))
            DUM2=TRMA(2)*R(IB(3))
            DUM3=TRMA(2)*R(IB(4))
            SQ1=1.0D0/DUM1**2
            SQ2=1.0D0/DUM2**2
            SQ3=1.0D0/DUM3**2
            TERM1=DA1(L)*DUM1
            TERM2=TRMA(1)*(DA2(L)*R(IB(2))+TRMA(2)*DFN2(L))
            DANG1(1)=DUMA(1)*SQ1*(TERM1-TERM2)-OUTPL(L,1)
            TERM1=DA1(L)*DUM2
            TERM2=TRMA(1)*(DA2(L)*R(IB(3))+TRMA(2)*DFN3(L))
            DANG1(2)=DUMA(2)*SQ2*(TERM1-TERM2)-OUTPL(L,2)
            TERM1=DA1(L)*DUM3
            TERM2=TRMA(1)*(DA2(L)*R(IB(4))+TRMA(2)*DFN4(L))
            DANG1(3)=DUMA(3)*SQ3*(TERM1-TERM2)-OUTPL(L,3)
            DUM1=TRMB(2)*R(IB(1))
            DUM2=TRMB(2)*R(IB(3))
            DUM3=TRMB(2)*R(IB(4))
            SQ1=1.0D0/DUM1**2
            SQ2=1.0D0/DUM2**2
            SQ3=1.0D0/DUM3**2
            TERM1=DB1(L)*DUM1
            TERM2=TRMB(1)*(DB2(L)*R(IB(1))+TRMB(2)*DFN1(L))
            DANG1(4)=DUMA(4)*SQ1*(TERM1-TERM2)-OUTPL(L,1)
            TERM1=DB1(L)*DUM2
            TERM2=TRMB(1)*(DB2(L)*R(IB(3))+TRMB(2)*DFN3(L))
            DANG1(5)=DUMA(5)*SQ2*(TERM1-TERM2)-OUTPL(L,4)
            TERM1=DB1(L)*DUM3
            TERM2=TRMB(1)*(DB2(L)*R(IB(4))+TRMB(2)*DFN4(L))
            DANG1(6)=DUMA(6)*SQ3*(TERM1-TERM2)-OUTPL(L,5)
            DUM1=TRMC(2)*R(IB(1))
            DUM2=TRMC(2)*R(IB(2))
            DUM3=TRMC(2)*R(IB(4))
            SQ1=1.0D0/DUM1**2
            SQ2=1.0D0/DUM2**2
            SQ3=1.0D0/DUM3**2
            TERM1=DC1(L)*DUM1
            TERM2=TRMC(1)*(DC2(L)*R(IB(1))+TRMC(2)*DFN1(L))
            DANG1(7)=DUMT(1)*SQ1*(TERM1-TERM2)-OUTPL(L,2)
            TERM1=DC1(L)*DUM2
            TERM2=TRMC(1)*(DC2(L)*R(IB(2))+TRMC(2)*DFN2(L))
            DANG1(8)=DUMT(2)*SQ2*(TERM1-TERM2)-OUTPL(L,4)
            TERM1=DC1(L)*DUM3
            TERM2=TRMC(1)*(DC2(L)*R(IB(4))+TRMC(2)*DFN4(L))
            DANG1(9)=DUMT(3)*SQ3*(TERM1-TERM2)-OUTPL(L,6)
            DUM1=TRMD(2)*R(IB(1))
            DUM2=TRMD(2)*R(IB(2))
            DUM3=TRMD(2)*R(IB(3))
            SQ1=1.0D0/DUM1**2
            SQ2=1.0D0/DUM2**2
            SQ3=1.0D0/DUM3**2
            TERM1=DD1(L)*DUM1
            TERM2=TRMD(1)*(DD2(L)*R(IB(1))+TRMD(2)*DFN1(L))
            DANG1(10)=DUMT(4)*SQ1*(TERM1-TERM2)-OUTPL(L,3)
            TERM1=DD1(L)*DUM2
            TERM2=TRMD(1)*(DD2(L)*R(IB(2))+TRMD(2)*DFN2(L))
            DANG1(11)=DUMT(5)*SQ2*(TERM1-TERM2)-OUTPL(L,5)
            TERM1=DD1(L)*DUM3
            TERM2=TRMD(1)*(DD2(L)*R(IB(3))+TRMD(2)*DFN3(L))
            DANG1(12)=DUMT(6)*SQ3*(TERM1-TERM2)-OUTPL(L,6)
C***********************************************************
C                                                          *
C   COMPUTE QUADRATIC TERM DERIVATIVES                     *
C                                                          *
C***********************************************************
            DUM1=FD1(NL,1)*DERV(J,1)*(DLTA(NL,2)+DLTA(NL,6)+DLTA(NL,10))
            DUM2=FD1(NL,2)*DERV(J,2)*(DLTA(NL,14)
     *           +DLTA(NL,18)+DLTA(NL,22))
            DUM3=FD1(NL,3)*DERV(J,3)*(DLTA(NL,26)
     *           +DLTA(NL,30)+DLTA(NL,34))
            DUM4=FD1(NL,4)*DERV(J,4)*(DLTA(NL,38)
     *           +DLTA(NL,42)+DLTA(NL,46))
            TERM1=DUM1+DUM2+DUM3+DUM4
            DUM1=FD(NL,1)*(DLTA(NL,1)*DANG1(1)+DLTA(NL,5)*DANG1(2)
     *           +DLTA(NL,9)*DANG1(3))
            DUM2=FD(NL,2)*(DLTA(NL,13)*DANG1(4)+DLTA(NL,17)*DANG1(5)
     *           +DLTA(NL,21)*DANG1(6))
            DUM3=FD(NL,3)*(DLTA(NL,25)*DANG1(7)+DLTA(NL,29)*DANG1(8)
     *           +DLTA(NL,33)*DANG1(9))
            DUM4=FD(NL,4)*(DLTA(NL,37)*DANG1(10)+DLTA(NL,41)*DANG1(11)
     *           +DLTA(NL,45)*DANG1(12))
            PDOT(IQ(L))=PDOT(IQ(L))+TERM1
     *                  +2.0D0*(DUM1+DUM2+DUM3+DUM4)
C***********************************************************
C                                                          *
C   COMPUTE QUARTIC TERM DERIVATIVES                       *
C                                                          *
C***********************************************************
            DUM1=HD1(NL,1)*DERV(J,1)*(DLTA(NL,4)+DLTA(NL,8)+DLTA(NL,12))
            DUM2=HD1(NL,2)*DERV(J,2)*(DLTA(NL,16)
     *           +DLTA(NL,20)+DLTA(NL,24))
            DUM3=HD1(NL,3)*DERV(J,3)*(DLTA(NL,28)
     *           +DLTA(NL,32)+DLTA(NL,36))
            DUM4=HD1(NL,4)*DERV(J,4)*(DLTA(NL,40)
     *           +DLTA(NL,44)+DLTA(NL,48))
            TERM1=DUM1+DUM2+DUM3+DUM4
            DUM1=HD(NL,1)*(DLTA(NL,3)*DANG1(1)+DLTA(NL,7)*DANG1(2)
     *           +DLTA(NL,11)*DANG1(3))
            DUM2=HD(NL,2)*(DLTA(NL,15)*DANG1(4)+DLTA(NL,19)*DANG1(5)
     *           +DLTA(NL,23)*DANG1(6))
            DUM3=HD(NL,3)*(DLTA(NL,27)*DANG1(7)+DLTA(NL,31)*DANG1(8)
     *           +DLTA(NL,35)*DANG1(9))
            DUM4=HD(NL,4)*(DLTA(NL,39)*DANG1(10)+DLTA(NL,43)*DANG1(11)
     *           +DLTA(NL,47)*DANG1(12))
            PDOT(IQ(L))=PDOT(IQ(L))+TERM1
     *                  +4.0D0*(DUM1+DUM2+DUM3+DUM4)
C***********************************************************
C                                                          *
C   CALL THE NON-DIAGONAL CUBIC SUBROUTINE TO COMPUTE      *
C   THE NON-DIAGONAL CUBIC DERIVATIVES                     *
C                                                          *
C***********************************************************
            CALL CUBEND(NL,I,J,DERV2)
            PDOT(IQ(L))=PDOT(IQ(L))+DERV2
  650    CONTINUE
  700 CONTINUE
      RETURN
      END
