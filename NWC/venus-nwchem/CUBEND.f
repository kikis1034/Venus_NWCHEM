      SUBROUTINE CUBEND(NL,N1,N2,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C***********************************************************
C                                                          *
C   THIS SUBROUTINE COMPUTES THE DERIVATIVES OF THE        *
C   NON-DIAGONAL CUBIC TERMS OF THE CH4 POTENTIAL          *
C                                                          *
C***********************************************************
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
      DIMENSION DGNDX(10)
C***********************************************************
C                                                          *
C   CALCULATE THE NON-DAIGONAL CUBIC FORCE CONSTANT        *
C                                                          *
C***********************************************************
      GN4(NL)=S3(1)*S3(2)*S3(3)*S3(4)*GN0(NL,4)
C***********************************************************
C                                                          *
C   CALCULATE DGN4/DXI                                     *
C                                                          *
C***********************************************************
      DGNDX(1)=DS3(1)*S3(2)*S3(3)*S3(4)*GN0(NL,4)
      DGNDX(2)=S3(1)*DS3(2)*S3(3)*S3(4)*GN0(NL,4)
      DGNDX(3)=S3(1)*S3(2)*DS3(3)*S3(4)*GN0(NL,4)
      DGNDX(4)=S3(1)*S3(2)*S3(3)*DS3(4)*GN0(NL,4)
      DGN4=0
      IF (N2.GT.4) THEN
         DO I=1,4
            DGN4=DGN4-DGNDX(I)
         ENDDO
      ELSE
         DGN4=DGNDX(N2)
      ENDIF
      DUM1=ANG1(NL,2,1)*ANG1(NL,4,1)+ANG1(NL,3,1)*ANG1(NL,5,1)
      DUM2=ANG1(NL,1,1)*ANG1(NL,4,1)+ANG1(NL,3,1)*ANG1(NL,6,1)
      DUM3=ANG1(NL,1,1)*ANG1(NL,5,1)+ANG1(NL,2,1)*ANG1(NL,6,1)
      DUM4=ANG1(NL,1,1)*ANG1(NL,2,1)+ANG1(NL,5,1)*ANG1(NL,6,1)
      DUM5=ANG1(NL,3,1)*ANG1(NL,1,1)+ANG1(NL,4,1)*ANG1(NL,6,1)
      DUM6=ANG1(NL,3,1)*ANG1(NL,2,1)+ANG1(NL,5,1)*ANG1(NL,4,1)
      DUM7=ANG1(NL,1,1)*ANG1(NL,2,1)*ANG1(NL,4,1)
      DUM8=ANG1(NL,3,1)*ANG1(NL,1,1)*ANG1(NL,5,1)
      DUM9=ANG1(NL,3,1)*ANG1(NL,2,1)*ANG1(NL,6,1)
      DUM10=ANG1(NL,5,1)*ANG1(NL,4,1)*ANG1(NL,6,1)
      L=N1+3*(N2-1)
      D=DGN4*(DUM7+DUM8+DUM9+DUM10)
     +    +GN4(NL)*(CBIC(L,1)*DUM1+CBIC(L,2)*DUM2
     +         +CBIC(L,3)*DUM3+CBIC(L,4)*DUM4
     +         +CBIC(L,5)*DUM5+CBIC(L,6)*DUM6)
      RETURN
      END
