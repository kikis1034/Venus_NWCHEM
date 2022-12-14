      SUBROUTINE NMODE(NATOM,NDIS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         DRIVER FOR NORMAL MODE ANALYSIS
C         DOES NOT ALTER EITHER COORDINATES NOR THE ENERGY GRADIENT.
C
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/TABLEB/TABLE(42*NDA)
      COMMON/RANCOM/RANLST(100),ISEED3(8),IBFCTR
      COMMON/GPATHB/WM(NDA3),TEMP(NDP),AI1D(5),AI(2),BI(2),SYMM(5),
     *SYMA,SYMB,GTEMP(NDP),NFLAG(NDP),N1DR,N2DR
      COMMON/TESTIN/VRELO,INTST
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      COMMON/VECTB/VI(4),OAMI(4),AMAI(4),AMBI(4),ETAI,ERAI,ETBI,ERBI
      DIMENSION X(NDA3),GZS(NDA3)
C
C         WRITE RELEVANT INFORMATION IN CHECKPOINT FILE
C         (NORMAL MODE ANALYSIS OR REACTION PATH FOLLOWING)
C
      IF (NSELT.LT.0) THEN
         OPEN(50,FORM='UNFORMATTED')
         REWIND(50)
         WRITE(50)Q,P,QDOT,PDOT,TABLE,VRELO,RANLST,GTEMP,NFLAG,
     *            ISEED0,ISEED3,NX,NC,NTZ,INTST,NAST,IBFCTR,
     *            VI,OAMI,AMAI,AMBI,ETAI,ERAI,ETBI,ERBI
         CLOSE(50)
      ENDIF
C
C         SAVE COORDINATES AND GRADIENT
C
      I3N=3*NATOM
      DO I=1,I3N
         X(I)=Q(I+3*NDIS)
         GZS(I)=PDOT(I+3*NDIS)
      ENDDO
C
      CALL FMTRX(NATOM,NDIS,I3N)
C
C         RESTORE COORDINATES AND GRADIENT
C
      DO I=1,I3N
         Q(I+3*NDIS)=X(I)
         PDOT(I+3*NDIS)=GZS(I)
      ENDDO
      RETURN
      END
