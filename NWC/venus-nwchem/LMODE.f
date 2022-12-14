      SUBROUTINE LMODE(I,ENU,EDELTU,ENL,EDELTL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         SELECT THE DESIRED BOND LENGTH (RWANT) FOR A DIATOM.
C         THE VALUE FOR RWANT IS SELECTED BY MILLER'S SEMICLASSICAL
C         EXPRESSION WITH J=0.
C
C         SUBROUTINE LMODE PROVIDES THE NECESSARY PARAMETERS FOR
C         INITIAL EXCITATION OF THE MORSE STRETCH.  FOR I = 0, IT
C         CALCULATES THE EXCITATION ENERGY AND FOR I = 1, RWANT IS
C         ADDED TO THE MORSE STRETCH.
C
      COMMON/LMODEB/ENON,EDELTA,RWANT,PWANT,NEXM,NLEV,JFLAG
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/MORSEB/RMZ(ND02),B(ND02),D(ND02),N2J(ND02),N2K(ND02),
     *CM1(ND02),CM2(ND02),CM3(ND02),CM4(ND02),CM5(ND02)
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      save umass
C
      DE=D(NEXM)
      BETA=B(NEXM)
      R0=RMZ(NEXM)
      IF (I.EQ.0) THEN
         WA=W(N2J(NEXM))
         WB=W(N2K(NEXM))
         WT=WA+WB
         UMASS=WA*WB/WT
         WE=SQRT(2.0D0*DE*BETA*BETA/UMASS)
     *        /2.9979246D0/TWOPI/3.49757D-2
         XE=WE*WE/DE/4.0D0
         XE=XE*C1
         DUM=DBLE(NLEV)+0.5D0
         ENON=DUM*WE-DUM*DUM*XE
         UDUM=DUM+1.0D0
         ENU=UDUM*WE-UDUM*UDUM*XE
         EDELTA=0.5D0*(ENU-ENON)
         EENU=(UDUM+1.0D0)*WE-(UDUM+1.0D0)*(UDUM+1.0D0)*XE
         ENL=(DUM-1.0D0)*WE-(DUM-1.0D0)*(DUM-1.0D0)*XE
         EDELTU=0.5D0*(EENU-ENU)
         EDELTL=0.5D0*(ENON-ENL)
      ELSEIF (I.EQ.1) THEN
         RAND=RAND0(ISEED)
         DUM=TWOPI*RAND
         SCOOR=SIN(DUM)
         EWANT=ENON*C1
         EFACT=SQRT(DE*EWANT)
         DELR=-LOG((DE-EWANT)/(DE+EFACT*SCOOR))
         DELR=DELR/BETA
         RWANT=R0+DELR
         VWANT=DE*(1.0D0-EXP(-BETA*(RWANT-R0)))**2
         TWANT=EWANT-VWANT
         PWANT=SQRT(2.0D0*UMASS*TWANT)
         RAND=RAND0(ISEED)
         IF (RAND.LT.0.5D0) PWANT=-PWANT
      ENDIF
      RETURN
      END
