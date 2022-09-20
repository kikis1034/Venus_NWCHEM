      SUBROUTINE SELECT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      parameter(nstep=5000)
C
C         SELECT INITIAL CONDITIONS FOR COORDINATES AND MOMENTA
C
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      COMMON/TRANSB/TRANS,NREL
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,L(NDA),NAM
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/FORCES/NATOMS,I3N,NST,NM,NBB,NAA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/HFIT/PSCALA,PSCALB,VZERO
      COMMON/PRFLAG/NFQP,NCOOR,NFR,NUMR,NFB,NUMB,NFA,NUMA,NFTAU,NUMTAU,
     *NFTET,NUMTET,NFDH,NUMDH,NFHT,NUMHT
      COMMON/FINALB/EROTA,EROTB,EA(3),EB(3),AMA(4),AMB(4),AN,AJ,BN,BJ,
     *OAM(4),EREL,ERELSQ,BF,SDA,SDB,DELH(NDP),ANG(NDG),NFINAL
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/CHEMAC/WWA(NDA3),CA(NDA3yf,NDA3yf),AI(3),ENMTA,
     *AMPA(NDA3),WWB(NDA3),CB(NDA3yf,NDA3yf),BI(3),ENMTB,
     *AMPB(NDA3),SEREL,S,BMAX,TROTA,TROTB,ANQA(NDA3),ANQB(NDA3),
     *TVIBA,TVIBB,NROTA,NROTB,NOB
      COMMON/DIATB/NA,JA,NB,JB
      COMMON/VECTB/VI(4),OAMI(4),AMAI(4),AMBI(4),ETAI,ERAI,ETBI,ERBI
      COMMON/ARRAYS/A(NDA3yf,NDA3yf),DAA(NDA3),
     *B(NDA3yf,NDA3yf),DBB(NDA3)
      COMMON/EIGVL/EIG(NDA3yf)
      COMMON/SADDLE/EBAR,TBAR,EZERO,NBAR 
      common/vrscal/nsel,nscale,nequal,thermotemp,nrgd
      common/thermobath/nthermb,nrscl,nthmid(nda)
      common/mdsp/mdflag
      dimension qtemp(3*NDA)
      DIMENSION QMAXA(NDA3),QMINA(NDA3),PMAXA(NDA),QMAXB(NDA3),
     *QMINB(NDA3),PMAXB(NDA)
      DIMENSION ENLOW(nda3),ENMOD(nda3),SUM(0:NSTEP)
c     SAVE NMBAR,NMA,NMB,RMINA,RMAXA,DH,RMASSA,ENJA,PTESTA,ALA,
c    *  RMINB,RMAXB,RMASSB,ENJB,PTESTB,ALB
      save
C
   45 FORMAT(//5X,'SELECT:NORMAL MODE QUANTUM NUMBERS')
   46 FORMAT(5X,10F10.2)
  100 FORMAT(5X,'DIATOM A FREQUENCY =',F7.1,' CM-1, AND ENERGY =',
     *F7.2,' KCAL/MOL'/)
  106 FORMAT(5X,'SELECT:   JXA,JYA,JZA=',1P3D13.5,' H-BAR'/)
  117 FORMAT(5X,'REACTANT A')
  124 FORMAT(/5X,'SELECT:   EROTA =',F7.3,' KCAL/MOL'/
     *15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
  128 FORMAT(5X,'DIATOM B FREQUENCY =',F7.1,' CM-1, AND ENERGY =',
     *F7.2,' KCAL/MOL'/)
  136 FORMAT(5X,'SELECT:   JXB,JYB,JZB=',1P3D13.5,' H-BAR'/)
  158 FORMAT(//5X,'REACTANT B')
  162 FORMAT(/15X,'IMPACT PARAMETER=',F7.3,' A')
  163 FORMAT(/5X,'SELECT:   EROTB =',F7.3,' KCAL/MOL'/
     *15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
  196 FORMAT(/5X,'RELATIVE TRANSLATIONAL ENERGY SELECTED: ',F7.2,
     *' KCAL/MOL'/)
  198 FORMAT(/5X,'CHOSEN:   LX,LY,LZ =',1P3D13.5,' H-BAR')
  200 FORMAT(/5X,'CHOSEN:   EROT =',F7.3,' KCAL/MOL'/
     *15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
 201  FORMAT('THE ZERO POINT ENERGY CANNOT BE LARGER THAN',
     &           ' THE AVAILABLE ENERGY')
 202  FORMAT('ZERO POINT ENERGY ',F7.2)
 203  FORMAT('AVAILABLE ENERGY ',F7.2)
 205  FORMAT('INCREASE THE VALUE OF PARAMETER NSTEP TO AT LEAST ',I8)
C
C         NSELT=-1 PROGRAM DOES NORMAL MODE ANALYSIS
C         NSELT=0  Q'S AND P'S ARE READ IN
C         NSELT=1  PROGRAM FINDS MINIMUM ENERGY GEOMETRY
C                  (Q'S AND P'S ARE READ IN)
C         NSELT=2  CHOOSE INITIAL CONDITIONS FOR ONE OR TWO MOLECULES
C         NSELT=3  CHOOSE INITIAL CONDITIONS FROM POTENTIAL BARRIER
C               NBAR=1  MICROCANONICAL SAMPLING
C                       WITH FIXED REACTION COORDINATE ENERGY
C               NBAR=2  THERMAL SAMPLING
C               NBAR=3  MICROCANONICAL QUASICLASSICAL SAMPLING
C                       INCLUDING REACTION COORDINATE ENERGY
C                   ENMTA is the total available energy
C                         + zero point energy at TS.
C
C         NACT=1  ACTIVATE WITH ORTHANT SAMPLING
C         NACT=2  ACTIVATE WITH MICROCANONICAL NORMAL MODE SAMPLING
C         NACT=3  ACTIVATE WITH NORMAL MODE SAMPLING
C         NACT=4  ACTIVATE WITH LOCAL MODE EXCITATION
C         NACT=5  ACTIVATE WITH BOLTZMANN VIBRATIONAL DISTRIBUTION
C         NACT=6   SAME AS NBAR=3 (CHANGES NACT TO 6 WHEN NBAR=3)
C         NACT=7   MOLECULAR DYNAMICS SAMPLING BY RESCALING VELOCITIES
c   rescale the velocities of the system according to given temp.
c   reference:  "MOLECULAR DYNAMICS SIMULATION" by JIM HAILE  p.458
c   Only fragment B can have this when it is surface.
C
c   mdflag is used to avoid the gaussian call while doing the md sampling.
c   mdflag=0: call DOGAUSS
c   mdflag=1: do not call DOGAUSS
c
      mdflag=0
      IF (NSELT.NE.2.AND.NSELT.NE.3) THEN
         READ(5,*)(Q(I),I=1,I3N)
         if(nselt.eq.0) then
           READ(5,*)(P(I),I=1,I3N)
         else
           do i=1,i3n
             p(i)=0.d0
           enddo
         endif
         CALL DVDQ
         CALL ENERGY
         GOTO 999
      ENDIF
C
C         IN THE MAIN PROGRAM THE COORDINATES Q ARE SET EQUAL TO QZ
C         AND THE MOMENTA P ARE SET EQUAL TO ZERO (FOR NSELT=2).
C
C         INITIAL CONDITIONS FOR REACTANT A
C
      IF (NATOMA(1).LE.1) THEN
C
C             IF A IS AN ATOM ZERO ITS Q, P ELEMENTS
C
         DO K=1,3
            J=3*LA(1,1)-3+K
            Q(K)=0.0D0
            P(K)=0.0D0
            QQ(K)=0.0D0
            PP(K)=0.0D0
         ENDDO
         GOTO 126
      ENDIF
c
C             DISPLACE REACTANT B BY 1000.0 ANGSTROMS
C
      WRITE(6,117)
      N=NATOMB(1)
      IF (N.NE.0) THEN
         DO I=1,N
            J3=3*LB(1,I)
            K3=3*I
            Q(J3)=QZB(1,K3)+1000.0D0
         ENDDO
      ENDIF
      
C
C             ENERGY REFERENCE FOR SEPARATED REACTANTS
C
      CALL DVDQ
      CALL ENERGY
      DH=V
C
C             SELECT INITIAL Q'S AND P'S FOR REACTANT A
C
      IF (NATOMA(1).GT.2.or.nlina.eq.0) GOTO 120
      L(1)=LA(1,1)
      L(2)=LA(1,2)
C
C          DIATOM A IS TREATED SEMICLASSICALLY
C
      IF (NSFLAG.NE.1) THEN
         N=NATOMA(1)
         CALL NMODE(N,0)
         DUM=EIG(6)
         ENJA=(DBLE(NA)+0.5D0)*DUM*CM2CAL*C1
         CALL INITEBK(NA,JA,RMINA,RMAXA,DH,RMASSA,ENJA,PTESTA,ALA)
         SDUM=ENJA/C1
         WRITE(6,100)DUM,SDUM
      ENDIF
C
C             SELECT INITIAL RELATIVE COORDINATE AND MOMENTUM.
C
      DUM=ALA**2/2.0D0/RMASSA
  102 RAND=RAND0(ISEED)
      R=RMINA+(RMAXA-RMINA)*RAND
      Q(3*L(1))=-0.5D0*R
      Q(3*L(1)-1)=0.0D0
      Q(3*L(1)-2)=0.0D0
      Q(3*L(2))= 0.5D0*R
      Q(3*L(2)-1)=0.0D0
      Q(3*L(2)-2)=0.0D0


      CALL DVDQ


      CALL ENERGY
      VDUM=(V-DH)*C1
      SUMM=ENJA-DUM/R**2-VDUM
      IF (SUMM.LE.0.0D0) THEN
         SUMM=0.0D0
         PR=0.0D0
      ELSE
         PR=SQRT(2.0D0*RMASSA*SUMM)
         SDUM=PTESTA/PR
         RAND=RAND0(ISEED)
         IF (SDUM.LT.RAND) GOTO 102
      ENDIF
      RAND=RAND0(ISEED)
      IF (RAND.LT.0.5D0) PR=-PR
C
C             CHOOSE INITIAL CARTESIAN COORDINATES AND MOMENTA, AND
C             ANGULAR MOMENTUM.  DIATOM LIES ALONG THE X-AXIS.
C             THEN RANDOMLY ROTATE THE CARTESIAN COORDINATES AND
C             MOMENTA IN THE CENTER OF MASS FRAME.
C
      CALL HOMOQP(R,PR,ALA,AMA,RMASSA,AI)
      AMAI(1)=AMA(1)/C7
      AMAI(2)=AMA(2)/C7
      AMAI(3)=AMA(3)/C7
      AMAI(4)=AMA(4)/C7
      WRITE(6,106)AMAI(1),AMAI(2),AMAI(3)
      GOTO 126
C
  120 CONTINUE
C
C             FRAGMENT A IS A POLYATOMIC
C
      IF (NACTA.EQ.1) GOTO 123
C
C             CALCULATE NORMAL MODE EIGENVALUES AND EIGENVECTORS
C
      IF(NSFLAG.EQ.1)THEN
        IF(NACTA.EQ.5)THEN
          GOTO 113
        ELSE
          GOTO 123
        ENDIF
      ENDIF
C     IF (NSFLAG.EQ.1.AND.NACTA.NE.5) GOTO 123
C     IF (NSFLAG.EQ.1.AND.NACTA.EQ.5) GOTO 113
      N=NATOMA(1)
      K=3*N
      M=6-NLINA
      NMA=K-M
C
C             TRANSITION STATE ONLY HAS 3N-7 NORMAL MODES FOR A 
C             NONLINEAR MOLECULE, AND 3N-6 MODES FOR A LINEAR ONE
C
      IF (NSELT.EQ.3) THEN
         NMBAR=NMA-1
         IBARR=1
      ELSE
         NMBAR=NMA
         IBARR=0
      ENDIF
C
      CALL NMODE(N,0)
      DO I=1,NMBAR
         WWA(I)=EIG(I+M+IBARR)*C6
         DO J=1,K
            CA(J,I)=A(I+M+IBARR,J)
         ENDDO
      ENDDO
      IF (NSELT.EQ.3) THEN
         WWA(NMA)=EIG(1)*C6
         DO J=1,K
            CA(J,NMA)=A(1,J)
         ENDDO
      ENDIF
C
C     THIS PART IS FOR MICROCANONICAL SAMPLING.
C     CHECK THAT THE TOTAL VIBRATIONAL ENERGY DOES NOT EXCEED THE
C     MAXIMUM TOTAL ENERGY THAT IS USED IN BEYER-SWINEHARDT COUNTING
C
      IF(NACTA.EQ.6) THEN
         NBAR=1                   ! FIXED ENERGY IN REACTION COORDINATE
         DUM=ENMTA*cal2cm
         ZPE=0.0D0
         DO I=1,NMBAR
            ENMOD(I)=WWA(I)/C6
            ZPE=ZPE+WWA(I)/C6
            ENLOW(I)=ENMOD(I)*0.5D0*CM2CAL*C1
         ENDDO
         ZPE=ZPE*0.5D0
         DUM=DUM-ZPE
         STEP=ENMOD(1)*0.1D0      ! ENMOD(1) IS THE LOWEST FREQUENCY
         IF(DUM/STEP.GE.DBLE(NSTEP)) THEN
            WRITE(6,205)IDNINT(DUM/STEP)+1
            STOP
         ENDIF
         IF(ZPE/cal2cm.GT.ENMTA) THEN
            WRITE(6,201)
            WRITE(6,202)ZPE*CM2CAL
            WRITE(6,203)ENMTA
            STOP
         ENDIF
         DUMNAC6=DUM
      ENDIF
C
      IF ((NACTA.EQ.2).OR.(NACTA.EQ.6)) GOTO 123
C
C             CHOOSE NORMAL MODE QUANTUM NUMBERS FROM A THERMAL
C             DISTRIBUTION IF NACT=5
C
  113 CONTINUE
      IF (NACTA.EQ.5) THEN
         CALL THRMAN(WWA,ANQA,TVIBA,NMBAR)
         WRITE(6,45)
         WRITE(6,46)(ANQA(I),I=1,NMBAR)
         WRITE(9,46)(ANQA(I),I=1,NMBAR)
         call flush(6)
      ENDIF      
C
C             CALCULATE NORMAL MODE ENERGIES AND AMPLITUDES
C
      ENMTA=0.0D0
      DO I=1,NMBAR
         DUM=(ANQA(I)+0.5D0)*WWA(I)/(C6*cal2cm)
         ENMTA=ENMTA+DUM
         DUM=DUM*C1
         AMPA(I)=SQRT(2.0D0*DUM)/WWA(I)
      ENDDO
C
C             CHOOSE THE ANGULAR MOMENTUM VECTOR
C
  123 CONTINUE
      CALL ROTEN(AMA,AI,TROTA,EROTA,NROTA,NLINA)
C
C             SAVE INITIAL ROTATIONAL ANGULAR MOMENTUM
C
      AMAI(1)=AMA(1)/C7
      AMAI(2)=AMA(2)/C7
      AMAI(3)=AMA(3)/C7
      AMAI(4)=SQRT(AMAI(1)**2+AMAI(2)**2+AMAI(3)**2)
      WRITE(6,124)EROTA,AMAI(1),AMAI(2),AMAI(3)
C
C             CHOOSE THE INITIAL COORDINATES AND MOMENTA
C
      N=NATOMA(1)
      DO I=1,N
         L(I)=LA(1,I)
         J3=3*L(I)
         J2=J3-1
         J1=J2-1
         K3=3*I
         K2=K3-1
         K1=K2-1
         QZ(J1)=QZA(1,K1)
         QZ(J2)=QZA(1,K2)
         QZ(J3)=QZA(1,K3)
      ENDDO
      DUM1=WTA(1)
      ETAI=EROTA+ENMTA
      IF (NACTA.EQ.1) THEN
         CALL ORTHAN(AMA,DUM1,ENMTA,ETAI,QMAXA,QMINA,PMAXA,
     *               PSCALA,ERAI,N)
         GOTO 126
      ENDIF
C
C             CALCULATE NORMAL MODE ENERGIES AND AMPLITUDES FOR A
C             MICROCANONICAL ENSEMBLE
C
      IF (NACTA.EQ.2) THEN
         DUM=ENMTA*C1
         NN=NMBAR-1
         DO I=1,NN
            RAND=RAND0(ISEED)
            SDUM=1.0D0/DBLE(NMBAR-I)
            SDUM=DUM*(1.0D0-RAND**SDUM)
            DUM=DUM-SDUM
            AMPA(I)=SQRT(2.0D0*SDUM)/WWA(I)
         ENDDO
         AMPA(NMBAR)=SQRT(2.0D0*DUM)/WWA(NMBAR)
      ENDIF
C
      IF (NACTA.EQ.6) THEN
         DUM=DUMNAC6     ! TOTAL ENERGY ABOVE ZPE IN 1/CM
         DO II=1,NMBAR   ! LOWEST TO HIGHEST FREQUENCY FOR EFFICIENCY
            IMAXST=IDINT(DUM/ENMOD(II)) ! MAXIMUM QUANTUM NO.
            IF(II.NE.NMBAR) THEN
               STEP=ENMOD(II+1)*0.1D0
C
C     THIS IS THE BEYER-SWINEHARDT ALGORITHM
C
               IEND=IDNINT(DUM/STEP)
               DO J=0,IEND
                  SUM(J)=1.0D0
               ENDDO
               DO J=II+1,NMBAR
                  ISTART=IDNINT(ENMOD(J)/STEP)
                  DO K=ISTART,IEND
                     SUM(K)=SUM(K)+SUM(K-ISTART)
                  ENDDO
               ENDDO
C     
C     
               PROB=0.0D0
               DO J=IMAXST,0,-1
                  ENDUM=DUM-DBLE(J)*ENMOD(II) ! ENERGY IN REST OF MOLECULE
                  IDDUM=IDNINT(ENDUM/STEP)
                  PROB=PROB+SUM(IDDUM)
               ENDDO
            ELSE
               STEP=ENMOD(II)
               PROB=0.0D0
               DO J=IMAXST,0,-1
                  ENDUM=DUM-DBLE(J)*ENMOD(II)
                  IDDUM=IDNINT(ENDUM/STEP)
                  SUM(IDDUM)=1.0D0
                  PROB=PROB+SUM(IDDUM)
               ENDDO
            ENDIF
            RAND=RAND0(ISEED)*PROB
            IF(RAND.GT.PROB)RAND=PROB
            PROB=0.0D0
            DO J=IMAXST,0,-1
               ENDUM=DUM-DBLE(J)*ENMOD(II)
               IDDUM=IDNINT(ENDUM/STEP)
               PROB=SUM(IDDUM)+PROB
               IF(RAND.LE.PROB) THEN
                  SDUM=ENMOD(II)*DBLE(J)
                  GOTO 1115
               ENDIF
            ENDDO
 1115       CONTINUE
            DUM=DUM-SDUM
            SDUM=SDUM*CM2CAL*C1
            AMPA(II)=SQRT(2.0D0*(SDUM+ENLOW(II)))/WWA(II)
         ENDDO
C     GIVE THE REMAINDER OF THE ENERGY TO THE REACTION COORDINATE
         EBAR=DUM*CM2CAL
         ETAI=ETAI-EBAR
      ENDIF
C
      CALL INITQP(WWA,AMPA,CA,AMA,DUM1,ETAI,EROTA,AI,ERAI,N,NMBAR)
C
C         INITIAL CONDITION FOR REACTION COORDINATE
C
      IF (NSELT.EQ.3) CALL BAREXC(DUM1,CA,AMA,ERAI,N,NMA)
C
C         INITIAL CONDITIONS FOR REACTANT B
C
  126 CONTINUE
      IF (NATOMB(1).EQ.0) THEN
         IF (NATOMA(1).NE.2) GOTO 999
         CALL DVDQ
         CALL ENERGY
         GOTO 999
      ENDIF
      IF (NATOMB(1).GT.1) GOTO 129
C
C             IF FRAGMENT B IS AN ATOM ZERO ITS Q ELEMENTS
C
      J3=3*LB(1,1)
      J2=J3-1
      J1=J2-1
      Q(J1)=0.0D0
      Q(J2)=0.0D0
      Q(J3)=0.0D0
      GOTO 160
C
C             SELECT INITIAL Q'S AND P'S FOR REACTANT B.
C
C             SAVE Q'S AND P'S OF A IN TEMPORARY STORAGE. (in fact,
C             this is unnecessory as we can store them from QQ,PP)
C             SET P ARRAY TO ZERO AND EQUATE THE Q AND QZ ARRAYS FOR A
C
  129 CONTINUE
      WRITE(6,158)
      N=NATOMA(1)
      DO I=1,N
         J3=3*LA(1,I)
         J2=J3-1
         J1=J2-1
         K3=3*I
         K2=K3-1
         K1=K2-1
         Q(J1)=QZA(1,K1)
         Q(J2)=QZA(1,K2)
         Q(J3)=QZA(1,K3)+1000.0D0
         P(J1)=0.0D0
         P(J2)=0.0D0
         P(J3)=0.0D0
      ENDDO
C
C             RESET Q ARRAY TO QZ ARRAY FOR B
C
      N=NATOMB(1)
      DO I=1,N
         J3=3*LB(1,I)
         K3=3*I
         Q(J3)=QZB(1,K3)
      ENDDO
C
      IF (NATOMB(1).GT.2.or.nlinb.eq.0) GOTO 150
      L(1)=LB(1,1)
      L(2)=LB(1,2)
C
C          DIATOM B IS TREATED SEMICLASSICALLY
C
      IF (NSFLAG.NE.1) THEN
         I=NATOMA(1)
         N=NATOMB(1)
         CALL NMODE(N,I)
         DUM=EIG(6)
         ENJB=(DBLE(NB)+0.5D0)*DUM*CM2CAL*C1
         CALL INITEBK(NB,JB,RMINB,RMAXB,DH,RMASSB,ENJB,PTESTB,ALB)
         SDUM=ENJB/C1
         WRITE(6,128)DUM,SDUM
      ENDIF
C
C             SELECT INITIAL RELATIVE COORDINATE AND MOMENTUM
C
      DUM=ALB**2/2.0D0/RMASSB
  132 RAND=RAND0(ISEED)
      R=RMINB+(RMAXB-RMINB)*RAND
      Q(3*L(1))=-0.5D0*R
      Q(3*L(1)-1)=0.0D0
      Q(3*L(1)-2)=0.0D0
      Q(3*L(2))= 0.5D0*R
      Q(3*L(2)-1)=0.0D0
      Q(3*L(2)-2)=0.0D0
      CALL DVDQ
      CALL ENERGY
      VDUM=(V-DH)*C1
      SUMM=ENJB-DUM/R**2-VDUM
      IF (SUMM.LE.0.0D0) THEN
         SUMM=0.0D0
         PR=0.0D0
      ELSE
         PR=SQRT(2.0D0*RMASSB*SUMM)
         SDUM=PTESTB/PR
         RAND=RAND0(ISEED)
         IF (SDUM.LT.RAND) GOTO 132
      ENDIF
      RAND=RAND0(ISEED)
      IF (RAND.LT.0.5D0) PR=-PR
C
C             CHOOSE INITIAL CARTESIAN COORDINATES AND MOMENTA, AND
C             ANGULAR MOMENTUM.  DIATOM LIES ALONG THE X-AXIS.
C             THEN RANDOMLY ROTATE THE CARTESIAN COORDINATES AND
C             MOMENTA IN THE CENTER OF MASS FRAME.
C
  142 CONTINUE
      CALL HOMOQP(R,PR,ALB,AMB,RMASSB,BI)
      AMBI(1)=AMB(1)/C7
      AMBI(2)=AMB(2)/C7
      AMBI(3)=AMB(3)/C7
      AMBI(4)=AMB(4)/C7
      WRITE(6,136)AMBI(1),AMBI(2),AMBI(3)
      GOTO 155
C
  150 CONTINUE
C
C             FRAGMENT B IS A POLYATOMIC
C
      IF (NACTB.EQ.1) GOTO 153
      if(nactb.eq.7) goto 111
C
C             CALCULATE NORMAL MODE EIGENVALUES AND EIGENVECTORS
C
      IF (NSFLAG.EQ.1.AND.NACTB.NE.5) GOTO 153
      IF (NSFLAG.NE.1.OR.NACTB.NE.5) THEN
         N=NATOMB(1)
         K=3*N
         M=6-NLINB
         NMB=K-M
         I=NATOMA(1)
         CALL NMODE(N,I)
         DO I=1,NMB
            WWB(I)=EIG(I+M)*C6
            DO J=1,K
               CB(J,I)=A(I+M,J)
            ENDDO
         ENDDO
      ENDIF
C
      IF (NACTB.EQ.2) GOTO 153
C
C             CHOOSE NORMAL MODE QUANTUM NUMBERS FROM A THERMAL
C             DISTRIBUTION IF NACT=5
C
      IF (NACTB.EQ.5) THEN
         CALL THRMAN(WWB,ANQB,TVIBB,NMB)
         WRITE(6,45)
         WRITE(6,46)(ANQB(I),I=1,NMB)
         WRITE(9,46)(ANQB(I),I=1,NMB)
      ENDIF
C             CALCULATE NORMAL MODE ENERGIES AND AMPLITUDES
C
      ENMTB=0.0D0
      DO I=1,NMB
         DUM=(ANQB(I)+0.5D0)*WWB(I)/(C6*cal2cm)
         ENMTB=ENMTB+DUM
         DUM=DUM*C1
         AMPB(I)=SQRT(2.0D0*DUM)/WWB(I)
      ENDDO
C
C             CHOOSE THE ANGULAR MOMUNTUM VECTOR
C
  153 CONTINUE
      CALL ROTEN(AMB,BI,TROTB,EROTB,NROTB,NLINB)
C
C             SAVE INITIAL ROTATIONAL ANGULAR MOMENTUM
C
      AMBI(1)=AMB(1)/C7
      AMBI(2)=AMB(2)/C7
      AMBI(3)=AMB(3)/C7
      AMBI(4)=SQRT(AMBI(1)**2+AMBI(2)**2+AMBI(3)**2)
      WRITE(6,163)EROTB,AMBI(1),AMBI(2),AMBI(3)
C
C             CHOOSE THE INITIAL COORDINATES AND MOMENTA
C
  111 continue
      N=NATOMB(1)
      DO I=1,N
         L(I)=LB(1,I)
         J3=3*L(I)
         J2=J3-1
         J1=J2-1
         K3=3*I
         K2=K3-1
         K1=K2-1
         QZ(J1)=QZB(1,K1)
         QZ(J2)=QZB(1,K2)
         QZ(J3)=QZB(1,K3)
      ENDDO
      if(nactb.eq.7) goto 155
      DUM1=WTB(1)
      ETBI=EROTB+ENMTB
      IF (NACTB.EQ.1) THEN
         CALL ORTHAN(AMB,DUM1,ENMTB,ETBI,QMAXB,QMINB,PMAXB,
     *               PSCALB,ERBI,N)
         GOTO 155
      ENDIF
C
C             CALCULATE NORMAL MODE ENERGIES AND AMPLITUDES FOR A
C             MICROCANONICAL ENSEMBLE
C
      IF (NACTB.EQ.2) THEN
         DUM=ENMTB*C1
         NN=NMB-1
         DO I=1,NN
            RAND=RAND0(ISEED)
            SDUM=1.0D0/DBLE(NMB-I)
            SDUM=DUM*(1.0D0-RAND**SDUM)
            DUM=DUM-SDUM
            AMPB(I)=SQRT(2.0D0*SDUM)/WWB(I)
         ENDDO
         AMPB(NMB)=SQRT(2.0D0*DUM)/WWB(NMB)
      ENDIF
      CALL INITQP(WWB,AMPB,CB,AMB,DUM1,ETBI,EROTB,BI,ERBI,N,NMB)
C
C         RESTORE Q'S AND P'S FROM TEMPORARY STORAGE
C         (in fact we can store from QQ, PP)
C
  155 CONTINUE
      call flush(6)
c--------------------------------------------------------------------
c   initialize a phase space point by randomly assigning velocities
c   at equalibrium coordinates
      if (nactb.eq.7) then
c
         mdflag=1
c
c       store the coordinates and momenta of reactant A first
c
         n=natoma(1)
c
         do i=1,n
            j3=3*la(1,i)
            j2=j3-1
            j1=j2-1
            qtemp(j3)=q(j3)
            qtemp(j2)=q(j2)
            qtemp(j1)=q(j1)
            p(j3) = 0.0d0
            p(j2) = 0.0d0
            p(j1) = 0.0d0
         end do
c
         ncoororg=ncoor
         ncoor=0
         n=natomb(1)-nrgd
c
c   randomly assigning velocities firstly then will be rescaled
c   to the temperature wanted. two steps
c      1. Generate Random Numbers: Generate 3N random numbers from a normal 
c         distribution with zero mean and unit variance. 
c      2. Generate momenta: Multiply these random numbers by sqrt(mkT) to obtain 
c         the x, y, and z components of momenta for the N particles
c   desket = sqrt(kb*T), kb - Boltzmann Const. T - temperature
c
	 desket=sqrt(0.00198717d0*thermotemp*C1) 
         do i = 1,n
            j3=3*lb(1,i)
            j2=j3-1
            j1=j2-1
	    j=lb(1,i)
	    p(j1)=gasdev_venus(ISEED)*desket*sqrt(w(j))
	    p(j2)=gasdev_venus(ISEED)*desket*sqrt(w(j))
	    p(j3)=gasdev_venus(ISEED)*desket*sqrt(w(j))
c            p(J1) = (0.5D0-RAND0(ISEED)) * W(lb(1,I))
c            p(j2) = (0.5D0-RAND0(ISEED)) * W(lb(1,I))
c            p(j3) = (0.5D0-RAND0(ISEED)) * W(lb(1,I))
c            pmag = sqrt( p(j1)**2 + p(j2)**2 + p(j3)**2 )
c            p(J1) = p(J1) / pmag 
c            p(J2) = p(J2) / pmag
c            p(J3) = p(J3) / pmag
         enddo
c
c   rescale the velocities of the system according to given temp.
c   reference:  "MOLECULAR DYNAMICS SIMULATION" by JIM HAILE  p.458
c   k = 1.38066 * 10(-23) J/K = 0.00198624 kcal/mol K
c
	  if (nscale .gt. 0) then
	     call thermo(0)
	  else
c
c		using thermo-bath to heat up the system
c
            if (nthermb.gt.0 .and. mod(nc,nrscl).eq.0) then
               call thermbath
	    endif
	  endif
c
c         perform the first six integration cycles required for adamsm
c
         nc=0
         call dvdq
         call energy
         call parti
c         if (nsel.eq.1) call gwrite
c
         do i=1,6
            nc=nc+1
            call rungek
            if (nscale .gt. 0) then
               call thermo(nc)
            else
               if (nthermb.gt.0 .and. mod(nc,nrscl).eq.0) then
                  call thermbath
	       endif
            endif
         enddo
c
c        run the trajectory with no A-B translational energy to
c        get to equilibrium
c
         nse=nscale+nequal-6  
         do i=1,nse
            nc=nc+1
            call adamsm
            call energy
             if (nc.lt.nscale) then
	       call thermo(nc)
	    else
               if (nthermb.gt.0 .and. mod(nc,nrscl).eq.0) then
                  call thermbath
	       endif
            endif
         enddo
c
c	 do i=1, n
c	    j3=3*lb(1,i)
c           j2=j3-1
c            j1=j2-1
c	    j=lb(1,i)
c            write(33,*)p(j1),p(j1)*p(j1)/w(j)/2/c1
c            write(33,*)p(j2),p(j2)*p(j2)/w(j)/2/c1
c	    write(33,*)p(j3),p(j3)*p(j3)/w(j)/2/c1
c	 enddo
c
c     set the coordinates and momenta for reactant A back to
c     the initial choice recall that the momenta were equated
c     to zero and stored in pp.  The coordinates have the z
c     component of fragrment b at 1000 a above the surface.
c     This will be corrected lated in the code.
c
        WT=WTB(NPATH)
        N=NATOMB(NPATH)
        CALL CENMAS(WT,QCMB,VCMB,N)
c
        n=natoma(1)
        do i=1,n
           j3=3*la(1,i)
           j2=j3-1
           j1=j2-1
           p(j3)=pp(j3)
           p(j2)=pp(j2)
           p(j1)=pp(j1)
           q(j3)=qtemp(j3)
           q(j2)=qtemp(j2)
           q(j1)=qtemp(j1)
        end do
c
        write(6,*)'EQUALIBRATION FOR NACTB EQ 7 IS NOW OVER'
        call flush(6)
c
        nfinal=0
        ncoor=ncoororg
        nc=0
        mdflag=0
c
      endif
c---------------------------------------------------------------
C
      N=NATOMA(1)
      DO I=1,N
         DO K=1,3
            J=3*LA(1,I)-3+K
            Q(J)=QQ(J)
            P(J)=PP(J)
         ENDDO 
      ENDDO
  160 CONTINUE
C
C         CHOOSE MOMENTA AND RELATIVE POSITIONS FOR REACTANTS.
C         USE THESE TO FIND REACTANTS' INITIAL P AND Q
C
C             GAS/SURFACE COLLISION
C
      IF (NSURF.NE.0) THEN
         CALL SURF(NSURF)
      ELSE
C
C             GAS PHASE COLLISION
C             SELECT IMPACT PARAMETER.  FIX POSITIONS OF A AND
C             B FOR CHOSEN IMPACT PARAMETER(B) AND SEPARATION(S).
C
         SB=BMAX
         IF (NOB.NE.1) THEN
            RAND=RAND0(ISEED)
            SB=BMAX*SQRT(RAND)
         ENDIF
         WRITE(6,162)SB
         DUM1=SQRT(S*S-SB*SB)
         N=NATOMB(1)
         DO 164 I=1,N
            J=3*LB(1,I)
            Q(J)=Q(J)+DUM1
            Q(J-1)=Q(J-1)+SB
  164    CONTINUE
C
C             IF NREL = 0 FOR GAS PHASE COLLISION
C             CHOOSE RELATIVE ENERGY FROM BOLTZMANN DISTRIBUTION
C
         IF (NREL.EQ.0) THEN
            DUM=GAMA(2,ISEED)
            SEREL=0.00198717D0*DUM*TRANS
            WRITE(6,196)SEREL
            SEREL=SEREL*C1
         ENDIF
C
C             ADD RELATIVE TRANSLATIONAL ENERGY FOR GAS PHASE
C             COLLISION
C
         WT=WTA(1)+WTB(1)
         SDUM=WTA(1)*WTB(1)/WT
         DUM=SQRT(2.0D0*SEREL/SDUM)
         VELA=DUM*WTB(1)/WT
         VELB=VELA-DUM
         N=NATOMA(1)
         DO I=1,N
            J=3*LA(1,I)
            P(J)=P(J)+VELA*W(LA(1,I))
         ENDDO
         N=NATOMB(1)
         DO I=1,N
            J=3*LB(1,I)
            P(J)=P(J)+VELB*W(LB(1,I))
         ENDDO
         CALL DVDQ
         CALL ENERGY
C
C             SAVE INITIAL RELATIVE VELOCITY AND ORBITAL ANGULAR
C             MOMENTUM FOR GAS PHASE COLLISION
C
         VI(1)=0.0D0
         VI(2)=0.0D0
         VI(3)=DUM
         VI(4)=DUM
         OAMI(1)=-SB*DUM*SDUM/C7
         OAMI(2)=0.0D0
         OAMI(3)=0.0D0
         OAMI(4)=ABS(OAMI(1))
         WRITE(6,198)(OAMI(I),I=1,3)
      ENDIF
C
  999 CONTINUE
C
C             SAVE CHOSEN INITIAL ROTATIONAL ANGULAR MOMENTUM
C
      IF (NATOMA(1).GT.1) THEN 
         AMAI(1)=AMA(1)/C7
         AMAI(2)=AMA(2)/C7
         AMAI(3)=AMA(3)/C7
         AMAI(4)=SQRT(AMAI(1)**2+AMAI(2)**2+AMAI(3)**2)
         WRITE(6,200)ERAI,AMAI(1),AMAI(2),AMAI(3)
      ENDIF
      IF (NATOMB(1).GT.1) THEN
         AMBI(1)=AMB(1)/C7
         AMBI(2)=AMB(2)/C7
         AMBI(3)=AMB(3)/C7
         AMBI(4)=SQRT(AMBI(1)**2+AMBI(2)**2+AMBI(3)**2)
         WRITE(6,200)ERBI,AMBI(1),AMBI(2),AMBI(3)
      ENDIF
C
      NSFLAG=1
      RETURN
      END
