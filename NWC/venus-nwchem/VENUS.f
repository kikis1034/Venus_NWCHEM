      PROGRAM VENUS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      character*80 title1,title2
C***********************************************************************
C                                                                      *
C                              VENUS05                                 *
C                                                                      *
C              A GENERAL CHEMICAL DYNAMICS COMPUTER PROGRAM            *
C                                                                      *
C                                 BY                                   *
C                                                                      *
C          William L. Hase, Kim Bolton, Pascal de Sainte Claire,       *
C        Ronald J. Duchovic, Xiche Hu,Andrew Komornicki,Guosheng Li,   *
C       Kieran F. Lim, Da-hong Lu,Gilles H. Peslherbe, Kihyung Song,   *
C         Kandadai N. Swamy, Scott R. Vande Linde,Antonio Varandas,    *
C                       Haobin Wang, and Ralph J. Wolf                 *
C                                                                      *
C                            December, 2004                                 *
C                                                                      *
C***********************************************************************
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/INTEGR/ATIME,NI,NID
      COMMON/TABLEB/TABLE(42*NDA)
      COMMON/PRFLAG/NFQP,NCOOR,NFR,NUMR,NFB,NUMB,NFA,NUMA,NFTAU,NUMTAU,
     *NFTET,NUMTET,NFDH,NUMDH,NFHT,NUMHT
      COMMON/PARRAY/KR(300),JR(300),KB(300),MB(300),IB(300),IA(300),
     *ITAU(300),ITET(300),IDH(300),IHT(300)
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      COMMON/TRANSB/TRANS,NREL
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/HFIT/PSCALA,PSCALB,VZERO
      COMMON/PSN2/PESN2,GA,RA,RB
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/CUBEB/S3(4),DS3(4),CBIC(15,6),ANG1(20,6,4),GN4(20)
      COMMON/TESTIN/VRELO,INTST
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      COMMON/TESTSN2/GAO,NSAD,NCBA,NCAB,IBAR
      COMMON/FINALB/EROTA,EROTB,EA(3),EB(3),AMA(4),AMB(4),AN,AJ,BN,BJ,
     *OAM(4),EREL,ERELSQ,BF,SDA,SDB,DELH(NDP),ANG(NDG),NFINAL
      COMMON/CHEMAC/WWA(NDA3),CA(NDA3yf,NDA3yf),AI(3),ENMTA,
     *AMPA(NDA3),WWB(NDA3),CB(NDA3yf,NDA3yf),BI(3),ENMTB,
     *AMPB(NDA3),SEREL,S,BMAX,TROTA,TROTB,ANQA(NDA3),ANQB(NDA3),
     *TVIBA,TVIBB,NROTA,NROTB,NOB
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,LL(NDA),NAM
      COMMON/DIATB/NNA,JA,NNB,JB
      COMMON/RSTART/HINC,NPTS
      COMMON/RKUTTA/RAA1,RA1,RA2,RA3,RB1,RB2,RB3,RC1,RC2
      COMMON/LMODEB/ENON,EDELTA,RWANT,PWANT,NEXM,NLEV,JFLAG
      COMMON/RANCOM/RANLST(100),ISEED3(8),IBFCTR
      COMMON/GPATHB/WM(NDA3),TEMP(NDP),AI1D(5),AAI(2),BBI(2),SYMM(5),
     *SYMA,SYMB,GTEMP(NDP),NFLAG(NDP),N1DR,N2DR
      COMMON/SYBB/TITLE1,TITLE2,SYBTI
      COMMON/SADDLE/EBAR,TBAR,EZERO,NBAR 
      COMMON/INERT/UXX,UXY,UXZ,UYY,UYZ,UZZ,AIXX,AIXY,AIXZ,AIYY,
     *AIYZ,AIZZ
      COMMON/EIGVL/EIG(NDA3yf)
      COMMON/VECTB/VI(4),OAMI(4),AMAI(4),AMBI(4),ETAI,ERAI,ETBI,ERBI
      COMMON/GAUSS/ICHRG,IMULTP,IAN(NDAyf),VGAUSS
      COMMON/SURFB/NN1,NN2,NN3,NN4,THTA,PHI,NCHI,CHI,RX0,RY0,RZ0,
     *NTHET,THET,NPHI1,PHI1,NPHI2,PHI2
      COMMON/VMAXB/QVMAX(NDA3),PVMAX(NDA3),VMAX,NCVMAX
      common/gdd/diatm(mxdiatom,4),gstop(mxpth,6),iatmp(mxpth,nda),
     *  idiatm(mxdiatom,2),iatsp(mxpth,3),isotope(nda),npth,ndiatom,
     *  nmem,ndisk,ns,nchkp,nhess,ncpu,nisotope,method
      common/vrscal/nsel,nscale,nequal,thermotemp,nrgd
      common/thermobath/nthermb,nrscl,nthmid(nda)
      common/mdsp/mdflag
      common/ghessb/FA(NDA3,NDA3),trad,rmin,dt,nstep,nhesup,
     * nhessf,nip
C
      DIMENSION ERAVA(500),ERAVB(500),QCM(3),VCM(3),AM(4)
      DIMENSION ETIM(2000),ESAV(2000),ESQ(2000),NEVIB(2000),
     *NEVIBU(2000),NEVIBL(2000),NEVIBN(2000)
      DIMENSION EBM(50),ENM(50),EBSAV(2,2000),EBSQ(2,2000),MNLM(50)
      DIMENSION QDUM(NDA3),TABDUM(42*NDA),GDUM(NDP),NFDUM(NDP)
      DOUBLE PRECISION, ALLOCATABLE:: ENSAV(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: ENSQ(:,:)
C
      CHARACTER*30 ADATE,method
      DATA RHO/1.0D-8/
C     SAVE TUSER1,TUSER2,TTOTAL1,TTOTAL2
      save lll,xl,niter

C
C      QM PARAMETERS
C

      character*2 nameatm(54)
      save nameatm
      data nameatm/'H','He','Li','Be','B','C','N','O','F','Ne','Na',
     * 'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
     * 'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
     * 'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
     * 'Sn','Sb','Te','I','Xe'/

      COMMON /QMINFO/ qm_choice

      logical initialize
      character*256 qm_ini_file
      character*256 theory
      character*256 basis
      integer natom
      integer stack,heap,global
      integer count
      integer n,i,j
      integer qm_choice
      double precision charge
      double precision coordinates (nda3)
      double precision QQCHEM(nda3)
      character*16 labels(nda)
      character*256 printName
      character*256 line_3
      character*256 chtmp(4)
      double precision qm_grad(nda3),qm_hessian(nda3,nda3)
      double precision qm_energy
      INTEGER MYID,IDOP, nnodes
      integer fnum, ferr
C
  802 FORMAT(/1X,A30/)
  809 FORMAT(3X,8F9.2)
  810 FORMAT(6F11.5)
  811 FORMAT('   NT=',I5,'  NS=',I10,'  NIP=',I10,'  NCROT=',I6)
  812 FORMAT('   ISEED=',I12,'  TIME=',F9.5)
  813 FORMAT('   ACTIVATE WITH ORTHANT SAMPLING'/)
  814 FORMAT(' NUMBER OF ATOMS=',I5/
     * ' NUMBER OF HARMONIC STRETCHES=',I5/
     * ' NUMBER OF MORSE STRETCHES=',I5/
     * ' NUMBER OF HARMONIC BENDS=',I5/
     * ' NUMBER OF ALPHA BENDS=',I5/
     * ' NUMBER OF LENNARD-JONES INTERACTIONS=',I5/
     * ' NUMBER OF TORSIONS=',I5/
     * ' NUMBER OF REPULSIONS=',I5/
     * ' NUMBER OF GHOST PAIRS=',I5/
     * ' NUMBER OF TETRAHEDRAL CENTERS=',I5/
     * ' NUMBER OF R-R COUPLINGS=',I5/
     * ' NUMBER OF R-THETA COUPLINGS=',I5/
     * ' NUMBER OF THETA-THETA COUPLINGS=',I5/
     * ' NUMBER OF DIHEDRAL ANGLES=',I5/
     * ' NUMBER OF AXILROD-TELLER THREE-BODY INTERACTIONS=',I5/
     * ' CHOICE OF SN2 POTENTIALS=',I5/
     * ' NUMBER OF RYDBERG POTENTIAL FUNCTIONS=',I5/
     * ' NUMBER OF HARTREE-FOCK DISPERSION FUNCTIONS=',I5/
     * ' NUMBER OF GENERALIZED LEPS(A) FUNCTIONS=',I5/
     * ' NUMBER OF GENERALIZED LEPS(B) FUNCTIONS=',I5/
     * ' CHOICE OF DMBE POTENTIAL ENERGY SURFACE=',I5/
     * ' NUMBER OF RELAXATION TERMS FOR H+DIAMOND INTERACTIONS=',I5/
     * ' NUMBER OF H---H NON-BONDED INTERACTION TERMS=',I5/
     * ' NUMBER OF ATOMS TREATED BY M.O. CALCULATIONS=',I5/
     * ' CHOICE OF CRCO6 POTENTIAL=', I5/
     * ' CHOICE OF QM_CODE', I5//)
  816 FORMAT('   NFQP=',I2,'  NCOOR= ',I2/)
  817 FORMAT('   NFR=',I2,6X,'NUMR=',I2)
  818 FORMAT(6X,'J-ATOM=',I4,4X,'K-ATOM=',I4)
  819 FORMAT('   NFB=',I2,6X,'NUMB=',I2)
  821 FORMAT(6X,20I4)
  822 FORMAT('   NFA=',I2,6X,'NUMA=',I2)
  823 FORMAT(15X,'INDICES OF ALPHA ANGLES TO BE PRINTED')
  824 FORMAT('   NFTAU=',I2,4X,'NUMTAU=',I2)
  825 FORMAT(15X,'INDICES OF TAU ANGLES TO BE PRINTED')
  826 FORMAT('   ACTIVATE WITH A BOLTZMANN DISTRIBUTION OF VIBRATIONAL',
     *' ENERGIES')
  827 FORMAT(28X,'NSELT = ',I2/)
  828 FORMAT('   MASSES OF ATOMS:',I5,' ATOMS'/)
  829 FORMAT(6X,'J-ATOM=',I4,4X,'M-ATOM=',I4,4X,'I-ATOM=',I4)
  835 FORMAT('   ACTIVATE WITH MICROCANONICAL NORMAL MODE SAMPLING'/)
  836 FORMAT(7X,F10.6,4X,F10.6,4X,F10.6)
  837 FORMAT('   VENUS CALCULATES AND SETS VZERO AND DELH VALUES SO',
     *'   THAT THE POTENTIAL ENERGY OF THE REACTANT(S) IS ZERO')
  838 FORMAT('   VZERO=',F14.3,' KCAL/MOLE')
  839 FORMAT('   PARAMETERS FOR PATH ',I1,':'/5X,'RMAX=',F6.2,'  RBAR='
     *,F6.2,'  NATOMA=',I3,'  NATOMB=',I3,'  DELH=',F6.2)
  840 FORMAT(5X,'INDICES FOR ATOMS OF FRAGMENT A:')
  841 FORMAT(5X,'INDICES FOR ATOMS OF FRAGMENT B:')
  842 FORMAT(5X,'DISTANCE BETWEEN THESE ATOMS DEFINES R.C. :',2I4,2X,'AN
     *D',2X,2I4)
  843 FORMAT('   NUMBER OF ADDITIONAL REACTION PATHS=',I2/)
  844 FORMAT('   PARAMETERS FOR REACTANT A')
  845 FORMAT(5X,'NORMAL MODE QUANTUM NUMBERS:')
  846 FORMAT('     REACTANT IS A DIATOMIC (TREATED SEMICLASSICALLY)')
  847 FORMAT(11X,'N:',I3,8X,'J:',I3)
  848 FORMAT(/'   INTERNUCLEAR PARAMETERS'/5X,'J-ATOM=',I3,4X,
     *'K-ATOM=',I3,4X,'RMAX=',F8.2,4X,'RBAR=',F8.2,4X,'DELH=',F8.3)
  851 FORMAT('   ACTIVATE WITH MICROCANONICAL QUASICLASSICAL NORMAL ',
     *     'MODE SAMPLING '/'      EBAR=',1PE13.5)
  852 FORMAT('   REACTANT B CANNOT BE INITIALIZED AT A BARRIER')
  853 FORMAT(/,5X,'L-ATOM=',I3,4X,'K-ATOM=',I3)
  890 FORMAT(/'   AN ERROR OCCURRED : THE CHECKPOINT FILE COULD',
     *' NOT BE READ')
  891 FORMAT(/'   THE OPTION NCHKP=-1 IS NOT AVAILABLE FOR ',
     *'TRAJECTORY CALCULATIONS'/)
  894 FORMAT(/5X,'COORDINATES AND MOMENTA ARE READ IN FROM ',
     *'CHECKPOINT FILE'/5X,'READING FROM UNIT 50'/)
  895 FORMAT(/'   CONTINUING RANDOM NUMBER SEQUENCE - READING FROM',
     *' UNIT 50'/)
  896 FORMAT(5X,'SAVE RANDOM NUMBER SEQUENCE IN UNIT 50',
     */,5X,12HNEXT SEED IS,1X,8I4,/)
  897 FORMAT(5X,47HCALCULATIONS ARE RESTARTED FROM CHECKPOINT FILE,
     */5X,'READING FROM UNIT 50'//)
  900 FORMAT(5X,'MOMENTS OF INERTIA IX, IY AND IZ(AMU-A**2):',3F9.3)
  901 FORMAT('   PARAMETERS FOR REACTANT B')
  902 FORMAT(5X,'RELATIVE ENERGY(KCAL):',F7.2,5X,
     *'INITIAL SEPARATION(A):',F6.2)
  903 FORMAT(5X,'NOB=',I2,5X,'BMAX(A)=',F5.1)
  904 FORMAT(5X,'NROT=',I2,5X,'TROT=',F9.2)
  906 FORMAT(5X,'EQUILIBRIUM COORDINATES FOR A:')
  907 FORMAT(5X,'EQUILIBRIUM COORDINATES FOR B:')
  909 FORMAT(6X,25I4)
  910 FORMAT(10X,'REACTION OCCURRED FOR PATH',I3)
  914 FORMAT('   ACTIVATE WITH NORMAL MODE SAMPLING'/)
  918 FORMAT(a)
  919 FORMAT(1X,a)
  920 FORMAT(5X,'HSCALE=',F8.3,'  PSCALE=',F5.2)
  921 FORMAT(5X,'HSCALE=',F8.3)
  925 FORMAT('   HINC=',F9.6,5X,'NPTS=',I2//)
  947 FORMAT('   NFTET=',I2,4X,'NUMTET=',I2)
  948 FORMAT(15X,'INDICES OF TETRAHEDRAL ANGLES TO BE PRINTED')
  953 FORMAT('   NFDH=',I2,5X,'NUMDH=',I2)
  954 FORMAT(15X,'INDICES OF DIHEDRAL ANGLES TO BE PRINTED')
  955 FORMAT('   ACTIVATE WITH LOCAL MODE EXCITATION'/)
  956 FORMAT('   PARAMETERS FOR LOCAL MODE EXCITATION'/6X,I3,
     *' THE MORSE STRETCH TO BE EXCITED INTO N =',I2,' LEVEL'/6X,
     *'LOCAL MODE ENERGY =',F8.3,6X,' EDELTA(BOXING)=',F6.3/)
  957 FORMAT(10X,'BOXING (N+1) AND (N-1) LEVELS'/
     *10X,'E(N+1)=',F10.3,'  EDEL=',F8.3/
     *10X,'E(N-1)=',F10.3,'  EDEL=',F8.3)
  958 FORMAT(5X,'MPLOT=',I2,'  NPLOT=',I3,'  NLM=',I3,'  MNTR=',I3)
  961 FORMAT('   THE INDICES FOR THE MORSE OSCILLATORS ',
     *'TO BE MONITORED ',4I3)
  962 FORMAT(5X,'VIBRATIONAL TEMPERATURE=',F7.1)
  963 FORMAT(/5X,'DIATOM CANNOT HAVE BARRIER EXCITATION!'/)
  964 FORMAT(13X,'START TRAJECTORIES AT A POTENTIAL BARRIER'//)
  965 FORMAT('   REACTION COORDINATE HAS FIXED ENERGY  ',F7.3,
     *'  KCAL/MOL'//)
  966 FORMAT('   REACTION COORDINATE HAS FIXED TEMPERATURE ',F7.1,
     *'  K'//)
  967 FORMAT(/8X,'ORTHANT SAMPLING CANNOT HAVE BARRIER EXCITATION!!!'/)
  968 FORMAT(14X,'INITIAL CONDITIONS ARE CHOSEN RANDOMLY'//)
  969 FORMAT(23X,'MINIMUM ENERGY SEARCH'//)
  970 FORMAT(19X,'INITIAL CONDITIONS ARE READ IN'//)
  971 FORMAT(23X,'NORMAL MODE ANALYSIS'//)
  972 FORMAT(22X,'REACTION PATH FOLLOWING'//)
  973 FORMAT(22X,'TRAJECTORY CALCULATIONS')
  974 FORMAT(/'   ACTIVATE WITH A BOLTZMAN DISTRIBUTION OF ',
     *       'TRANSLATIONAL ENERGIES'/5X,
     *       'TRANSLATIONAL TEMPERATURE = ',F10.1,' KELVIN',
     *       /5X,'INITIAL SEPARATION (A) = : ', F6.2/)
  982 FORMAT(/'   SURFACE PARAMETERS:'/5X,
     *       'ATOMS DEFINING THE REFERENCE PLANE IN THE SURFACE:',3I5)
  984 FORMAT(5X,'AIMING POINT ATOM:',I4/
     *       5X,'THTA=',F7.3,'  PHI=',F7.3,
     *       '  NCHI=',I2,'  CHI=',F7.3,' (DEG)')
  986 FORMAT(5X,'DISTANCES DEFINING THE ORIGIN:',
     *       ' RX0=',F8.3,' RY0=',F8.3,' RZ0=',F8.3,' (A)'/
     *       5X,'NTHET=',I2,'  THET=',F7.3,'   NPHI1=',I2,'  PHI1=',
     *       F7.3,'   NPHI2=',I2,'  PHI2=',F7.3,' (DEG)')
  990 FORMAT('   NFHT=',I2,5X,'NUMHT=',I2)
  992 FORMAT(6X,'INDICES OF ATOMS WHOSE HEIGHT ABOVE THE SURFACE TO BE',
     *' PRINTED:'/6X,'(HEIGHT VALUES ARE SAVED IN fort.17)')
  994 FORMAT(6X,20I4)
 1000 FORMAT(//'VZERO SET TO ',6E12.4,' KCAL/MOL SO THE REACTANT(S)',
     *'POTENTIAL ENERGY EQUALS ZERO'//)
C      INTIALIZE PARALLEL ENVIRNMENT
      call VENUS_parallel_init (myid,nnodes)
C
C         INITIALIZE ARRAYS AND PARAMETERS.  SPECIFIC VALUES FOR 
C         PARAMETERS MAY BE READ IN.
C
      stack=0
      heap=0
      global=0
      mplot=0
      nlm=0
      DO 11 I=1,2000
         NEVIB(I)=0
         NEVIBU(I)=0
         NEVIBL(I)=0
         NEVIBN(I)=0
         ESAV(I)=0.0D0
         ESQ(I)=0.0D0
   11 CONTINUE
      NMA=0
      NATOMA(1)=0
      NATOMB(1)=0
      NT=1
      NPATHS=0
      NABJ(1)=1
      NABK(1)=1
      RMAX(1)=100.D0
      RBAR(1)=100.D0
c
c   mdflag is used to avoid the gaussian call while doing the md sampling.
c   mdflag=0: call DOGAUSS
c   mdflag=1: do not call DOGAUSS
c
      mdflag=0
c
c   nhessf is a flag to control diagonalization in FMTRX
c     nhessf=1: not to diagonalize^M
c     nhessf=0: normal diagonalization for frequencies
c
      nhessf=0
c
C       CONSTANTS WITHIN THE COMPUTER PROGRAM
C
C           BOND ENERGIES FROM KCAL/MOLE:(C1)
C           HARMONIC STRETCH FORCE CONSTANT FROM MDYN/A:(C2)
C           HARMONIC BEND FORCE CONSTANT FROM MDYN-A/RAD**2:(C3)
C           EQUILIBRIUM ANGLES FROM DEGREES:(C4)
C           GAS LAW CONSTANT IN INTEGRATION UNITS:(C5)
C           FREQUENCIES(X TWOPI) FROM CM-1:(C6)
C           PLANCK'S CONSTANT DIVIDED BY TWOPI IN INTEGRATION UNITS:(C7)
C
      C1=0.04184D0
      C2=6.022045D0
      C3=6.022045D0
      C4=0.01745329D0
      C5=0.083144D-3
      C6=1.8836518D-3
      C7=0.063508D0
C
      RAA1=1.0D0/2.0D0
      RA1=1.0D0-SQRT(2.0D0)/2.0D0
      RA2=2.0D0*RA1
      RA3=2.0D0-3.0D0*SQRT(2.0D0)/2.0D0
      RB1=2.0D0-RA1
      RB2=2.0D0*RB1
      RB3=4.0D0-RA3
      RC1=1.0D0/6.0D0
      RC2=1.0D0/3.0D0
C
      PI=4.0D0*DATAN(1.0D0)
      HALFPI=PI/2.0D0
      TWOPI=2.0D0*PI
      GAO=0.0D0
      NSAD=0
      NCBA=0
      NCAB=0
      IBAR=0
      NPTS=2
      HINC=0.001D0

C YA: Read input file only on master node
      IF (MYID.EQ.0) THEN

C
C         READ AND WRITE TWO TITLE CARDS
C
      READ(5,918)TITLE1
      READ(5,918)TITLE2
      WRITE(6,919)TITLE1
      WRITE(6,919)TITLE2
C
C         THE FOLLOWING LINES ARE FOR WRITING THE DATE & TIME
C         IN THE OUTPUT, RIGHT AFTER THE TITLE. THIS, OF COURSE,
C         IS MACHINE AND OPERATING SYSTEM DEPENDENT.
C
C         SUBROUTINE vFDATE IS CURRENTLY FOR UNIX-SUPPORTED SYSTEMS.
C         FOR UNICOS-SUPPORTED SYSTEMS (CRAY COMPUTER SYSTEMS WITH
C            CFT77 COMPILER) COMMENT OUT THE CURRENT LINES AND UNCOMMENT
C            THE COMMENTED TWO LINES OF vFDATE.f.
C         
       CALL vFDATE(ADATE)
C
      WRITE(6,802)ADATE
C
C         READ # OF ATOMS AND # OF POTENTIAL TYPES
C
      READ(5,*)NATOMS,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,NVRR,NVRT,
     *NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,NRAX,NONB,NMO,
     *NCRCO6,qm_choice
      IF (NDMBE.NE.0) NDMBE=4
      WRITE(6,814)NATOMS,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,NVRR,
     *NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,NRAX,
     *NONB,NMO,NCRCO6,qm_choice
      I3N=3*NATOMS
C
C         READ POTENTIAL PARAMETERS
C
      CALL READPT

C YA: Read initialization information for NWChem from input file
      IF (NMO.NE.0) THEN
         if(nmo.gt.0) then
           natom=nmo
         else
           natom=natoms
         endif
C
C Read the QM paramenters from the Input file.
C
C     If line 3 starts with "read_qm_file" then all qm parameters
C     are read from qm file. The theory must be specified as set task: instead
C    of task. Print none is highly recommended.
c
            Read(5,*) qm_ini_file, basis, theory
            Read(5,*) (IAN(I),I=1,natom),charge
            write(6,*) 'The ab intio parameters including theory and
     * basis set are taken from', qm_ini_file
312     format('Basis Set : ',a16, ' Theory : ',a16,' Charge : ',f3.2)
313     format('The qm input is taken from ',a70)
        write (6,312) basis, theory, charge
         if (qm_ini_file.ne.'') then
           write (6,313) qm_ini_file
         end if
      END IF
C YA: End of reading input file only on master node
      END IF

C YA: Broadcast NMO to all the children nodes.
      call broadcast_data(NMO)

C It is not general. There are three choices: if VENUS is working with QM code or not.
C qm_choice defines choice of QM code (1 for NWChem). One needs to define variable that
C will check if code compiled with MPI or not.
      IF (NMO.NE.0) THEN

c All the nodes go to this function but only the master node continues
c and the rest of the nodes go into the sleep mode.

      idop=0
      if (myid.ne.0) call VENUS_parallel_loop (myid,idop,theory, 
     $  basis, natom,
     *  QQCHEM,labels,printName,qm_hessian,qm_grad,qm_energy,charge, 
     *  stack,heap,global,qm_ini_file)

       DO i=1,natom
            WRITE(6,'(4X,I3,8X,I4)')I,IAN(I)
       ENDDO

       do i=1,natom
       labels(i)= nameatm(ian(i))
       write(6,*)'The labels from the input file are ',labels(i)
       enddo

C To make sure that everything is intialized and avoid repetative
c Intialization
      initialize = .true.

C Print name can have the foollowing values
c High: Print a lot of Nwchem output
c Moderate : prints moderate amount of Nwchem Information
C Low : Prints only few parts of  Nwchem output
C None : don't print any Nwchem part
C Debug : Print every thing from Nwchem to depug
      printName = "none"

c Intializing the qm gradient parameter
      do i = 1,natom*3
       qm_grad(i) = 0.0d0
      enddo
      qm_energy = 0.0d0

C ENDIF for NMO
      endif


C         IF VZERO=999.0, VENUS DETERMINES VZERO SO THAT THE CLASSICAL
C         POTENTIAL ENERGY MINIMUM OF THE REACTANT(S) IS ZERO.  
C         THE POTENTIAL ENERGY MINIMA FOR THE PRODUCTS OF THE VARIOUS
C         REACTION PATHS ARE ALSO CALCULATED BY VENUS.
C
C 
      NVZERO=0
      READ(5,*)VZERO
      IF ((VZERO.GT.998.99d0).AND.(VZERO.LT.999.02d0)) THEN
         NVZERO=1
         VZERO=0.0D0
         WRITE(6,837)
      ELSE
         WRITE(6,838)VZERO
         VZERO=VZERO*C1
      ENDIF
C
C         INPUT FOR TRAJECTORIES
C
C            NSELT=-3 STATIONARY POINT SEARCH
C                     (OPTION NOT AVAILABLE IN the release VERSION yet)
C            NSELT=-2 PROGRAM FINDS REACTION PATH and CVTST rate constants
c              nmep=0 Reaction Path only (old venus way, MPATHO.f)
c              nmep=1 CVTST Rate Constants too.
c                      Use frequencies and rotational constants
c              nmep=2 Calculate frequencies and rotational constants
c                      Needs coordinates to do this
c              nmep>0 Needs extra input as specified in MPATH.f
C
C            NSELT=-1 PROGRAM DOES NORMAL MODE ANALYSIS
C            NSELT=0  INITIAL Q'S AND P'S ARE READ IN
C            NSELT=1  PROGRAM FINDS MINIMUM ENERGY GEOMETRY
C                     (INITIAL Q'S AND P'S ARE READ IN)
C            NSELT=2  CHOOSE INITIAL CONDITIONS FOR ONE OR TWO REACTANTS
C                     NACT ARE THE OPTIONS FOR CHOOSING INITIAL CONDITION
C            NSELT=3  CHOOSE INITIAL CONDITIONS FROM POTENTIAL BARRIER
C               NBAR=1  MICROCANONICAL SAMPLING
C                       WITH FIXED REACTION COORDINATE ENERGY
C               NBAR=2  THERMAL SAMPLING
C               NBAR=3  MICROCANONICAL QUASICLASSICAL SAMPLING
C                       INCLUDING REACTION COORDINATE ENERGY
C                   ENMTA is the total available energy
C                         + zero point energy at TS.
C                   IN SELECT, EBAR IS DETERMINED AND THEN NBAR IS SET TO 1
C                      BEFORE CALLING BAREXC
C
C            NACT=0   REACTANTS ARE ATOMS AND/OR DIATOMS.
C            NACT=1   ACTIVATE WITH ORTHANT SAMPLING
C            NACT=2   ACTIVATE WITH MICROCANONICAL NORMAL MODE SAMPLING
C            NACT=3   ACTIVATE WITH FIXED NORMAL MODE ENERGIES
C            NACT=4   LOCAL MODE EXCITATION
C            NACT=5   CHOOSE NORMAL MODE ENERGIES FROM A BOLTZMANN
C                     DISTRIBUTION
C            NACT=6   SAME AS NBAR=3 (CHANGES NACT TO 6 WHEN NBAR=3)
C            NACT=7   MOLECULAR DYNAMICS SAMPLING BY RESCALING VELOCITIES
c   rescale the velocities of the system according to given temp.
c   reference:  "MOLECULAR DYNAMICS SIMULATION" by JIM HAILE  p.458
c   Only fragment B can have this when it is surface.

C
      READ(5,*)NSELT,NSURF,NCHKP
  700 FORMAT('     NSURF=',I2/)
      WRITE(6,700)NSURF
      if(nselt.eq.-2)then
        read(5,*)nmep
        if(nmep.eq.0)then
          WRITE(6,'(a)')' Just calculate MEP'
        else
          WRITE(6,'(a)')' MEP with CVTST Rate constants'
        endif
      endif
      READ(5,*)(W(I),I=1,NATOMS)
      WRITE(6,*)
      WRITE(6,828)NATOMS
      WRITE(6,810)(W(I),I=1,NATOMS)
      WRITE(6,*)
      WRITE(6,827)NSELT
      IF (NSURF.EQ.1) THEN
         WRITE(6,'(14X,35HGAS/SURFACE COLLISION --- BEAM TYPE,/)')
      ELSEIF (NSURF.EQ.2) THEN
         WRITE(6,'(14X,38HGAS/SURFACE COLLISION --- THERMAL TYPE,/)')
      ENDIF
      IF (NSELT.LT.0) GOTO 87
      WRITE(6,973)
c
c  This is for integration method
c  INTEGRATOR=0: Runge-Kutta (set lll, niter=0 for this)
c             1: Gauss-Radau
c                 lll: converging digits
c                     6, 8: variable time step
c                     -6,-8: fixed time step
c                 niter=2 to check energy convergence, otherwise 6
c                    (2 takes only half time of 6)
c             2: Simplectic integrator
c         integrator              8th Line   9th Line  10th Line(timestep)
c   Case 1: Adams-Moulton:         0 0 0     1,600,10,0   0.01
c   Case 2: 6th-order Symplectic:  2 6 0     1,150, 2,0   0.04 
c   Case 3: 8th-order Symplectic:  2 8 0     1, 70, 1,0   0.085
c   Case 4: Gauss-Radau:           1 -8 2    1, 40, 1,0   0.15
c                         (fixed timestep, 29 Force each timestep)
c   Case 5: Gauss-Radau:           1 -8 6       1, 20, 1,0   0.30
c                         (fixed timestep, 57 Force each timestep)
c
c   This integrator choice is not necessary for nmo=-2 because of the
c     use of integrator in gaussian.
c
 333   format(5x,'INTEGRATOR=',i4,'   lll=',i4,'   niter=',i4)
      if(nselt.ne.1.and.nmo.ne.-2)then
        read(5,*)INTEGRATOR,lll,niter
      write(6,333)INTEGRATOR,lll,niter
        if(integrator.eq.1)then
          write(6,'(/a)')' Gauss-Radau integration'
          if(lll.lt.0)then
            write(6,'(a)')'   Fixed Time Step'
          else
            write(6,'(a)')'   Variable Time Step'
          endif
          write(6,'(a,i3)')'     Converging Digits:',abs(lll)
          if(niter.eq.2)then
            write(6,'(a)')'     Checking Energy Convergence'
          else
            write(6,'(a)')'     Not Checking Energy Convergence'
            niter=6
          endif
        elseif(integrator.eq.2)then
          if(lll.eq.4.or.lll.eq.6.or.lll.eq.8)then
            write(6,'(/i3,a)')lll,'th Order Symplectic Integration'
          else
            write(6,'(a)')' *** Wrong order for Symplectic Integration'
            write(6,'(a)')'     Must be 4, 6, or 8 for lll ***'
            stop
          endif
        elseif(integrator.eq.3)then
          if(lll.eq.1)then
            write(6,*)' Velocity Verlet Integration'
          else
            write(6,*)' Velocity-Corrected Verlet Integration'
          endif
        elseif(integrator.eq.4)then
                write(6,*)' Hessian based integrator'
                if(lll.eq.1)then
                write(6,*)' Rotated Coordinate interpolation'
          elseif(lll.eq.2)then
           write(6,*)' Rotated Coordinate interpolation using
     * Six polynomials'
         else
           write(6,*)' General Cartesian coordinate interpolation'
          endif
          read(5,*)trad,rmin,nhesup
          write(6,*)
          write(6,*)'  Number of Hessian updates : ',nhesup-1
c   trad = trust radius fluctuates around this value
c   rmin = tolerence for error in pred. and corr. trajectories (0.002)
c   nhesup = hessian updates: hessian updated every nhesup-th step.
        else
          write(6,'(/a)')' Adams-Moulton Integration (venus96 way)'
          integrator=0
        endif
      endif
      WRITE(6,*)
      READ(5,*)NT,NS,NIP,NCROT
      nstep = ns
      READ(5,*)TIME
c time step for hessian-based integrator
      dt = time
c
      SYBTI=10.0D0*TIME*DBLE(NIP)
      IF (NSELT.EQ.3) THEN
         WRITE(6,964)
         READ(5,*)NBAR,DUM
         IF (NBAR.EQ.1) THEN
            EBAR=DUM
            WRITE(6,965)EBAR
         ELSEIF (NBAR.EQ.2) THEN
            TBAR=DUM
            WRITE(6,966)TBAR
         ELSEIF(NBAR.EQ.3)THEN
c
c  ENMTA is the total available energy + zero point energy at TS.
c
            ENMTA=DUM
            WRITE(6,851)ENMTA
         ENDIF
      ENDIF
      IF (NSELT.EQ.2.OR.NSELT.EQ.3)then
        READ(5,*)NACTA,NACTB,ISEED
c--------------------------------------------------------------------
c
c  nactb = 7 for scaling velocities
c  nsel = 1: write out temperature for each step
c  nscale : # of time steps to rescale the velocities
c  nequal : # of time steps to anneal the system. the basic idea is
c           the velocities of the system might not be representing a
c           a boltamann distribution immediately after velocity re-
c           scaling steps. therefore, additional steps needed to let
c           the system anneal by itself with constant energy while it
c           acquires a boltzmann distribution. refer to
c           reference:  "MOLECULAR DYNAMICS SIMULATION" by JIM HAILE
c  nrgd: number of anchor atoms which need not be moved
c        These atoms should be at the end of the list.
c
c--------------------------------------------------------------------
        if(nactb.eq.7)then
          read(5,*)thermotemp,nscale,nequal,nsel,nrgd
           write(6,*)' Doing MD Sampling'
           write(6,*)'nsel =',nsel,', nscale = ',nscale,
     *       ', nequal =', nequal, ', thermotemp =',thermotemp,
     *       ', nrgd =',nrgd
          nthermb=0
          nrscl=0
        else
          nsel=0
          nscale=0
          nequal=0
        endif
      else
        iseed=0
      endif
      IF(NSELT.EQ.3.AND.NBAR.EQ.3)NACTA=6
      IF (NSELT.EQ.2) WRITE(6,968)
      IF (NSELT.EQ.1) WRITE(6,969)
      IF (NSELT.EQ.0) WRITE(6,970)
      WRITE(6,811)NT,NS,NIP,NCROT
      WRITE(6,812)ISEED,TIME
      IF (NSELT.EQ.3.AND.NACTA.EQ.1) THEN
         WRITE(6,967)
         STOP
      ENDIF
      IF (NSELT.EQ.0) GOTO 140
      WRITE(6,*) ''
      IF (NSELT.EQ.1) GOTO 145
      IF (NSELT.GE.2.AND.NACTA.LE.1.AND.NACTB.LE.1) GOTO 88
C
C         READ DISPLACEMENT INTERVAL, HINC.
C         READ NUMBER OF DISPLACEMENTS ABOUT CARTESIAN MINIMUM, NPTS
C
   87 READ(5,*)HINC,NPTS
      IF (NSELT.GT.0) WRITE(6,925)HINC,NPTS
      IF (NSELT.EQ.-1) THEN
         WRITE(6,971)
         WRITE(6,925)HINC,NPTS
         GOTO 145
      ELSEIF (NSELT.EQ.-2) THEN
         WRITE(6,972)
         WRITE(6,925)HINC,NPTS
C
C         STATIONARY POINT SEARCH
C         This version uses the routine adopted from Numerical Methods
C         book. If there is any legal problem, comment out this part and
C         do not include the STATPT.f file.
C
      ELSEIF (NSELT.EQ.-3) THEN
         READ(5,*)NT,NC,NX,TIME
         CALL STATPT
         WRITE(6,925)HINC,NPTS
         GOTO 145
      ENDIF
C
C         READ PARAMETERS FOR REACTANTS
C
C
C             PARAMETERS FOR REACTANT A
C
   88 READ(5,*)NATOMA(1),NLINA
      IF (NSELT.EQ.3)then
        if(NATOMA(1).EQ.2) THEN
          WRITE(6,963)
          STOP
        elseif(NATOMA(1).LT.2) THEN
          WRITE(6,*)"ERROR INPUT FOR BARRIER EXCITATION!!!"
          STOP
        endif
      ENDIF
C
C             NLINA=0, MOLECULE IS NONLINEAR
C             NLINA=1, MOLECULE IS LINEAR
C
      K=NATOMA(1)
      IF (NSELT.NE.-2) WRITE(6,844)
      IF (NACTA.EQ.1) WRITE(6,813)
      IF (NACTA.EQ.2) WRITE(6,835)
      IF (NACTA.EQ.3) WRITE(6,914)
      IF (NACTA.EQ.4) WRITE(6,955)
      IF (NACTA.EQ.5) WRITE(6,826)
      if (nacta.eq.6) write(6,*)' Fixed Energy Sampling including RC'
      if (nacta.eq.7) then
         write(6,*)' MD sampling for reactant A is not supported'
         stop
      endif
c
      WTA(1)=0.0D0
      DO J=1,K
         LA(1,J)=J
         LL(J)=LA(1,J)
         WTA(1)=WTA(1)+W(J)
      ENDDO
      IF (NSELT.EQ.-2) GOTO 107
      K=3*K
      READ(5,*)(QZA(1,J),J=1,K)
      WRITE(6,906)
      WRITE(6,836)(QZA(1,J),J=1,K)
c      if(natoma(1).eq.1)goto 106
c  changed to skip mode energies and mode populations part consistent with
c  the manual.   - Raj, Dec. 2007
      if(natoma(1).eq.1)goto 107
C
C             TRANSFORM QZA TO THE CENTER OF MASS / PRINCIPAL
C             AXES OF INERTIA FRAME
C
      DUM=0.0D0
      DO J=1,K/3
         DUM1=QZA(1,3*J)-QZA(1,3)
         DUM=DUM+DUM1*DUM1
      ENDDO
      IF (DUM.LE.1.0D-10) THEN
         DO J=1,K/3
            QZDUM=QZA(1,3*J)
            QZA(1,3*J)=QZA(1,3*J-1)
            QZA(1,3*J-1)=QZA(1,3*J-2)
            QZA(1,3*J-2)=QZDUM
         ENDDO
      ENDIF
      DO J=1,K
         Q(J)=QZA(1,J)
      ENDDO
      WT=WTA(1)
      N=NATOMA(1)
      CALL CENMAS(WT,QCM,VCM,N)
      CALL ROTN(AM,EROT,N)
      CA(1,1)=AIXX
      CA(2,1)=-AIXY
      CA(2,2)=AIYY
      CA(3,1)=-AIXZ
      CA(3,2)=-AIYZ
      CA(3,3)=AIZZ
      CALL EIGN(CA,CB,3,RHO)
      AI(1)=EIG(1)
      AI(2)=EIG(2)
      AI(3)=EIG(3)
C       
C        *** ADDED BY GILLES PESLHERBE TO CORRECT VENUS96 ***
C
C        TRANSPOSE MATRIX CB
C
      DO II=1,3
         DO JJ=1,II
            DUM=CB(II,JJ)
            CB(II,JJ)=CB(JJ,II)
            CB(JJ,II)=DUM
         END DO
      END DO
C
C        SET PHASE OF THE EIGENVECTOR
C
      DO I=1,3
         CBMAX=0.0D0
         DO J=1,3
            IF (ABS(CB(J,I)).GE.ABS(CBMAX)) CBMAX=CB(J,I)
         ENDDO
         IF (CBMAX.LT.0.0D0) THEN
            DO J=1,3
               CB(J,I)=-CB(J,I)
            ENDDO
         ENDIF
      ENDDO
C
C        CHECK FOR RIGHT-HANDED CARTESIAN SYSTEM
C
      DETCB1=CB(1,1)*(CB(2,2)*CB(3,3)-CB(2,3)*CB(3,2))
      DETCB2=CB(1,2)*(CB(2,3)*CB(3,1)-CB(2,1)*CB(3,3))
      DETCB3=CB(1,3)*(CB(2,1)*CB(3,2)-CB(2,2)*CB(3,1))
      DETCB=DETCB1+DETCB2+DETCB3
      IF (DETCB.LT.0.0) THEN
         DO J=1,3
            CB(J,1)=-CB(J,1)
         ENDDO
      END IF
C
C        *** END OF PESLHERBE CORRECTION***
C 
      DO J=1,K/3
         DO JL=1,3
            DUM=0.0D0
            DO L=1,3
               DUM=DUM+CB(L,JL)*QQ(3*J+L-3)
            ENDDO
            QZA(1,3*J+JL-3)=DUM
         ENDDO
      ENDDO
C
C             PARAMETERS FOR DIATOM
C             DIATOM IS TREATED SEMICLASSICALLY BY EBK QUANTIZATION
C
      IF (NATOMA(1).EQ.2) THEN
         READ(5,*)NNA,JA
         WRITE(6,846)
         WRITE(6,847)NNA,JA
         goto 107
C
C             REACTANT A IS A POLYATOMIC
C
      else
        IF (NACTA.EQ.1.OR.NACTA.EQ.2) THEN
          READ(5,*)ENMTA,PSCALA
          IF (NACTA.EQ.1) WRITE(6,920)ENMTA,PSCALA
          IF (NACTA.EQ.2) WRITE(6,921)ENMTA
          GOTO 105
        ELSEIF (NACTA.EQ.3 .OR. NACTA.EQ.4) THEN
          J=K-6+NLINA
          IF (NSELT.NE.3) THEN
            READ(5,*)(ANQA(I),I=1,J)
            WRITE(6,845)
            WRITE(6,809)(ANQA(I),I=1,J)
          ELSE
            READ(5,*)(ANQA(I),I=1,J-1)
            WRITE(6,845)
            WRITE(6,809)(ANQA(I),I=1,J-1)
          ENDIF
        ELSEIF(NACTA.EQ.5) THEN
          READ(5,*)TVIBA
          WRITE(6,962)TVIBA
        ENDIF
      endif
C
C             SET NMA AND INITIALIZE N.M. ARRAYS
C
      NMA=K-6+NLINA
C
      ALLOCATE(ENSAV(NMA,2000))
      ALLOCATE(ENSQ(NMA,2000))
      DO I=1,NMA
         DO J=1,2000
            ENSAV(I,J)=0.0D0
            ENSQ(I,J)=0.0D0
         ENDDO
      ENDDO
C
C             NROTA=0, SAMPLE ROTATIONAL ENERGY FROM A THERMAL
C                      DISTRIBUTION
C             NROTA=1, ROTATIONAL ENERGY OF EACH AXIS IS RT/2
C
  105 CONTINUE
      READ(5,*)NROTA,TROTA
      WRITE(6,904)NROTA,TROTA
      WRITE(6,900)(AI(I),I=1,3)
C
C             PARAMETERS FOR LOCAL MODE EXCITATION
C             NEXM    : THE INDEX OF THE MORSE OSCILLATOR TO BE EXCITED
C             NLEV   : THE QUANTUM NUMBER OF THE INITIAL EXCITATION
C             ENON    : THE INITIAL BOND ENERGY IN KCAL/MOLE
C             EDELTA  : BOXING LENGTH OF ENERGY
C
      IF (NACTA.EQ.4) THEN
         WRITE(6,*)
         READ(5,*)NEXM,NLEV
         CALL LMODE(0,ENU,EDELTU,ENL,EDELTL)
         WRITE(6,956)NEXM,NLEV,ENON,EDELTA
         WRITE(6,957)ENU,EDELTU,ENL,EDELTL
      ENDIF
C
C             PARAMETERS FOR PLOTTING MODE ENERGIES
C
C             NPLOT   : THE INCREMENT IN CYCLE COUNT FOR LISTING MODE
C                       ENERGIES
C             MPLOT   : FOR THE TIME AVERAGED PLOTTING
C                       MPLOT=0 : DO NOT PRINT ON FILE 12 THE AVERAGE
C                                 MODE ENERGIES AT CYCLES
C                       MPLOT=1 : PRINT ON FILE 12 THE AVERAGE MODE
C                                 ENERGIES AT CYCLES INCREMENTED BY
C                                 NPLOT
C             NLM     : NUMBER OF BOND MODES MONITORED IN ADDITION TO
C                       THE EXCITED ONE
C             MNTR    : INDEX FOR NORMAL MODE WHOSE POPULATION IS TO BE
C                       MONITORED
C
  106 READ(5,*)MPLOT,NPLOT,NLM,MNTR
      WRITE(6,958)MPLOT,NPLOT,NLM,MNTR
C
C             INITIALIZE THE LOCAL MODE ARRAYS EXCEPT THE EXCITED ONE
C
      IF (MPLOT.EQ.1.AND.NLM.NE.0) THEN
         IF (NLM.NE.0) THEN
            READ(5,*)(MNLM(I),I=1,NLM)
            WRITE(6,961)(MNLM(I),I=1,NLM)
            DO I=1,NLM
               DO J=1,2000
                  EBSAV(I,J)=0.0D0
                  EBSQ(I,J)=0.0D0
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      call flush(6)
C
C             PARAMETERS FOR REACTANT B
C
  107 READ(5,*)NATOMB(1),NLINB
      if(nselt.eq.-2)nacta=nmep
      IF(NSELT.EQ.3.AND.NATOMB(1).NE.0)THEN
         WRITE(6,*)"ERROR INPUT FOR BARRIER EXCITATION!!!"
         STOP
      ENDIF
      IF (NATOMB(1).EQ.0)THEN
         IF(NSELT.EQ.-2) GOTO 145
         GOTO 140
      ENDIF
      IF (NSELT.NE.-2) THEN
         WRITE(6,*)
         WRITE(6,901)
      ENDIF
      IF (NACTB.EQ.1) WRITE(6,813)
      IF (NACTB.EQ.2) WRITE(6,835)
      IF (NACTB.EQ.3) WRITE(6,914)
      IF (NACTB.EQ.4) WRITE(6,955)
      IF (NACTB.EQ.5) WRITE(6,826)
      if (nactb.eq.6.and.nselt.ne.3) then
         write(6,*)' Microcanonical Quasiclassical sampling for B?'
      endif
      if (nactb.eq.7)then
        if(nsurf.ne.0)then
          write(6,*)' MD Sampling for Surface'
        else
          write(6,*)' MD sampling for molecule is not supported'
          stop
        endif
      endif
c
      K=NATOMB(1)
      WTB(1)=0.0D0
      DO J=1,K
         M=J+NATOMA(1)
         LB(1,J)=M
         LL(J)=LB(1,J)
         WTB(1)=WTB(1)+W(M)
      ENDDO
      IF (NSELT.EQ.-2) GOTO 145
      K=3*K
      READ(5,*)(QZB(1,J),J=1,K)
      WRITE(6,907)
      WRITE(6,836)(QZB(1,J),J=1,K)

C         IF NVZERO=1, DETERMINE VZERO SO THE POTENTIAL ENERGY OF THE
C         REACTANTS IS ZERO.
C
      IF(NVZERO.EQ.1)THEN
         II=1
         CALL POTENZ(II)
         VZERO=-V
         WRITE(6,1000)VZERO
         VZERO=VZERO*C1
      ENDIF
C
C             TRANSFORM QZB TO THE CENTER OF MASS / PRINCIPAL
C             AXES OF INERTIA FRAME
C
      DUM=0.0D0
      DO J=1,K/3
         DUM1=QZB(1,3*J)-QZB(1,3)
         DUM=DUM+DUM1*DUM1
      ENDDO
      IF (DUM.LE.1.0D-10) THEN
         DO J=1,K/3
            QZDUM=QZB(1,3*J)
            QZB(1,3*J)=QZB(1,3*J-1)
            QZB(1,3*J-1)=QZB(1,3*J-2)
            QZB(1,3*J-2)=QZDUM
         ENDDO
      ENDIF
C
      I=3*NATOMA(1)
      DO J=1,K
         Q(J+I)=QZB(1,J)
      ENDDO
      WT=WTB(1)
      N=NATOMB(1)
      CALL CENMAS(WT,QCM,VCM,N)
      CALL ROTN(AM,EROT,N)
      CA(1,1)=AIXX
      CA(2,1)=-AIXY
      CA(2,2)=AIYY
      CA(3,1)=-AIXZ
      CA(3,2)=-AIYZ
      CA(3,3)=AIZZ
      CALL EIGN(CA,CB,3,RHO)
      BI(1)=EIG(1)
      BI(2)=EIG(2)
      BI(3)=EIG(3)
C       
C        *** ADDED BY GILLES PESLHERBE TO CORRECT VENUS96 ***
C
C
C        TRANSPOSE MATRIX CB
C
      DO II=1,3
         DO JJ=1,II
            DUM=CB(II,JJ)
            CB(II,JJ)=CB(JJ,II)
            CB(JJ,II)=DUM
         END DO
      END DO
C
C        SET PHASE OF THE EIGENVECTOR
C
      DO ICB=1,3
         CBMAX=0.0D0
         DO J=1,3
            IF (ABS(CB(J,ICB)).GE.ABS(CBMAX)) CBMAX=CB(J,ICB)
         ENDDO
         IF (CBMAX.LT.0.0D0) THEN
            DO J=1,3
               CB(J,ICB)=-CB(J,ICB)
            ENDDO
         ENDIF
      ENDDO
C
C        CHECK FOR RIGHT-HANDED CARTESIAN SYSTEM
C
      DETCB1=CB(1,1)*(CB(2,2)*CB(3,3)-CB(2,3)*CB(3,2))
      DETCB2=CB(1,2)*(CB(2,3)*CB(3,1)-CB(2,1)*CB(3,3))
      DETCB3=CB(1,3)*(CB(2,1)*CB(3,2)-CB(2,2)*CB(3,1))
      DETCB=DETCB1+DETCB2+DETCB3
      IF (DETCB.LT.0.0) THEN
         DO J=1,3
            CB(J,1)=-CB(J,1)
         ENDDO
      ENDIF
C
C        *** END OF PESLHERBE CORRECTION ***
C
      DO J=1,K/3
         DO JL=1,3
            DUM=0.0D0
            DO L=1,3
               DUM=DUM+CB(L,JL)*QQ(I+3*J+L-3)
            ENDDO
            QZB(1,3*J+JL-3)=DUM
         ENDDO
      ENDDO
C
C             PARAMETERS FOR DIATOM
C             DIATOM IS TREATED SEMICLASSICALLY BY EBK QUANTIZATION
C
      IF (NATOMB(1).LE.2) THEN
         IF (NATOMB(1).NE.1.and.nlinb.eq.1) THEN
            READ(5,*)NNB,JB
            WRITE(6,846)
            WRITE(6,847)NNB,JB
            WRITE(6,*)
         ENDIF
      ELSE
C
C             REACTANT B IS A POLYATOMIC
C
         IF (NACTB.EQ.1.OR.NACTB.EQ.2) THEN
            READ(5,*)ENMTB,PSCALB
            IF(NACTB.EQ.1)WRITE(6,920)ENMTB,PSCALB
            IF(NACTB.EQ.2)WRITE(6,921)ENMTB
         ELSEIF (NACTB.EQ.3.OR.NACTB.EQ.4) THEN
            J=K-6+NLINB
            READ(5,*)(ANQB(I),I=1,J)
            WRITE(6,845)
            WRITE(6,809)(ANQB(I),I=1,J)
         ELSEIF(NACTB.EQ.5)THEN
            READ(5,*)TVIBB
            WRITE(6,962)TVIBB
         ENDIF
C
C             NROTB=0, SAMPLE ROTATIONAL ENERGY FROM A THERMAL
C                      DISTRIBUTION
C             NROTB=1, ROTATIONAL ENERGY OF EACH AXIS IS RT/2
C
         if(nsurf.eq.0)then
           READ(5,*)NROTB,TROTB
           WRITE(6,904)NROTB,TROTB
           WRITE(6,900)(BI(I),I=1,3)
         else
           nrotb=1
           trotb=0.d0
         endif
      ENDIF

C
C        REACTANT A AND B INTERNUCLEAR PARAMETERS FOR THE FIRST
C        REACTION PATH.  TWO ATOM PAIRS MAY BE USED TO DEFINE THE 
C        REACTION PATH; NABJ-NABK AND NABL-NABM.  IF NABL=0 ONLY
C        THE FIRST ATOM PAIR IS USED.  IF NVZERO=1, VENUS CALCULATES
C        DELH.
C
      IF (NVZERO.EQ.1) THEN
         II=1
         CALL POTENZ(II)
         DELH(1)=V
         ENDIF
C
      READ(5,*)NABJ(1),NABK(1),NABL(1),NABM(1),RMAX(1),RBAR(1),DELH(1)
      WRITE(6,848)NABJ(1),NABK(1),RMAX(1),RBAR(1),DELH(1)
      WRITE(6,853)NABL(1),NABM(1)
C
C             NOB=0, RANDOMLY SAMPLE IMPACT PARAMETER BETWEEN 0 AND BMAX
C             NOB=1, IMPACT PARAMETER EQUALS BMAX
C
      READ(5,*)NREL,SEREL,S
      IF (NREL.EQ.1) THEN
         WRITE(6,902)SEREL,S
         SEREL=SEREL*C1
      ELSEIF (NREL.EQ.0) THEN
         TRANS=SEREL
         WRITE(6,974)TRANS,S
      ENDIF 
      READ(5,*)NOB,BMAX
      WRITE(6,903)NOB,BMAX
C
C         PARAMETERS FOR SURFACE
C
      IF (NSURF.EQ.1) THEN
          READ(5,*)NN1,NN2,NN3
          WRITE(6,982)NN1,NN2,NN3
          READ(5,*)NN4,THTA,PHI,NCHI,CHI
          WRITE(6,984)NN4,THTA,PHI,NCHI,CHI
          THTA=THTA*C4
          PHI=PHI*C4
          CHI=CHI*C4
       ELSEIF (NSURF.EQ.2) THEN
          READ(5,*)NN1,NN2,NN3
          WRITE(6,982)NN1,NN2,NN3
          READ(5,*)RX0,RY0,RZ0,NTHET,THET,NPHI1,PHI1,NPHI2,PHI2
          WRITE(6,986)RX0,RY0,RZ0,NTHET,THET,NPHI1,PHI1,NPHI2,PHI2
          THET=THET*C4
          PHI1=PHI1*C4
          PHI2=PHI2*C4
       ENDIF
C
C         PARAMETERS FOR CLASSIFYING REACTION EVENTS
C         NPATHS = NUMBER OF REACTION PATHS IN ADDITION TO THE
C         REACTION PATH 1, WHICH IS FOR REACTANTS A AND B.  TWO ATOM
C         PAIRS MAY BE USED TO DEFINE THE REACTION PATH; NABJ-NABK
C         AND NABL-NABM.  IF NABL=0 ONLY THE FIRST ATOM PAIR IS USED.
C
  140 CONTINUE
      WRITE(6,*)
      READ(5,*)NPATHS
c
c NPATHS < 0 added by ksong
c This is for SID only.
c Surface should not react
c Fragmentation of A will be monitored.
c Reaction will end when channel 1 occurred
c
      if(npaths.lt.0)then
        nsid=1
        npaths=-npaths
        write(6,*)' Surface does not react and fragmentation of A'
        write(6,*)' only will be considered'
      else
        nsid=0
      endif
      WRITE(6,843)NPATHS
      IF (NPATHS.NE.0) THEN
         M=NPATHS+1
         DO I=2,M
            READ(5,*)RMAX(I),RBAR(I),NATOMA(I),NATOMB(I),DELH(I)
            WRITE(6,839)I,RMAX(I),RBAR(I),NATOMA(I),NATOMB(I),DELH(I)
            READ(5,*)NABJ(I),NABK(I),NABL(I),NABM(I)
            WRITE(6,842)NABJ(I),NABK(I),NABL(I),NABM(I)
            K=NATOMA(I)
            READ(5,*)(LA(I,J),J=1,K)
            WTA(I)=0.0D0
            DO J=1,K
               WTA(I)=WTA(I)+W(LA(I,J))
            ENDDO
            WRITE(6,840)
            WRITE(6,909)(LA(I,J),J=1,K)
            K=3*NATOMA(I)
            READ(5,*)(QZA(I,J),J=1,K)
            WRITE(6,906)
            WRITE(6,836)(QZA(I,J),J=1,K)
            K=NATOMB(I)
            IF (K.NE.0) THEN
               READ(5,*)(LB(I,J),J=1,K)
               WTB(I)=0.0D0
               DO J=1,K
                  WTB(I)=WTB(I)+W(LB(I,J))
               ENDDO
               WRITE(6,841)
               WRITE(6,909)(LB(I,J),J=1,K)
               K=3*NATOMB(I)
               READ(5,*)(QZB(I,J),J=1,K)
               WRITE(6,907)
               WRITE(6,836)(QZB(I,J),J=1,K)
               WRITE(6,*)
            ENDIF
C
C           IF NVZERO=1, VENUS CALCULATES DELH.
C
            IF (NVZERO.EQ.1) THEN
            II=I
            CALL POTENZ(II)
            DELH(I) = V
            WRITE(6,*) DELH(I)
            ENDIF              
         ENDDO
      ENDIF
  145 CONTINUE
C
C         INFORMATION TO BE PRINTED
C
      READ(5,*)NFQP,NCOOR
      WRITE(6,816)NFQP,NCOOR
C
C         NFQP=0, DO NOT PRINT Q AND P ARRAYS
C         NCOOR=0, DO NOT WRITE COORDINATES INTO UNIT 8
C
      READ(5,*)NFR,NUMR
      WRITE(6,817)NFR,NUMR
      IF (NFR.NE.0) THEN
         READ(5,*)(JR(I),KR(I),I=1,NUMR)
         WRITE(6,*)
         DO I=1,NUMR
            WRITE(6,818)JR(I),KR(I)
         ENDDO
         WRITE(6,*)
      ENDIF
C
      READ(5,*)NFB,NUMB
      WRITE(6,819)NFB,NUMB
      IF (NFB.NE.0) THEN
         WRITE(6,*)
         READ(5,*)(KB(I),IB(I),MB(I),I=1,NUMB)
         WRITE(6,829)(KB(I),IB(I),MB(I),I=1,NUMB)
         WRITE(6,*)
      ENDIF
C
      READ(5,*)NFA,NUMA
      WRITE(6,822)NFA,NUMA
      IF (NFA.NE.0) THEN
         READ(5,*)(IA(I),I=1,NUMA)
         WRITE(6,823)
         WRITE(6,821)(IA(I),I=1,NUMA)
         WRITE(6,*)
      ENDIF
C
      READ(5,*)NFTAU,NUMTAU
      WRITE(6,824)NFTAU,NUMTAU
      IF (NFTAU.NE.0) THEN
         READ(5,*)(ITAU(I),I=1,NUMTAU)
         WRITE(6,825)
         WRITE(6,821)(ITAU(I),I=1,NUMTAU)
      ENDIF
C
      READ(5,*)NFTET,NUMTET
      WRITE(6,947)NFTET,NUMTET
      IF (NFTET.NE.0) THEN
         READ(5,*)(ITET(I),I=1,NUMTET)
         WRITE(6,948)
         WRITE(6,821)(ITET(I),I=1,NUMTET)
      ENDIF
C
      READ(5,*) NFDH,NUMDH
      WRITE(6,953) NFDH,NUMDH
      IF (NFDH.NE.0) THEN
         READ(5,*)(IDH(I),I=1,NUMDH)
         WRITE(6,954)
         WRITE(6,821)(IDH(I),I=1,NUMDH)
      ENDIF
      WRITE(6,*)
C
      READ(5,*)NFHT,NUMHT
      WRITE(6,990)NFHT,NUMHT
      IF (NFHT.NE.0) THEN
         READ(5,*)(IHT(I),I=1,NUMHT)
         WRITE(6,992)
         WRITE(6,994)(IHT(I),I=1,NUMHT)
      ENDIF
      WRITE(6,*)
      if(nselt.ge.0)then
        IF (INTEGRATOR.EQ.0) THEN
          ATIME=TIME/1440.0D0
        ELSE
          ATIME=TIME
        ENDIF
        IF (INTEGRATOR.EQ.1) THEN
c
c     if lll<0 then xl is the step size for fixed step size integration
c     if lll>0 then xl is the first step size for variable  step size 
c     integration and lll determines the energy conservation
            xl=atime
c     tf is the total time that the integration routine is called for
            if(lll.lt.0) then
               tf=xl
            else
               tf=atime*dble(ns)
            endif
        ENDIF
      endif
432   continue
C
C         SET FLAGS AND PARAMETERS
C
C             NSFLAG:  FOR CALCULATING QMIN,QMAX, AND PMAX IN 'INITQP'
C             NAM:     FOR CALCULATTING ANGULAR MOMENTUM IN 'ROTN',
C                      'INITQP', AND 'ORTHAN'
C             NAST:    FOR PASSING BARRIER IN 'MAIN'
C             NFINAL:  FOR CALCULATING EROT IN 'FINAL'.
C
      NSFLAG=0
      NAM=0
      NI=I3N
      NID=2*NI
      NTZ=0
C *** DETERMINE THE STARTING TIME.
c      CALL CPUSEC(TUSER1,TTOTAL1)
C     CALL CPUSEC(TUSER1,TTOTAL1,TCPU1)
      SECADD=0.0D0
C
C         IF NCHKP .NE. 0 THEN CALCULATIONS ARE RESTARTED FROM CHECKPOINT
C         FILE (ANY TYPE OF VENUS CALCULATION PRODUCES A CHECKPOINT FILE
C         IN UNIT 50, WHICH CONTAINS Q,P,QDOT,etc. AND THE RANDOM NUMBER
C         GENERATOR INFORMATION).
C
C         IF NCHKP IS 1, CALCULATIONS ARE RESTARTED WHERE LEFT, IN
C         CASE OF A COMPUTER CRASH, OR UNDER-ESTIMATION OF THE CYCLE
C         NUMBER IN REACTION PATH FOLLOWING, NUMBER OF TRAJECTORIES,
C         ETC. (THE CODE FOLLOWS INPUT UPDATES)
C
C         IF NCHKP IS -1, NEW CALCULATIONS ARE STARTED WITH CARTESIAN
C         COORDINATES AND MOMENTA READ FROM CHECKPOINT FILE. FOR EXAMPLE,
C         ONE CAN PERFORM A SEQUENCE OF CALCULATIONS WITH THIS OPTION :
C         MINIMUM ENERGY SEARCH, REFINEMENT BY REACTION PATH FOLLOWING,
C         NORMAL MODE ANALYSIS, WHERE Q AND P'S CAN BE PASSED THRU THE
C         CHECKPOINT FILE. THIS OPTION IS AVAILABLE FOR NSELT .NE. 2 ONLY.
C
C         IF NCHKP IS -2, THEN AN ARRAY FOR THE RANDOM NUMBER GENERATOR IS
C         READ IN FROM CHECKPOINT FILE IN UNIT 50.
C
      IF(NCHKP.EQ.0)GOTO 223

C
C  INITIALIZE COORDINATES FOR QM CALCULATIONS IN QQCHEM
C
c QM INIT
      IF (NMO.NE.0) THEN
         if(nmo.gt.0) then
           i3nn = natom*3
         else
C i3nn should be equal to natom*3 but it equals to i3n?
           i3nn = i3n
         endif
        DO count=1,I3NN
c         write(6,*)'The Initial co-ordinates are ',Q(count)
         QQCHEM(count) = Q(count)
c        write(6,*)'The Initial co-ordinates are ', QQCHEM(count)
        ENDDO
        call flush(6)
C Nwchem Intialization
       if(initialize.eqv..true.) then
        idop=1
        write (6,*) 'Before Init:stack=',stack,' heap=',heap,
     $              ' global=',global
        call flush(6)
        call VENUS_parallel_idop (myid,idop,theory, basis, natom,
     *  QQCHEM,labels,printName,qm_hessian,qm_grad,qm_energy,charge,
     *  stack,heap,global,qm_ini_file)
        initialize = .false.
c        write(6,*)'The initalize variable is ', initialize
       endif
C ENDIF for NMO
      endif


C
      IF (NCHKP.EQ.-2) THEN
         OPEN(50,FORM='UNFORMATTED')
         READ(50,ERR=222)QDUM,QDUM,QDUM,QDUM,TABDUM,VDUM,RANLST,GDUM,
     *        NFDUM,ISEED0,ISEED3,NDUM,NDUM,NDUM,NDUM,NDUM,IBFCTR,
     *        VI,OAMI,AMAI,AMBI,ETAI,ERAI,ETBI,ERBI
         CLOSE(50)
         WRITE(6,895)
         GOTO 451
      ENDIF
C
      OPEN(50,FORM='UNFORMATTED')
      READ(50,ERR=222)Q,P,QDOT,PDOT,TABLE,VRELO,RANLST,GTEMP,NFLAG,
     *ISEED0,ISEED3,NX,NC,NTZ,INTST,NAST,IBFCTR,
     *VI,OAMI,AMAI,AMBI,ETAI,ERAI,ETBI,ERBI
      CLOSE(50)
      nx=(nc/nip+1)*nip
      write(6,*) 'NX=',NX, 'NIP=', NIP
      KRE=1
      CALL DVDQ
      CALL ENERGY
      IF (NCHKP.EQ.-1) WRITE(6,894)
      IF (NCHKP.EQ.1) WRITE(6,897)
C
      IF (NSELT.LT.0) THEN
         IF (NSELT.NE.-3.OR.NCHKP.EQ.-1) NTZ=1
         NDUM=MAX0(1,NTZ)
         DO J=1,NDUM
            READ(5,*)(QDUM(I),I=1,I3N)
c           read(5,*)(q(i),i=1,i3n)
         ENDDO
         DO I=1,I3N
            P(I)=0.0D0
         ENDDO
         IF (NCHKP.EQ.-1) THEN
            NC=0
            NX=0
         ENDIF
c
cki Added mpatho to keep the old venus compatibility
c
         IF (NSELT.EQ.-2) then
            if(nmep.eq.0) then
              call mpatho
            else
              CALL MPATH
            endif
         endif
c
         CALL GWRITE
         IF (NSELT.EQ.-3) THEN
            CALL STATPT
            GOTO 451
         ENDIF
         I=0
         CALL NMODE(NATOMS,I)
      ENDIF
C
      IF (NCOOR.EQ.1) CALL SYBMOL
      IF (NSELT.EQ.0.OR.NSELT.EQ.1) THEN
         NDUM=MAX0(1,NTZ)
         IF (NC.EQ.NS) NDUM=NDUM-1
         DO J=1,NDUM
            READ(5,*)(QDUM(I),I=1,I3N)
            READ(5,*)(QDUM(I),I=1,I3N)
         ENDDO
         IF (NCHKP.EQ.-1) THEN
            NC=0
            NTZ=1
            NX=NIP
         ENDIF
      ENDIF
      call flush(6)
C
      IF (NC.EQ.NS.AND.NCHKP.EQ.1) THEN
         NTZ=NTZ-1
         GOTO 451
      ENDIF
      IF (NSELT.EQ.0.AND.NCHKP.EQ.-1) GOTO 425
      CALL GWRITE
      IF (NSELT.EQ.1) GOTO 447
      IF (NSELT.EQ.0) GOTO 400
      IF (NCHKP.EQ.-1) THEN
         WRITE(6,891)
         STOP
      ENDIF
      GOTO 400
  222 WRITE(6,890)
      STOP
C
C         INITIALIZE ARRAY OF RANDOM NUMBERS FOR ROUTINES RAND0 AND RAND1.
C
  223 IF (NSELT.EQ.2.OR.NSELT.EQ.3) CALL RANDST(ISEED)
C
C         INCREMENT NUMBER OF TRAJECTORIES NTZ  BY 1
C
  451 NTZ=NTZ+1
C *** DETERMINE THE COMPUTATION TIME FOR INTEGRATION OF TRAJECTORY.
C     CALL CPUSEC(TUSER2,TTOTAL2)
c      CALL CPUSEC(TUSER2,TTOTAL2,TCPU2)
C     WRITE(6,510) (SECADD+TUSER2-TUSER1),(SECADD+TTOTAL2-TTOTAL1)
C     TUSER1=TUSER2
C     TTOTAL1=TTOTAL2
C     TCPU1=TCPU2
C 510 FORMAT(//' USER AND REAL TIME FOR INTEGRATION:',
C    *       2F15.1,' SECONDS'//)
C     *       3F8.1,' SECONDS'//)
      call flush(6)
C
C          BARRIER CROSSING PARAMETERS INITIALIZATION
C
      IBAR=0
      GAO=0.0D0
C
C         WRITE RANDOM NUMBER ARRAY AND OTHER RELEVANT INFORMATION
C         TO CHECKPOINT FILE IN UNIT 50.
C
      OPEN(50,FORM='UNFORMATTED')
      REWIND(50)
      WRITE(50)Q,P,QDOT,PDOT,TABLE,VRELO,RANLST,GTEMP,NFLAG,ISEED0,
     *ISEED3,NX,NC,NTZ,INTST,NAST,IBFCTR,
     *VI,OAMI,AMAI,AMBI,ETAI,ERAI,ETBI,ERBI
      CLOSE(50)      
      IF (NSELT.EQ.2.OR.NSELT.EQ.3) WRITE(6,896)(ISEED3(9-I),I=1,8)
C
C         UPDATE PRINT ARRAYS FOR MODE POPULATIONS AFTER COMPLETING
C         NTZ TRAJECTORIES.
C
      IF (NTZ.GE.3 .AND. MPLOT.EQ.1) THEN
         NDUM=NTZ-1
C
C              NORMAL MODE POPULATIONS
C
         IF (NMA.NE.0) THEN
            CALL WENMOD(ETIM,ENSAV,ENSQ,NMA,ISAV,NDUM,NEVIBN,MNTR)
         ENDIF
C
C              THE EXCITED LOCAL MODE POPULATIONS
C
         IF (NACTA.EQ.4) THEN
            NOUT=0
            CALL WEBOND(ETIM,ESAV,ESQ,NEVIB,ISAV,NDUM,NOUT)
            NOUT=1
            CALL WEBOND(ETIM,ESAV,ESQ,NEVIBU,ISAV,NDUM,NOUT)
            NOUT=-1
            CALL WEBOND(ETIM,ESAV,ESQ,NEVIBL,ISAV,NDUM,NOUT)
C
C              REMAINING LOCAL MODE POPULATIONS
C
            IF (NLM.NE.0) THEN
               CALL WLBOND(ETIM,EBSAV,EBSQ,NLM,ISAV,NDUM)
            ENDIF
         ENDIF
      ENDIF
C
      IF (NTZ.GT.NT) then
C       CALL CPUSEC(TUSER2,TTOTAL2)
c        CALL CPUSEC(TUSER2,TTOTAL2,tcpu2)
C       WRITE(6,515) (SECADD+TUSER2-TUSER1),(SECADD+TTOTAL2-TTOTAL1)
C       TUSER1=TUSER2
C       TTOTAL1=TTOTAL2
C       TCPU1=TCPU2
  515 FORMAT(//' USER AND REAL TIME FOR VENUS:',
     *       2F15.1,' SECONDS'//)
c     *       3F8.1,' SECONDS'//)

C QM finalizing Block
c
       IF (NMO.NE.0) THEN
        idop=5
        call VENUS_parallel_idop (myid,idop,theory, basis, natom,
     *  QQCHEM,labels,printName,qm_hessian, qm_grad,qm_energy,charge,
     *  stack,heap,global,qm_ini_file)
        WRITE(6,*)'FINAL IS CALLED'
        call flush(6)
       ENDIF
       STOP
      endif

      NC=0
      DO I=1,8
         ISEED0(I)=ISEED3(I)
      ENDDO
      VRELO=0.0D0
      INTST=0
      NAST=2
      NFINAL=0
      KRE=1
C
C  INITIALIZE THE Q AND P ARRAYS FOR QM CODE
C
      J=3*NATOMA(1)
      DO I=1,J
         Q(I)=QZA(1,I)
         P(I)=0.0D0
      ENDDO
      K=3*NATOMB(1)
      DO I=1,K
         Q(J+I)=QZB(1,I)
         P(J+I)=0.0D0
      ENDDO

C
C  INITIALIZE COORDINATES FOR QM CALCULATIONS IN QQCHEM
C
c QM INIT
      IF (NMO.NE.0) THEN
         if(nmo.gt.0) then
           i3nn = natom*3
         else
C i3nn should be equal to natom*3 but it equals to i3n?
           i3nn = i3n
         endif
        DO count=1,I3NN
c         write(6,*)'The Initial co-ordinates are ',Q(count)
         QQCHEM(count) = Q(count)
c        write(6,*)'The Initial co-ordinates are ', QQCHEM(count)
        ENDDO
        call flush(6)
C Nwchem Intialization
       if(initialize.eqv..true.) then
        idop=1
        write (6,*) 'Before Init:stack=',stack,' heap=',heap,
     $              ' global=',global
        call flush(6)
        call VENUS_parallel_idop (myid,idop,theory, basis, natom,
     *  QQCHEM,labels,printName,qm_hessian,qm_grad,qm_energy,charge,
     *  stack,heap,global,qm_ini_file)
        initialize = .false.
c        write(6,*)'The initalize variable is ', initialize
       endif
C ENDIF for NMO
      endif

C
C         INITIALIZE VMAX FOR LOCATING MAXIMUM POTENTIAL ENERGY
C         DURING THE TRAJECTORY
C
      VMAX=-1.0D20
C
C         INITIALIZE THE Q AND P ARRAYS
C
      J=3*NATOMA(1)
      DO I=1,J
         Q(I)=QZA(1,I)
         P(I)=0.0D0
      ENDDO
      K=3*NATOMB(1)
      DO I=1,K
         Q(J+I)=QZB(1,I)
         P(J+I)=0.0D0
      ENDDO
C
C         SELECT INITIAL CONDITIONS
C
      CALL SELECT
C *** DETERMINE THE COMPUTATION TIME FOR SELECT
C     CALL CPUSEC(TUSER2,TTOTAL2)
C     WRITE(6,520) (SECADD+TUSER2-TUSER1),(SECADD+TTOTAL2-TTOTAL1)
C     TUSER1=TUSER2
C     TTOTAL1=TTOTAL2
  520 FORMAT(//' USER AND REAL TIME FOR SELECT:',
     *       2F15.1,' SECONDS'//)
c     *       3F8.1,' SECONDS'//)
      call flush(6)
c
      IF (NCOOR.EQ.1.AND.NSELT.NE.-2) CALL SYBMOL
C
C         SET ARRAYS FOR MODE POPULATIONS
C
      TI=0.0D0
      IF (MPLOT.NE.1) GOTO 425
      ISAV=1
      ETIM(ISAV)=TI
C
C              NORMAL MODE ARRAYS
C
      IF (NMA.NE.0) THEN
C
C              CALCULATE INITIAL ENERGY AND BIN WIDTH FOR THE NORMAL
C              MODE WHOSE POPULATION IS TO BE MONITORED.
C
         IF (MNTR.NE.0) THEN
            NEVIBN(ISAV)=NEVIBN(ISAV)+1
            I=MNTR
            WWA(I)=WWA(I)/C6*cm2cal
            ENN=(ANQA(I)+0.5D0)*WWA(I)
            EDELTN=0.5D0*WWA(I)
            WWA(I)=WWA(I)*C6/cm2cal
         ENDIF
C
         CALL ENMODE(ENM,NMA)
         DO I=1,NMA
            ENSAV(I,ISAV)=ENSAV(I,ISAV)+ENM(I)
            ENSQ(I,ISAV)=ENSQ(I,ISAV)+ENM(I)*ENM(I)
         ENDDO
      ENDIF
C
C              THE EXCITED LOCAL MODE ARRAYS
C
      IF (NACTA.EQ.4) THEN
         CALL EBOND(EBCH,EKCH,RCH,NEXM)
         NEVIB(ISAV)=NEVIB(ISAV)+1
         ESAV(ISAV)=ESAV(ISAV)+EBCH
         ESQ(ISAV)=ESQ(ISAV)+EBCH**2
C
C              REMAINING LOCAL MODE ARRAYS
C
         IF (NLM.NE.0) THEN
            DO I=1,NLM
               J=MNLM(I)
               CALL EBOND(EBM(I),EK,RCH,J)
               EBSAV(I,ISAV)=EBSAV(I,ISAV)+EBM(I)
               EBSQ(I,ISAV)=EBSQ(I,ISAV)+EBM(I)*EBM(I)
            ENDDO
         ENDIF
      ENDIF
C
  425 CONTINUE
      CALL GWRITE
C
C         SEARCH FOR STATIONARY POINTS IF NSELT=-3
C         (OPTION NOT AVAILABLE IN THIS VERSION)
C
C         PERFORM REACTION PATH FOLLOWING IF NSELT=-2
C         PERFORM NORMAL MODE ANALYSIS IF NSELT=-1
C
      I=0
      IF (NSELT.EQ.-1) THEN
         CALL NMODE(NATOMS,I)
         GOTO 451 
      ENDIF
c
cki  nmep added for keeping compatibility with old venus
c
      IF (NSELT.EQ.-2) then
         if(nmep.eq.0)then
           call mpatho
         else
           CALL MPATH
         endif
      endif
c
      NC=0
      NX=NIP
      IF (NSELT.EQ.-3) THEN
         CALL STATPT
         GOTO 451
      ENDIF
      NXPLOT=NPLOT
C
C         INTEGRATE TRAJECTORY
C
  402 NCZ=NC
C
      call flush(6)
      IF (INTEGRATOR.EQ.0) THEN
         CALL PARTI
         DO NC=1,6
            CALL RUNGEK
         ENDDO
         NC=NCZ+6
      ELSEIF (INTEGRATOR.EQ.2) THEN
         CALL DVDQ
      ENDIF
C
  400 NC=NC+1
C
c      write(6,*) 'ETot=', H 
      IF (INTEGRATOR.EQ.0) THEN
         CALL ADAMSM
      ELSEIF (INTEGRATOR.EQ.1) THEN
         if(lll.ge.0) then
            write(6,*)'DOING TIME VARIABLE INTEGRATION'
            call RADAU(tf,xl,lll,nip,niter)
            IF (NTZ.GT.NT) then
C             CALL CPUSEC(TUSER2,TTOTAL2)
C              CALL CPUSEC(TUSER2,TTOTAL2,tcpu2)
C             WRITE(6,515) (SECADD+TUSER2-TUSER1),(SECADD+TTOTAL2-
C    *          TTOTAL1)
C             TUSER1=TUSER2
C             TTOTAL1=TTOTAL2
              STOP
            endif
            goto 451
         else
            call RADAU(tf,xl,lll,nip,niter)
         endif
      ELSEIF (INTEGRATOR.EQ.2) THEN
         CALL SYMPLE(lll)
      elseif(integrator.eq.3)then
         call verlet(lll)
      elseif(integrator.eq.4)then
c   decide which integrator to use
c   lll = 1 - Projection method by Yu Zhuang
c   lll = 2 - Projection method similar to Gaussian
c   lll = 3 - Generalized coordinates projection method by Yu Zhuang
c   lll = 4 - Time based corrector by Yu Zhuang
c   options 3 and 4 not implimented yet !!!!
        call hessint(lll)
        if(ntz.gt.nt)stop
        goto 451
      ENDIF
C
C         DETERMINE MODE POPULATIONS
C
      IF (MPLOT.NE.1.or.NC.LT.NXPLOT) GOTO 426
      NXPLOT=NXPLOT+NPLOT
      TI=DBLE(NC)*TIME
      IF (ISAV.GE.2000) GOTO 426
      ISAV=ISAV+1
      ETIM(ISAV)=TI
C
C              NORMAL MODE ENERGIES
C
      IF (NMA.NE.0) THEN
         CALL ENMODE(ENM,NMA)
         IF (MNTR.NE.0) THEN
            I=MNTR
            IF (ABS(ENN-ENM(I)).LT.EDELTN) NEVIBN(ISAV)=NEVIBN(ISAV)+1
         ENDIF
         DO I=1,NMA
            ENSAV(I,ISAV)=ENSAV(I,ISAV)+ENM(I)
            ENSQ(I,ISAV)=ENSQ(I,ISAV)+ENM(I)*ENM(I)
         ENDDO
      ENDIF
C
C              THE EXCITED LOCAL MODE ENERGIES AND POPULATIONS
C
      IF (NACTA.EQ.4) THEN
         CALL EBOND(EBCH,EKCH,RCH,NEXM)
         IF (ABS(ENON-EBCH).LT.EDELTA) NEVIB(ISAV)=NEVIB(ISAV)+1
         IF (ABS(ENU-EBCH).LT.EDELTU) NEVIBU(ISAV)=NEVIBU(ISAV)+1
         IF (ABS(ENL-EBCH).LT.EDELTL) NEVIBL(ISAV)=NEVIBL(ISAV)+1
         ESAV(ISAV)=ESAV(ISAV)+EBCH
         ESQ(ISAV)=ESQ(ISAV)+EBCH*EBCH
C
C              REMAINING LOCAL MODE ENERGIES
C
         IF (NLM.NE.0) THEN
            DO I=1,NLM
               J=MNLM(I)
               CALL EBOND(EBM(I),EK,RCH,J)
               EBSAV(I,ISAV)=EBSAV(I,ISAV)+EBM(I)
               EBSQ(I,ISAV)=EBSQ(I,ISAV)+EBM(I)*EBM(I)
            ENDDO
         ENDIF 
      ENDIF
  426 CONTINUE
C
c      write (6,*) 'NC=', NC, 'NX=', NX
      IF (NFINAL.EQ.1) GOTO 414
      IF (NC.GE.NS) GOTO 450
      IF (NC.EQ.NX) GOTO 449
      IF (NSELT.EQ.1) GOTO 400
C
C         FIRST AND LAST TIME FOR A PARTICULAR ATOM THE SURFACE-ABOVE-
C         HEIGHT OF WHICH IS LESS THAN A PARTICULAR VALUE (GIVEN ATOM 
C         NUMBER AND HEIGHT VALUE). ALSO TO FIND MINIMUM HEIGHT VALUE 
C
c      IF (NSURF.EQ.1) CALL HTURN(1,13.0)
C
C         SUBROUTINE TEST CHECKS FOR EVENTS.
C
C             NTEST=0, A REACTION HAS NOT OCCURRED.
C             NTEST=1, REACTIVE COORDINATE EQUALS RBAR.
C             NTEST=2, REACTIVE COORDINATE EQUALS RMAX.
C
C         TO TERMINATE A TRAJECTORY IF NUMBER OF INNER TURNING POINTS
C         EXCEEDS 100:
C
C     IF (INTST.GT.99) THEN
C        NPATH=1+NPATHS
C        GOTO 410
C     ENDIF
C
      CALL VETEST
      IF (NTEST.EQ.1.AND.NAST.EQ.0) GOTO 406
      IF (NTEST.EQ.2.AND.NAST.EQ.1) GOTO 410
      IF (NTEST.EQ.2.AND.NAST.EQ.0) GOTO 410
      IF (NTEST.EQ.0.AND.NAST.NE.0) NAST=0
      IF (NTEST.EQ.1.AND.NAST.EQ.2) NAST=1
      call flush(6)
      GOTO 400
  406 NAST=1
      WRITE(6,*)
      WRITE(6,910)NPATH
      CALL ENERGY
      CALL GWRITE
      GOTO 400
  410 CONTINUE
      CALL DVDQ
      CALL ENERGY
      CALL GWRITE
  414 CONTINUE
      CALL FINAL
C
C         AVERAGE A AND B ROTATIONAL ENERGIES OVER NCROT CYCLES
C
      ERAVA(KRE)=EROTA
      ERAVB(KRE)=EROTB
      KRE=KRE+1
      IF (KRE.LE.NCROT) GOTO 400
      EROTA=0.0D0
      EROTB=0.0D0
      DO 412 I=1,NCROT
         EROTA=EROTA+ERAVA(I)
  412 EROTB=EROTB+ERAVB(I)
      EROTA=EROTA/DBLE(NCROT)
      EROTB=EROTB/DBLE(NCROT)
      SDA=0.0D0
      SDB=0.0D0
      DO 413 I=1,NCROT
         SDA=SDA+(EROTA-ERAVA(I))**2
         SDB=SDB+(EROTB-ERAVB(I))**2
  413 CONTINUE 
      SDA=SQRT(SDA/(NCROT-1))
      SDB=SQRT(SDB/(NCROT-1))
      CALL GFINAL
      GOTO 451
  447 DO 448 I=1,I3N
  448 P(I)=0.0D0
      GOTO 402
  449 CALL ENERGY
      REWIND 77
      NX=NX+NIP
      CALL GWRITE
      IF (NSELT.EQ.1) GOTO 447
      GOTO 400
  450 CALL ENERGY
      CALL GWRITE
      GOTO 451
      END

