      SUBROUTINE DVDQ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'

c QM Commom BLock
      COMMON /QMINFO/ qm_choice

C QM related variables
      INTEGER MYID
      character*256 theory
      character*256 basis
      integer natom,qm_choice,idop
      double precision charge
      double precision coordinates (nda3)
      character*16 labels(nda)
      character*256 qm_ini_file

C
C         CALCULATE POTENTIAL ENERGY PARTIAL DERIVATIVES WITH 
C         RESPECT TO COORDINATES (PDOT)
C
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/GAUSS/ICHRG,IMULTP,IAN(NDAyf),VGAUSS
      COMMON/GAUSS2/GAUHES(NDIHE)
      common/mdsp/mdflag
c
      common/ghessb/FA(NDA3,NDA3),trad,rmin,dt,nstep,
     *  nhesup,nhessf,nip
      DIMENSION AHESS(NDIHE)
      LOGICAL SAVFIL,OK
      data H2KCAL/627.5095D0/,B2A/0.529177249D0/
      save
C
C         ZERO PDOT'S
C
      DO I=1,I3N
         PDOT(I)=0.0D0
      ENDDO
C
C         CALCULATE PARTIALS OF HARMONIC STRETCHES
C
      IF (NST.NE.0) THEN
         DO I=1,NST
            CALL STRET(I)
         ENDDO
      ENDIF
C
C   code MODIFIED by navdeep on jun7,04.....STARTS
C   
C	CALCULATE PARTIALS OF CRCO6 POTENTIAL FUNCTION
C	Phys. Chem. Chem. Phys., 2001,3,2306-2314
C	
C	LOG: Ming add this feature on 11/13/2003
C
      IF (NCRCO6 .NE. 0) THEN
          CALL CRCO6
      ENDIF
C     
C
C
C         CALCULATE PARTIALS OF MORSE STRETCHES
C
        J=1
        IF (NCRCO6.NE.0) THEN 
              J=7
        ENDIF
        IF (NM.NE.0) THEN
            DO I=J,NM
               CALL MORSE(I)
            ENDDO
        ENDIF
C
C       code MODIFIED by navdeep on jun7,04.....ENDS
C              
C         CALCULATE PARTIALS OF HARMONIC BENDS
C
      IF (NB.NE.0) THEN
         DO I=1,NB
            CALL HBEND(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF ALPHA BENDS
C
      IF (NA.NE.0) THEN
         DO I=1,NA
            CALL HALPHA(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF LENNARD-JONES INTERACTIONS
C
      IF (NLJ.NE.0) THEN
         CALL LENJ(1,NLJ)
      ENDIF
C
C         CALCULATE PARTIALS OF TORSION POTENTIALS
C
      IF (NTAU.NE.0) THEN
         DO I=1,NTAU
            CALL HTAU(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF REPULSIONS AND ELECTROSTATIC
C         INTERACTIONS
C
      IF (NEXP.NE.0) THEN
         DO I=1,NEXP
            CALL HEXP(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF GHOST PAIRS
C
      IF (NGHOST.NE.0) THEN
         DO I=1,NGHOST
            CALL GHOST(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF TETRAHEDRAL CENTERS
C
      IF (NTET.NE.0) THEN
         DO I=1,NTET
            CALL TETRA(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF R-R COUPLING
C
      IF (NVRR.NE.0) THEN
         DO I=1,NVRR
            CALL VRR(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF R-THETA COUPLING
C
      IF (NVRT.NE.0) THEN
         DO I=1,NVRT
            CALL VRT(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF THETA-THETA COUPLING
C
      IF (NVTT.NE.0) THEN
         DO I=1,NVTT
            CALL VTT(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF DIHEDRAL ANGLES
C
      IF (NANG.NE.0) THEN
         DO I=1,NANG
            CALL DANGLE(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF AXILROD-TELLER FUNCTION
C
      IF (NAXT.NE.0) THEN
         DO I=1,NAXT
            CALL AXT(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF SN2 FUNCTION
C
      IF (NSN2.NE.0) THEN
         CALL VSN2
      ENDIF
C
C         CALCULATE PARTIALS OF RYDBERG POTENTIAL FUNCTION
C
      IF (NRYD.NE.0) THEN
         DO I=1,NRYD
            CALL RYDBG(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF HARTREE-FOCK DISPERSION FUNCTION
C
      IF (NHFD.NE.0) THEN
         CALL HFD(1,NHFD)
      ENDIF
C
C         CALCULATE PARTIALS OF LEPS(A) POTENTIAL FUNCTION
C
      IF (NLEPSA.NE.0) THEN
         DO I=1,NLEPSA
            CALL LEPS1(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF LEPS(B) POTENTIAL FUNCTION
C
      IF (NLEPSB.NE.0) THEN
         DO I=1,NLEPSB
            CALL LEPS2(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF DOUBLE MANY-BODY EXPANSION PES
C
      IF (NDMBE.NE.0) CALL DMBE
C
C         CALCULATE PARTIALS OF RELAX INTERACTIONS
C
      IF (NRAX.NE.0) THEN
         DO I=1,NRAX
            CALL RELAX(I)
         ENDDO
      ENDIF
C
C         CALCULATE PARTIALS OF H---H NONBONDED FUNCTIONS 
C
      IF (NONB.NE.0) THEN
         DO I=1,NONB 
            CALL HNONB(I)
         ENDDO
      ENDIF
C
C Compute QM gradient
c      if(qm_choice.gt.0) then
       IF (NMO.NE.0) THEN
         if(nmo.gt.0)then
           natom=nmo
           i3nn=natom*3
         else
           natom=natoms
           i3nn=i3n
         endif

c Call the the function which calculates Energy and Gradient for Nwchem
c        write(*,*)'charge ....',charge
         call QMCALC(0)
      endif
      DO I=1,I3N
         PDOT(I)=-PDOT(I)
      ENDDO
      RETURN
      END
