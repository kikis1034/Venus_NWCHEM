      SUBROUTINE MPATH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      character*80 title1,title2
C
C         PERFORM REACTION PATH FOLLOWING
C
C
C  In this version, rate constant for unimolecular dissociation will also
C   be performed using the reactant frequencies read in.
C  The routines from ART program are inserted.  They are marked by cART.
C  Rate constants are calculated only for C->A+B channel.
C  Usually the CVTST calculations are from A+B to C. So it is important
C  to enter correct frequencies for C.
C  Modified by Kihyung Song, July 2000.
C        nchn=1 for C -> A + B reaction (dissociation)
C        nchn=2 for A + B -> C reaction (recombination)
C        V0: energy which venus gives for C (combined) for nchn=1
C                                         A + B (far away) for nchn=2
C        NWM=# FREQUENCIES FOR THE MOLECULE.(CM-1)
C        NIRM=# INTERNAL ROTORS FOR THE MOLECULE.
C        WTM=WEIGHT OF MOLECULE (AMU)
C        SIGM=EXTERNAL SYMMETRY # FOR THE MOLECULE.
C        TXM,TYM,TZM,=MOMENTS OF INERTIA IN THE X,Y, AND Z DIRECTIONS
C        For linear molecules (2D-rotor) put tzm as 0.0 or less.
c        For Surfaces, put txm=0 to avoid the rotational d.o.f. in k
C        FOR THE MOLECULE.  IN AMU-A2.
C        WM(I)=FREQUENCIES FOR MOLECULE(CM-1)
C        TRM(I)=REDUCED MOMENTS OF INERTIA FOR INTERNAL ROTORS IN THE
C        MOLECULE(AMU-A2).
C        SYMM(I)=SYMMETRY NUMBERS FOR THE INTERNAL ROTORS IN THE
C        MOLECULE.
C        THE SYMBOLS ARE THE SAME FOR THE COMPLEX AND RADICALS EXCEPT C,
C        R1, AND R2 REPLACE M FOR THE COMPLEX, RADICAL1, AND RADICAL2,
C        RESPECTIVELY.
C        A-FACTOR FOR THE FORWARD REACTION IS IN UNITS OF SEC-1 AND
C        A-FACTOR FOR REVERSE REACTION IS IN LT./MOLE-SEC.
C        akust(i): Miller's UST rate constant at T(i)
C
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
      COMMON/EIGVL/EIG(NDA3yf)
      COMMON/RKUTTA/RAA1,RA1,RA2,RA3,RB1,RB2,RB3,RC1,RC2
      COMMON/GPATHB/WM(NDA3),TEMP(NDP),AI1D(5),AAI(2),BBI(2),SYMM(5),
     *SYMA,SYMB,GTEMP(NDP),NFLAG(NDP),N1DR,N2DR
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/SYBB/TITLE1,TITLE2,SYBTI
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      common/frqc/wc(nda3),wco(nda3),ezp,ezpo,ezm,nrq,nvq
      common/frqr/wr1(nda3),wr2(nda3),txm,tym,tzm,txr1,tyr1,tzr1,
     * txr2,tyr2,tzr2,wtm,wtr1,wtr2,nchn,nwm,nwr1,nwr2,nirm,nirr1,nirr2
      DIMENSION BQ(nda3),PJ(nda3yf,nda3yf),FAP(nda3yf,nda3yf),
     * QM(nda3),GZ(nda3),aa(nda3yf,nda3yf),qtemp(nda3),
     * sm(ndp),hm(ndp),akust(ndp),akmin(ndp),QCM(3),VCM(3),trm(5),
     * trr1(5),trr2(5),symr1(5),trc(5),symc(5),symr2(5),AM(4)
      DATA RHO/1.0D-08/
      DATA TFACT/2.81837552D05/
C
  150 FORMAT(' XXXXXXXXXXXXXXXXXXXXX REACTION PATH  NC = ',I7,
     *' XXXXXXXXXXXXXXXXXXXXXXX')
  151 FORMAT(//10X,'*****CALCULATE REACTION PATH*****'//
     *4X,'NS=',I10,'  NIP=',I10,'  DS=',F13.10/)
  152 FORMAT(/4X,'NC=',I8,5X,'POTENTIAL ENERGY(KCAL/MOLE)=',F9.4)
  153 FORMAT(10X,'CARTESIAN COORDINATES')
  154 FORMAT(3X,1P3D18.10,3X,3D18.10)
  155 FORMAT(10X,'EIGENVALUES')
  156 FORMAT(3X,1P7D14.6)
  157 FORMAT(10X,'PRINCIPAL MOMENTS OF INERTIA: IX,IY,IZ')
  158 FORMAT(4X,' ATOMS  ',5X,'BOND LENGTH(A)')
  159 FORMAT(2X,2I4,10X,1PE15.8)
  160 FORMAT(1X,3F11.7)
  161 FORMAT(10X,'REACTANT VIBRATIONAL FREQUENCIES(CM-1)')
  162 FORMAT(4X,'TEMPERATURES=',10F9.2)
  164 FORMAT(4X,'NUMBER OF 1-DIMENSIONAL ROTORS =',I2)
  165 FORMAT(10X,'SYMA=',F4.1)
  166 FORMAT(10X,'SYMB=',F4.1)
  167 FORMAT(10X,'AI =',F8.3,'  SYMM =',F5.1)
  168 FORMAT(4X,'NUMBER OF 2-DIMENSIONAL ROTORS =',I2)
  171 FORMAT(10X,'MOMENTS OF INERTIA FOR FRAGMENT A: IX,IY,IZ')
  172 FORMAT(10X,'REDUCED MOMENTS OF INERTIA FOR 2-D ROTOR:',
     *'  IX=',1PE14.6,'  IZ=',E14.6)
  173 FORMAT(5X,'ATOMS',7X,'ANGLE (DEGREES)')
  174 FORMAT(2X,3I3,8X,F8.3)
  175 FORMAT(10X,'MOMENTS OF INERTIA FOR FRAGMENT B: IX,IY,IZ')
6     FORMAT(' NUMBER OF FREQUENCIES=',I5,6X,'NUMBER OF ROTORS=',I5//)
8     FORMAT(' MASS=',F7.2/' SIGMA=',F7.2/
     *  ' ELECTRONIC PARTITION FUNCTION=',F5.2//)
9     FORMAT(' IXX=',1pe11.3,6X,'IYY=',e11.3,6X,'IZZ=',e11.3//)
13    FORMAT(' REDUCED MOMENT OF INERTIA=',F7.2,3X,'SYMMETRY NUMBER=',
     *F5.1)
C
C  Initialize moment of inertia for complex
C  It will be changed to actual values if necessary
C
      txco=0.d0
      tyco=0.d0
      tzco=0.d0
      txm=0.d0
      tym=0.d0
      tzm=0.d0
      txr1=0.d0
      tyr1=0.d0
      tzr1=0.d0
      txr2=0.d0
      tyr2=0.d0
      tzr2=0.d0
C
C         READ CYCLE COUNTS AND INTEGRATION PARAMETERS
C
      READ(5,*)NS,NIP,DS
      WRITE(6,151)NS,NIP,DS
C
C         READ TEMPERATURES(K) FOR CALCULATING THE FREE ENERGY
C         ALONG THE REACTION PATH
C
      READ(5,*)NTEMP
      READ(5,*)(TEMP(I),I=1,NTEMP)
      WRITE(6,162)(TEMP(I),I=1,NTEMP)
C
C         CALCULATE TOTAL MASS
C
      WTM=0.0D0
      DO I=1,NATOMS
         WTM=WTM+W(I)
      enddo
c
c  Additional data to calculate rate constant k(T)
c  enter 1 for C -> A+B dissociation reaction
c  enter 2 for A+B -> C recombination reaction
c
      nrq=1
      read(5,*)nchn,v0
      nwm=I3N-6
C
C         READ REACTANT VIBRATIONAL FREQUENCIES(CM-1)
C
      if(nchn.eq.1)then
        write(6,'(a)')' C -> A + B channel rate constants'
        write(6,'(a)')'  V0 should be the energy for C'
        write(6,'(a,1p3d10.3)')'   V0=',v0
        write(6,*)' **** Input data for Molecule (C) ****'
C
C         READ PARAMETERS FOR CALCULATING 1- AND 2-DIMENSIONAL INTERNAL
C         ROTOR PARTITION FUNCTIONS.  N2DR MAY EQUAL 0, 1, OR 2.
C         In this implication, n2dr is not used. (Ki Song, July 2000)
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
        nirm=n1dr
        nwm=nwm-nirm
        if(nacta.eq.1)then
          READ(5,*)(WM(I),I=1,nwm)
          read(5,*)txm,tym,tzm
        else
          call cfreq(0)
        endif
        WRITE(6,161)
        WRITE(6,156)(WM(I),I=1,nwm)
        WRITE(6,6)nwm,nirm
        read(5,*)sigm,qelm
        WRITE(6,8)wtm,sigm,qelm
        if(txm.eq.0.)nrq=0
        WRITE(6,9)txm,tym,tzm
        if(nirm.gt.0)then
          WRITE(6,*)'  Data for Internal Rotor'
          do i=1,nirm
            trm(i)=ai1d(i)
            WRITE(6,13)trm(i),symm(i)
          enddo
        endif
      else if(nchn.eq.2) then
        write(6,'(a)')' A + B -> C channel rate constants'
        write(6,'(a)')'  V0 should be the energy for A+B'
        write(6,'(a,1p3d10.3)')'   V0=',v0
c   Input for Fragment 1 (A)
c
c   nlina,nlinb=1 for linear, 0 for nonlinear
c
        write(6,*)' **** Input data for Fragment 1 (A) ****'
        read(5,*)natoma(2),nlina,nirr1
        nwr1=max(0,3*natoma(2)-6+nlina-nirr1)
        WRITE(6,6)nwr1,nirr1
        k=natoma(2)
        read(5,*)(la(2,j),j=1,k)
        wtr1=0.d0
        do j=1,natoma(2)
          wtr1=wtr1+w(la(2,j))
        enddo
        if(nwr1.gt.0)then
          if(nacta.eq.1)then
            read(5,*)(wr1(i),i=1,nwr1)
            read(5,*)txr1,tyr1,tzr1
          else
            call cfreq(1)
          endif
          write(6,161)
          write(6,156)(wr1(i),i=1,nwr1)
          if(txr1.eq.0.)nrq=0
          WRITE(6,9)txr1,tyr1,tzr1
          if(nirr1.gt.0)then
            WRITE(6,*)'  Data for Internal Rotor'
            do i=1,nirr1
              read(5,*)trr1(i),symr1(i)
              WRITE(6,13)trr1(i),symr1(i)
            enddo
          endif
        else
          do i=1,3
            qza(2,i)=0.d0
          enddo
        endif
        read(5,*)sigr1,qelr1
        WRITE(6,8)wtr1,sigr1,qelr1
c   Input for Fragment 2 (B)
        write(6,*)' **** Input data for Fragment 2 (B) ****'
        read(5,*)natomb(2),nlinb,nirr2
        nwr2=max(0,3*natomb(2)-6+nlinb-nirr2)
        WRITE(6,6)nwr2,nirr2
        k=natomb(2)
        read(5,*)(lb(2,j),j=1,k)
        wtr2=0.d0
        do j=1,natomb(2)
          wtr2=wtr2+w(lb(2,j))
        enddo
        if(nwr2.gt.0)then
          if(nacta.eq.1)then
            read(5,*)(wr2(i),i=1,nwr2)
            read(5,*)txr2,tyr2,tzr2
          else
            call cfreq(2)
          endif
          write(6,161)
          write(6,156)(wr2(i),i=1,nwr2)
          if(txr2.eq.0.)nrq=0
          WRITE(6,9)txr2,tyr2,tzr2
          if(nirr2.gt.0)then
            WRITE(6,*)'  Data for Internal Rotor'
            do i=1,nirr2
              read(5,*)trr2(i),symr2(i)
              WRITE(6,13)trr2(i),symr2(i)
            enddo
          endif
        else
          do i=1,3
            qzb(2,i)=0.d0
          enddo
        endif
        read(5,*)sigr2,qelr2
        WRITE(6,8)wtr2,sigr2,qelr2
        nwm=nwm-nirr1-nirr2
      else
        write(6,'(a)')' nchn should be 1 or 2'
        stop
      endif
c   Input for Complex
c   The moment of inertia are either set to be 0 for surface
c   or calculated on the way.
c   The only thing needed is sigma and qel (electronic multiplicity).
      write(6,*)' **** Input data for Complex (TS) ****'
      nwc=nwm-1
      if(nchn.eq.1)then
        nirc=nirm
        if(nirc.gt.0)then
          do i=1,nirc
            trc(i)=trm(i)
            symc(i)=symm(i)
          enddo
        endif
      else
        nirc=nirr1+nirr2
        if(nirr1.gt.0)then
          do i=1,nirr1
            trc(i)=trr1(i)
            symc(i)=symr1(i)
          enddo
        endif
        if(nirr2.gt.0)then
          do i=1,nirr2
            trc(i+nirr1)=trr2(i)
            symc(i+nirr1)=symr2(i)
          enddo
        endif
      endif
      wtc=wtm
      WRITE(6,6)nwc,nirc
      read(5,*)sigc,qelc
      WRITE(6,8)wtc,sigc,qelc
C
C         CALCULATE ZERO-POINT ENERGY OF REACTANT
C
      EZM=0.0D0
      IF (nchn.eq.1) THEN
         DO I=1,nwm
            EZM=EZM+WM(I)
         enddo
      else
        if(nwr1.ge.1) then
          do i=1,nwr1
            ezm=ezm+wr1(i)
          enddo
        endif
        if(nwr2.ge.1) then
          do i=1,nwr2
            ezm=ezm+wr2(i)
          enddo
        endif
      ENDIF
      ezm=ezm*0.5d0
c
c  Calculate V0
c
      if(nacta.eq.2)then
        do i=1,i3n
          qtemp(i)=q(i)
        enddo
        if(nchn.eq.2)then
          do i=1,natoma(2)
            do j=0,2
              q(la(2,i)*3-j)=qza(2,i*3-j)
              if(j.eq.0)q(la(2,i)*3)=q(la(2,i)*3)+1000.d0
            enddo
          enddo
          do i=1,natomb(2)
            do j=0,2
              q(lb(2,i)*3-j)=qzb(2,i*3-j)
            enddo
          enddo
        else
          do i=1,i3n
            q(i)=qza(2,i)
          enddo
        endif
        call dvdq
        call energy
        v0=v
        write(6,*)' V0=',v0
        do i=1,i3n
          q(i)=qtemp(i)
        enddo
      endif
c
c set nacta=0 to avoid eigenvalue calculation during the path
c
      nacta=0
C
C  Calculate S and H for molecule
C  nvq is set to 0 to use harmonic oscillator partition function
C
      nvq=0
      do i=1,ntemp
        tk=temp(i)
        akust(i)=0.d0
        akmin(i)=1.d50
        write(6,'(a,f10.3,a)')' Temperature =',tk,' K'
        if(nchn.eq.1)then
          write(6,*)' Enthanlpy and Entropy for Molecule (C)'
          call sandh(tk,wtm,txm,tym,tzm,sigm,qelm,nwm,wm,nirm,
     *             trm,symm,sm(i),hm(i),i)
        else
          write(6,*)' Enthalpy and Entropy for Fragment A'
          call sandh(tk,wtr1,txr1,tyr1,tzr1,sigr1,qelr1,nwr1,wr1,
     *             nirr1,trr1,symr1,sr1,hr1,i)
          write(6,*)' Enthanlpy and Entropy for Fragment B'
          call sandh(tk,wtr2,txr2,tyr2,tzr2,sigr2,qelr2,nwr2,wr2,
     *             nirr2,trr2,symr2,sr2,hr2,i)
          sm(i)=sr1+sr2
          hm(i)=hr1+hr2
        endif
      enddo
      call flush(6)
C
C         CALCULATE DIM ARRAY USED FOR MASS-WEIGHTING
C
      K=0
      DO I=1,NATOMS
         DO J=1,3
            K=K+1
            DIM(K)=1.D0/SQRT(W(I))
         enddo
      enddo
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
      call dvdq
      call energy
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
      REWIND 77
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
      if(natoms.le.20) then
        iprn=6
      else
        iprn=19
        write(iprn,152)nc,v
      endif
      WRITE(iprn,153)
      WRITE(iprn,154)(Q(I),I=1,I3N)
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
      WRITE(iprn,155)
      WRITE(iprn,156)(EIG(I),I=1,I3N)
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
C  Diagonalize the inertia tensor and store it to txc,tyc,tzc
C    Only if nrq=1
C
      if(nrq.eq.1)then
        if(nc.gt.1)then
          txco=txc
          tyco=tyc
          tzco=tzc
          txc=EIG(1)
          tyc=EIG(2)
          tzc=EIG(3)
        else
          txco=EIG(1)
          tyco=EIG(2)
          tzco=EIG(3)
        endif
      endif
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
C    Additional routines are added to incorporate ART routines.
C    The frequencies passed to sandh is the ones from one step before.
C
      if(nc.gt.0) then
        vzo=vz
        ezpo=ezp
      endif
      vz=v-v0
      DO I=1,NTEMP
        if(nc.gt.1)nflago=nflag(i)
        CALL GPATH(G,vz,I,NC)
        if(nc.gt.1.and.nflag(i).ne.nflago)then
          tk=temp(i)
          Write(6,*)' Enthalpy and Entropy for Complex'
          call sandh(tk,wtc,txco,tyco,tzco,sigc,qelc,nwc,wco,nirc,trc,
     *             symc,sc,hc,i)
          if(nflag(i).eq.-1)then
            write(6,*)'  $$$ Minimum rate constant'
          else if(nflag(i).eq.1)then
            write(6,*)'  %%% Maximum rate constant'
          endif
          dels=sc-sm(i)
          delh=hc-hm(i)
          ez0=vzo+(ezpo-ezm)/349.755D0
          ak=ratec(ez0,tk,dels,delh,nchn)
          if(nchn.eq.1)then
            write(6,293)' T, k(T)=',tk,ak,' sec-1'
          else
            write(6,293)' T, k(T)=',tk,ak,' Lt/mol-sec'
          endif
          if(nflag(i).eq.-1)then
            akust(i)=akust(i)+1.d0/ak
            if(ak.lt.akmin(i)) akmin(i)=ak
          else if(nflag(i).eq.1)then
            akust(i)=akust(i)-1.d0/ak
          endif
        endif
      enddo
293   format(a,f8.1,1pe13.5,a/)
C
C         SOLVE FOR REACTION PATH
C
   47 DO I=1,I3N
         RDUM=RAA1*DS
         QM(I)=QM(I)+RDUM*GZ(I)
         BQ(I)=GZ(I)
      ENDDO
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
      IF (NC.GT.NS) then
        write(6,'(///a)')' %&%&%& Final rate constants'
        do i=1,ntemp
          if(nchn.eq.1)then
            write(6,*)' Minimum rate constant'
            write(6,293)'  T, k(T)=',temp(i),akmin(i),' sec-1'
            write(6,*)' Universal statistical theory rate constant'
            write(6,293)'  T, k(T)=',temp(i),1.d0/akust(i),' sec-1'
          else
            write(6,*)' Minimum rate constants'
            write(6,293)'  T, k(T)=',temp(i),akmin(i),' Lt/mol-sec'
            write(6,*)' Universal statistical theory rate constant'
            write(6,293)'  T, k(T)=',temp(i),1.d0/akust(i),' Lt/mol-sec'
          endif
        enddo
        STOP
      endif
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
c   This program calculates the enthalpies and entropies for
c   translational, rotational, vibrational, electronic, and internal
c   rotors, if present.
c   The units are in cal/mol/K for entropy and kcal/mol for enthalpy.
      subroutine sandh(temp,wt,tx,ty,tz,sigma,qelec,nw,w,nir,tr,sym,
     *  s,h,kt)
      implicit DOUBLE PRECISION (a-h,o-z)
      INCLUDE 'SIZES'
      parameter(c0=1.9873d0,c1=0.5*c0,c2=1.5*c0,c3=2.5*c0,c4=7.2836,
     *  c5=0.014837d0,c6=0.041228d0,c7=0.3599d0,c8=1.438769d0)
      common/frqc/wc(nda3),wco(nda3),ezp,ezpo,ezm,nrq,nvq
      DIMENSION W(*),TR(*),SYM(*)
c
14    FORMAT('-TRANSLATIONAL S=',F9.3/' ROTATIONAL S=',F9.3/
     *' VIBRATIONAL S=',F9.3/' INTERNAL ROTATIONAL S=',F9.3/
     *' ELECTRONIC S=',F9.3/' TOTAL S=',F9.3//)
15    FORMAT('-TRANSLATIONAL H=',F9.3/' ROTATIONAL H=',F9.3/
     *' VIBRATIONAL H=',F9.3/' INTERNAL ROTATIONAL H=',F9.3/
     *' TOTAL H=',F9.3/)
c
      STRANS=c3+c2*LOG(WT)+c3*LOG(TEMP)-c4
      HTRANS=c3*TEMP*1.d-3
      SVIB=0.0
      HVIB=0.0
      SIR=0.0
      HIR=0.0
      SELEC=c0*log(QELEC)
      if(nw.le.0) then
        srot=0.d0
        hrot=0.d0
      else
        if(tx.gt.0.d0) then
          IF(TZ.gt.0.d0) then
            QROT=c5*sqrt(TX*TY*TZ*TEMP**3)/SIGMA
          else
            QROT=c6*sqrt(TX*TY)*TEMP/SIGMA
          endif
          SROT=c0*log(QROT)+c2
          HROT=c2*TEMP*1.d-3
        endif
        DO I=1,NW
          U=c8*W(I)/TEMP
          cu=1.d0/(EXP(U)-1.d0)
          stemp=c0*(U*cu-log(1.d0-EXP(-U)))
          SVIB=SVIB+stemp
          htemp=c0*U*TEMP*cu
          HVIB=HVIB+htemp
        enddo
        hvib=hvib*1.d-3
        IF(NIR.gt.0) then
          DO I=1,NIR
            QIR=c7*sqrt(TR(I)*TEMP)/SYM(I)
            SIR=SIR+c0*log(QIR)+c1
            HIR=HIR+c1*TEMP
          enddo
          hir=hir*1.d-3
        endif
      endif
      h=htrans+hrot+hvib+hir
      s=strans+srot+svib+sir+selec
      WRITE(6,14)strans,srot,svib,sir,selec,s
      WRITE(6,15)htrans,hrot,hvib,hir,h
      return
      end
c   This program calculates the rate constant using delta-S and delta-H
c   for the reaction.  The activation energy and temperature are also
c   needed.
c   n is number of reactants, for A+B->C: n=2, for C->A+B: n=1
      DOUBLE PRECISION function ratec(ez0,temp,ds,dh,n)
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter(c0=1.9873d-3,c1=1.d0/1.9873d0)
20    FORMAT(' DS+=',1pE12.5/' DH+=',e12.5,//)
21    FORMAT(' A-FACTOR=',1pE12.5/' ACTIVATION ENERGY=',e12.5/
     *  ' EZERO=',e12.5//)
      if(n.eq.2)then
        ca=0.1263d11
      else
        ca=0.5664d11
      endif
      A=ca*TEMP**n*EXP(DS*c1)
      EA=ez0+DH+dble(n)*c0*TEMP
      ratec=A*EXP(-EA/(c0*TEMP))
      WRITE(6,20)DS,DH
      WRITE(6,21)A,EA,ez0
      return
      end
c
c This program calculates vibrational frequencies and rotational constants
c from equilibrium coordinates read in for molecule (C) or fragments (A and B)
c
      subroutine cfreq(nfrg)
      implicit DOUBLE PRECISION (a-h,o-z)
      include 'SIZES'
      common/frqr/wr1(nda3),wr2(nda3),txm,tym,tzm,txr1,tyr1,tzr1,
     * txr2,tyr2,tzr2,wtm,wtr1,wtr2,nchn,nwm,nwr1,nwr2,nirm,nirr1,nirr2
      COMMON/GPATHB/WM(NDA3),TEMP(NDP),AI1D(5),AAI(2),BBI(2),SYMM(5),
     *SYMA,SYMB,GTEMP(NDP),NFLAG(NDP),N1DR,N2DR
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,LL(NDA),NAM
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/INERT/UXX,UXY,UXZ,UYY,UYZ,UZZ,AIXX,AIXY,AIXZ,AIYY,AIYZ,
     *AIZZ
      COMMON/CHEMAC/WWA(NDA3),CA(NDA3yf,NDA3yf),AI(3),ENMTA,
     *AMPA(NDA3),WWB(NDA3),CB(NDA3yf,NDA3yf),BI(3),ENMTB,
     *AMPB(NDA3),SEREL,S,BMAX,TROTA,TROTB,ANQA(NDA3),ANQB(NDA3),
     *TVIBA,TVIBB,NROTA,NROTB,NOB
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      COMMON/EIGVL/EIG(NDA3yf)
      COMMON/RSTART/HINC,NPTS
      dimension qtemp(nda3),qqtemp(nda3),am(4),QCM(3),VCM(3)
      DATA RHO/1.0D-8/,tfact/2.81837552D05/
      save rho,tfact
c
c Store coordinate
c
      do i=1,i3n
        qtemp(i)=q(i)
        qqtemp(i)=qq(i)
      enddo
c
c fragment A
c
        if(nfrg.eq.0)then
          wt=wtm
          n=natoms
          do j=1,n
            ll(j)=j
          enddo
          k=n*3
          read(5,*)(qza(2,j),j=1,k)
          do j=1,i3n
            q(j)=qza(2,j)
          enddo
          m=6+nirm
          nnm=k-m
          nns=0
          if(nsurf.eq.0)then
            ncen=1
          else
            txm=0.d0
            tym=0.d0
            tzm=0.d0
            ncen=0
          endif
        elseif(nfrg.eq.1)then
      	  if(natoma(2).le.1) return
          n=natoma(2)
          DO J=1,n
            ll(j)=la(2,j)
          ENDDO
          K=3*n
          READ(5,*)(QZA(2,J),J=1,K)
          WT=wtr1
          DO J=1,K
            Q(J)=QZA(2,J)
          ENDDO
          m=k-nwr1
          nnm=nwr1
          nns=0
          ncen=1
        else
          if(natomb(2).le.1) return
          n=natomb(2)
          do j=1,n
            ll(j)=lb(2,j)
          enddo
          k=3*n
          read(5,*)(qzb(2,j),j=1,k)
          wt=wtr2
          i=3*natoma(2)
          do j=1,k
            q(j+i)=qzb(2,j)
          enddo
          m=k-nwr2
          nnm=nwr2
          nns=natoma(2)
          if(nsurf.eq.0)then
            ncen=1
          else
            txr2=0.d0
            tyr2=0.d0
            tzr2=0.d0
            ncen=0
          endif
        endif
c
c moment of inertia
c
        if(ncen.eq.1)then
          CALL CENMAS(WT,QCM,VCM,N)
          CALL ROTN(AM,EROT,N)
          CA(1,1)=AIXX
          CA(2,1)=-AIXY
          CA(2,2)=AIYY
          CA(3,1)=-AIXZ
          CA(3,2)=-AIYZ
          CA(3,3)=AIZZ
          CALL EIGN(CA,CB,3,RHO)
          if(n.eq.2)then
            eig(1)=eig(3)
            eig(2)=eig(3)
            eig(3)=0.d0
          endif
          if(nfrg.eq.0)then
            txm=eig(1)
            tym=eig(2)
            tzm=eig(3)
          elseif(nfrg.eq.1)then
            txr1=EIG(1)
            tyr1=EIG(2)
            tzr1=EIG(3)
          else
            txr2=eig(1)
            tyr2=eig(2)
            tzr2=eig(3)
          endif
        endif
c
c vibrational frequencies
c
        call nmode(n,nns)
        do i=1,nnm
          if(nfrg.eq.0)then
            wm(i)=eig(i+m)
          elseif(nfrg.eq.1)then
            wr1(i)=eig(i+m)
          else
            wr2(i)=eig(i+m)
          endif
        enddo
c
c restore coordinates
c
      do i=1,i3n
        q(i)=qtemp(i)
        qq(i)=qqtemp(i)
      enddo
      return
      end

