      SUBROUTINE GPATH(g,VZ,K,NCC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CALCULATE FREE ENERGY ALONG THE REACTION PATH
C
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,L(NDA),NAM
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/EIGVL/EIG(NDA3yf)
      COMMON/GPATHB/WM(NDA3),TEMP(NDP),AI1D(5),AAI(2),BBI(2),SYMM(5),
     *SYMA,SYMB,GTEMP(NDP),NFLAG(NDP),N1DR,N2DR
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      common/frqc/wc(nda3),wco(nda3),ezp,ezpo,ezm,nrq,nvq
      DIMENSION QTEM(nda3)
C
   50 FORMAT(5X,'A MAXIMUM IN FREE ENERGY,  TEMPERATURE=',F8.2,
     *',  FREE ENERGY=',1PE12.5,',  VQ(KCAL/MOL)=',0PF8.3,'*****')
   51 FORMAT(5X,'A MINIMUM IN FREE ENERGY,  TEMPERATURE=',F8.2,
     *',  FREE ENERGY=',1PE12.5,',  VQ(KCAL/MOL)=',0PF8.3,'*****')
   52 FORMAT(32X,'TEMPERATURE=',F8.2,',  FREE ENERGY=',1PE12.5,
     *',  VQ(KCAL/MOL)=',0PF8.3)
C
      DUM=-1.98721587169D-3*TEMP(K)
      G=0.0D0
C
C         CALCULATE FREE ENERGY
C
C             VIBRATIONAL FREE ENERGY.  IT IS ASSUMED THAT THE LOWEST
C             FREQUENCIES ARE FOR THE 1- AND 2-DIMENSIONAL INT. ROTORS.
C
      ezp=0.d0
      J=8+N1DR+2*N2DR
      ratio=1.d0
      ii=0
      DO I=J,I3N
        ii=ii+1
        if(k.eq.1)then
          if(ncc.gt.0)wco(ii)=wc(ii)
          wc(ii)=eig(i)
        endif
        QPF=1.0D0/(1.0D0-EXP(-1.438769D0*EIG(I)/TEMP(K)))
        ezp=ezp+eig(i)
        G=G+DUM*LOG(QPF)
      ENDDO
      ezp=ezp*0.5d0
C
C             ROTATIONAL FREE ENERGY
C             Calculate only if nrq=1 (molecule, not surface)
C
      if(nrq.eq.1)then
        QPF=0.014837D0*SQRT(EIG(4)*EIG(5)*EIG(6))*TEMP(K)*SQRT(TEMP(K))
        G=G+DUM*LOG(QPF)
      endif
C
C             INTERNAL ROTATION FREE ENERGY(1-DIMENSIONAL)
C
      IF (N1DR.NE.0) THEN
         DO I=1,N1DR
            QPF=0.35990D0*SQRT(AI1D(I)*TEMP(K))/SYMM(I)
            G=G+DUM*LOG(QPF)
         ENDDO
      ENDIF
C
C             INTERNAL ROTATION FREE ENERGY(2-DIMENSIONAL)
C             N2DR IS THE NUMBER OF 2-D ROTORS(N2DR = 0, 1, OR 2).
C
      IF (N2DR.NE.0) THEN
C
C             PUT CURRENT COORDINATES(Q) IN TEMPORARY STORAGE.
C
         DO I=1,I3N
            QTEM(I)=Q(I)
         ENDDO
C
C             FOR A POLYATOMIC-ATOMIC REACTANT PAIR WITH ONE 2-D ROTOR,
C             FRAGMENT A MUST BE THE POLYATOMIC.
C
         N=NATOMA(1)
         WT=WTA(1)
         DO I=1,N
            L(I)=LA(1,I)
         ENDDO
         CALL GINROT(AAI,SYMA,VZ,WT,QPF,N,K)
         G=G+DUM*LOG(QPF)
C
C                  RESET Q FOR FRAGMENT A
C
         NAA=3*N
         DO I=1,NAA
            Q(I)=QTEM(I)
         ENDDO
C
         IF (N2DR.NE.1) THEN
            N=NATOMB(1)
            WT=WTB(1)
            DO I=1,N
               L(I)=LB(1,I)
            ENDDO
            CALL GINROT(BBI,SYMB,VZ,WT,QPF,N,K)
            G=G+DUM*LOG(QPF)
C
C                  RESET Q FOR FRAGMENT B
C
            II=NAA+1
            DO I=II,I3N
               Q(I)=QTEM(I)
            ENDDO
C
C                  RESET L ARRAY
C
            DO I=1,NATOMS
               L(I)=I
            ENDDO
         ENDIF
      ENDIF
C
C             TOTAL FREE ENERGY
C
      VQ=VZ+(EZP-EZM)*cm2cal
      G=G+VQ
C
C         FREE ENERGY TESTS
C
      IF (NCC.NE.0) THEN
         IF (NCC.LE.1) THEN
            IF (G.GT.GTEMP(K)) NFLAG(K)=1
            IF (G.LT.GTEMP(K)) NFLAG(K)=-1
            GOTO 8
         ENDIF
         IF (NFLAG(K).GE.1) THEN
            IF (G.GT.GTEMP(K)) GOTO 8
            NFLAG(K)=-1
            WRITE(6,50)TEMP(K),G,VQ
            GOTO 9
         ENDIF
         IF (G.GE.GTEMP(K)) THEN
            NFLAG(K)=1
            WRITE(6,51)TEMP(K),G,VQ
            GOTO 9
         ENDIF
      ENDIF
    8 WRITE(6,52)TEMP(K),G,VQ
    9 GTEMP(K)=G
      RETURN
      END
