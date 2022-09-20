      SUBROUTINE RADAU(TF,XL,LL,NIP,NITER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C******************************************************************
C         A 15TH ORDER INTEGRATOR THAT USES GAUSS-RADAU SPACINGS.
C         THIS MEANS THAT EACH SEQUENCE IS SPACED INTO 7 SUBSTEPS,
C         THE POSITION AND VELOCITY ARE ACCURATE TO THE 15TH ORDER.
C         THE FIRST SEQUENCE IS ITERNATED 6 TIMES AND LATER SEQUENCES
C         ARE ITERNATED 2 TIMES. 
C         IN OTHER WORDS, IT DOES 7 FORCE CALCULATIONS FOR EACH 
C         ITERATION OF SEQUENCE AND ON THE FIRST SEQUENCE (6 
C         ITERATIONS) THROUGH IT DOES 1+6*7, WHILE ON THE LATER 
C         SEQUENCES (2 ITERATIONS) THROUGH IT DOES 1+2*7.
C         NEW SEQUENCE WILL BE CONTINOUSLY GENERATED AS LONG AS
C         TOTAL INTEGRATION LENGTH IS LESS THAN THE GIVEN TIME 
C         INTERVAL (TF).
C 
C         ADAPTED FROM SUBROUTINE RA15 IN REF:
C             EDGAR EVERHART, 'AN EFFICIENT INTEGRATOR THAT USES
C             GAUSS-RADAU SPACINGS'. IN 'DYNAMICS OF COMETS:
C             THEIR ORIGIN AND EVOLUTION', 185-202,
C             A. CARUSI AND G.B. VALSECCI (EDS.),
C             D. REIDEL PUBLISHING COMPANY 1985.
C         WITH SLIGHT MODIFICATION SEQUENTLY BY N. MARKOVIC(1992), 
C         K. BOLTON(1997) AND G. LI(1999)
C
C         NOTE THAT NAMES OF SOME VARIABLES IN RA15 HAVE BEEN CHANGED 
C         BECAUSE THEY ARE ALREADY DEFINED AS OTHER VARIABLES IN VENUS.
C         THESE VARIABLES AND THEIR CHANGES ARE: NI->NITER, NV->NI, 
C         Q->QT, X->Q, W->WC, WW->WTMP, H->HH AND T->TT.        
C
C         TF IS THE LENGTH OF INTEGRATION (NEGATIVE WHEN INTEGRATING 
C         BACKWARD). LL CONTROLS ACCURACY. THUS IF SS=10.**(-LL) 
C         CONTROLS THE SIZE OF THE LAST TERM IN A SERIES. TRY LL=8 AND 
C         WORK UP OR DOWN FROM THERE. HOWEVER, IF LL.LT.0, THEN XL IS 
C         THE CONSTANT SEQUENCE SIZE USED. A NON-ZERO XL SETS THE SIZE 
C         OF THE FIRST SEQUENCE REGARDLESS OF LL'S SIGN. NIP IS THE 
C         NUMBER OF SEQUENCES BETWEEN INTERMEDIATE PRINTOUT.
C         TYPE OF DIFFERENTIAL EQUATION:
C                 NCLASS=1  :  Y' = F(Y,T)
C                 NCLASS=-2 :  Y''= F(Y,T)
C                 NCLASS=2  :  Y''= F(Y',Y,T)
C************************************************************************
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/INTEGR/ATIME,NI,NID
C
      DIMENSION C(21),D(21),R(21),HH(8),WC(7),U(7),NW(8)
      DIMENSION V(NDA3),F1(NDA3),FJ(NDA3),Y(NDA3),Z(NDA3),
     2          B(7,NDA3),G(7,NDA3),E(7,NDA3),BD(7,NDA3)
      DIMENSION QSAV(NDA3),WR(NDA3)
C
      LOGICAL FIRST,NPQ,NSF,NPER,NCL,NES
      SAVE FIRST,NCL,NPQ,NES,PW,WC,U,C,D,R,WR
      SAVE V,BD,B,E,SS
      DATA FIRST/.TRUE./
      DATA NW/0,0,1,3,6,10,15,21/
      DATA ZERO,HALF,ONE,SR/0.0D0,0.5D0,1.0D0,1.4D0/
      DATA HH/       0.0D0, .05626256053692215D0, .18024069173689236D0,
     *.35262471711316964D0, .54715362633055538D0, .73421017721541053D0,
     *.88532094683909577D0, .97752061356128750D0/
C
 100  FORMAT('  CURRENT SEQUENCE NO. :',I6,
     *       '         TOTAL NO. OF FORCE CALLING:',I6/
     *       '  CURRENT SEQUENCE SIZE:',F12.5,3X,'(fs)'/
     *       '  TOTAL TIME FINISHED  :',F12.5,
     *       '   TOTAL TIME REQURIED:',F12.5,2X,'(fs)')
 110  FORMAT('  NOR,NCOUNT,TT,TP:',2I3,1P2D18.10)
C
C        INITIALIZATION OF FLAGS, CONSTANTS, VELOCIETIES AND SO ON 
C
C             NPER= .TRUE.   ON THE LAST SEQUENCE OF THE INTEGRATION,
C                   .FALSE.  OTHERWISE;
C             NSF = .FALSE.  ON STARTING SEQUENCE, 
C                   .TRUE.   OTHERWISE;
C             NCL = .TRUE.   IF  NCLASS=1,
C                   .FALSE.  IF  NCLASS=-2 OR NCLASS=2. 
C             NPQ = .TRUE.   IF  NCLASS=1  OR  NCLASS=-2,
C                   .FALSE.  IF  NCLASS=2;
C             NES = .TRUE.   IF LL IS NEGATIVE 
C                            THE SEQUENCE SIZE IS THEN FIXED TO BE XL,
C                   .FALSE.  IF LL IS POSITIVE  
C                            VARIABLE TIMESTEP (SEQUENCE SIZE) INTEGRATION.
C
      DIR=ONE
      IF(TF.LT.ZERO) DIR=-ONE
      XL=ABS(XL)*DIR
      NPER=.FALSE.
      NSF=.FALSE.
C             
C             THE FOLLOWING FLAGS AND CONSTANTS ARE THE SAME IN ALL 
C             INTEGRATION STEPS AND ARE CALCULATED ONLY ONCE
C
      IF (FIRST) THEN
         NCLASS=-2
         NCL=NCLASS.EQ.1
         NPQ=NCLASS.LT.2
         NES=LL.LT.0
         PW=1.D0/9.D0
C
C             CONSTANTS IN THE W-,U-,C-,D- AND R-VECTORS
C
         DO N=2,8
            WTMP=N+N*N
            IF(NCL) WTMP=N
            WC(N-1)=ONE/WTMP
            WTMP=N
            U(N-1)=ONE/WTMP
         ENDDO
         W1=HALF
         IF(NCL) W1=ONE
         C(1)=-HH(2)
         D(1)=HH(2)
         R(1)=ONE/(HH(3)-HH(2))
         LA=1
         LC=1
         DO K=3,7
            LB=LA
            LA=LC+1
            LC=NW(K+1)
            C(LA)=-HH(K)*C(LB)
            C(LC)=C(LA-1)-HH(K)
            D(LA)=HH(2)*D(LB)
            D(LC)=-C(LC)
            R(LA)=ONE/(HH(K+1)-HH(2))
            R(LC)=ONE/(HH(K+1)-HH(K))
            IF(K.NE.3) THEN
               DO L=4,K
                  LD=LA+L-3
                  LE=LB+L-4
                  C(LD)=C(LE)-HH(K)*C(LE+1)
                  D(LD)=D(LE)+HH(L-1)*D(LE+1)
                  R(LD)=ONE/(HH(K+1)-HH(L-1))
               ENDDO
            ENDIF
         ENDDO
C
C             SAVE THE RECIPROCAL OF MASS
C
         DO I=1,NI
            KN=1+(I-1)/3
            WR(I)=1.D0/W(KN)
         ENDDO
C
C             RA15 USES POSITION(Q)-VELOCITY(V) VECTOR, TRANSFORM 
C             THE MOMENTA (P) IN VENUS TO V
C 
         DO I=1,NI
            V(I)=P(I)*WR(I)
         ENDDO
C
C             ZERO INITIAL BD-, B- AND E-VALUES. VELOCITY VALUES 
C             ARE ALSO ZEROED IF NCL=.TRUE. (i.e. Y'=F(Y,T))
C 
         DO K=1,NI
            IF(NCL) V(K)=ZERO
            DO L=1,7
               BD(L,K)=ZERO
               B(L,K)=ZERO
               E(L,K)=ZERO
            ENDDO
         ENDDO
C
C             ERROR MAGNITUDE CONTROL (SIZE OF LAST TERM)
C             
         SS=10.**(-LL)
C
      ENDIF
C
      FIRST=.FALSE.
C
C        END OF INITIALIZATION SECTION.
C
C        NEXT ESTIMATE SEQUENCE SIZE (TP, SAME SIGN AS DIR), STARTING
C        WITH INPUT VALUE. NCOUNT IS TO COUNT TOTAL TIMES THAT TP IS 
C        SCALED BY 0.8 IN THE FIRST SEQUENCE
C
      TP=0.1D0*DIR
      IF(XL.NE.ZERO) TP=XL
C     LINE BELOW IN THE ORIGINAL CODE IS REPLACED BY THE LINE FOLLOWING
C     IF (TP/TF.GT.HALF) TP=HALF*TF      
      IF(NES) TP=XL
      NCOUNT=0
C
C        4000 IS THE STARTING PLACE OF THE FIRST SEQUENCE
C        NS: NO. OF SEQUENCES 
C        NF: TOTAL TIMES TO CALL FORCE CALCULATION
C        TM: LENGTH OF INTEGRATION (SUM OVER SEQUENCES)
C
 4000 NS=0
      NF=0
c      NITER=6
      TM=ZERO
      CALL DVDQ
      DO  I=1,NI
         F1(I)=PDOT(I)*WR(I)
      ENDDO
      NF=NF+1
C
C        LINE 722 BEGINS EVERY SEQUENCE AFTER THE FIRST. FIRST FIND
C        NEW G-VALUES FROM THE PREDICATED B-VALUES, EQS.(7) IN REF 
C
  722 DO K=1,NI
        G(1,K)=B(1,K)+D(1)*B(2,K)+
     *    D(2)*B(3,K)+D(4)*B(4,K)+D( 7)*B(5,K)+D(11)*B(6,K)+D(16)*B(7,K)
        G(2,K)=            B(2,K)+
     *    D(3)*B(3,K)+D(5)*B(4,K)+D( 8)*B(5,K)+D(12)*B(6,K)+D(17)*B(7,K)
        G(3,K)=B(3,K)+D(6)*B(4,K)+D( 9)*B(5,K)+D(13)*B(6,K)+D(18)*B(7,K)
        G(4,K)=            B(4,K)+D(10)*B(5,K)+D(14)*B(6,K)+D(19)*B(7,K)
        G(5,K)=                         B(5,K)+D(15)*B(6,K)+D(20)*B(7,K)
        G(6,K)=                                      B(6,K)+D(21)*B(7,K)
        G(7,K)=                                                   B(7,K)
      ENDDO
      TT=TP
      T2=TT*TT
      IF(NCL) T2=TT
      TVAL=ABS(TT)
C
C        LOOP 175 IS 6 ITERATIONS ON FIRST SEQUENCE AND 2 ITERATIONS  
C        THERAFTER, LOOP 174 IS FOR EACH SUBSTEP WITHIN A SEQUENCE
C
      DO M=1,NITER
        DO J=2,8
          JD=J-1
C         JDM=J-2
          S=HH(J)
          QT=S
          IF(NCL) QT=ONE
C
C        POSITION PREDICTOR. EQ.(9) IN REF. EXPRESSION BROKEN INTO TWO
C    
          DO K=1,NI
            A=WC(3)*B(3,K)+S*(WC(4)*B(4,K)+S*(WC(5)*B(5,K)+
     *         S*(WC(6)*B(6,K)+S*WC(7)*B(7,K))))
            Y(K)=Q(K)+QT*(TT*V(K)+T2*S*(F1(K)*W1+S*(WC(1)*B(1,K)+
     *         S*(WC(2)*B(2,K)+S*A))))
            IF(.not.NPQ) then
C
C        VELOCITY PREDICTOR. EQ.(10) IN REF. EXPRESSION BROKEN INTO TWO
C
              A=U(3)*B(3,K)+S*(U(4)*B(4,K)+S*(U(5)*B(5,K)+
     *          S*(U(6)*B(6,K)+S*U(7)*B(7,K))))
              Z(K)=V(K)+S*TT*(F1(K)+S*(U(1)*B(1,K)+S*(U(2)*B(2,K)+S*A)))
            endif
          ENDDO
C
C        FIND FORCES AT EACH SUBSTEP
C
          DO I=1,NI
            QSAV(I)=Q(I)
            Q(I)=Y(I)
          ENDDO
          CALL DVDQ
          DO I=1,NI
            Q(I)=QSAV(I)
            FJ(I)=PDOT(I)*WR(I)
          ENDDO
          NF=NF+1
C
          DO K=1,NI
C
C        FIND G-VALUES FROM THE FORCE FJ FOUND AT THE CURRENT SUBSTEP.
C        EQS.(4) IN REF.
C
            TEMP=G(JD,K)
            GOTO (102,102,103,104,105,106,107,108),J
  102       G(1,K)=GK
            GOTO 160
  103       G(2,K)=(GK-G(1,K))*R(1)
            GOTO 160
  104       G(3,K)=((GK-G(1,K))*R(2)-G(2,K))*R(3)
            GOTO 160
  105       G(4,K)=(((GK-G(1,K))*R(4)-G(2,K))*R(5)-G(3,K))*R(6)
            GOTO 160
  106       G(5,K)=((((GK-G(1,K))*R(7)-G(2,K))*R(8)-G(3,K))*R(9)-
     *           G(4,K))*R(10)
            GOTO 160
  107       G(6,K)=(((((GK-G(1,K))*R(11)-G(2,K))*R(12)-G(3,K))*R(13)-
     *           G(4,K))*R(14)-G(5,K))*R(15)
            GOTO 160
  108       G(7,K)=((((((GK-G(1,K))*R(16)-G(2,K))*R(17)-G(3,K))*R(18)-
     *         G(4,K))*R(19)-G(5,K))*R(20)-G(6,K))*R(21)
C
C        UPGRADE B-VALUES. EQS.(5) IN REF.
C        TEMP IS THE IMPROVEMENT ON G(JD,K) OVER ITS FORMER VALUE.
C
  160       TEMP=G(JD,K)-TEMP
            B(JD,K)=B(JD,K)+TEMP
            GOTO (171,171,203,204,205,206,207,208),J
  203       B(1,K)=B(1,K)+C(1)*TEMP
            GOTO 171
  204       B(1,K)=B(1,K)+C(2)*TEMP
            B(2,K)=B(2,K)+C(3)*TEMP
            GOTO 171
  205       B(1,K)=B(1,K)+C(4)*TEMP
	    B(2,K)=B(2,K)+C(5)*TEMP
	    B(3,K)=B(3,K)+C(6)*TEMP
	    GOTO 171
  206       B(1,K)=B(1,K)+C(7)*TEMP
            B(2,K)=B(2,K)+C(8)*TEMP
            B(3,K)=B(3,K)+C(9)*TEMP
            B(4,K)=B(4,K)+C(10)*TEMP
            GOTO 171
  207       B(1,K)=B(1,K)+C(11)*TEMP
            B(2,K)=B(2,K)+C(12)*TEMP
            B(3,K)=B(3,K)+C(13)*TEMP
            B(4,K)=B(4,K)+C(14)*TEMP
            B(5,K)=B(5,K)+C(15)*TEMP
            GOTO 171
  208       B(1,K)=B(1,K)+C(16)*TEMP
            B(2,K)=B(2,K)+C(17)*TEMP
            B(3,K)=B(3,K)+C(18)*TEMP
            B(4,K)=B(4,K)+C(19)*TEMP
            B(5,K)=B(5,K)+C(20)*TEMP
            B(6,K)=B(6,K)+C(21)*TEMP
  171       CONTINUE
          enddo
        enddo
C
        IF(NES.OR.M.ge.NITER) then
C
C        INTEGRATION OF SEQUENCE IS OVER. 
C
C        NEXT IS SEQUENCE SIZE (TP) CONTROL (ONLY FOR THE FIRST SEQUENCE). 
C        RESTART WITH TP=0.8*TT IF NEW TP IS SMALLER THAN ORIGINAL TT 
C        ON THE FIRST SEQUENCE
C
          HV=ZERO
          DO K=1,NI
            HV=MAX(HV,ABS(B(7,K)))
          ENDDO
          HV=HV*WC(7)/TVAL**7
        endif
      enddo
      IF(.not.NSF) then
        if(nes)then
          TP=XL
        else
          TP=(SS/HV)**PW*DIR
          IF(TP/TT.le.ONE) then
            TP=.8D0*TP
            NCOUNT=NCOUNT+1
            IF (NCOUNT.GT.10) THEN
              WRITE(6,110)NOR,NCOUNT,TT,TP
              write(6,*) 'RADAU: RESTARTED SEQUENCE SIZE (TP) ',
     *             'BECOMES TOO SMALL !'
              STOP
            ENDIF
            GOTO 4000
          endif
        endif
        NSF=.TRUE.
      endif
C
C        POSITION AND VELOCITY CORRECTOR (NEW Q AND V VALUES AT THE END
C        OF SEQUENCE). EQS.(11) AND (12) IN REF.
C        TRANSFORMATION OF VELOCITIES (V) BACK TO MOMENTA (P) IN VENUS
C
      DO K=1,NI
        Q(K)=Q(K)+V(K)*TT+T2*(F1(K)*W1+B(1,K)*WC(1)+B(2,K)*WC(2)+B(3,K)*
     *       WC(3)+B(4,K)*WC(4)+B(5,K)*WC(5)+B(6,K)*WC(6)+B(7,K)*WC(7))
        IF(.not.NCL) then
          V(K)=V(K)+TT*(F1(K)+B(1,K)*U(1)+B(2,K)*U(2)+B(3,K)*U(3)+
     *         B(4,K)*U(4)+B(5,K)*U(5)+B(6,K)*U(6)+B(7,K)*U(7))
          KN=1+(K-1)/3
          P(K)=V(K)*W(KN)
        endif
      ENDDO
C
      TM=TM+TT
      NS=NS+1
C
C        PRINT OUT EVERY NIPth SEQUENCE FOR VARIABLE TIMESTEP INTEGRATION
C
      IF (.NOT.NES) THEN 
         CALL VETEST
         IF (NS/NIP*NIP.EQ.NS) THEN 
            CALL ENERGY
            CALL GWRITE
            WRITE(6,100)NS,NF,TP*10.,TM*10.,TF*10. 
         ENDIF
      ENDIF
C
C        RETURN IF DONE, OTHERWISE CONTINUE WITH NEXT SEQUENCE
C
      IF(NPER) then
C  
        IF (.NOT.NES) THEN
          CALL ENERGY
          CALL GWRITE
          WRITE(6,100)NS,NF,TP*10.,TM*10.,TF*10.
        ENDIF
        RETURN
      endif
C
C        CONTROL ON SIZE OF NEXT SEQUENCE AND ADJUST LAST SEQUENCE TO 
C        EXACTLY COVER THE INTEGRATION SPAN. SET NPER=.TRUE. ON LAST 
C        SEQUENCE.
C
      CALL DVDQ
      DO I=1,NI
         F1(I)=PDOT(I)*WR(I)
      ENDDO
      NF=NF+1
C
      IF(.not.NES) then
        TP=DIR*(SS/HV)**PW
        IF(TP/TT.GT.SR) TP=TT*SR
      endif
      TP=XL
C     THE LINE BELOW IS A BIT DIFFERENT TO ORIGINAL CODE
      IF(DIR*(TF-TM-TP).le.XL) then
        TP=TF-TM
        NPER=.TRUE.
      endif
C
C        NOW PREDICT B-VALUES FOR THE NEXT STEP. EQ.(13) IN REF.
C        VALUES FROM THE PRECEDING SEQUENCE WERE SAVED IN THE E-MATRIX.
C        THE CORRECTION BD IS APPLIED IN LOOP 39 TO REFINE B
C 
      QT=TP/TT
      DO K=1,NI
         IF(NS.NE.1) THEN
           DO J=1,7
             BD(J,K)=B(J,K)-E(J,K)
           ENDDO
         ENDIF
         E(1,K)=    QT*(B(1,K)+ 2.D0*B(2,K)+ 3.D0*B(3,K)+
     *            4.D0*B(4,K)+ 5.D0*B(5,K)+ 6.D0*B(6,K)+ 7.D0*B(7,K))
         E(2,K)=              QT**2*(B(2,K)+ 3.D0*B(3,K)+
     *            6.D0*B(4,K)+10.D0*B(5,K)+15.D0*B(6,K)+21.D0*B(7,K))
         E(3,K)=                           QT**3*(B(3,K)+
     *            4.D0*B(4,K)+10.D0*B(5,K)+20.D0*B(6,K)+35.D0*B(7,K))
         E(4,K)= QT**4*(B(4,K)+ 5.D0*B(5,K)+15.D0*B(6,K)+35.D0*B(7,K))
         E(5,K)=              QT**5*(B(5,K)+ 6.D0*B(6,K)+21.D0*B(7,K))
         E(6,K)=                           QT**6*(B(6,K)+ 7.D0*B(7,K))
         E(7,K)=                                         QT**7*B(7,K)
         DO L=1,7
            B(L,K)=E(L,K)+BD(L,K)
         ENDDO
      ENDDO
C
C        TWO ITERATIONS FOR EVERY SEQUENCE AFTER FIRST SEQUENCE
C
      NITER=2
      GOTO 722
C
      RETURN
      END
