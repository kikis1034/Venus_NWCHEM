      SUBROUTINE LENJ(INL,LNL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CALCULATE LENNARD-JONES POTENTIAL ENERGY DERIVATIVES
C
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/LENJB/ALJ(ND05),BLJ(ND05),CLJ(ND05),N5J(ND05),N5K(ND05),
     *NREP(ND05),MREP(ND05),LREP(ND05)
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      DIMENSION JKA(ND05),RNA(ND05),RMB(ND05),RLC(ND05)
      LOGICAL FIRST,DIFFN,DIFFM,DIFFL
      DATA FIRST,DIFFN,DIFFM,DIFFL/.TRUE.,.FALSE.,.FALSE.,
     *          .FALSE./
      save
C
      IF (FIRST) THEN
         DO NL=INL,LNL
            JKA(NL)=ISHFT((N5J(NL)-1)*(2*NATOMS-N5J(NL)),-1)+N5K(NL)
     *             -N5J(NL)
            RNA(NL)=-NREP(NL)*ALJ(NL)
            RMB(NL)=-MREP(NL)*BLJ(NL)
            RLC(NL)=-LREP(NL)*CLJ(NL)
         ENDDO
C
         NNEXP=NREP(INL)
         NTEXP=NNEXP+2
         DO NL=INL+1,LNL
            IF (NREP(NL).NE.NNEXP) THEN
               DIFFN=.TRUE.
               GOTO 8005
            ENDIF
         ENDDO
 8005    CONTINUE
         MEXP=MREP(INL)
         MTEXP=MEXP+2
         DO NL=INL+1,LNL
            IF (MREP(NL).NE.MEXP) THEN
               DIFFM=.TRUE.
               GOTO 8015
            ENDIF
         ENDDO
 8015    CONTINUE
         LEXP=LREP(INL)
         LTEXP=LEXP+2
         DO NL=INL+1,LNL
           IF (LREP(NL).NE.LEXP) THEN
              DIFFL=.TRUE.
              GOTO 8025
           ENDIF
         ENDDO
 8025    CONTINUE
         FIRST=.FALSE.
      ENDIF
C
C         CODE FOR GENERAL LENNARD-JONES
C
      IF (DIFFN.OR.DIFFM.OR.DIFFL) THEN
         DO NL=INL,LNL
            J3=3*N5J(NL)
            J2=J3-1
            J1=J2-1
            K3=3*N5K(NL)
            K2=K3-1
            K1=K2-1
            JK=JKA(NL)
            T1=Q(K1)-Q(J1)
            T2=Q(K2)-Q(J2)
            T3=Q(K3)-Q(J3)
            R(JK)=SQRT(T1*T1+T2*T2+T3*T3)
            RRJK=1.0/R(JK)
            DUM1=RNA(NL)*RRJK**(2+NREP(NL))
            DUM1=DUM1+RMB(NL)*RRJK**(MREP(NL)+2)
            DUM1=DUM1+RLC(NL)*RRJK**(LREP(NL)+2)
            TDUM1=DUM1*T1
            TDUM2=DUM1*T2
            TDUM3=DUM1*T3
            PDOT(K1)=PDOT(K1)+TDUM1
            PDOT(K2)=PDOT(K2)+TDUM2
            PDOT(K3)=PDOT(K3)+TDUM3
            PDOT(J1)=PDOT(J1)-TDUM1
            PDOT(J2)=PDOT(J2)-TDUM2
            PDOT(J3)=PDOT(J3)-TDUM3
         ENDDO
C
C         CODE FOR 12-6-1 LENNARD-JONES
C
      ELSE IF (NNEXP.EQ.12.AND.MEXP.EQ.6.AND.LEXP.EQ.1) THEN
         DO NL=INL,LNL
            J3=3*N5J(NL)
            J2=J3-1
            J1=J2-1
            K3=3*N5K(NL)
            K2=K3-1
            K1=K2-1
            JK=JKA(NL)
            T1=Q(K1)-Q(J1)
            T2=Q(K2)-Q(J2)
            T3=Q(K3)-Q(J3)
            R(JK)=SQRT(T1*T1+T2*T2+T3*T3)
            RRJK=1.0/R(JK)
            DUM1=RNA(NL)*RRJK**8*RRJK**6 
            DUM1=DUM1+RMB(NL)*RRJK**8
            DUM1=DUM1+RLC(NL)*RRJK**3
            TDUM1=DUM1*T1
            TDUM2=DUM1*T2
            TDUM3=DUM1*T3
            PDOT(K1)=PDOT(K1)+TDUM1
            PDOT(K2)=PDOT(K2)+TDUM2
            PDOT(K3)=PDOT(K3)+TDUM3
            PDOT(J1)=PDOT(J1)-TDUM1
            PDOT(J2)=PDOT(J2)-TDUM2
            PDOT(J3)=PDOT(J3)-TDUM3
         ENDDO
C
C         CODE FOR 12-6-0 LENNARD-JONES
C
      ELSE IF (NNEXP.EQ.12.AND.MEXP.EQ.6.AND.LEXP.EQ.0) THEN
         DO NL=INL,LNL
            J3=3*N5J(NL)
            J2=J3-1
            J1=J2-1
            K3=3*N5K(NL)
            K2=K3-1
            K1=K2-1
            JK=JKA(NL)
            T1=Q(K1)-Q(J1)
            T2=Q(K2)-Q(J2)
            T3=Q(K3)-Q(J3)
            R(JK)=SQRT(T1*T1+T2*T2+T3*T3)
            RRJK=1.0/R(JK)
            DUM1=RNA(NL)*RRJK**8*RRJK**6 
            DUM1=DUM1+RMB(NL)*RRJK**8
            TDUM1=DUM1*T1
            TDUM2=DUM1*T2
            TDUM3=DUM1*T3
            PDOT(K1)=PDOT(K1)+TDUM1
            PDOT(K2)=PDOT(K2)+TDUM2
            PDOT(K3)=PDOT(K3)+TDUM3
            PDOT(J1)=PDOT(J1)-TDUM1
            PDOT(J2)=PDOT(J2)-TDUM2
            PDOT(J3)=PDOT(J3)-TDUM3
         ENDDO
      ELSE
C
C         CODE FOR LENNARD-JONES WITH OTHER CONSTANT NREP, MREP AND LREP
C
         DO NL=INL,LNL
            J3=3*N5J(NL)
            J2=J3-1
            J1=J2-1
            K3=3*N5K(NL)
            K2=K3-1
            K1=K2-1
            JK=JKA(NL)
            T1=Q(K1)-Q(J1)
            T2=Q(K2)-Q(J2)
            T3=Q(K3)-Q(J3)
            R(JK)=SQRT(T1*T1+T2*T2+T3*T3)
            RRJK=1.0/R(JK)
            DUM1=RNA(NL)*RRJK**NTEXP
            DUM1=DUM1+RMB(NL)*RRJK**MTEXP
            DUM1=DUM1+RLC(NL)*RRJK**LTEXP
            TDUM1=DUM1*T1
            TDUM2=DUM1*T2
            TDUM3=DUM1*T3
            PDOT(K1)=PDOT(K1)+TDUM1
            PDOT(K2)=PDOT(K2)+TDUM2
            PDOT(K3)=PDOT(K3)+TDUM3
            PDOT(J1)=PDOT(J1)-TDUM1
            PDOT(J2)=PDOT(J2)-TDUM2
            PDOT(J3)=PDOT(J3)-TDUM3
         ENDDO
      ENDIF
      RETURN
      END
