      SUBROUTINE EBOND(H,T,R,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CALCULATE LOCAL MODE (MORSE OSCILLATOR) ENERGIES 
C
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/MORSEB/RMZ(ND02),B(ND02),D(ND02),N2J(ND02),N2K(ND02),
     *CM1(ND02),CM2(ND02),CM3(ND02),CM4(ND02),CM5(ND02)
C
      DE=D(I)
      BETA=B(I)
      R0=RMZ(I)
      WA=W(N2J(I))
      WB=W(N2K(I))
      WT=WA+WB
      WJ=WB/WT
      WI=WA/WT
      UMASS=WA*WB/WT
      IZ=3*N2J(I)
      JZ=3*N2K(I)
      IY=IZ-1
      IX=IY-1
      JY=JZ-1
      JX=JY-1
      R=SQRT((Q(IX)-Q(JX))**2 + (Q(IY)-Q(JY))**2 + (Q(IZ)-Q(JZ))**2)
      T=((WJ*P(IX)-WI*P(JX))**2 + (WJ*P(IY)-WI*P(JY))**2 +
     *                            (WJ*P(IZ)-WI*P(JZ))**2)/2.0D0/UMASS
      V=EXP(-BETA*(R-R0))
      V=DE*(V**2-2*V+1.0D0)
      T=T/C1
      V=V/C1
      H= V + T
      RETURN
      END
