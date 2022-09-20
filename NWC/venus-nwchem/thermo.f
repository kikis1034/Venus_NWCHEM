      SUBROUTINE thermo(ipr)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      COMMON/TRANSB/TRANS,NREL
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,L(NDA),NAM
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,3*NDA),QZB(NDP,3*NDA),NATOMA(NDP),NATOMB(NDP)
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      common/vrscal/nsel,nscale,nequal,thermotemp,nrgd
c
      svelsq=0.00d0
C
C    make sure that the coordinates are not written during minization
C
      N=NATOMB(1)-nrgd
C	
C    calculate velsq.  Recall vnew = vold *sqrt(aheat/velsq)
C
	do i=1,n
	  j3 = 3 * lb(1,i)
	  j2 = j3 - 1
	  j1 = j2 - 1
          j = lb(1,i)
          svelsq =svelsq+(p(j1)**2+p(j2)**2+p(j3)**2)/w(j)
	end do
C
C   reference:  "MOLECULAR DYNAMICS SIMULATION" by JIM HAILE  p.458
C   k = 1.38066 * 10(-23) J/K = 0.00198624 kcal/mol K
C   factemp is aheat above
C
	factemp =3.0d0*dble(n)*0.00198717d0*thermotemp*c1
	factor = sqrt( factemp / svelsq )
        tempinits=svelsq/(3.0 * dble(n) * 0.00198717d0 * c1)
        if (mod(ipr,100).eq.0.and.nsel.eq.1) then
           write(6,'(1x,a,F12.6)')'system temp',
     *             tempinits
        endif
	do i=1,n
	  j3 = 3 * lb(1,i)
	  j2 = j3 - 1
	  j1 = j2 - 1
	  p( j1 ) = p( j1 ) * factor
	  p( j2 ) = p( j2 ) * factor
	  p( j3 ) = p( j3 ) * factor
	end do
C
C    make sure that the total linear momentum equals zero
C
	sumx = 0.00
	sumy = 0.00
	sumz = 0.00
	do i=1,n
	  j3 = 3 * i + 3*natoma(1)
	  j2 = j3 - 1
	  j1 = j2 - 1
	  sumx = sumx + p( j1 )
	  sumy = sumy + p( j2 )
	  sumz = sumz + p( j3 )
	end do
	sumx = sumx / real( n )
	sumy = sumy / real( n )
	sumz = sumz / real( n )
	do i=1,n
	  j3 = 3 * i + 3*natoma(1)
	  j2 = j3 - 1
	  j1 = j2 - 1
	  p(j1) = p(j1) - sumx
	  p(j2) = p(j2) - sumy
	  p(j3) = p(j3) - sumz
	end do
C
C     calculate the temperature
C
      svelsq=0.00
      do i=1,n
         j3 = 3 * lb(1,i)
         j2 = j3 - 1
         j1 = j2 - 1                                                                                       
         j = lb(1,i)
	 svelsq = svelsq+(p(j1)**2+p(j2)**2+p(j3)**2)/w(j)
      end do
C
C  comput the temperature after the rescaling. 
C
      tempinits=svelsq/(3.0 * dble(n) * 0.00198717d0 * c1)
      if (mod(ipr,100).eq.0.and.nsel.eq.1) then
         write(6,'(1x,a,F12.6)')'system temp after RESCALING',
     *        tempinits
      endif
 999  RETURN
      END
