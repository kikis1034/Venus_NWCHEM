      SUBROUTINE thermbath
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
      common/thermobath/nthermb,nrscl,nthmid(nda)
c
      svelsq=0.00d0
      n=nthermb
C	
C    calculate velsq.  Recall vnew = vold *sqrt(aheat/velsq)
C
	do i=1,n
	  j3 = 3 * nthmid(i)
	  j2 = j3 - 1
	  j1 = j2 - 1
	  j = nthmid(i)
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
        if (nsel.eq.1 .and. mod(nc,1000).eq.0) then
           write(6,'(1x,a,F12.6)')'initial thermo-bath temp',
     *             tempinits
        endif
c
c	rescale thermo-bath temperature
c
	do i=1,n
	  j3 = 3 * nthmid(i)
	  j2 = j3 - 1
	  j1 = j2 - 1
	  p( j1 ) = p( j1 ) * factor
	  p( j2 ) = p( j2 ) * factor
	  p( j3 ) = p( j3 ) * factor
	end do
C
C     calculate the temperature
C
      svelsq=0.00
      do i=1,n
         j3 = 3 * nthmid(i)
         j2 = j3 - 1
         j1 = j2 - 1                                                                                       
	 j = nthmid(i)
	 svelsq = svelsq+(p(j1)**2+p(j2)**2+p(j3)**2)/w(j)
      end do
C
C  comput the temperature after the rescaling. 
C
      tempinits=svelsq/(3.0 * dble(n) * 0.00198717d0 * c1)
      if (nsel.eq.1 .and. mod(nc,1000).eq.0) then
         write(6,'(1x,a,F12.6)')'thermo-bath temp after RESCALING',
     *        tempinits
      endif
 999  RETURN
      END




