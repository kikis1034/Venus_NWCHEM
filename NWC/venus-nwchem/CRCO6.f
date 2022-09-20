      SUBROUTINE CRCO6      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DO I=1,6
         CALL MORSECRCO6(I)
      ENDDO
      CALL DDMORSECRCO6
      
      RETURN
      END
C
C     SUBROUTINE FROM OSAMMA MORSECRCO6.f
C
      SUBROUTINE MORSECRCO6(NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C      PARAMETER(ND1=2100)
      COMMON/QPDOT/Q(3*NDA),PDOT(3*NDA)
      COMMON/PQDOT/P(3*NDA),QDOT(3*NDA),W(NDA)
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/MORSEB/RMZ(ND02),B(ND02),D(ND02),N2J(ND02),N2K(ND02),
     *CM1(ND02),CM2(ND02),CM3(ND02),CM4(ND02),CM5(ND02)
      COMMON/CRCO6B/a1,a2,a3,a4,a5,a6,a7,a8,dcrco6,rtmp(3*NDA),
     *dfixed(100),deltad
      COMMON/FORCES/N,I3N,NS,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,
     *NTET,NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI       

      sum1j=0.0D0
      sum2j=0.0D0
      sum3j=0.0D0
      sum=0.0D0      

      DO i=1,6

	 nj=n2j(i)
	 nk=n2k(i)
	 j3 = 3 * nj
         j2 = j3 - 1
	 j1 = j2 - 1
	 k3 = 3 * nk
	 k2 = k3 - 1
	 k1 = k2 - 1

	 d1 = q(k1) - q(j1)
	 d2 = q(k2) - q(j2)
	 d3 = q(k3) - q(j3)
	
	 rtmp(i) = sqrt(d1*d1+d2*d2+d3*d3)

	 sum = sum + (1.0D0 / rtmp(i)) 

C
C        compute dX/dq for j
C

	 rr = rtmp(i)*rtmp(i)
	 rrr = rr * rtmp(i)
	 sum1j = sum1j + (d1 / rrr)
	 sum2j = sum2j + (d2 / rrr)
	 sum3j = sum3j + (d3 / rrr)

      ENDDO

C
C     CALCULATE INDICES FOR COORDINATES
C

      J3=3*N2J(NL)
      J2=J3-1
      J1=J2-1
      K3=3*N2K(NL)
      K2=K3-1
      K1=K2-1

C
C     CALCULATE INDEX FOR R
C

      JK=(N2J(NL)-1)*(2*N-N2J(NL))/2+N2K(NL)-N2J(NL)

C
C     CALCULATE RELATIVE COORDINATES AND R
C

      T1=Q(K1)-Q(J1)
      T2=Q(K2)-Q(J2)
      T3=Q(K3)-Q(J3)

      R(JK)=DSQRT(T1*T1+T2*T2+T3*T3)

C
C     CALCULATE (DV/DQ)'S
C
C     compute dx/dr for k atom
C

      rs = r(jk)*r(jk)
      rcubed = rs * r(jk)
      sum1k = (1.0D0 * t1) / rcubed
      sum2k = (1.0D0 * t2) / rcubed
      sum3k = (1.0D0 * t3) / rcubed

C
C     compute D. Note use sqrt(c1) instead of c1 to convert to 
C     venus units of energy.
C

      DIFF = sum - a6
      CC = DIFF * DIFF 
      CCC = CC * DIFF
      CCCC=((CCC-a7)/a8)*((CCC-a7)/a8)
      d_1=a1+sqrt(c1)*dtanh(a2*sum-a3)
      d_2=a4-a5*dexp(-1.0D0*CCCC)
      dcrco6 = d_1 * d_2

C     write(6,*)'from morsecrco6',sum,dcrco6/c1
C
C     compute dD/dR
C

      d1=6.0D0*a5*(CCC-a7)*CC
      d2=exp(-CCCC)
      d4=d1*d2
      d5=a8*a8
      dd_2=d4/d5
	
      dd_1=a2*sqrt(c1)*(1.0D0/dcosh(a2*sum-a3))**2
      d3=d_1*dd_2+d_2*dd_1

C
C     compute second term of dvdq  D * dVmorse/dq
C

      DUM0=DEXP(-B(NL)*(R(JK)-RMZ(NL))) 
      DUM4=d3*(1.0D0-DUM0)*(1.0D0-DUM0)
      DUML2=0.0D0

      DO i=1,6
	 DUML1=DEXP(-B(NL)*(rtmp(i)-RMZ(NL)))
	 DUML2=DUML2+(1.0D0-DUML1)*(1.0D0-DUML1)
      ENDDO

      DUM3=d3*DUML2

      DUM1=2.0D0*B(NL)*dcrco6*DUM0*(1.0D0-DUM0)/R(JK)

      DUM2=T1*DUM1
      T1J= (-1.0D0*DUM2) + DUM4*sum1j
      T1K= DUM2 - DUM3*sum1k
      PDOT(K1)=PDOT(K1) + T1K
      PDOT(J1)=PDOT(J1) + T1J
      
      DUM2=T2*DUM1
      T2J= (-1.0D0*DUM2) + DUM4*sum2j
      T2K= DUM2 - DUM3*sum2k
      PDOT(K2)=PDOT(K2) + T2K
      PDOT(J2)=PDOT(J2) + T2J
         
      DUM2=T3*DUM1
      T3J= (-1.0D0*DUM2) + DUM4*sum3j
      T3K= DUM2 - DUM3*sum3k
      PDOT(K3)=PDOT(K3) + T3K
      PDOT(J3)=PDOT(J3) + T3J
      
      RETURN
      END
      
C
C     SUBROUTINE FROM OSAMMA DDMORSECRCO6.f
C
      SUBROUTINE DDMORSECRCO6
C      PARAMETER(ND1=2100)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON/QPDOT/Q(3*NDA),PDOT(3*NDA)
      COMMON/PQDOT/P(3*NDA),QDOT(3*NDA),W(NDA)
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/MORSEB/RMZ(ND02),B(ND02),D(ND02),N2J(ND02),N2K(ND02),
     *CM1(ND02),CM2(ND02),CM3(ND02),CM4(ND02),CM5(ND02)
      COMMON/CRCO6B/a1,a2,a3,a4,a5,a6,a7,a8,dcrco6,rtmp(3*NDA),
     *dfixed(100),deltad
      COMMON/FORCES/N,I3N,NS,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,
     *NTET,NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI  
	dimension d1(100),d2(100),d3(100)
	dimension sum1k(100),sum2k(100),sum3k(100)
	dimension indx(100)
	dimension fd(100),dfd(100),swe(100),dswe(100)
	dimension dum(100,100)
		
      sum1j=0.0D0
      sum2j=0.0D0
      sum3j=0.0D0
      sum=0.0D0

      DO i=1,6

	 nj=n2j(i)
	 nk=n2k(i)
	 j3 = 3 * nj
         j2 = j3 - 1
	 j1 = j2 - 1
	 k3 = 3 * nk
	 k2 = k3 - 1
	 k1 = k2 - 1

	 d1(i) = q(k1) - q(j1)
	 d2(i) = q(k2) - q(j2)
	 d3(i) = q(k3) - q(j3)
	
	 d1s = d1(i) * d1(i)
	 d2s = d2(i) * d2(i)
	 d3s = d3(i) * d3(i)

	 rtmp(i) = sqrt(d1s+d2s+d3s)
	 indx(i) = i

	 sum = sum + (1.0D0 / rtmp(i)) 

         rr = rtmp(i)*rtmp(i)
         rrr = rr * rtmp(i)

         sum1j = sum1j + (d1(i) / rrr)
         sum2j = sum2j + (d2(i) / rrr)
         sum3j = sum3j + (d3(i) / rrr)

	 sum1k(i) = -d1(i) / rrr
	 sum2k(i) = -d2(i) / rrr
	 sum3k(i) = -d3(i) / rrr

      ENDDO

C
C	order the bonds in increasing order
C

	call indexx(6,rtmp,indx)

C
C       order bonds in decreasing order
C

	DO i=1,6
	  indx(6+i) = indx(6-i+1)
	ENDDO

	do i=1,6
	  indx(i) = indx(6+i)
	end do	

C
C calculate the derivatives
C

C
C       define some variables first
C

	xx=(sum-a6)*(sum-a6)	
	xxx=xx*(sum-a6)
	yy=((xxx-a7)/a8)**2
	yyy=a2*sum-a3

C
C  calculate the derivatives
C

	do i=1,6

	drj=rtmp(indx(i))-rmz(indx(i))
	
	dum1=a1+sqrt(c1)*dtanh(yyy)
	dum2=a4-a5*dexp(-yy)
	dum12=dum1*dum2
	fd(i)=dfixed(i)-dum12

	ddum1=a2*sqrt(c1)*(1.0D0/dcosh(yyy))**2
	ddum_1=6.0D0*a5/a8**2
	ddum_2=xx*(xxx-a7)*dexp(-yy)
	ddum2=ddum_1*ddum_2


	dfd(i)=-1.0d0*(dum1*ddum2+dum2*ddum1)

C
C  note multiply DDOMEGA by dsum/dq.
C

	if(drj.gt.0.00) then

	 dum1=1.0D0-dexp(-b(indx(i))*drj)
	 swe(i)=dum1*dum1
	 dswe(i)=2.0D0*b(indx(i))*(1.0D0-dum1)*dum1

	else
	  swe(i)=0.0D0
	  dswe(i)=0.0D0
	end if

	end do

C
C  calculate deltad
C

	deltad=0.0D0

	do i=1,6
	dummy=fd(i)

	 do j=1,i
	  dummy=dummy*swe(j)
	 end do

	deltad=deltad+dummy

	end do

C
C  calculate the derivatives with respect to atom 1
C

      dum(1,1)=dfd(1)*swe(1)  
      dum(1,2)=fd(1)*dswe(1)  

      dum(2,1)=dfd(2)*swe(1)*swe(2) 
      dum(2,2)=fd(2)*dswe(1)*swe(2) 
      dum(2,3)=fd(2)*swe(1)*dswe(2) 

      dum(3,1)=dfd(3)*swe(1)*swe(2)*swe(3) 
      dum(3,2)=fd(3)*dswe(1)*swe(2)*swe(3) 
      dum(3,3)=fd(3)*swe(1)*dswe(2)*swe(3)  
      dum(3,4)=fd(3)*swe(1)*swe(2)*dswe(3) 

      dum(4,1)=dfd(4)*swe(1)*swe(2)*swe(3)*swe(4)
      dum(4,2)=fd(4)*dswe(1)*swe(2)*swe(3)*swe(4) 
      dum(4,3)=fd(4)*swe(1)*dswe(2)*swe(3)*swe(4) 
      dum(4,4)=fd(4)*swe(1)*swe(2)*dswe(3)*swe(4) 
      dum(4,5)=fd(4)*swe(1)*swe(2)*swe(3)*dswe(4) 

      dum(5,1)=dfd(5)*swe(1)*swe(2)*swe(3)*swe(4)*swe(5) 
      dum(5,2)=fd(5)*dswe(1)*swe(2)*swe(3)*swe(4)*swe(5)
      dum(5,3)=fd(5)*swe(1)*dswe(2)*swe(3)*swe(4)*swe(5) 
      dum(5,4)=fd(5)*swe(1)*swe(2)*dswe(3)*swe(4)*swe(5) 
      dum(5,5)=fd(5)*swe(1)*swe(2)*swe(3)*dswe(4)*swe(5) 
      dum(5,6)=fd(5)*swe(1)*swe(2)*swe(3)*swe(4)*dswe(5)

      dum(6,1)=dfd(6)*swe(1)*swe(2)*swe(3)*swe(4)*swe(5)*swe(6) 
      dum(6,2)=fd(6)*dswe(1)*swe(2)*swe(3)*swe(4)*swe(5)*swe(6) 
      dum(6,3)=fd(6)*swe(1)*dswe(2)*swe(3)*swe(4)*swe(5)*swe(6)
      dum(6,4)=fd(6)*swe(1)*swe(2)*dswe(3)*swe(4)*swe(5)*swe(6) 
      dum(6,5)=fd(6)*swe(1)*swe(2)*swe(3)*dswe(4)*swe(5)*swe(6)
      dum(6,6)=fd(6)*swe(1)*swe(2)*swe(3)*swe(4)*dswe(5)*swe(6)  
      dum(6,7)=fd(6)*swe(1)*swe(2)*swe(3)*swe(4)*swe(5)*dswe(6) 
	
	ninc=0
	dumx=0.0D0
	dumy=0.0D0
	dumz=0.0D0

	do i=1,6

          j3 = 3 * n2j(indx(i))
          j2 = j3 - 1
          j1 = j2 - 1

         dumx = dum(i,1) * sum1j
         dumy = dum(i,1) * sum2j
         dumz = dum(i,1) * sum3j

	    pdot(j1) = pdot(j1) + dumx 
	    pdot(j2) = pdot(j2) + dumy
	    pdot(j3) = pdot(j3) + dumz 


	 do j=1,6

c
c add derivative of first term of ddeltad to atom i
c
          j3 = 3 * n2j(indx(j))
          j2 = j3 - 1
          j1 = j2 - 1

          k3 = 3 * n2k(indx(j))
          k2 = k3 - 1
          k1 = k2 - 1

	dum1 = d1(indx(j))/rtmp(indx(j))
	dum2 = d2(indx(j))/rtmp(indx(j))
	dum3 = d3(indx(j))/rtmp(indx(j))

	if(j.le.i) then
C
C  derivatives with respect to atom i
C


 	 dumx = -dum(i,j+1) * dum1
 	 dumy = -dum(i,j+1) * dum2
 	 dumz = -dum(i,j+1) * dum3

	    pdot(j1) = pdot(j1) + dumx 
	    pdot(j2) = pdot(j2) + dumy
	    pdot(j3) = pdot(j3) + dumz 


C
C  derivatives with respect to atom k
C

	     dumx = dum(i,1) * sum1k(indx(j))
	     dumy = dum(i,1) * sum2k(indx(j))
	     dumz = dum(i,1) * sum3k(indx(j))

         dumx = dumx + dum(i,j+1) * dum1 
         dumy = dumy + dum(i,j+1) * dum2
         dumz = dumz + dum(i,j+1) * dum3

	     pdot(k1) = pdot(k1) + dumx 
	     pdot(k2) = pdot(k2) + dumy 
 	     pdot(k3) = pdot(k3) + dumz 

	else

	 dumx = dum(i,1) * sum1k(indx(j))
	 dumy = dum(i,1) * sum2k(indx(j))
	 dumz = dum(i,1) * sum3k(indx(j))

	     pdot(k1) = pdot(k1) + dumx 
	     pdot(k2) = pdot(k2) + dumy 
 	     pdot(k3) = pdot(k3) + dumz 

	end if	

        end do
	dumx=0.0D0
	dumy=0.0D0
	dumz=0.0D0
      end do
      
      RETURN
      END

C
C     SUBROUTINE FROM OSAMMA INDEXX.f
C

      SUBROUTINE indexx(n,arr,indx)
      IMPLICIT NONE
      INTEGER n,indx(n),M,NSTACK
      double precision arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      double precision a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
