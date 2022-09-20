************************************************************************
       subroutine hessint(lll)
c   Program: MD simulation using predictor-corrector algorithm
c     lll=1: using rotated coordinate projection method correction
c     lll=2: using redeveloped rotated coordinate correction (gaussian way)
c     lll=3: using general caretesian coordinate correction
c   Prediction is made from q1 to q2. true potential is calculated
c   at r2. then correction is done from q1 to q2 using a 5th order
c   interpolated potential along the prediction step and using a
c   2nd order potential along the other directions.  the prediction
c   step is valid locally around q1 called trust region.
c
       implicit double precision(a-h,o-z)
      INCLUDE 'SIZES'
      parameter(cc1=4.184d-4,cc2=1.d0/cc1,cc3=0.04184d0)
      COMMON/FORCES/NAT,nat3,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
c   common blocks
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/INTEGR/ATIME,NI,NID
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      COMMON/GAUSS/ICHRG,IMULTP,IAN(NDAyf),VGAUSS
      common/ghessb/fa(NDA3,NDA3),trad,rmin,dt,nstep,
     *  nhesup,nhessf,nip
      dimension vin(nat3),vfin(nat3),g1(nat3),g2(nat3),q1(nat3),
     *  q2(nat3),q3(nat3),h1(nat3,nat3),h2(nat3,nat3),fmod1(nat3),
     *  fmod(nat3),fmod3(nat3),gg1(nat3),fmodfin(nat3)
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
c
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

230      format(2x,f8.4,3x,3(f15.10,3x),2f8.5)

c   nhessf is the flag to control FMTRX diagonalization
c     nhessf=1: not to diagonalize
c     nhessf=0: normal diagonalization for frequencies
c
      nhessf=1
      ihc = 1
c   trust region parameters....
       rmax = 2.d0 - rmin
       dscale = 0.5d0
       tol = rmin
       trad1 = trad
c  threshold for ratio and p-c-distance
       thrs1 = 0.5d0
       thrs2 = 0.002d0
c   weights for rho function
       ww1 = 1.d0
c
c   Initialize some parameters
       dtmax = 0.d0
       avener  = 0.d0
       ntstep = 0
       irf = 0
       tmstart=0.d0
       tt = 0.d0
       dr = 0.0d0
       nv = 0
       ic = 1
       avdt = 0.d0
c
c   Convert dt (time interval) to fs.
       dt = dt*10.d0
c
       write(43,*)'Time step within Pred. and Corr. steps(fs)=',dt
       write(43,*)'Mass of the atoms ',(w(i),i=1,nat)
c
c   energy, gradient and hessian at q1(initial coordinates)
c
       call dvdq
       call energy
       call QMCALC(1)
c
c   Set up potentials, coordinates, momentum, gradients and 
c   Hessians for the integration
       v1=v
       vp = v1
       do i = 1, nat3
          q1(i) = q(i)
          q3(i) = q(i)
          g1(i)=-pdot(i)/0.04184d0
          vin(i)=0.1d0*p(i)/w((i+2)/3)
          do j=1,nat3
            h1(i,j)=fa(i,j)/0.04184d0
          enddo
       enddo
c   Get model potential, gradient, pdot
        call ener_model(nat3,q1,q,v1,g1,h1,emod)       
        call dvdq_model(1,nat3,q1,q,g1,h1,fmod)
        call dvdq_model(0,nat3,q1,q,g1,h1,fmod3)
c
       tke = 0.d0
       do i = 1, nat3
          fmod1(i) = fmod(i)
          gg1(i) = g1(i)
              tke = tke + 0.01d0*p(i)*p(i)/2.d0/w((i+2)/3)
       enddo
         write(23,230)tt,tke*cc2,v1,(tke*cc2+v1),dr,dt*nv
c
c     Start the integration
150    continue
        write(43,*)
        write(43,*)'------NSTEP------',nc
c   prediction step
c
          write(43,*)'Begin Prediction '
          call predict(nat3,trad,dt,fmod1,q3,vin,gg1,h1,w,q2,dr,nv)   
c
          write(43,*)'End of prediction  '
          write(43,*)'Integration steps for prediction  ',nv
c   Calculate energy, gradient and Hessian at q2.
c   Hessians at q2 are either exact or updated values
c
          do i=1,nat3
            q(i)=q2(i)
          enddo
          call dvdq
          do i = 1, nat3
            g2(i)=-pdot(i)/0.04184d0
          enddo
          call energy
          v2=v
          write(43,*)'Ab initio potential after prediction ',v2
          call ener_model(nat3,q3,q2,vp,gg1,h1,emod)
          write(43,*)'Model potential after prediction  ',emod
          delE = abs(emod - v2)
          erat = (v1 - v2)/(v1 - emod)
          erat1 = delE/abs(v1 - emod)
          write(43,*)'Energy error ratio, delE  :',erat,delE
c  check for the error : if error is in the opposite directions,
c  redo the prediction with a smaller trust radius
c
          if ( erat .le. 0.d0 ) then
              write(43,*)'eratio in check ',erat,irf
               trad = max(trad*dscale,0.05d0)
             if ( irf .eq. 0 ) then
                irf = 1
                write(43,*)'Redo the step:irf ',irf
                goto 150
             else
           write(43,*)'Redo did not work after 1st interation:cont.'
                irf = 0
                trad = 0.01d0
             endif
          endif
c
c    Calculate Hessian
          write(43,*)'Calculate Hessian '
          write(43,*)'ihc,nhesup = ',ihc,nhesup
          call hess2(nat3,q1,q2,g1,g2,h1,h2,ihc)
c
c   corrector step
c
          write(43,*)
          write(43,*)'Begin correction '
          if (lll.eq.1) then
c
c   lll=1 : Rotated coordinate corrector
          call correct1(nat3,fmod1,q1,v1,g1,h1,q2,v2,g2,h2,w,
     \ vin,q3,vfin,fmodfin,vp,tt,dt,nv,nc,ic,pcdist,q3p)
c
c   lll=2 : Rotated coordinate corrector - Gaussian way
          elseif (lll.eq.2) then
          call correct2(nat3,fmod1,q1,v1,g1,h1,q2,v2,g2,h2,w,
     \ vin,q3,vfin,fmodfin,vp,tt,dt,nv,nc,ic,pcdist,q3p)
c 
c   lll=3 : General Cartesian coodinate corrector
          else
          call correct3(nat3,fmod1,q1,v1,g1,h1,q2,v2,g2,h2,w,
     \ vin,q3,vfin,fmodfin,vp,tt,dt,nv,nc,ic,pcdist,q3p)
c
          endif       
c
           write(43,*)'End of Correction '
           write(43,*)'TR step: count(nc) ',nc
c   Get out of the integrator, if the nc > nsteps
c
           if ( nc.gt.nstep ) return
c
c   update trust radius using gaussian function of rho
c
          rho = (1.d0-ww1)*erat1/thrs1 + ww1*q3p/thrs2
c
          if ( rho .gt. 0.d0 .and. rho .lt. 2.d0 ) then
            fac = 0.5d0 + 1.d0*exp(-(rho*rho)/1.442695d0)
            trad = fac*trad1
          write(43,*)'Error, Factor, TR',rho,fac,trad
          else
            trad = trad1*dscale
          write(43,*)'Error, Factor, TR ',rho,dscale,trad
          endif
c
c   write erat,erat1, delE,rho, pcdist, trad
          write(43,*)'erat,erat1, delE,rho, pcdist, trad'
          write(43,110)erat,erat1,(v2-emod),rho,pcdist,q3p,trad
 110      format(2x,7f12.8)
c
c   update coordinates, energy, gradient, forces and hessian
          tke = 0.d0
          do i = 1, nat3
              q1(i) = q2(i)
              vin(i) = vfin(i)
              p(i)=vin(i)*w((i+2)/3)
              g1(i) = g2(i)
              fmod1(i) = fmodfin(i)
              gg1(i) = -fmodfin(i)*cc2
              tke = tke + p(i)*p(i)/(2.d0*w((i+2)/3))
          enddo
c
              v1 = v2
              tv = vp
c
          do i = 1, nat3
            do j = 1, nat3
                h1(i,j) = h2(i,j)
            enddo
          enddo
c
c   for testing the reaction criteria, call radout here
c
         tm=tt*0.1d0
         write(43,*)'  Trajectory time =',tt,' fs'
         write(23,230)tt,tke*cc2,tv,(tke*cc2+tv),dr,dt*nv
          avener = avener + (tke*cc2+tv)
          tstep = dt*nv
          ntstep = ntstep + 1
          if (dtmax .le. tstep ) dtmax=tstep
          IF ((NTEST.EQ.2.AND.(NAST.EQ.0.OR.NAST.EQ.1)).OR.
     \ (NTEST.EQ.1.AND.NAST.EQ.0)) then
          write(43,*)'Out of corrector'
        write(43,*)'Average ener ',avener/ntstep
        write(43,*)'dtmax = ',dtmax
        write(43,*)'No. of t.r. step = ',ntstep
          return
          endif
c
        if ( nc.le.nstep) goto 150
c
        write(43,*)'Average ener ',avener/ntstep
        write(43,*)'dtmax = ',dtmax
        write(43,*)'No. of TR step = ',ntstep
       nhessf=0
       return
2010       end
c
c************************************************************************
       subroutine matmult(m1,n1,mm,nn,a,b,c,iflag)
       implicit real*8(a-h,o-z)
       dimension a(m1,n1),b(mm,nn),c(m1,nn)
c      this subroutine does matrix multiplication
c      iflag = 0/1 - for matrix order i.e. a or a^t
c      iflag = 0
c      c(i,j) = sum(k) a(i,k)*b(k,j)
c      iflag = 1
c      c(i,j) = sum(k) a(k,i)*b(k,j)
c
       if (iflag.eq.0) then
         do i=1,m1
           do j=1,nn
             s2=0.0d0
             do k=1,mm
             s2=s2+a(i,k)*b(k,j)
             enddo
             c(i,j)=s2
           enddo
         enddo
       else
         do i=1,n1
           do j=1,nn
             s2=0.0d0
             do k=1,mm
             s2=s2+a(i,k)*b(j,k)
             enddo
             c(i,j)=s2
           enddo
         enddo
       endif
       return
       end

************************************************************************
c
       subroutine matvec(iflag,m1,n1,a,x,c)
       implicit real*8(a-h,o-z)
c     this subroutine does matrix-vector multiplication
c      c = A*x , where A is (m,n) matrix and x(n) is the vector
c
       dimension a(m1,n1), x(n1), c(n1)
c
      if (iflag .eq. 0) then
         do i = 1,m1
           s1 = 0.d0
           do j = 1,n1
              s1 = s1 + a(i,j)*x(j)
           enddo
              c(i) = s1
         enddo
      else
         do i = 1,n1
           s1 = 0.d0
           do j = 1,m1
              s1 = s1 + a(j,i)*x(j)
           enddo
              c(i) = s1
         enddo
      endif
c
      return
      end
c
************************************************************************
       subroutine ener_model(nat3,q0,q,v0,dvdq,d2vdq,emod)
       implicit real*8(a-h,o-z)
c    this subroutine finds the model potential at a given set of 
c    coordinates using quadratic model
c
       dimension q0(*),q(*),dvdq(*),d2vdq(nat3,*),dq(nat3),temp(nat3)
c    get del_q
c
       do i = 1,nat3
         dq(i) = (q(i)-q0(i))
       enddo
c    get dvdq*del_q
c
         first = 0.d0
       do i = 1,nat3
         first = first + dvdq(i)*dq(i)
       enddo
c    get del_q*d2vdq*del_q
c
       call matvec(0,nat3,nat3,d2vdq,dq,temp)
       second = 0.d0
       do i = 1,nat3
         second = second + dq(i)*temp(i)
       enddo
c   get model potential
c
       emod = v0 + first + second*0.5d0
c
      return
      end
************************************************************************
       subroutine dvdq_model(iflag,nat3,q0,q,dvdq,d2vdq,dvdq_mod)
       implicit real*8(a-h,o-z)
c     this subroutine returns dVdq for the quadratic model potential
c     iflag = 0   only dVdq
c     iflag = 1   force = -dVdq for integration of trajectory
c
       dimension q0(*),q(*),dvdq(nat3),d2vdq(nat3,nat3),temp(nat3),
     \ dq(nat3),dvdq_mod(nat3)
       cc1 = 4.184d-4
c     get del_q
       do i = 1,nat3
         dq(i) = (q(i)-q0(i))
       enddo
c
c     get d2vdq*del_q
       call matvec(0,nat3,nat3,d2vdq,dq,temp)
       do i = 1,nat3
         if (iflag.eq.0) then
           dvdq_mod(i) = dvdq(i)+temp(i)
         else
           dvdq_mod(i) = - (dvdq(i)+temp(i))*cc1
         endif
       enddo
c 
       return
       end
c
************************************************************************
c
      subroutine predict(nat3,trad,dt0,fmod1,q0,vin,dvdq,d2vdq,
     \  am,q,dr,nv)
      implicit real*8(a-h,o-z)
c     integrate EOM using velocity-verlet algorithm
c
      dimension q(nat3),q0(nat3),fmod(nat3),v0(nat3),dvdq(nat3),am(*)
      dimension v(nat3),qtemp(nat3),vin(nat3),dq(nat3),fmodt(nat3),
     \  fmod1(nat3)
c
        cc1 = 4.184d-4
        qfac = 1.d-3
        itmax = 100000
        dt = dt0
        dth = dt*0.5d0
        nv = 0
c
           write(43,*)'Coord. and Vel. while entering Pred.'
         do i = 1,nat3
           q(i) = q0(i)
           qtemp(i) = q(i)
           fmod(i) = fmod1(i)
           v(i) = vin(i)
           write(43,*)'   qt..',q0(i),v(i)
         enddo
c
      do 20 k = 1, itmax
c  velocity-verlet
          do i = 1,nat3
            kk = (i+2)/3
            fmodt(i) = fmod(i)
            v0(i) = v(i)
            q(i) = q(i) + dt*(v0(i) + dth*fmodt(i)/am(kk))
          enddo
c
            call dvdq_model(1,nat3,q0,q,dvdq,d2vdq,fmod)
c
          do i = 1,nat3
            kk = (i+2)/3
            v(i) = v0(i) + dth*(fmodt(i) + fmod(i))/am(kk)
          enddo
c  end of verlet ..
c
         nv = nv + 1
c  check conditions  : the end of trust radius
           drsq = 0.d0       
         do i = 1,nat3
           dq(i) = abs(q(i) - q0(i))
           drsq = drsq + dq(i)*dq(i)
         enddo
         dr = dsqrt(drsq)
         rratio = 100.d0*abs(dr - trad)/trad
         if ( dr .gt. trad .or. rratio .le. 1.d0 ) then
         write(43,*)'rratio .., adjusted trust radius', rratio,dr
         return
         endif
c
 20   continue
      return
      end
c
************************************************************************
       subroutine correct1(nat3,fmod1,q1,v1,g1,h1,q2,v2,g2,h2,am,vin,
     \ q3,v,fmodfin,vp,t,dt,nv,nc,ic,pcdist,q3p)
c    this subroutine corrects the predictor step in multidimensions(3N).
c    the system is integrated from q to q3 using the quadratic potential.
c    the coordinates(q3) after integration are projected to a hyperplanes
c    that contain q1 and q2 and are perpendicular to the vector connecting
c    q1 and q2.  the projected points q4 and q5 are then fitted to a
c    5th order polynomial using gauss-elimination by using V, dV and d2V
c    at q4 and q5. then the corrected potential is obtained at q3 using
c    the fitted potential.  the projection of q3 on the hyperplanes is
c    obtained by gauss-elimination.
c
       implicit real*8(a-h,o-z)
      INCLUDE 'SIZES'
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      common/ghessb/fa(NDA3,NDA3),trad,rmin,ddt,nstep,
     *  nhesup,nhessf,nip

c
      dimension q1(nat3),q11(nat3),q2(nat3),q3(nat3),u(nat3,nat3),
     \ q33(nat3),q22(nat3),g11(nat3),g22(nat3),h11(nat3,nat3),
     \ h22(nat3,nat3),fmod1(nat3),am(*),nj(3),nk(3),qq3(nat3)
      dimension v0(nat3),v(nat3),vin(nat3),fmodt(nat3),fmod(nat3),
     \ g1(nat3),g2(nat3),h1(nat3,*),h2(nat3,*),q4(nat3),q5(nat3),
     \ fmodfin(nat3)
c
c
       cc1 = 4.184d-4
       dth = dt*0.5d0
       ical = 0
       dtol = 0.01d0
c
       write(43,*)'   CORRECTOR STEP : Rotated coordinate method '
c
       do i = 1, nat3
          v(i) = vin(i)
          fmodfin(i)=fmod1(i)
       enddo
c
c   translate q1 and q2 and get unitary matrix for rotation
       write(43,*)'Do translation and rotation to obtain unitary matrix'
          call utrans(nat3,q1,q2,q11,q22,u)
c
          call rotate1(nat3,g1,g2,h1,h2,u,g11,g22,h11,h22,1)
c
c   velocity-verlet
      write(43,*)'nv inside corrector', nv
c
       do j = 1, nv
c
         if ( j.gt.1) vpm1 = vp
c
      write(43,*)'-------  nv = ',j
c   integrate from q1 to q3 by dt
          do i = 1,nat3
            kk = (i+2)/3
            fmodt(i) = fmodfin(i)
            v0(i) = v(i)
c  remember old coord.
            qq3(i) = q3(i)
            q3(i) = q3(i) + dt*(v0(i) + dth*fmodt(i)/am(kk))
          enddo
c
c    translate q3 to new coord. and rotate q1, q2, q3 from (q1,q2..q3n)
c    to (r,0,0.....0) coordinate system
          do i = 1, nat3
            q3(i) = q3(i)-q1(i)
          enddo
          call matvec(0,nat3,nat3,u,q3,q33)
c
c   get max(q_|_) and mag(q_|_)
        if ( j .eq. nv) then
          q3m = abs(q33(2))
          q3p2 = 0.d0
          do i =2,nat3
            q3m = max(q3m,abs(q33(i)))
            q3p2 = q3p2 + q33(i)*q33(i)
          enddo
            q3p = sqrt(q3p2)
          write(44,440) q3m,q3p,q33(1)
 440    format(2x,3f8.5)
            write(43,*) 'Perpendicular distance between pred. and 
     * corr. ',q3p
        endif
c
c
c   project q3 on hyperplanes to get q4 and q5
          write(43,*)'Projecting q3 to q4 and q5 '
          call project(nat3,q22,q33,q4,q5)
c
c   polynomial fit between q4 and q5 and gradient at q3 
          call polyfit1(nat3,q11,v1,g11,h11,q22,v2,g22,h22,q33,q4,q5,
     \ vp,fmod,issflag,ical)
c
       if ( issflag .eq. 1 ) then
           write(43,*)'issflag  ',issflag
           write(43,*)'coming our corrector'
           do i = 1, nat3
              q3(i) = qq3(i)
              fmodfin(i) = fmodt(i)
           enddo
           vp = vpm1
c  distance between pred. and corr.
             pcdist = 0.d0
             do i = 1 ,nat3
               pcdist = pcdist + (q2(i)-q3(i))*(q2(i)-q3(i))
             enddo
               pcdist = sqrt(pcdist)
            write(43,*) 'distance between pred. and corr. ',pcdist
c calculate perpendicular part
          q3p2 = 0.d0
          do i =2,nat3
            q3p2 = q3p2 + q33(i)*q33(i)
          enddo
            q3p = sqrt(q3p2)
            write(43,*) 'Perpendicular distance between pred. and 
     * corr. ',q3p
c
         return
       endif
c  
c   rotate coordinates and gradients back to cartesian coordinates for q3
          call rotate2(nat3,u,q1,q33,q3,fmod,fmodfin)
c
       do i = 1,nat3
         fmodfin(i) = - fmodfin(i)*cc1
       enddo
c
          t = t + dt
c
          do i = 1,nat3
            kk = (i+2)/3
            v(i) = v0(i) + dth*(fmodt(i) + fmodfin(i))/am(kk)
c          write(43,*)'v0(i),fmodt(i),fmod(i),v(i) ',i,v0(i),fmodt(i),
c     \ fmod(i),v(i)
          enddo
c
c     check the distance between the predicted and the corrected step
          if ( j .eq. nv ) then
             pcdist = 0.d0
             do i = 1 ,nat3
               pcdist = pcdist + (q2(i)-q3(i))*(q2(i)-q3(i))
             enddo
               pcdist = sqrt(pcdist)
             write(43,*) 'distance between pred. and corr. ',pcdist
c calculate perpendicular part
          q3p2 = 0.d0
          do i =2,nat3
            q3p2 = q3p2 + q33(i)*q33(i)
          enddo
            q3p = sqrt(q3p2)
            write(43,*) 'Perpendicular distance between pred. and 
     * corr. ',q3p
c
          endif
            
c     print results
          call hessout(nat3,dt,vp,v,am,q3,ic)
          if ( nc.gt.nstep) return
          IF ((NTEST.EQ.2.AND.(NAST.EQ.0.OR.NAST.EQ.1)).OR.
     \ (NTEST.EQ.1.AND.NAST.EQ.0)) then
          write(*,*)'going out of corr.'
          return
          endif
c
       enddo
c  end of verlet ..
c
       return
       end
c
************************************************************************
c  
      subroutine polyfit1(nat3,q11,v1,g11,h11,q22,v2,g22,h22,q33,q4,q5,
     \ vp,delv5,issflag,ical)
      implicit real*8(a-h,o-z)
c   this subroutine converts from q3,q4,q5 variables to 's' variable
c   such that r = s(q5-q4)-q4
c   for s=0 ; r = q4   and s=1 ; r = q5
c   then it fits the variable s to a 5th order polynomial
c
      dimension q11(nat3),q22(nat3),q4(nat3),q5(nat3),g11(nat3),
     \ g22(nat3),h11(nat3,*),h22(nat3,*),dvdq_mod(nat3),
     \ d2vdq(nat3,nat3),temp(nat3),delv5(nat3),q33(nat3),d(6),
     \ ain(nat3,nat3)
c
       s0 = 0.d0
       s1 = 1.d0
c    get potentials at q4 and q5
c
       write(43,*)'Fifth-order polyonmial fitting'
       call ener_model(nat3,q11,q4,v1,g11,h11,emod_q4)
       call ener_model(nat3,q22,q5,v2,g22,h22,emod_q5)
       write(43,*)'emodq4, emodq5  :',emod_q4,emod_q5
c
c   get the dVm and d2Vm with 's' variable at q4 and q5
c    dVm/ds = dVm/dq * (q5-q4) or dVm/dq * (q2-q1)
c
       call dvdq_model(0,nat3,q11,q4,g11,h11,dvdq_mod)
c
         dvds_q4 =  dvdq_mod(1)*q22(1)
c
       call dvdq_model(0,nat3,q22,q5,g22,h22,dvdq_mod)
c
         dvds_q5 = dvdq_mod(1)*q22(1)
       write(43,*)'dvds-q4, dvds-q5  : ',dvds_q4,dvds_q5
c 
c    d2Vm/ds = (q5-q4) * d2Vm/dq * (q5-q4)        
c    d2Vm/dq = d2V/dq ; again (q5-q4) = (q2-q1)
c    (q2-q1) is actually q2 as q1 = 0
          d2vds_q4 = q22(1)*h11(1,1)*q22(1)
          d2vds_q5 = q22(1)*h22(1,1)*q22(1)
c
       write(43,*)'d2vds-q4, d2vds-q5  : ',d2vds_q4,d2vds_q5
c
c    get coefficients of the polynomial: Ad = b
c
       call coeff(emod_q4,emod_q5,dvds_q4,dvds_q5,d2vds_q4,d2vds_q5,d)
c
c    get the inverse of matrix A : A^(-1) 
c
       call ainverse(6,ain)  
c  
c    calculate potential and gradient at q3 as functions of q_|_ and s
c
        s = (q33(1)-q11(1))/abs((q22(1)-q11(1)))
       rr = abs((q22(1)-q11(1)))
c
c   get max of x _|_
c
       write(43,*)'ical  ',ical
       write(43,*)'s ...: ',s
      if ( s.ge.1.d0 .and. ical.eq.0) then
        qp = 0.d0
       do i = 2, nat3
          qp = qp + q33(i)*q33(i)
       enddo
       sthres = 0.5d0*dsqrt(qp)/rr
       ical = 1
       write(43,*)'qp, rs :',qp,rr
       write(43,*)'qp/rs : ',sthres
      endif
c
         write(43,*)'1+sthres :',1.d0+sthres,sthres
       if ( s .gt. (1.d0 + sthres)  ) then
         issflag = 1
         write(43,*)'s :',s,' s > 1:coming out of poly'
         sthres = 0.d0
         return
       else
         issflag = 0
       endif

       if ( s .lt. 0.d0 ) then
         call ener_model(nat3,q11,q33,v1,g11,h11,vp)
         call dvdq_model(0,nat3,q11,q33,g11,h11,delv5)
         write(43,*)'s,vp  2nd order :', s,vp
       elseif ( s .gt. 1.d0 ) then
         call ener_model(nat3,q22,q33,v2,g22,h22,vp)
         call dvdq_model(0,nat3,q22,q33,g22,h22,delv5)
         write(43,*)'s,vp  2nd order :', s,vp
       else
       call vpoly(s,d,vp)
       write(43,*)'s,vp 5th order  :', s,vp
       call gradient(nat3,d,ain,s,q22,q33,q4,g11,g22,h11,h22,delv5)
       endif
c
       return 
       end
c
************************************************************************
c 
      subroutine vpoly(s,ai,vp)
      implicit real*8(a-h,o-z)
c     this subroutine returns the polynomial fitted potential at q3
c     as a function of 's'. the coefficients for the fifth order fit
c     were obtained from 'subroutine coeff'
c
      dimension ai(6)
c
c
       vp = ai(1) + ai(2)*s + ai(3)*s**2 + ai(4)*s**3 + ai(5)*s**4
     \ + ai(6)*s**5
c
       return
       end
c
************************************************************************
c
      subroutine gradient(n,d,b,s,q22,q33,q4,g11,g22,h11,h22,grad)
      implicit real*8(a-h,o-z)
c     this subroutine gives the gradient of the polynomial fitted
c     surface. the gradient has 2 terms, 1st term is the derivative 
c     w.r.t. 's' and the 2nd term is the derivative w.r.t. q_|_
c      b : A_inv
c      d : coeff of the polynomial 
c
      dimension b(6,6),g11(n),g22(n),h11(n,n),h22(n,n),grad(n),
     \ dedq(6,n),q22(n),q33(n),q4(n),dq(n),first(n),second(n),
     \ temp(n),d(6)
c
       cc1 = 4.184d-4
c      
c     1st term :  (dvds)*(q2)/(abs(q2))**2 ; since q2-q1 = q2
c
c   dvds : Sum of ( k*d(k+1)*s**(k-1) ) , k = 1-5
c       s = q3(1)
       dvds = 0.d0
       do k = 1, 5
          dvds = dvds + k*d(k+1)*s**(k-1)
       enddo
c
         first(1) = dvds/q22(1)
C
       do i = 2, n
         first(i) = 0.d0
       enddo
c
c    2nd term :  Sum over k,j ( b(k,j)*(de(k)dq)*s(k-1) ), k,j = 1-6
c    dedq(k,n)
       do i = 1, n
          dq(i) = q4(i)
       enddo
c    k =1
       call matvec(0,n,n,h11,dq,temp)
        do i = 1, n
          dedq(1,i) = g11(i) + temp(i)
        enddo
c    k = 2
c       call matvec(0,n,n,h11,q22,temp)
        do i = 1, n
          dedq(2,i) = h11(1,i)*q22(1)
        enddo
c    k = 3, 6
        do i = 1, n
          dedq(3,i) = 0.d0
          dedq(6,i) = 0.d0
        enddo
c    k = 4
       call matvec(0,n,n,h22,dq,temp)
        do i = 1, n
          dedq(4,i) = g22(i) + temp(i)
        enddo
c    k = 5
c       call matvec(0,n,n,h22,q22,temp)
        do i = 1, n
          dedq(5,i) = h22(1,i)*q22(1)
        enddo
c   
c   A_inv * dedq(k,n) * s**(k-1) : n components
       do i = 1, n
         sum3 = 0.d0
         do k = 1, 6
           sum4 = 0.d0
           do j = 1, 6
               sum4 = sum4 + b(k,j)*dedq(j,i)
           enddo
           sum3 = sum3 + sum4*s**(k-1)
         enddo
         second(i) = sum3
       enddo
         second(1) = 0.d0
c 
       do i = 1, n
          grad(i) = (first(i) + second(i))
       enddo
c
       return
       end
c    
************************************************************************
c
      subroutine coeff(e1,e2,gr1,gr2,he1,he2,a)
      implicit real*8(a-h,o-z)
c     this subroutine returns the coefficients of the 5th order polynomial
c     in the variable 's'.
c
      dimension a(6)
c
      a(1) = e1
      a(2) = gr1
      a(3) = he1/2.d0
c
      tp1 = a(2)+ a(3)
      p1 = e2 - ( tp1 + a(1))
      tp2 = gr2 - 2.d0*a(3)
      p2 = tp2 - a(2)
      p3 = he2 - a(3) - a(3) 
c
      a(4) = 10.d0*p1 - 4.d0*p2 + 0.5d0*p3
      a(5) = -15.d0*p1 + 7.d0*p2 - p3
      a(6) =   6.d0*p1 - 3.d0*p2 + 0.5d0*p3
c
      return
      end
c                
************************************************************************
c     
       subroutine ainverse(m,ain)
       implicit real*8(a-h,o-z)
c     this subroutine returns the inverse matrix elements of A
c
       dimension ain(m,m)
c
       do i = 1, m
         do j = 1, m
            ain(i,j) = 0.d0
         enddo
       enddo
c
       ain(1,1) =   1.d0
       ain(2,2) =   1.d0
       ain(3,3) =   1.d0
       ain(4,1) = -10.d0
       ain(5,1) =  15.d0
       ain(6,1) = - 6.d0
       ain(4,2) = - 6.d0
       ain(5,2) =   8.d0
       ain(6,2) = - 3.d0
       ain(4,3) = - 3.d0
       ain(5,3) =   3.d0
       ain(6,3) = - 1.d0
       ain(4,4) =  10.d0
       ain(5,4) = -15.d0
       ain(6,4) =   6.d0
       ain(4,5) = - 4.d0
       ain(5,5) =   7.d0
       ain(6,5) = - 3.d0
       ain(4,6) =   0.5d0
       ain(5,6) = - 1.d0
       ain(6,6) =   0.5d0
c
       return
       end
c
************************************************************************
c 
      subroutine project(n,q2,q3,q4,q5)
      implicit real*8(a-h,o-z)
c    this subroutine gets the projection of q3 on the hyperplanes to get
c    q4 and q5.  since the coordinates are translated and rotated,
c    q4 has the same coordinates as q3 with q4(1) projection as 0. q5 has the
c    same coordinates as q3 with q5(1) projection as q2(1).
c
      dimension q2(*),q3(*),q4(*),q5(*),b4(n),b5(n)
c
c     get q4
c       write(43,*)'q3(i)  ',(q3(i),i=1,n)
           q4(1) = 0.d0
        do i = 2, n
           q4(i) = q3(i)
        enddo
c     get q5
           q5(1) = q2(1)
         do i = 2, n
           q5(i) = q3(i)
         enddo
c    
      return
      end
c
************************************************************************
c   
      subroutine utrans(n,q1,q2,q11,q22,u)
c     this subroutine translates the coordinates the cartesian origin
c     to q1 and then generates the U^t matrix needed for rotation
c
      implicit real*8(a-h,o-z)
      dimension q1(n),q2(n),u(n,n),q11(n),q22(n)
c
c
c     translate the coordinate origin to (q1,q2....q3n)
c
         do i = 1,n
           q2(i) = q2(i)-q1(i)
           q11(i) = 0.d0
         enddo
c
c     get U^t matrix for rotation
c     note: our U matrix is already a transposed matrix!!!!
         call umat(n,q2,u)
         call matvec(0,n,n,u,q2,q22)
         do i = 1,n
           q2(i) = q2(i)+q1(i)
         enddo
c
       return
       end
c
************************************************************************
c
      subroutine rotate1(n,g1,g2,h1,h2,u,grad1,grad2,hess1,hess2,igrad)
c     this subroutine rotates the coordinate system from 
c     (q1,q2,q3,....q3n) to (r,0,0,,,,,0), such that the direction of
c     the predictor vector is parallel to r and other components are  0.
c     this is done by a rotation matrix U. this rotation matrix is used
c     the transform the gradients and the hessians to the new coordinates
c        igrad = 0,  transform gradient only
c        igrad = 1,  transform gradient and hessian 
c     the coordinates were translated before rotation such that origin 
c     lies on q1.
c
      implicit real*8(a-h,o-z)
      dimension q1(n),q2(n),u(n,n),g1(n),h1(n,n),g2(n),h2(n,n)
      dimension grad1(n),grad2(n),hess1(n,n),hess2(n,n),temp(n,n)
c
c     transform from (q1,q2...q3n) to (r,0,...0) coord. system
c     rotate gradients only
         if ( igrad .eq. 0) then
            call matvec(0,n,n,u,g1,grad1)
            call matvec(0,n,n,u,g2,grad2)
         else
c    rotate both gradient and hessians
            call matvec(0,n,n,u,g1,grad1)
            call matvec(0,n,n,u,g2,grad2)
            call matmult(n,n,n,n,h1,u,temp,1)
            call matmult(n,n,n,n,u,temp,hess1,0)
            call matmult(n,n,n,n,h2,u,temp,1)
            call matmult(n,n,n,n,u,temp,hess2,0)
         endif
c
      return 
      end
c 
************************************************************************
c
       subroutine rotate2(n,u,q1,q,qq,g,gg)
c      this subroutine rotates the coordinates and the gradients from 
c      (r,0,..0) to (q1,q1....q3n) coord. system
c
       implicit real*8(a-h,o-z)
       dimension u(n,n),g(n),q1(n),q(n),qq(n),gg(n)
c
c    tranform from (r,0,0...0) to (q1,q2,.....q3n)
         call matvec(1,n,n,u,q,qq)
         do i = 1,n
           qq(i) = qq(i)+q1(i)
         enddo
c    gradient transformation
       call matvec(1,n,n,u,g,gg)
c
       return
       end
c
************************************************************************
c      
       subroutine umat(n,q,u)
c     this subroutine sets up the unitary matrix (U) for the transformation 
c
       implicit real*8(a-h,o-z)
       dimension q(n),r(n-1),u(n,n)
c
       do i = 1,n
         do j = 1,n
           u(i,j) = 0.d0
         enddo
       enddo
c
c       write(43,*)'q2 inside umat ', (q(i),i=1,n)
c
          sum6 = q(n)*q(n)
       do m = 1,n-1
           ii = n-m
             sum6 = (sum6 + q(ii)*q(ii))
           r(m) = dsqrt(sum6)
       enddo
c
         do j = 1, n
           u(1,j) = q(j)/r(n-1)
         enddo
c        
        do i = 2, n-1
          do j = i, n
             u(i,j) = q(i-1)*q(j)/(r(n-i)*r(n-i+1))
          enddo
        enddo
c      
        do i = 2, n-1
            j = i-1
             u(i,j) = - r(n-i)/r(n-i+1)
        enddo
c
        u(n,n-1) = - q(n)/r(1)
        u(n,n) = q(n-1)/r(1)
c
c  r(1) = 0
c
       if (r(1) .eq. 0.d0) then
         u(n-1,n-1) = 1.d0
         u(n-1,n) = 0.d0
         u(n,n-1) = 0.d0
         u(n,n)   = 1.d0
       endif
        return
        end
************************************************************************
c
       subroutine hess2(n,q1,q2,g1,g2,h1,h2,ihc)
c
c  this subroutine gives the hessian at q2. the hessians are obtained
c  either accurately or by Bofill's updating scheme depending on the
c  interation number.
c      ihc - hessian update counter
c      nhup - no. of hessian updates
c
       implicit real*8(a-h,o-z)
      include 'SIZES'
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/GAUSS/ICHRG,IMULTP,IAN(NDAyf),VGAUSS
      COMMON/FORCES/NAT,nat3,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
c
c QM Commom BLock
      COMMON /QMINFO/ qm_choice


      common/ghessb/FA(NDA3,NDA3),trad,rmin,dtver,nstep,nhup,
     *  nhessf,nip
c
       dimension q1(n),q2(n),dq(n),h1(n,n),h2(n,n),g1(n),g2(n)
c
        if ( ihc .eq. nhup ) then
c      get hessian accurately
c
        write(43,*)'getting exact Hessian '
          do i = 1,n
             q(i) = q2(i)
          enddo
           call QMCALC(1)
          do i=1,n
            do j=1,n
              h2(i,j)=fa(i,j)/0.04184d0
            enddo
          enddo
          ihc = 1
c
        else
c      get hessian by Bofill's update
c
             do i = 1, n
                dq(i) = q2(i) - q1(i)
             enddo
             call hessupd(n,dq,g1,g2,h1,h2)
             ihc = ihc + 1
             write(43,*)'Hessian using Bofill update'
          endif
c
          return
          end
c  ***********************************************************************
c
        subroutine hessupd(n,dq,g1,g2,h1,h2)
c
        implicit real*8(a-h,o-z)
c
c    this subroutine updates the hessians by using Boffil's
c    method. Ref: JCC 15 (1994) 1.  this is a combination of
c    Murtagh and Sargent(MS) and Powell-symmetric-Broyden(PSB)
c    methods.
c    Bofill's update:
c    H2 = H1 + phi*del_H_ms + (1 - phi)*del_H_psb
c
c    Murtagh-Sargent part:
c               [(g2(i) - g1(i) - H1*dq)*dq][(g2(i) - g1(i) - H1*dq)*dq]
c    del_H_ms = --------------------------------------------------------
c                          (g2(i) - g1(i) - H1*dq)*dq
c
c    Powell-symmetric-Broyden part:
c                (g2(i) - g1(i) - H1*dq)*dq + dq*(g2(i) - g1(i) - H1*dq)
c    del_H_psb = -------------------------------------------------------
c                                       dq*dq
c                     dq*(g2(i) - g1(i) - H1*dq)*dq*dq
c                  -  --------------------------------
c                                (dq*dq)^2
c
c
c           [(g2(i) - g1(i) - H1*dq)*dq]^2
c    Phi = --------------------------------
c          (dq^2)*(g2(i) - g1(i) - H1*dq)^2
c
c
        dimension dq(n),g1(n),g2(n),h1(n,n),h2(n,n),temp(n),xi(n)
c
c    dq*dq
        dqdq = 0.d0
        do i = 1, n
           dqdq = dqdq + dq(i)*dq(i)
        enddo
c
c    H1(i,j)*dq(j)
        call matvec(0,n,n,h1,dq,temp)
c
c    g2(i) - g1(i) - H1*dq
        do i = 1, n
          xi(i) = g2(i) - g1(i) - temp(i)
        enddo
c
c    (g2(i) - g1(i) - H1*dq)(g2(i) - g1(i) - H1*dq)
        xixi = 0.d0
        do i = 1, n
          xixi = xixi + xi(i)*xi(i)
        enddo
c
c    (g2(i) - g1(i) - H1*dq)*dq
        xidq = 0.d0
        do i = 1, n
          xidq = xidq + xi(i)*dq(i)
        enddo
c
        write(43,*)'xidq  ',xidq
c     Phi
c
        cphi = 1.d0 - (xidq*xidq)/(dqdq*xixi)
        rphi = 1.d0 - cphi
c
c    Bofill's update:
c    H2 = H1 + phi*del_H_ms + (1 - phi)*del_H_psb
c
        do i = 1, n
           do j = 1, i
             delh_ms =  h1(i,j) + xi(i)*xi(j)/xidq
             delh_psb = h1(i,j) + (xi(i)*dq(j) + dq(i)*xi(j))/dqdq
     \        - dq(i)*dq(j)*xidq/(dqdq)**2
             h2(i,j) =  rphi*delh_ms + cphi*delh_psb
             h2(j,i) = h2(i,j)
           enddo
        enddo
c
        return
        end
c
************************************************************************       
       subroutine correct2(nat3,fmod1,q1,v1,g1,h1,q2,v2,g2,h2,am,vin,
     \ q3,v,fmodfin,vp,t,dt,nv,nc,ic,pcdist,q3p)
c    this subroutine corrects the predictor step in multidimensions(3N).
c    the system is integrated from q1 to q3 using the quadratic potential.
c    the coordinates(q3) after integration are projected to a hyperplanes
c    that contain q1 and q2 and are perpendicular to the vector connecting
c    q1 and q2.  the projected points q4 and q5 are then fitted to a
c    5th order polynomial using V, dV and d2V at q4 and q5. then the 
c    corrected potential and the gradients is obtained at q3 using 
c    the fitted potential.  
c                    
       implicit real*8(a-h,o-z)
      INCLUDE 'SIZES'
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      common/ghessb/fa(NDA3,NDA3),trad,rmin,ddt,nstep,
     *  nhesup,nhessf,nip

c
      dimension q1(nat3),q11(nat3),q2(nat3),q3(nat3),u(nat3,nat3),
     \ q33(nat3),q22(nat3),g11(nat3),g22(nat3),h11(nat3,nat3),
     \ h22(nat3,nat3),fmod1(nat3),am(*),nj(3),nk(3),qq3(nat3)
      dimension v0(nat3),v(nat3),vin(nat3),fmodt(nat3),fmod(nat3),
     \ g1(nat3),g2(nat3),h1(nat3,*),h2(nat3,*),q4(nat3),q5(nat3),
     \ fmodfin(nat3)
c
c
       cc1 = 4.184d-4
       dth = dt*0.5d0
       ical = 0
       dtol = 0.01d0
c
       write(43,*)'   CORRECTOR STEP : GAUSSIAN LIKE PROJECTION METHOD '
c
       do i = 1, nat3
          v(i) = vin(i)
          fmodfin(i) = fmod1(i)
       enddo
c
c   translate q1 and q2 and get unitary matrix for rotation
          call utrans(nat3,q1,q2,q11,q22,u)
c
          call rotate1(nat3,g1,g2,h1,h2,u,g11,g22,h11,h22,1)
c
c   velocity-verlet
       do j = 1, nv
c
         if ( j.gt.1) vpm1 = vp
c
c   integrate from q1 to q3 by dt
         write(43,*)'----------- nv =  ',j
          do i = 1,nat3
            kk = (i+2)/3
            fmodt(i) = fmodfin(i)
            v0(i) = v(i)
c  remember old coord.
            qq3(i) = q3(i)
            q3(i) = q3(i) + dt*(v0(i) + dth*fmodt(i)/am(kk))
          enddo
c
c    translate q3 to new coord. and rotate q1, q2, q3 from (q1,q2..q3n)
c    to (r,0,0.....0) coordinate system
          do i = 1, nat3
            q3(i) = q3(i)-q1(i)
          enddo
          call matvec(0,nat3,nat3,u,q3,q33)
c
c   project q3 on hyperplanes to get q4 and q5
          write(43,*)'Projecting q3 to q4 and q5'
          call project(nat3,q22,q33,q4,q5)
c
c   polynomial fit between q4 and q5 and gradient at q3
          call polyfit2(nat3,q11,v1,g11,h11,q22,v2,g22,h22,q33,q4,q5,
     \ vp,fmod,issflag,ical)
c
       if ( issflag .eq. 1 ) then
           write(43,*)'issflag  ',issflag
           write(43,*)'coming our corrector'
           do i = 1, nat3
              q3(i) = qq3(i)
              fmodfin(i) = fmodt(i)
           enddo
           vp = vpm1
c  distance between pred. and corr.
             pcdist = 0.d0
             do i = 1 ,nat3
               pcdist = pcdist + (q2(i)-q3(i))*(q2(i)-q3(i))
             enddo
               pcdist = sqrt(pcdist)
c             write(*,*) 'distance between pred. and corr. ',pcdist
c calculate perpendicular part
          q3p2 = 0.d0
          do i =2,nat3
            q3p2 = q3p2 + q33(i)*q33(i)
          enddo
            q3p = sqrt(q3p2)
            write(43,*) 'Perp. distance between pred. and corr. ',q3p
c
         return
       endif
c
c   rotate coordinates and gradients back to cartesian coordinates for q3
          call rotate2(nat3,u,q1,q33,q3,fmod,fmodfin)
c
       do i = 1,nat3
         fmodfin(i) = -fmodfin(i)*cc1
       enddo
c
          t = t + dt
c
          do i = 1,nat3
            kk = (i+2)/3
            v(i) = v0(i) + dth*(fmodt(i) + fmodfin(i))/am(kk)
          enddo
c
c     check the distance between the predicted and the corrected step
          if ( j .eq. nv ) then
             pcdist = 0.d0
             do i = 1 ,nat3
               pcdist = pcdist + (q2(i)-q3(i))*(q2(i)-q3(i))
             enddo
               pcdist = sqrt(pcdist)
             write(43,*) 'distance between pred. and corr. ',pcdist
c calculate perpendicular part
          q3p2 = 0.d0
          do i =2,nat3
            q3p2 = q3p2 + q33(i)*q33(i)
          enddo
            q3p = sqrt(q3p2)
            write(43,*) 'Perp. distance between pred. and corr. ',q3p
c
          endif

c     print results
          call hessout(nat3,dt,vp,v,am,q3,ic)
          if ( nc.gt.nstep) return
          IF ((NTEST.EQ.2.AND.(NAST.EQ.0.OR.NAST.EQ.1)).OR.
     \ (NTEST.EQ.1.AND.NAST.EQ.0)) then
          write(43,*)'it is going out of corr.'
          return
          endif
c
       enddo
c  end of verlet ..
c
       return
       end
c
************************************************************************
c  
      subroutine polyfit2(nat3,q11,v1,g11,h11,q22,v2,g22,h22,q33,q4,q5,
     \ vp,delv5,issflag,ical)
      implicit real*8(a-h,o-z)
c   this subroutine converts from q3,q4,q5 variables to 's' variable
c   such that r = s(q5-q4)-q4
c   for s=0 ; r = q4   and s=1 ; r = q5
c   then it fits the variable s to a 5th order polynomial
c
      dimension q11(nat3),q22(nat3),q4(nat3),q5(nat3),g11(nat3),
     \ g22(nat3),h11(nat3,*),h22(nat3,*),dvdq_mod(nat3),
     \ d2vdq(nat3,nat3),temp(nat3),delv5(nat3),q33(nat3),d(6),
     \ ain(nat3,nat3)
c
c    get potentials at q4 and q5
c
       write(43,*)'5th order polyomial fit ...'
       call ener_model(nat3,q11,q4,v1,g11,h11,emod_q4)
       call ener_model(nat3,q22,q5,v2,g22,h22,emod_q5)
       write(43,*)'emodq4, emodq5  :',emod_q4,emod_q5
c
c  get the dVm and d2Vm at q4 and q5 w. r. t. parallel component.
c   
        hdq = 0.d0
        do i = 2, nat3
          hdq = hdq + h11(1,i)*q4(i)
        enddo
         dvds_q4 = g11(1) + hdq
c       
        hdq = 0.d0
        do i = 2, nat3
          hdq = hdq + h22(1,i)*q4(i)
        enddo
         dvds_q5 = g22(1) + hdq
c       
       write(43,*)'dvds-q4, dvds-q5  : ',dvds_q4,dvds_q5
c 
c    d2Vm with respect to parallel component at q4 and q5
c
          d2vds_q4 = h11(1,1)
          d2vds_q5 = h22(1,1)
       write(43,*)'d2vds-q4, d2vds-q5  : ',d2vds_q4,d2vds_q5
c
c    get coefficients of the polynomial: Ad = b
c
       s = (q33(1)-q11(1))/abs((q22(1)-q11(1)))
       rr = abs((q22(1)-q11(1)))
c
c   get max of x _|_
c
       write(43,*)'ical  ',ical
       write(43,*)'s ...: ',s
      if ( s.ge.1.d0 .and. ical.eq.0) then
        qp = 0.d0
       do i = 2, nat3
          qp = qp + q33(i)*q33(i)
       enddo
       sthres = 0.5d0*dsqrt(qp)/rr
       ical = 1
       write(43,*)'qp, rs :',qp,rr
       write(43,*)'qp/rs : ',sthres
      endif
c
         write(43,*)'1+sthres :',1.d0+sthres,sthres
       if ( s .gt. (1.d0 + sthres)  ) then
         issflag = 1
         write(43,*)'s :',s,' s > 1:coming out of poly'
         sthres = 0.d0
         return
       else
         issflag = 0
       endif
c
       if ( s .lt. 0.d0 ) then
         call ener_model(nat3,q11,q33,v1,g11,h11,vp)
         call dvdq_model(0,nat3,q11,q33,g11,h11,delv5)
         write(43,*)'s,vp  2nd order :', s,vp
       elseif ( s .gt. 1.d0 ) then
         call ener_model(nat3,q22,q33,v2,g22,h22,vp)
         call dvdq_model(0,nat3,q22,q33,g22,h22,delv5)
         write(43,*)'s,vp  2nd order :', s,vp
       else
       call poly(nat3,s,rr,q22,q33,q4,emod_q4,emod_q5,g11,g22,
     \ h11,h22,dvds_q4,d2vds_q4,dvds_q5,d2vds_q5,vp,delv5)
       write(43,*)'s,vp 5th order  :', s,vp
       endif
c
       return 
       end
c
************************************************************************
c
      subroutine poly(n,s,rr,q22,q33,q4,emod_q4,emod_q5,g11,g22,
     \ h11,h22,dvds_q4,d2vds_q4,dvds_q5,d2vds_q5,vp,grad)
      implicit real*8(a-h,o-z)
c     this subroutine gives the gradient of the polynomial fitted
c     surface. the gradient has 2 terms, 1st term is the derivative 
c     w.r.t. ||(parallel) and the 2nd term is the derivative w.r.t. q_|_
c
c     here subscripts a and b are actually 3 and 4
c     potential:
c      V = Va*p1 + ga*p2 + ha*p3 + Vb*p4 + gb*p5 + hb*p6
c     where: 
c       p1 = 1.d0 - 10.d0*s**3 + 15.d0*s**4 - 6*s**5
c       p2 = rr*(s - 6.d0*s**3 + 8.d0*s**4 - 3.d0*s**5)
c       p3 = 0.5d0*rr*rr*(s**2 - 3.d0*s**3 + 3.d0*s**4 - s**5)
c       p4 = 10.d0*s**3 - 15.d0*s**4 + 6.d0*s**5
c       p5 = rr*(-4*s**3 + 7.d0*s**4 - 3.d0*s**5)
c       p6 = 0.5d0*rr*rr*(s**3 -2.d0*s**4 + s**5)
c
c       polynomials p1,p2,p3,p4,p5,p6 are obtained using the following
c       conditions:
c     p1(0)=1,p1(1)=0,p1'(0)=0,p1'(1)=0,p1"(0)=0,p1"(1)=0
c     p2(0)=0,p2(1)=1,p2'(0)=0,p2'(1)=0,p2"(0)=0,p2"(1)=0
c     p3(0)=0,p3(1)=0,p3'(0)=1,p3'(1)=0,p3"(0)=0,p3"(1)=0
c     p4(0)=0,p4(1)=0,p4'(0)=0,p4'(1)=1,p4"(0)=0,p4"(1)=0
c     p5(0)=0,p5(1)=0,p5'(0)=0,p5'(1)=0,p5"(0)=1,p5"(1)=0
c     p6(0)=0,p6(1)=0,p6'(0)=0,p6'(1)=0,p6"(0)=0,p6"(1)=1
c       s = del_q/abs( q1   - q2  )
c                        ||     ||
c       rr = abs( q1  - q2  ) 
c                   ||    ||
c
c     gradient: has parallel and perpendicular part
c   dVdq|| = dVds * dsdq = 1/rr * dVds
c          = 1/rr * (Va*dp1 + ga*dp2 + ha*dp3 + Vb*dp4 + gb*dp5 + hb*dp6  
c   where:   
c       dp1 = -30.d0*s**2 + 60.d0*s**3 - 30*s**4
c       dp2 = rr*(1 - 18.d0*s**2 + 32.d0*s**3 - 15.d0*s**4)
c       dp3 = 0.5d0*rr*rr*(2*s - 9.d0*s**2 + 12*s**3 - 5.d0*s**4)
c       dp4 = 30.d0*s**2 - 60*s**3 + 30.d0*s**4
c       dp5 = rr*(-12.d0*s**2 + 28.d0*s**3 -15.d0*s**4)
c       dp6 = 0.5d0*rr*rr*(3.d0*s**2 - 8.d0*s**3 + 5.d0*s**4)
c     
c   dvdq    = dVadq   *p1 + dgadq   *p2 + dVbdq   *p4 + dgbdq   *p5
c       _|_        _|_           _|_           _|_           _|_          
c
c      dVadq    = g    + h       *dq   
c           _|_    _|_    _|_,_|_   _|_
c      dgadq    = h      
c           _|_    ||,_|_
c
      dimension g11(n),g22(n),h11(n,n),h22(n,n),grad(n),
     \ q22(n),q33(n),q4(n),dq(n),first(n),second(n),
     \ temp(n),dv4_dq(n),dv5_dq(n),dg4_dq(n),dg5_dq(n)
c
       cc1 = 4.184d-4
c      
c    Fifth order potential
c
       p1 = 1.d0 - 10.d0*s**3 + 15.d0*s**4 - 6*s**5
       p2 = rr*(s - 6.d0*s**3 + 8.d0*s**4 - 3.d0*s**5)
       p3 = 0.5d0*rr*rr*(s**2 - 3.d0*s**3 + 3.d0*s**4 - s**5)
       p4 = 10.d0*s**3 - 15.d0*s**4 + 6.d0*s**5
       p5 = rr*(-4*s**3 + 7.d0*s**4 - 3.d0*s**5)
       p6 = 0.5d0*rr*rr*(s**3 -2.d0*s**4 + s**5)
       vp = emod_q4*p1 + dvds_q4*p2 + d2vds_q4*p3 + emod_q5*p4 +
     \ dvds_q5*p5 + d2vds_q5*p6
c
c     gradient
c     parallel component:      
c
       dp1 = -30.d0*s**2 + 60.d0*s**3 - 30*s**4
       dp2 = rr*(1 - 18.d0*s**2 + 32.d0*s**3 - 15.d0*s**4)
       dp3 = 0.5d0*rr*rr*(2*s - 9.d0*s**2 + 12*s**3 - 5.d0*s**4)
       dp4 = 30.d0*s**2 - 60*s**3 + 30.d0*s**4
       dp5 = rr*(-12.d0*s**2 + 28.d0*s**3 -15.d0*s**4)
       dp6 = 0.5d0*rr*rr*(3.d0*s**2 - 8.d0*s**3 + 5.d0*s**4)
c
       grad(1) = (emod_q4*dp1 + dvds_q4*dp2 + d2vds_q4*dp3 + emod_q5*dp4
     \  + dvds_q5*dp5 + d2vds_q5*dp6)/rr
c
c      perpendicular part:
c       dv4/dq_|_
        call matvec(0,n,n,h11,q4,temp)
        do i = 2, n
          dv4_dq(i) = g11(i) + temp(i)
        enddo
          dv4_dq(1) = 0.d0
c       dg4/dq_|_
        do i = 2, n
          dg4_dq(i) = h11(1,i)
        enddo
          dg4_dq(1) = 0.d0
c      dv5/dq_|_
        call matvec(0,n,n,h22,q4,temp)
        do i = 2, n
          dv5_dq(i) = g22(i) + temp(i)
        enddo
          dv5_dq(1) = 0.d0
c       dg4/dq_|_
        do i = 2, n
          dg5_dq(i) = h22(1,i)
        enddo
          dg5_dq(1) = 0.d0
c
        do i = 2, n
          grad(i) = dv4_dq(i)*p1 + dg4_dq(i)*p2 + dv5_dq(i)*p4 + 
     \ dg5_dq(i)*p5 
        enddo
c
c
       return
       end
c    
************************************************************************
      SUBROUTINE HESSOUT(nat3,tt,vp,vel,am,q3,ic)
      implicit double precision (a-h,o-z)
      INCLUDE 'SIZES'
      parameter (cc1=0.04184D0)
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      COMMON/QPDOT/q(nda3),PDOT(nda3)
      COMMON/PQDOT/P(nda3),QDOT(nda3),VW(NDA)
      COMMON/INTEGR/ATIME,nv,NID
      common/ghessb/FA(NDA3,NDA3),trad,rmin,dtver,nstep,nhup,
     *  nhessf,nip
c
      dimension am(nat3),vel(nat3),q3(nat3)
c
c     this is where output will be written when using the variable
c     time step integrator
c
c     modified to use venus GWRITE and FINAL 
c
      time=tt*0.1d0
          t = 0.d0
          do i = 1, nat3
              q(i) = q3(i)
              p(i)=vel(i)*am((i+2)/3)*10.d0
              t = t + p(i)*p(i)/(2.d0*am((i+2)/3))
          enddo
          v = vp
          t = t/cc1
          h = t + v
      if(nc.ge.0)then
      call test2(v)
      IF (NTEST.EQ.1.AND.NAST.EQ.0) GOTO 406
      IF (NTEST.EQ.2.AND.NAST.EQ.1) GOTO 410
      IF (NTEST.EQ.2.AND.NAST.EQ.0) GOTO 410
      IF (NTEST.EQ.0.AND.NAST.NE.0) NAST=0
      IF (NTEST.EQ.1.AND.NAST.EQ.2) NAST=1
      ncounter = nc/(nip*ic)
      if (ncounter.eq.1) then
      call gwrite
      ic = ic + 1
      endif
      nc=nc+1
      RETURN
  406 NAST=1
      WRITE(6,*)
  910 FORMAT(10X,'REACTION OCCURRED FOR PATH',I3)
      WRITE(6,910)NPATH
      CALL GWRITE
      RETURN
  410 CONTINUE
      CALL GWRITE
  414 CONTINUE
      CALL FINAL
      CALL GFINAL
      endif
      return
      END
************************************************************************
      SUBROUTINE TEST2(vv)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CHECK FOR INTERMEDIATE AND FINAL EVENTS
C
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/PSN2/PESN2,GA,RA,RB
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,L(NDA),NAM
      COMMON/TESTIN/VRELO,INTST
      COMMON/TESTSN2/GAO,NSAD,NCBA,NCAB,IBAR
      COMMON/FINALB/EROTA,EROTB,EA(3),EB(3),AMA(4),AMB(4),AN,AJ,BN,BJ,
     *OAM(4),EREL,ERELSQ,ETCM,BF,SDA,SDB,DELH(NDP),ANG(NDG),NFINAL
      COMMON/VMAXB/QVMAX(NDA3),PVMAX(NDA3),VMAX,NCVMAX      
c      COMMON/SELTB/QZ(nda3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      DIMENSION QCMA(3),VCMA(3),QCMB(3),VCMB(3),QR(3),VR(3)
      CHARACTER*10 TYPE
      CHARACTER*15 COMP
C
 900  FORMAT(4X,' TURNING POINT #  ','  CYCLE ','  RCM(A)',
     *       '    EA     ','    EB     ','    EROTA  ','    EROTB  ',
     *       '    JA     ','    JB     ','    L      ')
 903  FORMAT(8X,A10,I4,I8,F8.3,1P7D11.4)
 905  FORMAT(/5X,'$$$$BARRIER CROSSING NUMBER$$$$ ',I6,
     &       '  AT CYCLE',I8)
 906  FORMAT(5X,'$$$$BARRIER CROSSING FROM B TO A $$$$')
 907  FORMAT(5X,'$$$$BARRIER CROSSING FROM A TO B $$$$')
 935  FORMAT(7X,'RA= ',F7.3,3X,'RB= ',F7.3,3X,'GA= ',F7.3)
 910  FORMAT(4X,' TURNING POINT #  ',5X,' COMPLEX ',6X,
     *'  CYCLE ','  RCM(A)',
     *       '    EA     ','    EB     ','    EROTA  ','    EROTB  ',
     *       '    JA     ','    JB     ','    L      ')
 913  FORMAT(8X,A10,I4,2X,A17,I8,F8.3,1P7D11.4)
C
C         TEST FOR MAXIMUM IN POTENTIAL ENERGY
C
c      CALL ENERGY      
       v = vv
c       t = tke
c       h = v + t
      IF (V.GT.VMAX) THEN
	 VMAX=V
	 NCVMAX=NC
	 DO I=1,NDA3
	 QVMAX(I)=Q(I)
	 PVMAX(I)=P(I)
	 ENDDO
      ENDIF
C
C         COUNT INNER TURNING POINTS BETWEEN THE TWO REACTANT'S
C         CENTERS OF MASS.
C
      IF (NSN2.EQ.0) THEN
         IF (NABJ(1).EQ.NABK(1)) GOTO 1010
C
C         CALCULATE CENTER OF MASS COORDINATES AND VELOCITIES FOR A
C
         WT=WTA(1)
         NN=NATOMA(1)
         DO I=1,NN
            L(I)=LA(1,I)
         ENDDO
         CALL CENMAS(WT,QCMA,VCMA,NN)
C
C         CALCULATE CENTER OF MASS COORDINATES AND VELOCITIES FOR B
C
         WT=WTB(1)
         NN=NATOMB(1)
         DO I=1,NN
            L(I)=LB(1,I)
         ENDDO
         CALL CENMAS(WT,QCMB,VCMB,NN)
C
C         CALCULATE INTERNUCLEAR SEPARATION BETWEEN THE CENTERS-OF-
C         MASS, AND THE RELATIVE VELOCITY ALONG THE INTERNUCLEAR AXIS.
C
         DO I=1,3
            QR(I)=QCMA(I)-QCMB(I)
            VR(I)=VCMA(I)-VCMB(I)
         ENDDO
         RCM=0.0D0
         VREL=0.0D0
         DO I=1,3
            RCM=RCM+QR(I)**2
            VREL=VREL+QR(I)*VR(I)
         ENDDO
         RCM=SQRT(RCM)
         VREL=VREL/RCM
C
         SGN=VRELO*VREL
         IF (SGN.LT.0.0D0) THEN
            INTST=INTST+1
            IF (INTST.EQ.1) WRITE(6,900)
C
C         WRITE ENERGIES AND ANGULAR MOMENTUM FOR PATH=1 AT THE
C         TURNING POINT
C
            NPATH=1
            CALL FINAL
            IF (VRELO.LT.0.0D0) TYPE=' **INNER** '
            IF (VRELO.GT.0.0D0) TYPE=' **OUTER** '
            WRITE(6,903)TYPE,INTST,NC,RCM,EA(3),EB(3),EROTA,
     *                  EROTB,AMA(4),AMB(4),OAM(4)
            NFINAL=0
         ENDIF
         VRELO=VREL
C 
C         CODE FOR SN2 DYNAMICS
C
      ELSE
         IF (NSN2.LE.6) GATS=0.0D0
         IF (NSN2.EQ.7.OR.NSN2.EQ.8) GATS=0.008D0
         CBAR=GAO*(GA-GATS)
         IF (CBAR.LT.0.0D0) THEN
            IBAR=IBAR+1
            WRITE(6,905)IBAR,NC
            WRITE(6,935)RA,RB,GA
            IF (GAO.LT.0.0D0) THEN
               NCBA=NCBA+1
               WRITE(6,906)
            ENDIF
            IF (GAO.GT.0.0D0) THEN
               NCAB=NCAB+1
               WRITE(6,907)
            ENDIF
         ENDIF 
C
         GAO=GA-GATS
         IF (NABJ(1).EQ.NABK(1)) THEN
            IF (GA.GE.GATS) NCPL=2
            IF (GA.LT.GATS) NCPL=3
         ELSE
            IF (GA.GE.GATS) NCPL=1
            IF (GA.LT.GATS) NCPL=2
         ENDIF
         WT=WTA(NCPL)
         NN=NATOMA(NCPL)
         DO I=1,NN
            L(I)=LA(NCPL,I)
         ENDDO
         CALL CENMAS(WT,QCMA,VCMA,NN)
C
C         CALCULATE CENTER OF MASS COORDINATES AND VELOCITIES FOR B
C
         WT=WTB(NCPL)
         NN=NATOMB(NCPL)
         DO I=1,NN
            L(I)=LB(NCPL,I)
         ENDDO
         CALL CENMAS(WT,QCMB,VCMB,NN)
C
C         CALCULATE INTERNUCLEAR SEPARATION BETWEEN THE CENTERS OF 
C         MASS, AND THE RELATIVE VELOCITY ALONG THE INTERNUCLEAR AXIS.
C
         DO I=1,3
            QR(I)=QCMA(I)-QCMB(I)
            VR(I)=VCMA(I)-VCMB(I)
         ENDDO
         RCM=0.0D0
         VREL=0.0D0
         DO I=1,3
            RCM=RCM+QR(I)*QR(I)
            VREL=VREL+QR(I)*VR(I)
         ENDDO
         RCM=SQRT(RCM)
         VREL=VREL/RCM
C
         SGN=VRELO*VREL
         IF (SGN.LT.0.0D0) THEN
            INTST=INTST+1
            IF (INTST.EQ.1) WRITE(6,910)
C
C         WRITE ENERGIES AND ANGULAR MOMENTUM FOR PATH=NCPL AT THE
C         TURNING POINT
C
            NPATH=NCPL
            CALL FINAL
            IF (NABJ(1).EQ.NABK(1)) THEN
               IF (NCPL.EQ.2) COMP='## X---CH3Y ## '
               IF (NCPL.EQ.3) COMP='## XCH3---Y ## '
            ELSE
               IF (NCPL.EQ.1) COMP='## X---CH3Y ## '
               IF (NCPL.EQ.2) COMP='## XCH3---Y ## '
            ENDIF
            IF (VRELO.LT.0.0D0) TYPE=' **INNER** '
            IF (VRELO.GT.0.0D0) TYPE=' **OUTER** '
            WRITE(6,913)TYPE,INTST,COMP,NC,RCM,EA(3),EB(3),
     *                  EROTA,EROTB,AMA(4),AMB(4),OAM(4)
            NFINAL=0
         ENDIF
         VRELO=VREL
      ENDIF
C
C         END OF CODE WHICH TESTS FOR INNER TURNING POINTS
C
 1010 CONTINUE
C
C         DEALING WITH TURNING POINT OF HEIGHT ABOVE SURFACE 
C
c      IF (NSURF.NE.0) CALL HTURN(1,13.0)
c
c         TEST FOR REACHING RBAR(I) OR RMAX(I)
c
      NTEST=0
      M=NPATHS+1
      DO I=1,M
         NPATH=I
         J3=3*NABJ(I)
         J2=J3-1
         J1=J2-1
         K3=3*NABK(I)
         IF (J3.NE.K3) THEN
            K2=K3-1
            K1=K2-1
            JK=(NABJ(I)-1)*(2*NATOMS-NABJ(I))/2+NABK(I)-NABJ(I)
            T1=Q(K1)-Q(J1)
            T2=Q(K2)-Q(J2)
            T3=Q(K3)-Q(J3)
            R(JK)=SQRT(T1*T1+T2*T2+T3*T3)
            IF (NABL(I).GT.0)THEN
                J3=3*NABL(I)
                J2=J3-1
                J1=J2-1
                K3=3*NABM(I)
                K2=K3-1
                K1=K2-1
                LM=(NABL(I)-1)*(2*NATOMS-NABL(I))/2+NABM(I)-NABL(I)
                T1=Q(K1)-Q(J1)
                T2=Q(K2)-Q(J2)
                T3=Q(K3)-Q(J3)
                R(LM)=SQRT(T1*T1+T2*T2+T3*T3)
                IF((R(JK).GT.RBAR(I)).AND.(R(LM).GT.RBAR(I))) NTEST =1
                IF((R(JK).GT.RMAX(I)).AND.(R(LM).GT.RMAX(I))) NTEST =2 
                IF (NTEST.GT.0) GOTO 3
            ELSE
                IF (R(JK).GT.RBAR(I)) NTEST=1
                IF (R(JK).GT.RMAX(I)) NTEST=2
                IF (NTEST.GT.0) GOTO 3
            ENDIF
             
         ENDIF
      ENDDO
      NPATH=1
    3 RETURN
      END
************************************************************************
       subroutine correct3(nat3,fmod1,q1,v1,g1,h1,q2,v2,g2,h2,am,vin,
     \ q3,v,fmodfin,vp,t,dt,nv,nc,ic,pcdist,q3p)
c    this subroutine corrects the predictor step in multidimensions(3N).
c    the system is integrated from q1 to q3 using the quadratic potential.
c    the coordinates(q3) after integration are projected to a hyperplanes
c    that contain q1 and q2 and are perpendicular to the vector connecting
c    q1 and q2.  the projected points q4 and q5 are then fitted to a
c    5th order polynomial using gauss-elimination by using V, dV and d2V
c    at q4 and q5. then the corrected potential is obtained at q3 using
c    the fitted potential.  the projection of q3 on the hyperplanes is
c    obtained by gauss-elimination.
c
       implicit real*8(a-h,o-z)
      INCLUDE 'SIZES'
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      common/ghessb/fa(NDA3,NDA3),trad,rmin,ddt,nstep,
     *  nhesup,nhessf,nip

      dimension q1(nat3),q2(nat3),q3(nat3),fmod1(nat3),am(*)
      dimension v0(nat3),v(nat3),vin(nat3),fmodt(nat3),fmod(nat3),
     \ g1(nat3),g2(nat3),h1(nat3,*),h2(nat3,*),q4(nat3),q5(nat3),
     \ fmodfin(nat3),qq3(nat3)
c                    

      write(43,*)' General Cartesian Coordinate corrector '
c
       cc1 = 4.184d-4
       dth = dt*0.5d0
       ical = 0
       dtol = 0.01d0
c
       do i = 1, nat3
          v(i) = vin(i)
          fmodfin(i)=fmod1(i)
       enddo
c   integrate from q to q3 by dt
c   velocity-verlet
       do j = 1, nv
c
         if ( j.gt.1) vpm1 = vp
c
c   integrate from q to q3 by dt
         write(43,*)'-----------------nv =  ',j
          do i = 1,nat3
            kk = (i+2)/3
            fmodt(i) = fmodfin(i)
            v0(i) = v(i)
c  remember old coord.
            qq3(i) = q3(i)
            q3(i) = q3(i) + dt*(v0(i) + dth*fmodt(i)/am(kk))
          enddo
c
c   project q3 on hyperplanes to get q4 and q5
c
       call project3(nat3,q1,q2,q3,q4,q5)
c
c   polynomial fit between q4 and q5. get potential and gradient at q3
          call polyfit3(nat3,q1,v1,g1,h1,q2,v2,g2,h2,q3,q4,q5,
     \ vp,fmod,issflag,ical)
c
c
       if ( issflag .eq. 1 ) then
           write(43,*)'issflag  ',issflag
           write(43,*)'coming out of corrector'
           do i = 1, nat3
              q3(i) = qq3(i)
              fmodfin(i) = fmodt(i)
           enddo
           vp = vpm1
c  distance between pred. and corr.
             pcdist = 0.d0
             do i = 1 ,nat3
               pcdist = pcdist + (q2(i)-q3(i))*(q2(i)-q3(i))
             enddo
               pcdist = sqrt(pcdist)
             write(43,*) 'distance between pred. and corr. ',pcdist
c calculate perpendicular part
          q3p2 = 0.d0
          do i =2,nat3
            q3p2 = q3p2 + (q2(i)-q3(i))*(q2(i)-q3(i))
          enddo
            q3p = sqrt(q3p2)
            write(43,*)'Perp. distance between pred. and corr. ',q3p
c
         return
       endif
c
       do i = 1,nat3
         fmodfin(i) = -fmod(i)*cc1
       enddo
c
          t = t + dt
c
          do i = 1,nat3
            kk = (i+2)/3
            v(i) = v0(i) + dth*(fmodt(i) + fmodfin(i))/am(kk)
          enddo
c
c     check the distance between the predicted and the corrected step
          if ( j .eq. nv ) then
             pcdist = 0.d0
             do i = 1 ,nat3
               pcdist = pcdist + (q2(i)-q3(i))*(q2(i)-q3(i))
             enddo
               pcdist = sqrt(pcdist)
             write(43,*) 'distance between pred. and corr. ',pcdist
c calculate perpendicular part
          q3p2 = 0.d0
          do i =2,nat3
            q3p2 = q3p2 + (q2(i)-q3(i))*(q2(i)-q3(i))
          enddo
            q3p = sqrt(q3p2)
            write(43,*)'Perp. distance between pred. and corr. ',q3p
c
          endif
c  print results
          call hessout(nat3,dt,vp,v,am,q3,ic)
          if ( nc.gt.nstep) return
          IF ((NTEST.EQ.2.AND.(NAST.EQ.0.OR.NAST.EQ.1)).OR.
     \ (NTEST.EQ.1.AND.NAST.EQ.0)) then
          write(43,*)'it is going out of corr.'
          return
          endif
c
       enddo
c  end of verlet ..
c
       return
       end
c          
************************************************************************
c
      subroutine project3(n,q1,q2,q3,q4,q5)
      implicit real*8(a-h,o-z)
c    this subroutine gets the projection of q3 on the hyperplanes to get
c    q4 and q5.  this is done by setting up 3n equations and equating them
c    to the equation with the maximum of q2-q1. once the corresponding q4 
c    and q5 are calculated, the rest of the q4s and q5s are obtained.
c
c    xi4 - xi3    xi04 - xi03
c    ---------  = -----------   where, ai = q2(i) - q1(i)
c       ai            ai0 
c
c    xi04 = Sum_over [ xi03*ai^2/ai0 + ai*xi1 - ai*xi3 ]/ Sum_over[ai^2/ai0]
c
c    xi4 = ai/ai0 * ( xi04 - xi03 ) + xi3
c
      dimension q1(*),q2(*),q3(*),q4(*),q5(*),ai(n)
c
       do i = 1,n
         ai(i) = q2(i)-q1(i)
       enddo
c
c    get max of a(i) to equate
          ai0 = abs(ai(1))
          k = 1
       do i = 2,n
          if ( abs(ai(i)) .ge. ai0 ) then
          ai0 = abs(ai(i))
          k = i
          endif
       enddo
          ai0 = ai(k)
c    assign coord. corresponding to ai0
       qi03 = q3(k)
c    get qi04, qi05
         xnum1 = 0.d0
         xnum2 = 0.d0
         den = 0.d0
       do i = 1, n
         xnum1=xnum1 + (ai(i)*ai(i)*qi03/ai0) + ai(i)*q1(i) -ai(i)*q3(i)
         xnum2=xnum2 + (ai(i)*ai(i)*qi03/ai0) + ai(i)*q2(i) -ai(i)*q3(i)
         den = den + ai(i)*ai(i)
       enddo
          den = den/ai0
          qi04 = xnum1/den
          qi05 = xnum2/den
c   get q4 and q5
       do i = 1, n
          q4(i) = ai(i)*(qi04 - qi03)/ai0 + q3(i)
          q5(i) = ai(i)*(qi05 - qi03)/ai0 + q3(i)
       enddo
c  
      return
      end
c
************************************************************************
c
      subroutine polyfit3(nat3,q1,v1,g1,h1,q2,v2,g2,h2,q3,q4,q5,
     \ vp,delv5,issflag,ical)
      implicit real*8(a-h,o-z)
c   this subroutine converts from q3,q4,q5 variables to 's' variable
c   such that r = s(q5-q4)-q4
c   for s=0 ; r = q4   and s=1 ; r = q5
c   then it fits the variable s to a 5th order polynomial
c
      dimension q1(nat3),q2(nat3),q4(nat3),q5(nat3),g1(nat3),
     \ g2(nat3),h1(nat3,*),h2(nat3,*),dvdq_mod(nat3),dq(nat3),
     \ d2vdq(nat3,nat3),temp(nat3),delv5(nat3),q3(nat3),d(6),
     \ ain(nat3,nat3)
c
c    get potentials at q4 and q5
c
       write(43,*)'5th order polyomial fit ...'
       call ener_model(nat3,q1,q4,v1,g1,h1,emod_q4)
       call ener_model(nat3,q2,q5,v2,g2,h2,emod_q5)
       write(43,*)'emodq4, emodq5  :',emod_q4,emod_q5
c
       do i = 1, nat3
          dq(i) = q2(i) - q1(i)
       enddo
c
c   get the dVm and d2Vm with 's' variable at q4 and q5
c    dVm/ds = dVm/dq * (q5-q4) or dVm/dq * (q2-q1)
c
       call dvdq_model(0,nat3,q1,q4,g1,h1,dvdq_mod)
       dvds_q4 = 0.d0
       do i = 1,nat3
         dvds_q4 = dvds_q4 + dvdq_mod(i)*dq(i)
       enddo
       call dvdq_model(0,nat3,q2,q5,g2,h2,dvdq_mod)
       dvds_q5 = 0.d0
       do i = 1,nat3
         dvds_q5 = dvds_q5 + dvdq_mod(i)*dq(i)
       enddo
       write(43,*)'dvds-q4, dvds-q5  : ',dvds_q4,dvds_q5
c
c    d2Vm/ds = (q5-q4) * d2Vm/dq * (q5-q4)
c    d2Vm/dq = d2V/dq ; again (q5-q4) = (q2-q1)
c    
       call matvec(0,nat3,nat3,h1,dq,temp)
       d2vds_q4 = 0.d0
       do i = 1,nat3
          d2vds_q4 = d2vds_q4 + dq(i)*temp(i)
       enddo
       call matvec(0,nat3,nat3,h2,dq,temp)
       d2vds_q5 = 0.d0
       do i = 1,nat3
          d2vds_q5 = d2vds_q5 + dq(i)*temp(i)
       enddo
c
       write(43,*)'d2vds-q4, d2vds-q5  : ',d2vds_q4,d2vds_q5
c
c    get coefficients of the polynomial: Ad = b
c
       call coeff(emod_q4,emod_q5,dvds_q4,dvds_q5,d2vds_q4,d2vds_q5,d)
c
c    get the inverse of matrix A : A^(-1)
c
       call ainverse(6,ain)
c
c    calculate potential and gradient at q3 as functions of q_|_ and s
c
        abs_q = 0.d0
        s = 0.d0
        do i = 1, nat3
           abs_q = abs_q + (q2(i)-q1(i))*(q2(i)-q1(i))
           s = s + (q3(i)-q4(i))*(q2(i)-q1(i))
        enddo
           s = s/abs_q
       write(43,*)'ical  ',ical
       write(43,*)'s ...: ',s
c
      if ( s.ge.1.d0 .and. ical.eq.0) then
        qp = 0.d0
       do i = 2, nat3
          qp = qp + q3(i)*q3(i)
       enddo
       sthres = 0.5d0*dsqrt(qp)/rr
       ical = 1
       write(43,*)'qp, rs :',qp,rr
       write(43,*)'qp/rs : ',sthres
      endif
c
         write(43,*)'1+sthres :',1.d0+sthres,sthres
       if ( s .gt. (1.d0 + sthres)  ) then
         issflag = 1
         write(43,*)'s :',s,' s > 1:coming out of poly'
         sthres = 0.d0
         return
       else
         issflag = 0
       endif

       if ( s .lt. 0.d0 ) then
       call ener_model(nat3,q1,q3,v1,g1,h1,vp)
       call dvdq_model(0,nat3,q1,q3,g1,h1,delv5)
         write(43,*)'s,vp  2nd order :', s,vp
       elseif ( s .gt. 1.d0 ) then
       call ener_model(nat3,q2,q3,v2,g2,h2,vp)
       call dvdq_model(0,nat3,q2,q3,g2,h2,delv5)
         write(43,*)'s,vp  2nd order :', s,vp
       else
       call vpoly(s,d,vp)
       call gradient3(nat3,d,ain,s,q1,q2,q3,q4,q5,g1,g2,h1,h2,delv5)
       write(43,*)'s,vp 5th order  :', s,vp
       endif
c
       return
       end
c
************************************************************************
c
      subroutine gradient3(n,d,b,s,q1,q2,q3,q4,q5,g1,g2,h1,h2,grad)
      implicit real*8(a-h,o-z)
c     this subroutine gives the gradient of the polynomial fitted
c     surface. the gradient has 2 terms, 1st term is the derivative
c     w.r.t. 's' and the 2nd term is the derivative w.r.t. q_|_
c      b : A_inv
c      d : coeff of the polynomial
c
      dimension b(6,6),g1(n),g2(n),h1(n,n),h2(n,n),q1(n),
     \ dedq(6,n),q2(n),q3(n),q4(n),dq(n),first(n),second(n),q5(n),
     \ temp(n),d(6),dq1(n),grad(n)
c
       cc1 = 4.184d-4
       do i = 1, n
          dq1(i) = q2(i)-q1(i)
       enddo
c
c     1st term :  (dvds)*(q2-q1)/(abs(q2))**2
c
       sumq2 = 0.d0
       do i = 1, n
         sumq2 = sumq2 + dq1(i)*dq1(i)
       enddo
         abs_q2_sq = sumq2
c
c   dvds : Sum of ( k*d(k+1)*s**(k-1) ) , k = 1-5
c       s = q3(1)
       dvds = 0.d0
       do k = 1, 5
          dvds = dvds + k*d(k+1)*s**(k-1)
       enddo
c
       do i = 1, n
         first(i) = dq1(i)*dvds/abs_q2_sq
       enddo
c
c    2nd term :  Sum over k,j ( b(k,j)*(de(k)dq)*s(k-1) ), k,j = 1-6
c    dedq(k,n)
       do i = 1, n
          dq(i) = q4(i)-q1(i)
       enddo
       write(43,*)'dq ',(dq(i),i=1,n)
c    k =1
       call matvec(0,n,n,h1,dq,temp)
        do i = 1, n
          dedq(1,i) = g1(i) + temp(i)
        enddo
       write(43,*)'temp  :',(g1(i),temp(i),i=1,n)
c    k = 2
       call matvec(0,n,n,h1,dq1,temp)
        do i = 1, n
          dedq(2,i) = temp(i)
        enddo
c    k = 3, 6
        do i = 1, n
          dedq(3,i) = 0.d0
          dedq(6,i) = 0.d0
        enddo
c    k = 4
       do i = 1, n
          dq(i) = q5(i)-q2(i)
       enddo
       call matvec(0,n,n,h2,dq,temp)
        do i = 1, n
          dedq(4,i) = g2(i) + temp(i)
        enddo
c    k = 5
       call matvec(0,n,n,h2,dq1,temp)
        do i = 1, n
          dedq(5,i) = temp(i)
        enddo
c
c   A_inv * dedq(k,n) * s**(k-1) : n components
       do i = 1, n
         sum3 = 0.d0
         do k = 1, 6
           sum4 = 0.d0
           do j = 1, 6
               sum4 = sum4 + b(k,j)*dedq(j,i)
           enddo
           sum3 = sum3 + sum4*s**(k-1)
         enddo
         second(i) = sum3
       enddo
       call project_grad(n,dq1,g1,g2,second,grad)
c
c    gradient = first + second
       do i = 1, n
c          grad(i) = (first(i) + second(i))*cc1
c          grad(i) = (first(i) + second(i))
          grad(i) = (first(i) + grad(i))
       enddo
c
       return
       end
c
************************************************************************
      subroutine project_grad(n,ai,g1,g2,g3,g4)
      implicit real*8(a-h,o-z)
c    this subroutine gets the projection of gradients at q3 on the 
c    hyperplanes that passes through the origin. it is done in a similar
c    way as the projection of q3 to get q4 and q5
c
      dimension g1(*),g2(*),g3(*),g4(*),ai(n),q1(n),q2(n)
c
c       do i = 1,n
c         ai(i) = q2(i)-q1(i)
c       enddo
c
c    get max of a(i) to equate
          ai0 = abs(ai(1))
          k = 1
       do i = 2,n
          if ( abs(ai(i)) .ge. ai0 ) then
          ai0 = abs(ai(i))
          k = i
          endif
       enddo
          ai0 = ai(k)
c    assign grad. corresponding to ai0
       gi03 = g3(k)
c    get qi04, qi05
         xnum1 = 0.d0
         den = 0.d0
       do i = 1, n
         xnum1 = xnum1 + (ai(i)*ai(i)*gi03/ai0) - ai(i)*g3(i)
         den = den + ai(i)*ai(i)
       enddo
          den = den/ai0
          gi04 = xnum1/den
c   get q4 and q5
       do i = 1, n
          g4(i) = ai(i)*(gi04 - gi03)/ai0 + g3(i)
       enddo
c
      return
      end
c
************************************************************************

