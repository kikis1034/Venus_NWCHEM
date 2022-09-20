c
c    Velocity verlet (1) and Velocity-corrected verlet (2) Integrator
c    Coded by Ki Song, 5/23/05
c
c    Velocity Verlet
c      x(t+dt)=x(t)+v(t)dt+f(t)/2m*dt**2
c      v(t+dt)=v(t)+(f(t+dt)+ft)/2m*dt
c
c    Velocity-corrected verlet
c      x(t+dt)=2x(t)-x(t-dt)+f(t)/2m*dt**2
c      v(t)=(8(x(t+dt)-x(t-dt)-(x(t+2dt)-x(t-2dt)))/12dt
c
c    f(t)=-dv/dq=pdot
c    venus put pdot like this. SYMPLE.f uses this expression.
c
      subroutine verlet(nverlet)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/INTEGR/ATIME,NI,NID
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      dimension ft(nda3),qtp1(nda3),qtm1(nda3),qtp2(nda3),qtm2(nda3)
      dimension qc(nda3)
      save qc,qtm1,qtm2,qtp1
c
c   Velocity Verlet
c
      if(nverlet.eq.1) then
        do k=1,ni
          kk=(k+2)/3
          ft(k)=pdot(k)
          q(k)=q(k)+p(k)*atime/w(kk)+pdot(k)/(2.d0*w(kk))*atime**2
        enddo
        call dvdq
c
c   p=mv so p(t+dt)=p(t)+(f(t+dt)+f(t))*dt/2
c
        do k=1,ni
          p(k)=p(k)+(ft(k)+pdot(k))*atime*0.5d0
        enddo
c*********************************************************************
c   Velocity-corrected Verlet
c*********************************************************************
      elseif(nverlet.eq.2)then
c   if nc=1 then do normal verlet and store q to qtm1
        if(nc.eq.1)then
          do k=1,ni
            kk=(k+2)/3
            qtm2(k)=q(k)
            ft(k)=pdot(k)
            q(k)=q(k)+p(k)*atime/w(kk)+pdot(k)/(2.d0*w(kk))*atime**2
            qtm1(k)=q(k)
          enddo
          call dvdq
          do k=1,ni
            p(k)=p(k)+(pdot(k)+ft(k))*atime*0.5d0
          enddo
          return
c   if nc=2 then do normal verlet and store the q in qc
        elseif(nc.eq.2)then
          do k=1,ni
            kk=(k+2)/3
            ft(k)=pdot(k)
            q(k)=q(k)+p(k)*atime/w(kk)+pdot(k)/(2.d0*w(kk))*atime**2
            qc(k)=q(k)
          enddo
          call dvdq
          do k=1,ni
            p(k)=p(k)+(ft(k)+pdot(k))*atime*0.5d0
          enddo
c    another normal verlet to get qtp1
c    the ft(k) serves as pdot(t)
          do k=1,ni
            kk=(k+2)/3
            ft(k)=pdot(k)
            q(k)=q(k)+p(k)*atime/w(kk)+pdot(k)/(2.d0*w(kk))*atime**2
            qtp1(k)=q(k)
          enddo
          call dvdq
          do k=1,ni
            p(k)=p(k)+(ft(k)+pdot(k))*atime*0.5d0
          enddo
        else
c
c   for nc>2 need to evaluate pdot at qtp1 to calculate qtp2
c
          do k=1,ni
            q(k)=qtp1(k)
          enddo
          call dvdq
        endif
c   for nc>2, general velocity-corrected verlet
c   now we have qtm2,qtm1,qtp1. Need to calculate qtp2
        do k=1,ni
          kk=(k+2)/3
          qtp2(k)=2.d0*qtp1(k)-qc(k)+pdot(k)/w(kk)*atime**2
        enddo
c
c   now we are ready to calculate corrected velocyt at time t
c   after calculating corrected momenta=v*m, shift everything for next step
c
        do k=1,ni
          kk=(k+2)/3
          p(k)=8.d0*(qtp1(k)-qtm1(k))-(qtp2(k)-qtm2(k))
          p(k)=p(k)*w(kk)/(12.d0*atime)
          q(k)=qc(k)
          qtm2(k)=qtm1(k)
          qtm1(k)=qc(k)
          qc(k)=qtp1(k)
          qtp1(k)=qtp2(k)
        enddo
c  Is this call to dvdq necessary?
c  tried to put previous pdot but didn't work
        call dvdq
c
c   need to calculate pdot for next step using current q(k) which will be
c   qtp1 at next step using normal verlet
c
      endif
      return
      end
