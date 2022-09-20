C A subroutine in which all the children stay in wait state
      Subroutine VENUS_parallel_loop (myid,idop,theory, basis,
     *  natom, coordinates,labels,printName,qm_hessian,qm_grad,
     * qm_energy,charge,stack,heap,global,qm_ini_file)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      integer myid, idop
      character*256 theory
      character*256 basis
      character*256 printName
      integer natom,i,j
      double precision coordinates (nda3)
      character*16 labels(nda)
      double precision qm_grad (nda3),qm_hessian(nda3,nda3)
      double precision  qm_energy,charge
      integer stack,heap,global
      character*256 qm_ini_file

      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     * NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     * NRAX,NONB,NMO,NCRCO6

c      natom = 0
c      theory = " "
c      basis = " "
c      qm_ini_file = " "
c      stack = 0
c      heap = 0
c      global =0
c      do i = 1,nda3
c         coordinates(i) = 0.0d0
c      enddo
c      do i = 1,I3N
c           labels(i) = " "
c      enddo
c      charge = 0.D0
      do while (.true.)
       call VENUS_parallel_idop (myid,idop,theory, basis,
     *  natom, coordinates,labels,printName,qm_hessian,qm_grad,
     *  qm_energy,charge,stack,heap,global,qm_ini_file)
      end do
      return
      end

      Subroutine VENUS_parallel_idop (myid,idop,theory, basis,
     *  natom, coordinates,labels,printName,qm_hessian,qm_grad, 
     *  qm_energy,charge, stack,heap,global,qm_ini_file)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'

      integer myid, idop
      character*256 theory
      character*256 basis
      character*256 printName
      integer natom
      double precision coordinates (nda3)
      character*16 labels(nda)
      double precision qm_grad (nda3)
      double precision qm_hessian (nda3,nda3)
      double precision  qm_energy,charge
      integer stack,heap,global
      character*256 qm_ini_file

      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     * NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     * NRAX,NONB,NMO,NCRCO6

       call broadcast_data(idop)
c       write (6,*) 'myid=',myid,' Before VENUS_QM_init stack=',stack,
c     $             ' heap=',heap,' global=',global
       call flush(6)
       call flush(6)
       if (idop.eq.1) then
c        write(*,*) '_idop call to _init, theory, basis', theory,basis
        call VENUS_QM_init(theory, basis, natom, coordinates,
     *                 charge,labels, printName,stack,heap,global,
     *  qm_ini_file)
c        write(*,*) '_idop ret, theory, basis', theory,basis
       end if
c
       if (idop.eq.2) then
         call VENUS_QM_gradient(coordinates, qm_grad)
       end if
c
       if (idop.eq.3) then
         call VENUS_QM_gradient_energy(coordinates, qm_energy,
     *       qm_grad)
       end if
c
       if (idop.eq.4) then
         call VENUS_QM_energy(coordinates, qm_energy)
       end if
c
       if (idop.eq.6) then
c         write(*,*) 'test..'
         call VENUS_QM_hessian(coordinates, qm_hessian)
c       write(*,*) 'test...'
       end if
c
       if (idop.eq.5) then
        call VENUS_parallel_final()
        if (myid.ne.0) then
c         write (6,*) 'myid=',myid,' is killed'
         stop
        end if
       end if

      return
      end

C This Subroutine make calls to approriate functions of
C  VENUS_NWCHEM.f to calculate energu and gradient
C
      SUBROUTINE QMCALC(iqmflag)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'SIZES'
      COMMON/QPDOT/Q(3*NDA),PDOT(3*NDA)
      COMMON/PQDOT/P(3*NDA),QDOT(3*NDA),W(NDA)
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON /KEYWRD/ KEYWRD
      COMMON /LAST  / LAST
      COMMON/GAUSS/ICHRG,IMULTP,IAN(NDAyf),VGAUSS
c
      common/ghessb/hess(NDA3,NDA3),trad,rmin,dt,nstep,
     *  nhesup,nhessf,nip
      data H2KCAL/627.5095D0/,B2A/0.529177249D0/

csbo----------------------------------------------------------
      COMMON/QMMDAT/NQMMD,NLINKS,NLINKH(15),NBOUNMD(15),
     *NBOUNQM(15)
csbo----------------------------------------------------------

      COMMON /QMINFO/ qm_choice

      character*256 theory
      character*256 basis
      integer natom
      integer qm_choice
      double precision charge
      double precision coordinates (nda3)
      character*16 labels(nda)
      character*256 printName
      character*256 qm_ini_file
      integer stack,heap,global

      double precision QCHEM(nda3)
      dimension grad(nda3)

c       if (nmo.gt.0) then
c        i3nn=nmo*3
c       else
c        i3nn=natoms*3
c       endif
         if(nmo.gt.0)then
           natom=nmo
           i3nn=natom*3
         else
           natom=natoms
c           i3nn=i3n
           i3nn=natoms*3
         endif


       DO I=1,I3NN
         grad(I)= 0.0d0
       ENDDO
C
C Init QM coordinates
C    
       DO I=1,I3NN
         QCHEM(I)=Q(I)
c         WRITE (6,*) 'QCHEM=',QCHEM(I)
       ENDDO
c
       if ( iqmflag .eq. 0 ) then
C
C Compute QM ENERGY and GRADIENT together
C
       idop=3
c        write(*,*)'natom in grad  ',natom
c        write(*,*)'qchem in grad :',(qchem(i),i=1,i3nn)
       call VENUS_parallel_idop (myid,idop,theory, basis, natom,
     * QCHEM,labels,printName,hess,grad,escf,charge,stack,heap,global,
     * qm_ini_file)
c       write(6,*) 'ESCF=', escf
c       call flush(6)
       VGAUSS=ESCF*C1*h2kcal

c Making the Value of PDot zero when all atoms are to be treated by Nwchem
c so that it doesn't add ay value from Venus
        DO I=1,I3NN
           PDOT(I) = 0.d0
        ENDDO

        DO I=1,I3NN
            PDOT(I)=PDOT(I)+GRAD(I)*C1*H2KCAL/B2A
c            WRITE (6,*) 'PDOT=',PDOT(I)
        ENDDO
        else
       idop=6
c        write(*,*)'natom in hess ',natom
       call VENUS_parallel_idop (myid,idop,theory, basis, natom,
     * QCHEM,labels,printName,hess,grad,escf,charge,stack,heap,global,
     * qm_ini_file)
c        write(*,*)'out of hessian'
        do i=1,i3nn
          do j=1,i3nn
            hess(i,j)=hess(i,j)*c1*h2kcal/(b2a*b2a)
          enddo
        enddo
        endif
C
      return
      END
