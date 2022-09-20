c A subroutine used to Intialize the parallel environment for NWChem
      subroutine VENUS_parallel_init (myid,nnodes)
      implicit none
      integer myid, nnodes
      call net_init(myid, nnodes)
      write (6,*) 'Parallel env init with myid=',myid,
     *            ' nnodes1=',nnodes
      return
      end

c A subroutine to finalize the Nwchem and parallel enviornment
      Subroutine VENUS_parallel_final()
      implicit none
      call net_final()
      write(6,*)'Parallel venus finalized'
      return
      end

c A subroutine to broadcast data
      subroutine broadcast_data(cdata)
      implicit none
      integer cdata
      call net_broadcast(cdata)
      return
      end

C A subroutine to intialize nwchem with the parameters.
      Subroutine VENUS_QM_init(theory, basis, natom, coordinates,
     *                    charge, labels,printName,stack,heap,global,
     *                    nwchem_ini_file)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      character*256 theory
      character*256 basis
      integer natom
      double precision coordinates (nda3)
      double precision charge
      character*16 labels(nda)
      character*256 printName
      integer stack,heap,global
      character*256 nwchem_ini_file 
      
      write (6,*) 'VENUS_QM_init:stack=',stack,' heap=',heap,
     $            ' global=',global 

      call initialize_qm(theory, basis, natom, coordinates, labels
     *                          ,charge, printName,stack,heap,global,
     *                          nwchem_ini_file) 
      return
      end


C A subroutine to obtain energy value from Nwchem.
      Subroutine VENUS_QM_energy(coordinates,nw_energy)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      external evaluateObjective 
      logical  evaluateObjective
      double precision coordinates (nda3)
      double precision nw_grad (nda3)
      double precision  nw_energy
      logical ans
      ans = .false.
      ans = evaluateObjective(coordinates, nw_energy)
      if (.not.ans) write (6,*) 'VENUS_QM_energy failed'
      call flush(6)
      return
      end

c A subroutine to obtain the value of gradient from Nwchem
      Subroutine VENUS_QM_gradient(coordinates, nw_grad)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      external evaluateGradient
      logical  evaluateGradient
      double precision coordinates (nda3)
      double precision nw_grad (nda3)
      logical ans
      ans = .false.
      ans = evaluateGradient(coordinates,nw_grad)
      if (.not.ans) write (6,*) 'VENUS_QM_gradient failed'
      return
      end

C A subroutine to obtaine both the value of gradient 
c and energy from Nwchem.
      Subroutine VENUS_QM_gradient_energy(coordinates,nw_energy,
     *                                           nw_grad)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      external evaluateObjectiveAndGradient
      logical evaluateObjectiveAndGradient
      double precision coordinates (nda3)
      double precision nw_grad (nda3)
      logical ans
      ans = .false.
      ans = evaluateObjectiveAndGradient(coordinates,nw_energy ,nw_grad)
      if (.not.ans) write (6,*) 'VENUS_QM_gradient_energy failed'
      return
      end

      Subroutine VENUS_QM_hessian(coordinates, nw_hessian) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      external evaluateHessian
      logical evaluateHessian
      COMMON/FORCES/NATOMS,I3N,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
c
      double precision coordinates (nda3),hessian(nda3*nda3)
      double precision nw_hessian (nda3,nda3)
      logical ans
c
       if (nmo.gt.0) then
        i3nn=nmo*3
       else
        i3nn=natoms*3
       endif
c
      ans = .false.
c      ans = evaluateHessian(coordinates,nw_hessian)
      ans = evaluateHessian(coordinates,hessian)
c  modified - raj
        write(53,*)'nw link hessian full triangle'
        write(53,*)(hessian(i),i=1,i3nn*i3nn)
        k = 0
        do i = 1,i3nn
        do j = 1,i3nn
            k=k+1
            nw_hessian(j,i)=hessian(k)
        enddo
        enddo
      if (.not.ans) write (6,*) 'VENUS_QM_hessian failed'
      return
      end

