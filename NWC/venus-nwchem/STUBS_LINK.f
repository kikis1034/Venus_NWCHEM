c A subroutine used to Intialize the parallel enviornment for NWChem
      subroutine QM_parallel_init (myid)
      implicit none
      integer myid, nnodes1
      return
      end

c A subroutine to finalize the Nwchem and parallel enviornment
      Subroutine QM_parallel_final()
      implicit none
      return
      end

c A subroutine to broadcast data
      subroutine broadcast_data(cdata)
      implicit none
      integer cdata
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
      return
      end

      Subroutine VENUS_QM_hessian(coordinates, nw_hessian) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      external evaluateHessian
      logical evaluateHessian
      double precision coordinates (nda3)
      double precision nw_hessian (nda3,nda3)
      logical ans
      return
      end

