      SUBROUTINE FMTRX(NATOM,NDIS,I3N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      logical ngflag
C
C         EVALUATE THE FORCE CONSTANT MATRIX BY DIFFERENCING
C         THE GRADIENT OF THE POTENTIAL ENERGY FUNCTION.
C
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/ARRAYS/A(NDA3yf,NDA3yf),DA(NDA3),B(NDA3yf,NDA3yf),DB(NDA3)
      COMMON/EIGVL/EIG(NDA3yf)
      COMMON/RSTART/HINC,NPTS
      COMMON/FR2/DG(NDA3,2),DIM(NDA3)
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      COMMON/FORCES/NATOMS,I3NS,NST,NM,NB,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/GAUSS2/GAUHES(NDIHE)
      COMMON/GAUSS/ICHRG,IMULTP,IAN(NDAyf),VGAUSS
      COMMON/FREQFLAG/IARB
c
      common/ghessb/FA(NDA3,NDA3),trad,rmin,dt,nstep,
     *  nhesup,nhessf,nip
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
c
c      DIMENSION DIST(2),GZ(NDA3),FA(NDA3yf,NDA3yf)
      DIMENSION DIST(2),GZ(NDA3)
      DATA PT5/0.5D0/
  100 FORMAT(10X,'***** CALCULATION OF FORCE CONSTANT  *****')
  102 FORMAT(///)
  103 FORMAT(15X,' FORCE CONSTANT MATRIX'/)
  104 FORMAT(/10X,'SYMMETRIZED CARTESIAN FORCE CONSTANT MATRIX'/)
c
c   This is to determine whether the gaussian will calculate the force
c   constant matrix or use numerical second derivative
c
c       Uncomment for G98
      if(nmo.lt.0.or.(nmo.gt.0.and.ndis.eq.0.and.natom.eq.nmo))then
        ngflag=.true.
      else
        ngflag=.false.
      endif
C
C         SET PRINT FILE FOR EIGOUT
C
      IP=77
      DIST(1)=HINC
      DIST(2)=-HINC
      DO I=1,I3N
         GZ(I)=-PDOT(I+3*NDIS)
      ENDDO
C
C         INITIALIZE SOME ARRAYS
C
      WRITE(ip,100)
C
C              GETTING FORCE MATRIX FROM GAUSSIAN 
c  For NMO < 0, this is entirely from Gaussian Hessian matrix.
c  Otherwise use the difference method to get the second derivative.
C
      IF (ngflag) THEN
c         WRITE(6,*)'            GETTING FORCE MATRIX FROM NWCHEM'
          call qmcalc(1)
c          write(*,*)'outside hessian matrix from nwchem'
c          do i = 1, i3n
c             write(*,*)(fa(i,j),j=1,i3n)
c          enddo
      ELSE
C
C             GETTING FORCE MATRIX FROM ANALYTICAL POTENTIALS
C
         IF (NPTS.GT.2.or.npts.le.0) NPTS=2
         DO I=1,I3N
            DO J=1,NPTS
               Q(I+3*NDIS)=Q(I+3*NDIS)+DIST(J)
C
C              DISPLACE COORDINATE AND CALCULATE GRADIENT INTO DG
C
               CALL DVDQ
               Q(I+3*NDIS)=Q(I+3*NDIS)-DIST(J)
               DO K=1,I3N
                  DG(K,J)=-PDOT(K+3*NDIS)
               ENDDO
            ENDDO
C
C              TWO POINT DIFFERENCE FORMULA
C
            IF (NPTS.NE.1) THEN
               DO K=1,I3N
                  FA(K,I)=(DG(K,1)-DG(K,2))*PT5/HINC
               ENDDO
            ELSE
C
C              SIMPLE ONE POINT FORMULA
C
               DO K=1,I3N
                  FA(K,I)=(DG(K,1)-GZ(K))/HINC
               ENDDO
            ENDIF
         ENDDO
C
      ENDIF
c
c  FA(I,J) is the hessian matrix necessary for hessian-based integrator.
c  so it stops here when hessian based integrator calls this routine.
c
c      write(*,*) 'nhessf 1 ',nhessf
      if(nhessf.eq.1) return
c      write(*,*) 'nhessf 2 ',nhessf
C
      DO I=1,I3N
         EIG(I)=0.0D0
      ENDDO
      WRITE(ip,102)
      WRITE(ip,103)
      CALL EIGOUT(FA,I3N,IP)
      DO I=1,I3N
         DO J=1,I
            IF (ngflag) THEN
               DU=FA(I,J)
               DL=0.0D0
            ELSE
               DU=PT5*(FA(I,J)+FA(J,I))
               DL=PT5*(FA(I,J)-FA(J,I))
            ENDIF
            A(I,J)=DU
            A(J,I)=DL
         ENDDO
      ENDDO
      DO I=1,I3N
         A(I,I)=FA(I,I)
      ENDDO
      IF (NSELT.EQ.-2.and.nacta.ne.2) RETURN
      WRITE(ip,102)
      WRITE(ip,104)
      CALL EIGOUT(A,I3N,IP)
C
C         CALCULATE ARRAY DIM(150) USED FOR MASS-WEIGHTING
C         CONVERT TO MASS WEIGHTED COORDINATES AND CALCULATE THE NORMAL
C         MODES AND THE SPECTROSCOPIC FREQUENCIES.
C
      K=0
      DO I=1,NATOM
         DO J=1,3
            K=K+1
            DIM(K)=1.D0/SQRT(W(I+NDIS))
         ENDDO
      ENDDO
      CALL FGMTRX(I3N)
      IF (NSELT.EQ.-1) STOP
      RETURN
223   continue
      stop
      END
