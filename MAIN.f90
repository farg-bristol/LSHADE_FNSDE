      PROGRAM main

!     NEIGHBOURHOOD-BASED SPECIATION DIFFERNETIAL EVOLUTION MULTIMODAL
!     OPTIMIZER WITH FEASIBLE SELECTION FOR CONSTRAINTS.
!
!     Copyright 2018 Daniel Poole
!
!     Permission is hereby granted, free of charge, to any person 
!     obtaining a copy of this software and associated documentation 
!     files (the "Software"), to deal in the Software without restriction, 
!     including without limitation the rights to use, copy, modify, merge, 
!     publish, distribute, sublicense, and/or sell copies of the Software, 
!     and to permit persons to whom the Software is furnished to do so, 
!     subject to the following conditions:
!
!     The above copyright notice and this permission notice shall be 
!     included in all copies or substantial portions of the Software.
!
!     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
!     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
!     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
!     BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
!     AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
!     IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
!     THE SOFTWARE.

      IMPLICIT NONE

      integer::dims,i,func,nruns,irun,nfuncs,cons,fmax,&
     &   nglobalreal,nglobal,j
      double precision::accuracy(5)
      double precision, dimension(:,:),allocatable:: bound
      double precision, dimension(:,:,:),allocatable:: peakratio, &
     &  successrate

      nruns=1!50
      nfuncs=18
      fmax=400000

      accuracy(1)=1.0E-1
      accuracy(2)=1.0E-2
      accuracy(3)=1.0E-3
      accuracy(4)=1.0E-4
      accuracy(5)=1.0E-5

      ALLOCATE(peakratio(nruns,nfuncs,5))
      ALLOCATE(successrate(nruns,nfuncs,5))
      PRINT*,''
      PRINT'(A)','---------------------------------------------------------'
      PRINT'(A)','Fun  Run  PR(E-1)   PR(E-2)   PR(E-3)   PR(E-4)   PR(E-5)'

      DO func=1,nfuncs

      PRINT'(A)','---------------------------------------------------------'

        DO irun=1,nruns

          peakratio(irun,func,:)=0.0d0
          successrate(irun,func,:)=0.0d0

!         PROBLEM DEFINITION
          OPEN(1,FILE='nichingcons_cases.dat')
          DO i=1,func
            READ(1,*)
          END DO
          READ(1,*)i,dims,cons
          CLOSE(1)

!         ALLOCATE DATA
          allocate(bound(dims,2))

!         DESIGN SPACE BOUNDARY
          CALL niching_func_bound_cons(func, dims, bound)

!         OPTIMIZER CALL
!          fmax=INT(2000*dims*SQRT(DBLE(numberopts(func,consy))))
          CALL de(func,dims,bound(:,1),bound(:,2),cons,fmax)

!         SOLUTION PERFORMANCE ANALYSIS
          DO j=1,5
           CALL perfanal(func,dims,bound,nglobalreal,nglobal,&
     &        accuracy(j))
           peakratio(irun,func,j)=DBLE(nglobal)/DBLE(nglobalreal)
           IF(nglobal.GE.nglobalreal)THEN
            successrate(irun,func,j)=1.0d0
           END IF
          END DO
          DEALLOCATE(bound)

!         WRITE INDIIVIDUAL RUN PR
          WRITE(6,240)func,irun,peakratio(irun,func,1:5)

!         CLEAR UP RUN
          CALL SYSTEM('rm solution')

        END DO

      END DO

      PRINT'(A)','---------------------------------------------------------'
      PRINT*,''

240   FORMAT(I3,I4,5F10.5)

      END PROGRAM
