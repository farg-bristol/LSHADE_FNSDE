      SUBROUTINE perfanal(func,dims,bound,nglobal,nsols,accuracy)

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

      INTEGER, INTENT(IN):: func, dims
      DOUBLE PRECISION, INTENT(IN):: accuracy,bound(dims,2)
      INTEGER, INTENT(OUT):: nglobal, nsols

      INTEGER:: ii, junk, particles, bestid, nn, k, nfeas, dd
      INTEGER:: numberopts
      DOUBLE PRECISION:: rad, gopt, dist, consaccuracy, avobj
      DOUBLE PRECISION:: maxspan, span, influencerad, optimal
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: fitness, &
     &   fsolutions, constraints
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: location, &
     &   xsolutions, sorting, fparticle
      LOGICAL:: found

!     GET INFLUENCE RADIUS AND ACCURACY
      rad=influencerad(func)
      consaccuracy=1.0E-4
      gopt=optimal(func)
      nglobal=numberopts(func)

!     GET MAXIMUM SPAN
      maxspan=-10000000.0d0
      DO dd=1,dims
        span=bound(dd,2)-bound(dd,1)
        IF(span.GE.maxspan)maxspan=span
      END DO

!     READ DATA
      OPEN(5,FILE='solution')
      READ(5,*)particles,junk
      ALLOCATE(location(particles,dims))
      ALLOCATE(fitness(particles))
      ALLOCATE(constraints(particles))
      ALLOCATE(fparticle(particles,2))
      avobj=0.0d0
      nfeas=0
      DO ii=1,particles
        READ(5,*)location(ii,:),fitness(ii),constraints(ii)
        fparticle(ii,1)=fitness(ii)
        fparticle(ii,2)=constraints(ii)
      END DO
      CLOSE(5)

!     SORT DATA
      ALLOCATE(sorting(particles,dims+1))
      sorting(:,1:dims)=location(:,1:dims)
      sorting(:,1+dims)=constraints(:)
      CALL sort(fitness, sorting, particles, dims+1)
      location(:,1:dims)=sorting(:,1:dims)
      constraints(:)=sorting(:,1+dims)
      DEALLOCATE(sorting)

!     GET NUMBER OF SOLUTIONS BY ALGORITHM 1 IN CEC2013
      nsols=0
      ALLOCATE(xsolutions(particles,dims))
      ALLOCATE(fsolutions(particles))
      DO ii=1,particles
        bestid=MINLOC(fitness,1)
        found=.FALSE.
        IF(ABS(gopt-fitness(ii)).LE.accuracy)THEN
         IF(constraints(ii).LE.consaccuracy)THEN
           IF(ii.NE.1)THEN
            DO nn=1,nsols
             dist=0.0d0
             DO k=1,dims
               dist=dist+(xsolutions(nn,k)-location(ii,k))**2
             END DO
             dist=SQRT(dist)
             IF(dist.LE.rad)THEN
               found=.TRUE.
               EXIT
             END IF
            END DO
           END IF
           IF(.NOT.found)THEN
             nsols=nsols+1
             xsolutions(nsols,:)=location(ii,:)
             fsolutions(nsols)=fitness(ii)
           END IF
         END IF
        END IF
      END DO
      nsols=MIN(nsols,nglobal)

      END SUBROUTINE

      DOUBLE PRECISION FUNCTION influencerad(func)
      IMPLICIT NONE
      INTEGER:: func
      SELECT CASE(func)
      CASE (10)
        influencerad=0.05d0
      CASE (11,12)
        influencerad=0.05d0
      CASE (13,14,15,16,17,18)
        influencerad=0.05d0
      CASE DEFAULT
        influencerad=0.5d0
      END SELECT
      END FUNCTION

      DOUBLE PRECISION FUNCTION optimal(func)
      IMPLICIT NONE
      INTEGER:: func,i,cons,dims
      DOUBLE PRECISION:: n
      SELECT CASE(func)
      CASE(1)
        optimal=1.0d0
      CASE(2)
        optimal=1.0d0
      CASE(3)
        optimal=1.6d0
      CASE(4)
        optimal=1.0d0
      CASE(5)
        optimal=1.729843561973525d0
      CASE(6)
        optimal=2.27272727272727d0
      CASE(7)
        optimal=1.0d0
      CASE(8)
        optimal=1.393589488353574d0
      CASE(9)
        optimal=1.826836992013024d0
      CASE(10)
        optimal=0.875d0
      CASE(11,12)
        optimal=1.0d0
      CASE(13,14)
        optimal=(10.0d0*1.0d0)-(5.0d0*1.0d0*SQRT(2.0d0))
      CASE(15,16)
        optimal=(10.0d0*2.0d0)-(5.0d0*2.0d0*SQRT(2.0d0))
      CASE(17)
        optimal=(10.0d0*3.0d0)-(5.0d0*3.0d0*SQRT(2.0d0))
      CASE(18)
        optimal=(10.0d0*5.0d0)-(5.0d0*5.0d0*SQRT(2.0d0))
      CASE(30:39)
        optimal=1.0d0
      CASE(40:49)
        OPEN(1,FILE='nichingcons_cases.dat')
        DO i=1,func
          READ(1,*)
        END DO
        READ(1,*)i,dims,cons
        CLOSE(1)
        n=DBLE(dims)
        optimal=n*n*(n*n+2*n-2)/((n*n+n-1)*(n*n-n+1))
      CASE(50:59)
        OPEN(1,FILE='nichingcons_cases.dat')
        DO i=1,func
          READ(1,*)
        END DO
        READ(1,*)i,dims,cons
        CLOSE(1)
        n=DBLE(dims)
        optimal=n**3*(n**3+4*n*n-2*n-8)/(n**6-2*n**4+4*n**3+7*n**2-20*n+8)
      CASE(60:69)
        OPEN(1,FILE='nichingcons_cases.dat')
        DO i=1,func
          READ(1,*)
        END DO
        READ(1,*)i,dims,cons
        CLOSE(1)
        n=DBLE(dims)
        optimal=n**3*(n**4+6*n**3-24*n-24)/(n**7-3*n**5+6*n**4+34*n**3-80*n*n-60*n+96)
      END SELECT
      END FUNCTION

      INTEGER FUNCTION numberopts(func)
      IMPLICIT NONE
      INTEGER, INTENT(IN):: func
      SELECT CASE(func)
      CASE(1)
        numberopts=2
      CASE(2)
        numberopts=2
      CASE(3)
        numberopts=4
      CASE(4)
        numberopts=2
      CASE(5)
        numberopts=8
      CASE(6)
        numberopts=32
      CASE(7)
        numberopts=2
      CASE(8)
        numberopts=8
      CASE(9)
        numberopts=32
      CASE(10)
        numberopts=10
      CASE(11,12)
        numberopts=4
      CASE(13)
        numberopts=2
      CASE(14)
        numberopts=10
      CASE(15)
        numberopts=8
      CASE(16)
        numberopts=24
      CASE(17)
        numberopts=16
      CASE(18)
        numberopts=64
      CASE(30:39)
        numberopts=2
      CASE(40:49)
        numberopts=4
      CASE(50:59)
        numberopts=8
      CASE(60:69)
        numberopts=16
      END SELECT
      END FUNCTION


      SUBROUTINE optimalx(func,nopt,dims,xopt)
      IMPLICIT NONE
      INTEGER, INTENT(IN):: func,nopt,dims
      DOUBLE PRECISION, INTENT(OUT)::xopt(nopt,dims)
      INTEGER::i
      DOUBLE PRECISION:: n, x1, x2
      SELECT CASE(func)
      CASE(1)
        xopt(1,1)=1.0d0
        xopt(2,1)=-1.0d0
      CASE(2)
        xopt(:,1)=0.0d0
        xopt(1,2)=1.0d0
        xopt(2,2)=-1.0d0
      CASE(3)
        xopt(1,1)=SQRT(0.8d0)
        xopt(1,2)=SQRT(0.8d0)
        xopt(2,1)=-SQRT(0.8d0)
        xopt(2,2)=SQRT(0.8d0)
        xopt(3,1)=SQRT(0.8d0)
        xopt(3,2)=-SQRT(0.8d0)
        xopt(4,1)=-SQRT(0.8d0)
        xopt(4,2)=-SQRT(0.8d0)
      CASE(4)
        xopt(1,1:4)=0.0d0
        xopt(1,5)=1.0d0
        xopt(2,1:4)=0.0d0
        xopt(2,5)=-1.0d0
      CASE(5)
        xopt(:,3)=0.0d0
        xopt(:,4)=0.0d0
        xopt(1,1)= 0.629371992938506d0
        xopt(1,2)= 0.645108721639740d0
        xopt(1,5)= 0.957898321192014d0
        xopt(2,1)=-0.629371992938506d0
        xopt(2,2)= 0.645108721639740d0
        xopt(2,5)= 0.957898321192014d0
        xopt(3,1)= 0.629371992938506d0
        xopt(3,2)=-0.645108721639740d0
        xopt(3,5)= 0.957898321192014d0
        xopt(4,1)=-0.629371992938506d0
        xopt(4,2)=-0.645108721639740d0
        xopt(4,5)= 0.957898321192014d0
        xopt(5,1)= 0.629371992938506d0
        xopt(5,2)= 0.645108721639740d0
        xopt(5,5)=-0.957898321192014d0
        xopt(6,1)=-0.629371992938506d0
        xopt(6,2)= 0.645108721639740d0
        xopt(6,5)=-0.957898321192014d0
        xopt(7,1)= 0.629371992938506d0
        xopt(7,2)=-0.645108721639740d0
        xopt(7,5)=-0.957898321192014d0
        xopt(8,1)=-0.629371992938506d0
        xopt(8,2)=-0.645108721639740d0
        xopt(8,5)=-0.957898321192014d0
      CASE(6)
        xopt(1,1)=SQRT(0.4545454545454545454545d0)
        xopt(1,2)=SQRT(0.4545454545454545454545d0)
        xopt(1,3)=SQRT(0.4545454545454545454545d0)
        xopt(1,4)=SQRT(0.4545454545454545454545d0)
        xopt(1,5)=SQRT(0.4545454545454545454545d0)
        xopt(2,:)=xopt(1,:)
        xopt(2,1)=-xopt(1,1)
        xopt(3:4,:)=xopt(1:2,:)
        xopt(3:4,2)=-xopt(1:2,2)
        xopt(5:8,:)=xopt(1:4,:)
        xopt(5:8,3)=-xopt(1:4,3)
        xopt(9:16,:)=xopt(1:8,:)
        xopt(9:16,4)=-xopt(1:8,4)
        xopt(17:32,:)=xopt(1:16,:)
        xopt(17:32,5)=-xopt(1:16,5)
      CASE(7)
        xopt(1:2,1:9)=0.0d0
        xopt(1,10)=1.0d0
        xopt(2,10)=-1.0d0
      CASE(8)
        xopt(1,1)=0.442990013929914d0
        xopt(1,2)=0.455649419140285d0
        xopt(1,3:9)=0.0d0
        xopt(1,10)=0.994853226737024d0
        xopt(2,:)=xopt(1,:)
        xopt(2,1)=-xopt(1,1)
        xopt(3:4,:)=xopt(1:2,:)
        xopt(3:4,2)=-xopt(1:2,2)
        xopt(5:8,:)=xopt(1:4,:)
        xopt(5:8,10)=-xopt(1:4,10)
      CASE(9)
        xopt(1,1)=0.462013517275990d0
        xopt(1,2)=0.468626768691906d0
        xopt(1,3)=0.476441515992551d0
        xopt(1,4)=0.485653252797998d0
        xopt(1,5:9)=0.0d0
        xopt(1,10)=0.964838770685610d0
        xopt(2,:)=xopt(1,:)
        xopt(2,1)=-xopt(1,1)
        xopt(3:4,:)=xopt(1:2,:)
        xopt(3:4,2)=-xopt(1:2,2)
        xopt(5:8,:)=xopt(1:4,:)
        xopt(5:8,3)=-xopt(1:4,3)
        xopt(9:16,:)=xopt(1:8,:)
        xopt(9:16,4)=-xopt(1:8,4)
        xopt(17:32,:)=xopt(1:16,:)
        xopt(17:32,10)=-xopt(1:16,10)
      CASE(10)
        DO i=1,10
          xopt(i,1)=0.05d0+DBLE(i-1)*0.1d0
        END DO
      CASE(11,12)
        xopt(1,1)=3.0d0
        xopt(1,2)=2.0d0
        xopt(2,1)=-2.805118d0
        xopt(2,2)=3.131312d0
        xopt(3,1)=-3.779310d0
        xopt(3,2)=-3.283186d0
        xopt(4,1)=3.584428d0
        xopt(4,2)=-1.848126d0
      CASE(13)
        xopt(1,1)=0.375d0
        xopt(2,1)=0.625d0
      CASE(14)
        xopt(1,1)=0.075d0
        xopt(2,1)=0.125d0
        xopt(3,1)=0.275d0
        xopt(4,1)=0.325d0
        xopt(5,1)=0.475d0
        xopt(6,1)=0.525d0
        xopt(7,1)=0.675d0
        xopt(8,1)=0.725d0
        xopt(9,1)=0.875d0
        xopt(10,1)=0.925d0
      CASE(15)
        xopt(1:4,1)=0.375d0
        xopt(1,2)=0.1875d0
        xopt(2,2)=0.3125d0
        xopt(3,2)=0.6875d0
        xopt(4,2)=0.8125d0
        xopt(5:8,1)=0.625d0
        xopt(5,2)=0.1875d0
        xopt(6,2)=0.3125d0
        xopt(7,2)=0.6875d0
        xopt(8,2)=0.8125d0
      CASE(16)
        xopt(1:6,1)=0.1875d0
        xopt(7:12,1)=0.3125d0
        xopt(13:18,1)=0.6875d0
        xopt(19:24,1)=0.8125d0
        xopt(1,2)=0.125000d0
        xopt(2,2)=5.0d0/24.0d0
        xopt(3,2)=11.0d0/24.0d0
        xopt(4,2)=13.0d0/24.0d0
        xopt(5,2)=19.0d0/24.0d0
        xopt(6,2)=0.875000d0
        xopt(7,2)=0.125000d0
        xopt(8,2)=5.0d0/24.0d0
        xopt(9,2)=11.0d0/24.0d0
        xopt(10,2)=13.0d0/24.0d0
        xopt(11,2)=19.0d0/24.0d0
        xopt(12,2)=0.875000d0
        xopt(13,2)=0.125000d0
        xopt(14,2)=5.0d0/24.0d0
        xopt(15,2)=11.0d0/24.0d0
        xopt(16,2)=13.0d0/24.0d0
        xopt(17,2)=19.0d0/24.0d0
        xopt(18,2)=0.875000d0
        xopt(19,2)=0.125000d0
        xopt(20,2)=5.0d0/24.0d0
        xopt(21,2)=11.0d0/24.0d0
        xopt(22,2)=13.0d0/24.0d0
        xopt(23,2)=19.0d0/24.0d0
        xopt(24,2)=0.875000d0
      CASE(17)
        xopt(1:8,1)=0.375d0
        xopt(9:16,1)=0.625d0
        xopt(1:4,2)=0.375d0
        xopt(5:8,2)=0.625d0
        xopt(9:12,2)=0.375d0
        xopt(13:16,2)=0.625d0
        xopt(1,3)=0.1875d0
        xopt(5,3)=0.1875d0
        xopt(9,3)=0.1875d0
        xopt(13,3)=0.1875d0
        xopt(2,3)=0.3125d0
        xopt(6,3)=0.3125d0
        xopt(10,3)=0.3125d0
        xopt(14,3)=0.3125d0
        xopt(3,3)=0.6875d0
        xopt(7,3)=0.6875d0
        xopt(11,3)=0.6875d0
        xopt(15,3)=0.6875d0
        xopt(4,3)=0.8125d0
        xopt(8,3)=0.8125d0
        xopt(12,3)=0.8125d0
        xopt(16,3)=0.8125d0
      CASE(18)
        xopt(1:32,1)=0.375d0
        xopt(33:64,1)=0.625d0
        xopt(1:16,2)=0.375d0
        xopt(17:32,2)=0.625d0
        xopt(33:48,2)=0.375d0
        xopt(49:64,2)=0.625d0
        xopt(1:8,3)=0.375d0
        xopt(9:16,3)=0.625d0
        xopt(17:24,3)=0.375d0
        xopt(25:32,3)=0.625d0
        xopt(33:40,3)=0.375d0
        xopt(41:48,3)=0.625d0
        xopt(49:56,3)=0.375d0
        xopt(57:64,3)=0.625d0
        xopt(1:4,4)=0.375d0
        xopt(5:8,4)=0.625d0
        xopt(9:12,4)=0.375d0
        xopt(13:16,4)=0.625d0
        xopt(17:20,4)=0.375d0
        xopt(21:24,4)=0.625d0
        xopt(25:28,4)=0.375d0
        xopt(29:32,4)=0.625d0
        xopt(33:36,4)=0.375d0
        xopt(37:40,4)=0.625d0
        xopt(41:44,4)=0.375d0
        xopt(45:48,4)=0.625d0
        xopt(49:52,4)=0.375d0
        xopt(53:56,4)=0.625d0
        xopt(57:60,4)=0.375d0
        xopt(61:64,4)=0.625d0
        do i=1,16
          xopt(1+(i-1)*4,5)=0.1875d0
          xopt(2+(i-1)*4,5)=0.3125d0
          xopt(3+(i-1)*4,5)=0.6875d0
          xopt(4+(i-1)*4,5)=0.8125d0
        end do
      CASE(30:39)
        xopt(:,:)=0.0d0
        xopt(1,dims)=1.0d0
        xopt(2,dims)=-1.0d0
      CASE(40:49)
        n=DBLE(dims)
        xopt(:,:)=0.0d0
        x1=SQRT(n*n*(2*n-1)/((n*n+n-1)*(n*n-n+1)))
        x2=-x1
        xopt(1,1)=x1
        xopt(2,1)=x1
        xopt(3,1)=x2
        xopt(4,1)=x2
        x1=SQRT(n*n*(n*n-1)/((n*n+n-1)*(n*n-n+1)))
        x2=-x1
        xopt(1,dims)=x1
        xopt(2,dims)=x2
        xopt(3,dims)=x1
        xopt(4,dims)=x2
      CASE(50:59)
        n=DBLE(dims)
        xopt(:,:)=0.0d0
        x1=SQRT(n*n*(2*n**3-n**2+4*n-8)/(n**6-2*n**4+4*n**3+7*n**2-20*n+8))
        x2=-x1
        xopt(1,1)=x1
        xopt(2,1)=x1
        xopt(3,1)=x1
        xopt(4,1)=x1
        xopt(5,1)=x2
        xopt(6,1)=x2
        xopt(7,1)=x2
        xopt(8,1)=x2
        x1=SQRT(n*n*(2*n**3+n**2-6*n+4)/(n**6-2*n**4+4*n**3+7*n**2-20*n+8))
        x2=-x1
        xopt(1,2)=x1
        xopt(2,2)=x1
        xopt(3,2)=x2
        xopt(4,2)=x2
        xopt(5,2)=x1
        xopt(6,2)=x1
        xopt(7,2)=x2
        xopt(8,2)=x2
        x1=SQRT(n*n*(n**4-2*n**2-6*n+4)/(n**6-2*n**4+4*n**3+7*n**2-20*n+8))
        x2=-x1
        xopt(1,dims)=x1
        xopt(2,dims)=x2
        xopt(3,dims)=x1
        xopt(4,dims)=x2
        xopt(5,dims)=x1
        xopt(6,dims)=x2
        xopt(7,dims)=x1
        xopt(8,dims)=x2
      CASE(60:69)
        n=DBLE(dims)
        xopt(:,:)=0.0d0
        x1=SQRT(n*n*(2*n**4-n**3+12*n**2-24*n-24)/&
     &         (n**7-3*n**5+6*n**4+34*n**3-80*n**2-60*n+96))
        x2=-x1
        xopt(1,1)=x1
        xopt(2,1)=x1
        xopt(3,1)=x1
        xopt(4,1)=x1
        xopt(5,1)=x1
        xopt(6,1)=x1
        xopt(7,1)=x1
        xopt(8,1)=x1
        xopt(9,1)=x2
        xopt(10,1)=x2
        xopt(11,1)=x2
        xopt(12,1)=x2
        xopt(13,1)=x2
        xopt(14,1)=x2
        xopt(15,1)=x2
        xopt(16,1)=x2
        x1=SQRT(n*n*(2*n**4+n**3-2*n**2-24)/&
     &         (n**7-3*n**5+6*n**4+34*n**3-80*n**2-60*n+96))
        x2=-x1
        xopt(1,2)=x1
        xopt(2,2)=x1
        xopt(3,2)=x1
        xopt(4,2)=x1
        xopt(5,2)=x2
        xopt(6,2)=x2
        xopt(7,2)=x2
        xopt(8,2)=x2
        xopt(9,2)=x1
        xopt(10,2)=x1
        xopt(11,2)=x1
        xopt(12,2)=x1
        xopt(13,2)=x2
        xopt(14,2)=x2
        xopt(15,2)=x2
        xopt(16,2)=x2
        x1=SQRT(n*n*(2*n**4+3*n**3-12*n**2-4*n+24)/&
     &         (n**7-3*n**5+6*n**4+34*n**3-80*n**2-60*n+96))
        x2=-x1
        xopt(1,3)=x1
        xopt(2,3)=x1
        xopt(3,3)=x2
        xopt(4,3)=x2
        xopt(5,3)=x1
        xopt(6,3)=x1
        xopt(7,3)=x2
        xopt(8,3)=x2
        xopt(9,3)=x1
        xopt(10,3)=x1
        xopt(11,3)=x2
        xopt(12,3)=x2
        xopt(13,3)=x1
        xopt(14,3)=x1
        xopt(15,3)=x2
        xopt(16,3)=x2
        x1=SQRT(n*n*(n**5-3*n**3-22*n**2+4*n+24)/&
     &         (n**7-3*n**5+6*n**4+34*n**3-80*n**2-60*n+96))
        x2=-x1
        xopt(1,dims)=x1
        xopt(2,dims)=x2
        xopt(3,dims)=x1
        xopt(4,dims)=x2
        xopt(5,dims)=x1
        xopt(6,dims)=x2
        xopt(7,dims)=x1
        xopt(8,dims)=x2
        xopt(9,dims)=x1
        xopt(10,dims)=x2
        xopt(11,dims)=x1
        xopt(12,dims)=x2
        xopt(13,dims)=x1
        xopt(14,dims)=x2
        xopt(15,dims)=x1
        xopt(16,dims)=x2
      END SELECT
      END SUBROUTINE


      INTEGER FUNCTION numberoflines(fid)
      IMPLICIT NONE
      INTEGER, INTENT(IN):: fid
      numberoflines = 0 
      DO 
          READ (fid,*, END=10) 
          numberoflines = numberoflines + 1 
      END DO 
10    CONTINUE
      END FUNCTION

















      SUBROUTINE sort(sorting_vector, arrange_matrix, N, D)

!      Bubble sort subroutine that sorts based on the 'sorting_vector'.
!      Moves the rows of the two arrange matrices to shift with 
!      their corresponding entry in vector.
!      N is sorting vector size and height of matrices. D is width of matrices

      IMPLICIT NONE

      INTEGER, INTENT(IN):: N,D
      DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT):: sorting_vector
      DOUBLE PRECISION, DIMENSION(N,D), INTENT(INOUT):: arrange_matrix

      INTEGER::sorted,i, krow
      DOUBLE PRECISION::value_temp
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: a_row

      ALLOCATE(a_row(D))
      goto 25
      DO
      sorted=0
      DO i=2,N
      IF(sorting_vector(i)<sorting_vector(i-1)) THEN
      value_temp=sorting_vector(i-1)
      a_row(:)=arrange_matrix(i-1,:)

      sorting_vector(i-1)=sorting_vector(i)
      arrange_matrix(i-1,:)=arrange_matrix(i,:)
      sorting_vector(i)=value_temp
      arrange_matrix(i,:)=a_row(:)
      sorted=1
      END IF
      END DO

      IF(sorted==0) THEN
      EXIT
      END IF
      END DO

25    continue

      DO i=1,N
      krow=MINLOC(sorting_vector(i:N), dim=1 ) + i - 1
      value_temp=sorting_vector(i)
      a_row(:)=arrange_matrix(i,:)

      sorting_vector(i) = sorting_vector(krow)
      arrange_matrix(i,:) = arrange_matrix(krow,:)

      sorting_vector(krow) = value_temp
      arrange_matrix(krow,:) = a_row(:)

      END DO


      DEALLOCATE(a_row)

      END SUBROUTINE sort
