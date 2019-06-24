      SUBROUTINE niching_func_cons(func,x,dims,cons,fitness,consfit)

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

      INTEGER, INTENT(IN):: func, dims, cons
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(dims):: x
      DOUBLE PRECISION, INTENT(OUT):: fitness
      DOUBLE PRECISION, INTENT(OUT):: consfit(cons)

      INTEGER:: i, heavi
      DOUBLE PRECISION:: pi, obj1, obj2
      DOUBLE PRECISION:: x1,x2,x3,x4,y1,y2,y3,y4,conseps
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: k

      pi=4.0d0*DATAN(1.0d0)
      conseps=1.0E-4

      SELECT CASE(func)
      CASE(10)

        fitness=-SIN(5.0d0*pi*x(1))**6+1.0d0
        consfit(1)=-COS(10.0d0*pi*x(1))

      CASE(11)

!       HIMMELBLAU CONSTRAINED BY QUADRATICS SUCH THAT ONE CONSTRAINT
!       TOUCHES TWO OPTIMAL POINTS, THEREFORE AT EACH OPTIMAL POINT,
!       TWO CONSTRAINTS ARE ACTIVE

        fitness=(x(1)**2+x(2)-11.0d0)**2+(x(1)+x(2)**2-7.0d0)**2+1.0d0

        x1=3.0d0
        y1=2.0d0
        x2=-2.805118d0
        y2=3.131312d0
        x3=-3.779310d0
        y3=-3.283186d0
        x4=3.584428d0
        y4=-1.848126d0

        consfit(1)=(x1*x2)-(x1*x(1))-(x2*x(1))+(y1*y2)-(y1*x(2))-(y2*x(2))+(x(1)**2)+(x(2)**2)
        consfit(2)=(x2*x3)-(x2*x(1))-(x3*x(1))+(y2*y3)-(y2*x(2))-(y3*x(2))+(x(1)**2)+(x(2)**2)
        consfit(3)=(x3*x4)-(x3*x(1))-(x4*x(1))+(y3*y4)-(y3*x(2))-(y4*x(2))+(x(1)**2)+(x(2)**2)
        consfit(4)=(x4*x1)-(x4*x(1))-(x1*x(1))+(y4*y1)-(y4*x(2))-(y1*x(2))+(x(1)**2)+(x(2)**2)

        consfit(:)=-consfit(:)

      CASE(12)

!       HIMMELBLAU CONSTRAINED BY QUADRATICS AND PLANES SUCH THAT
!       AT EACH OPTIMAL POINT THREE CONSTRAINTS ARE ACTIVE.
!       OPTIMAL POINTS AT SINGLE POINT WHERE CONS ARE ACTIVE

        fitness=(x(1)**2+x(2)-11.0d0)**2+(x(1)+x(2)**2-7.0d0)**2+1.0d0

        x1=3.0d0
        y1=2.0d0
        x2=-2.805118d0
        y2=3.131312d0
        x3=-3.779310d0
        y3=-3.283186d0
        x4=3.584428d0
        y4=-1.848126d0

        consfit(1)=(x1*x2)-(x1*x(1))-(x2*x(1))+(y1*y2)-(y1*x(2))-(y2*x(2))+(x(1)**2)+(x(2)**2)
        consfit(2)=(x2*x3)-(x2*x(1))-(x3*x(1))+(y2*y3)-(y2*x(2))-(y3*x(2))+(x(1)**2)+(x(2)**2)
        consfit(3)=(x3*x4)-(x3*x(1))-(x4*x(1))+(y3*y4)-(y3*x(2))-(y4*x(2))+(x(1)**2)+(x(2)**2)
        consfit(4)=(x4*x1)-(x4*x(1))-(x1*x(1))+(y4*y1)-(y4*x(2))-(y1*x(2))+(x(1)**2)+(x(2)**2)
        consfit(1:4)=-consfit(1:4)

        consfit(5)=(x(1)-x1)+(x(2)-y1)
        consfit(6)=-(x(1)-x2)+(x(2)-y2)
        consfit(7)=-(x(1)-x3)-(x(2)-y3)
        consfit(8)=(x(1)-x4)-(x(2)-y4)

      CASE(13,14,15,16,17,18)

!       MODULO-MODIFIED RASTRIGIN: DEB'S MODIFIED RASTRIGIN BUT FOR X<1, NO
!       X^2 UNDERLYING TREND - CAUSES PROBLEM TO HAVE GLOBAL AND LOCAL
!       MINIMA. CREATED BY D. POOLE.

        ALLOCATE(k(dims))

        SELECT CASE(func)
        CASE(13)
          k(1)=1.0d0
        CASE(14)
          k(1)=5.0d0
        CASE(15)
          k(1)=1.0d0
          k(2)=2.0d0
        CASE(16)
          k(1)=2.0d0
          k(2)=3.0d0
        CASE(17)
          k(1)=1.0d0
          k(2)=1.0d0
          k(3)=2.0d0
        CASE(18)
          k(1)=1.0d0
          k(2)=1.0d0
          k(3)=1.0d0
          k(4)=1.0d0
          k(5)=2.0d0
        END SELECT

        obj1=0.0d0
        obj2=0.0d0
        DO i=1,dims
          heavi=heaviside(x(i)-1.0d0)
          obj1=obj1+(10.0d0*(1.0d0+COS(2.0d0*pi*k(i)*x(i)))+&
     &               (2.0d0*k(i)*(x(i)-1.0d0)**2*heavi))

          obj2=obj2+(20.0d0*COS(4.0d0*pi*k(i)*x(i)))
        END DO
        fitness=obj1
        consfit(1)=obj2

      CASE(20)

        obj1=0.0d0
        obj2=0.0d0
        DO i=1,dims
          obj1=obj1+(10.0d0*COS(0.1d0*x(i)))
          obj2=obj2+(x(i)*COS(SQRT(ABS(x(i)))))
        END DO
        fitness=obj1
        consfit(1)=obj2
        consfit(1)=ABS(consfit(1))-conseps

      CASE(21)

        x1=x(1)
        x2=x(2)
        fitness=(x1**2+x1+x2**2+2.1d0*x2)+10.0d0*(1.0d0-COS(2.0d0*pi*x1))+&
     &                                    10.0d0*(1.0d0-COS(2.0d0*pi*x2))
        fitness=10.0d0*(1.0d0-COS(2.0d0*pi*x1))+&
     &                                    10.0d0*(1.0d0-COS(2.0d0*pi*x2))
        consfit(1)=-0.0d0

      CASE DEFAULT
        CALL CMMP(dims, cons, x, fitness, consfit)
      END SELECT


      CONTAINS

      INTEGER FUNCTION heaviside(x)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)::x
      IF(x.GE.0.0d0)THEN
        heaviside=1
      ELSE
        heaviside=0
      END IF
      END FUNCTION


      END SUBROUTINE




      SUBROUTINE niching_func_bound_cons(func, dims, bounds)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: func, dims
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(dims,2):: bounds

      SELECT CASE(func)
      CASE(10)
       bounds(:,1)=0.0d0
       bounds(:,2)=1.0d0
      CASE(11,12)
       bounds(:,1)=-6.0d0
       bounds(:,2)=6.0d0
      CASE(13,14,15,16,17,18)
       bounds(:,1)=0.0d0
       bounds(:,2)=2.0d0
      CASE(20)
       bounds(:,1)=-100.0d0
       bounds(:,2)=100.0d0
      CASE(21)
       bounds(:,1)=0.5d0
       bounds(:,2)=4.5d0
      CASE(30:39)
       bounds(:,1)=-DBLE(dims+1)
       bounds(:,2)=DBLE(dims+1)
      CASE DEFAULT
       bounds(:,1)=-DBLE(dims+1)
       bounds(:,2)=DBLE(dims+1)
      END SELECT

      END SUBROUTINE



      SUBROUTINE CMMP(n, J, x, f, g)
      IMPLICIT NONE
!     FROM DEB AND SAHA, GECCO2010
      INTEGER, INTENT(IN):: n, J
      DOUBLE PRECISION, INTENT(IN):: x(n)
      DOUBLE PRECISION, INTENT(OUT):: f, g(J)
!     NOTE: L not currently used
      INTEGER:: i, jj, ii, C
      INTEGER, ALLOCATABLE, DIMENSION(:):: Cvec

!     OBJECTIVE
      f=0.0d0
      DO i=1,n
        f=f+x(i)**2
      END DO

!     CONSTRAINTS
      ALLOCATE(Cvec(J*n))
      ii=0
      DO jj=1,J
       DO i=1,n
         ii=ii+1
         Cvec(ii)=i
       END DO
      END DO
      ii=1
      DO jj=1,J
        ii=ii-1
        g(jj)=0.0d0
        DO i=1,n
            ii=ii+1
            C=Cvec(ii)
            g(jj)=g(jj)+(DBLE(C)**2 * x(i)**2)
        END DO
        g(jj)=n**2 - g(jj)
      END DO

      END SUBROUTINE
