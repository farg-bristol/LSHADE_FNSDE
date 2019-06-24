      SUBROUTINE de(func,dims,xlower,xupper,cons,fmax)

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

!***********************************************************************
!     VARIABLE DECLARATIONS
!***********************************************************************

      INTEGER, INTENT(IN):: dims,func,cons,fmax
      DOUBLE PRECISION, INTENT(IN):: xlower(dims), xupper(dims)

      INTEGER:: t, i, np, tmax, nb, nbmax, m, p, k, feval, nbpmin, &
     &  nbp_used, tc, j, nsuc_cross, nsuc_mut, ninit, &
     &  ksuccess, numberopts, narchive, archive_yn, nmin, npcur
      INTEGER, ALLOCATABLE, DIMENSION(:):: nbp, isort
      INTEGER, ALLOCATABLE, DIMENSION(:,:):: ineighbour

      DOUBLE PRECISION:: seedradius, minimum, cross_pool, mut_pool
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: mutvector, &
     &   pparticle, mem_cross, mem_mut, success_cross, success_mut, &
     &   cur_cross, cur_mut
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: xparticle, &
     &   fparticle, cparticle, trialvector, euclidean, &
     &   xsort, csort, fsort, xarchive, fold, ftrial, ctrial
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:):: xneighbour, &
     &   fneighbour, cneighbour

!***********************************************************************
!     SET-UP SOLVER
!***********************************************************************

!     SEED
      CALL init_random_seed()

!     CONSTANTS
      CALL de_config(np, nbpmin, seedradius, cross_pool, mut_pool)
      np=INT(40*SQRT(DBLE(dims*numberopts(func,.TRUE.))))

!     MAXIMUM ITERATIONS
      ninit=np
      nmin=np/2
      feval=np
      npcur=np
      t=0
      DO
        t=t+1
        npcur=NINT((DBLE(nmin-ninit)/DBLE(fmax))*feval+ninit)
        feval=feval+npcur
        IF(feval.GE.fmax)EXIT
      END DO
      tmax=2*t
      tc=INT(FLOOR(0.4d0*tmax))
      nmin=nmin+1

!     SEED RADIUS FRACTION
      minimum=1.0E10
      DO k=1,dims
        IF((xupper(k)-xlower(k)).LT.minimum)THEN
          minimum=xupper(k)-xlower(k)
        END IF
      END DO
      seedradius=seedradius*minimum

!     MEMORY ALLOCATION
      ALLOCATE(xparticle(np,dims))
      ALLOCATE(fparticle(np,2))
      ALLOCATE(cparticle(np,cons))
      ALLOCATE(pparticle(np))
      ALLOCATE(fold(np,2))
      ALLOCATE(mutvector(dims))
      ALLOCATE(trialvector(np,dims))
      ALLOCATE(euclidean(np,np))
      ALLOCATE(ftrial(np,2))
      ALLOCATE(ctrial(np,cons))
      ALLOCATE(xsort(np,dims))
      ALLOCATE(fsort(np,2))
      ALLOCATE(csort(np,cons))
      ALLOCATE(isort(np))
      ALLOCATE(mem_cross(np))
      ALLOCATE(mem_mut(np))
      ALLOCATE(success_cross(np))
      ALLOCATE(success_mut(np))
      ALLOCATE(cur_cross(np))
      ALLOCATE(cur_mut(np))
      ALLOCATE(xarchive(np,dims))
      ALLOCATE(nbp(1))
      narchive=0
      feval=0
      archive_yn=0

!***********************************************************************
!     INITIALISE PARTICLES
!***********************************************************************

!     PARTICLE INITIALISATION
      CALL particleinitial(dims, np, xlower(:), xupper(:), &
     &     xparticle(:,:))

!     INITIAL OBJECTIVE AND CONSTRAINTS
      DO i=1,np
        CALL particleobjective(func,xparticle(i,:),dims,cons,&
     &     fparticle(i,:),cparticle(i,:),feval)
      END DO

!     NEED TO SET INITIAL VALUES FOR MEMORY
      mem_cross(:)=0.5d0
      mem_mut(:)=0.5d0
      ksuccess=1

!***********************************************************************
!     MAIN OPTIMISER LOOP
!***********************************************************************

!     START OPTIMISER LOOP
      t=0
      DO

!       ITERATION COUNTER
        t=t+1

!       POPULATION SIZE REDUCTION
        npcur=np
        np=NINT((DBLE(nmin-ninit)/DBLE(fmax))*feval+ninit)
        IF(np.NE.npcur)THEN
          CALL SDE_sort(npcur,dims,cons,xparticle(1:npcur,:),&
     &         fparticle(1:npcur,:),cparticle(1:npcur,:),&
     &         xsort(1:npcur,:),fsort(1:npcur,:),csort(1:npcur,:),&
     &         isort(1:npcur))
          xparticle(1:np,:)=xsort(1:np,:)
          fparticle(1:np,:)=fsort(1:np,:)
          cparticle(1:np,:)=csort(1:np,:)
        END IF

!       OLD VALUES
        fold(1:np,:)=fparticle(1:np,:)

!       SHADE, RESET SUCCESSFUL CROSSOVER AND MUTATION
        nsuc_cross=0
        nsuc_mut=0
        success_cross(:)=0.0d0
        success_mut(:)=0.0d0

!       NUMBER OF PARTICLES/POPULATIONS
!       NEIGHBOURHOOD SPECIES DE
        DEALLOCATE(nbp)
        nb=np/nbpmin
        ALLOCATE(nbp(nb))
        nbp(:)=nbpmin
        nbmax=nb

!***********************************************************************
!       GET EUCLIDEAN
!***********************************************************************

        CALL xeuc(np,dims,xparticle(1:np,:),euclidean(1:np,1:np))

!***********************************************************************
!       GENERATE SUB-POPULATIONS/NEIGHBOURHOODS
!***********************************************************************

!       ARRAYS FOR NEIGHBOURHOODS
        ALLOCATE(xneighbour(nb,MAXVAL(nbp),dims))
        ALLOCATE(fneighbour(nb,MAXVAL(nbp),2))
        ALLOCATE(cneighbour(nb,MAXVAL(nbp),cons))
        ALLOCATE(ineighbour(nb,MAXVAL(nbp)))

!       NEIGHBOURHOOD ASSIGNMENTS
        CALL neighbourhoods(np,dims,cons,xparticle(1:np,:),&
     &     fparticle(1:np,:),cparticle(1:np,:),nb,nbp,xneighbour, &
     &     fneighbour,cneighbour,ineighbour,euclidean(1:np,1:np))

        DO m=1,nb
          nbp_used=nbp(m)
          DO p=1,nbp_used
            i=ineighbour(m,p)

!***********************************************************************
!           L-SHADE - GET CROSSOVER AND MUTATION
!***********************************************************************

            CALL lshade_F_CR(np,mem_cross,mem_mut,cross_pool,mut_pool)
            cur_cross(i)=cross_pool
            cur_mut(i)=mut_pool

!***********************************************************************
!           GENERATE MUTATION VECTOR
!***********************************************************************

            CALL mutationvec(nbp(m),dims,cons,xneighbour(m,:,:),&
     &        fneighbour(m,:,:),cneighbour(m,:,:),p,mutvector,&
     &        mut_pool,narchive,xarchive,archive_yn,np)

!***********************************************************************
!           GENERATE TRIAL VECTORS
!***********************************************************************

!           TRIAL VECTOR GENERATION
            CALL trialvec(dims,xneighbour(m,p,:),cross_pool,&
     &                    mutvector,trialvector(i,:))

!           MAKE SURE TRIAL IS IN RANGE
            CALL boundaries(dims,xupper,xlower,trialvector(i,:))

!***********************************************************************
!           FUNCTION EVALUATIONS
!***********************************************************************

!           TRIAL VECTOR EVALUATION
            CALL particleobjective(func,trialvector(i,:),dims,cons,&
     &        ftrial(i,:),ctrial(i,:),feval)

!           FOR NSDE, REPLACE CHILD IF FITNESS IS SAME
            IF(ftrial(i,1).EQ.fparticle(i,1))THEN
                CALL mutationvec(np,dims,cons,xparticle(:,:),&
     &            fparticle(:,:),cparticle(:,:),&
     &            p,mutvector,mut_pool,narchive,xarchive,archive_yn,np)
                CALL trialvec(dims,xparticle(i,:),cross_pool,&
     &                    mutvector,trialvector(i,:))
                CALL boundaries(dims,xupper,xlower,trialvector(i,:))
                CALL particleobjective(func,trialvector(i,:),dims,cons,&
     &            ftrial(i,:),ctrial(i,:),feval)
             END IF 

         END DO
        END DO

!***********************************************************************
!       SELECTION
!***********************************************************************

        DO m=1,nb
          nbp_used=nbp(m)
          DO p=1,nbp_used
            i=ineighbour(m,p)

!           SELECTION
            CALL selection(i,dims,cons,np,&
     &        trialvector(i,:),ftrial(i,:),ctrial(i,:),xparticle(1:np,:),&
     &        fparticle(1:np,:),cparticle(1:np,:),cur_cross(i),cur_mut(i),&
     &        nsuc_cross,nsuc_mut,success_cross,success_mut,narchive,xarchive)

!           UPDATE NEIGHBOURHOODS WITH NEW POSITIONS
            xneighbour(m,p,:)=xparticle(i,:)
            fneighbour(m,p,:)=fparticle(i,:)
            cneighbour(m,p,:)=cparticle(i,:)

          END DO
        END DO

!       L-SHADE - MEMORY UPDATE
        CALL lshade_memupdate(np,fparticle(1:np,:),ftrial(1:np,:),&
     &    fold(1:np,:),nsuc_cross,nsuc_mut,success_cross,success_mut,&
     &    ksuccess,mem_cross(ksuccess),mem_mut(ksuccess))

        DEALLOCATE(xneighbour,fneighbour,cneighbour,ineighbour)

!***********************************************************************
!       STOPPING CRITERIA
!***********************************************************************

!       MAXIMUM FUNCTION EVALUATIONS CHECK
        IF(feval.GE.fmax)EXIT

!     END OPTIMISER LOOP
      END DO

!     FINAL SOLUTION
      OPEN(5,FILE='solution')
      WRITE(5,*)np,dims
      DO i=1,np
        j=0
        CALL particleobjective(func,xparticle(i,:),dims,cons,&
     &     fparticle(i,:),cparticle(i,:),j)
        WRITE(5,*)xparticle(i,:),fparticle(i,:)
      END DO
      CLOSE(5)


!
!***********************************************************************
      CONTAINS
!***********************************************************************
!

      INTEGER FUNCTION intrand(a,b)
!     Generate random integer between a and b
      INTEGER, INTENT(IN):: a,b
      DOUBLE PRECISION:: rand
      CALL RANDOM_NUMBER(rand)
      intrand=INT(rand*DBLE(b+1-a))+a
      END FUNCTION

      DOUBLE PRECISION FUNCTION dblerand(a,b)
!     Generate random number between a and b
      DOUBLE PRECISION, INTENT(IN):: a,b
      DOUBLE PRECISION:: rand
      CALL RANDOM_NUMBER(rand)
      dblerand=a+(b-a)*rand
      END FUNCTION



      END SUBROUTINE





!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE init_random_seed()

      IMPLICIT NONE

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  
      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))
  
      CALL SYSTEM_CLOCK(COUNT=clock)
  
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)
  
      DEALLOCATE(seed)

      END SUBROUTINE
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE particleinitial(dims,samples,xl,xu,mat)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: dims, samples
      DOUBLE PRECISION, DIMENSION(dims), INTENT(IN):: xl, xu
      DOUBLE PRECISION, DIMENSION(samples, dims):: mat

      INTEGER:: jj, ii
      DOUBLE PRECISION:: rand

!     RANDOM
      DO jj=1,dims
        DO ii=1,samples
          CALL RANDOM_NUMBER(rand)
          mat(ii,jj)=xl(jj)+(xu(jj)-xl(jj))*rand
        END DO
      END DO

      END SUBROUTINE particleinitial
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE particleobjective(func,x,n,c,obj,cons,feval)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: func, n, c
      INTEGER, INTENT(INOUT):: feval
      DOUBLE PRECISION, INTENT(IN):: x(n)
      DOUBLE PRECISION, INTENT(OUT):: obj(2), cons(c)

      INTEGER:: ii

      CALL niching_func_cons(func,x,n,c,obj(1),cons)
      feval=feval+1

      obj(2)=0.0d0
      DO ii=1,c
        obj(2)=obj(2)+MAX(0.0d0,cons(ii))
      END DO


      END SUBROUTINE
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE xeuc(n,d,x,euclidean)
      IMPLICIT NONE
      INTEGER, INTENT(IN):: n, d
      DOUBLE PRECISION, INTENT(IN):: x(n,d)
      DOUBLE PRECISION, INTENT(OUT):: euclidean(n,n)
      INTEGER:: i, j, ii, k
      DOUBLE PRECISION:: euc
      DO j=1,n
        ii=0
        DO i=1,n
          euc=0.0d0
          DO k=1,d
            euc=euc+(x(i,k)-x(j,k))**2
          END DO
          euc=SQRT(euc)
          euclidean(i,j)=euc
        END DO
      END DO
      END SUBROUTINE
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE SDE_sort(np,dims,cons,xparticle,fparticle,cparticle,&
     &           xsort,fsort,csort,isort)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: np, dims, cons
      DOUBLE PRECISION, INTENT(IN)::xparticle(np,dims),fparticle(np,2),&
     &   cparticle(np,cons)
      INTEGER, INTENT(OUT)::isort(np)
      DOUBLE PRECISION, INTENT(OUT)::xsort(np,dims),fsort(np,2),&
     &   csort(np,cons)

      INTEGER:: i
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: fitsort
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: sorting

      ALLOCATE(sorting(np,dims+2+cons+1))
      ALLOCATE(fitsort(np))

!     SORT PARTICLES WITH FITNESS FIRST, THEN CONSTARIN VIOLATION SECOND
      DO i=1,np
            sorting(i,1:dims)=xparticle(i,1:dims)
            sorting(i,1+dims:2+dims)=fparticle(i,1:2)
            sorting(i,3+dims:2+dims+cons)=cparticle(i,1:cons)
            sorting(i,3+dims+cons)=DBLE(i)
            IF(fparticle(i,2).LE.0.0d0)THEN
              fitsort(i)=fparticle(i,1)
            ELSE
              fitsort(i)=1.0E10+fparticle(i,2)
            END IF
      END DO
      CALL sort(fitsort, sorting, np, 3+dims+cons)
      xsort(:,1:dims)=sorting(:,1:dims)
      fsort(:,1:2)=sorting(:,1+dims:2+dims)
      csort(:,1:cons)=sorting(:,3+dims:2+dims+cons)
      isort(:)=INT(sorting(:,3+dims+cons))

      END SUBROUTINE
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE neighbourhoods(np,dims,cons,xparticle,fparticle,&
     &   cparticle,nb,nbp,xneighbour,fneighbour,cneighbour,&
     &   ineighbour,euclidean)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: np,dims,cons,nb,nbp(nb)
      DOUBLE PRECISION, INTENT(IN):: xparticle(np,dims), &
     &   fparticle(np,2),cparticle(np,cons),euclidean(np,np)
      DOUBLE PRECISION, INTENT(OUT):: xneighbour(nb,MAXVAL(nbp),dims), &
     &   fneighbour(nb,MAXVAL(nbp),2),cneighbour(nb,MAXVAL(nbp),cons)
      INTEGER, INTENT(OUT):: ineighbour(nb,MAXVAL(nbp))

      INTEGER:: i, m, nignore, icur, j, icount, istar, k, iwin, bts
      INTEGER, ALLOCATABLE, DIMENSION(:):: iignore, ieuc
      DOUBLE PRECISION:: minf, minc
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: euctemp
      LOGICAL:: ignore

!     NSDE HAS EQUAL NUMBER OF NEIGHBOURHOODS DETERMINED BY
!     CLOSEST EUCLIDEAN TO BEST FUNCTION VALUE
      ALLOCATE(iignore(np))
      ALLOCATE(euctemp(np))
      ALLOCATE(ieuc(np))
      nignore=0
      DO m=1,nb

!         FIND BEST FITNESS OF PARTICLES NOT ALREADY USED TO USE AS SEED
          minf=1.0E10
          minc=1.0E10
          DO j=1,np
             ignore=.FALSE.
             IF(nignore.GT.0)THEN
              DO i=1,nignore
                IF(iignore(i).EQ.j)THEN
                  ignore=.TRUE.
                  EXIT
                END IF
              END DO
             END IF
             IF(.NOT.ignore)THEN
              iwin=BTS(fparticle(j,1),fparticle(j,2),minf,minc)
              IF(iwin.EQ.1)THEN
                  minf=fparticle(j,1)
                  minc=fparticle(j,2)
                  istar=j
              END IF
             END IF
          END DO
          xneighbour(m,1,1:dims)=xparticle(istar,1:dims)
          fneighbour(m,1,1:2)=fparticle(istar,1:2)
          cneighbour(m,1,1:cons)=cparticle(istar,1:cons)
          ineighbour(m,1)=istar
          nignore=nignore+1
          iignore(nignore)=istar

!         FIND NEAREST PARTICLES TO FORM NEIGHBOURHOOD
          DO k=2,nbp(m)
!           EUCLIDEAN FROM LEAD PARTICLE TO ALL OTHERS NOT ALREADY USED
            icount=0
            DO j=1,np
              ignore=.FALSE.
              DO i=1,nignore
                IF(iignore(i).EQ.j)THEN
                  ignore=.TRUE.
                  EXIT
                END IF
              END DO
              IF(.NOT.ignore)THEN
                icount=icount+1
                euctemp(icount)=euclidean(istar,j)
                ieuc(icount)=j
              END IF
            END DO
!           FIND SMALLEST OF THOSE AND ADD TO NEIGHBOURHOOD
            icur=MINLOC(euctemp(1:icount),1)
            icur=ieuc(icur)
            xneighbour(m,k,1:dims)=xparticle(icur,1:dims)
            fneighbour(m,k,1:2)=fparticle(icur,1:2)
            cneighbour(m,k,1:cons)=cparticle(icur,1:cons)
            ineighbour(m,k)=icur
            nignore=nignore+1
            iignore(nignore)=icur
          END DO

      END DO

      END SUBROUTINE
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE mutationvec(nbp,dims,cons,xneighbour,fneighbour,&
     &        cneighbour,curvec,mutvector,&
     &        mut_pool,narchive,xarchive,archive_yn,np)
      IMPLICIT NONE

      INTEGER, INTENT(IN):: nbp,dims,cons,curvec,&
     &   narchive,archive_yn,np
      DOUBLE PRECISION, INTENT(IN):: xneighbour(nbp,dims), &
     &   fneighbour(nbp,2), cneighbour(nbp,cons), &
     &   mut_pool, xarchive(np,dims)
      DOUBLE PRECISION, INTENT(OUT):: mutvector(dims)

      INTEGER:: r1, r2, r3
      INTEGER, ALLOCATABLE, DIMENSION(:):: isort
      DOUBLE PRECISION:: mutation, pval, rand
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)::xsort,&
     &     fsortout,csort

!     MUTATION CONSTANT FROM POOL
      mutation=mut_pool

!     FIND RANDOM BEST PARTICLE BY BINARY TOURNAMENT SELECTION 
!     FOR DE/CURRENT-TO-PBEST/1/BIN
      CALL RANDOM_NUMBER(rand)
      pval=(2.0d0/DBLE(nbp))+(0.2d0-(2.0d0/DBLE(nbp)))*rand
      ALLOCATE(xsort(nbp,dims))
      ALLOCATE(fsortout(nbp,2))
      ALLOCATE(csort(nbp,cons))
      ALLOCATE(isort(nbp))
      CALL SDE_sort(nbp,dims,cons,xneighbour,fneighbour,cneighbour,&
     &           xsort,fsortout,csort,isort)
      r1=isort(intrand(1,INT(pval*nbp)))

!     GET SECOND RANDOM INDIVIDUAL
      DO 
       r2=intrand(1,nbp)
       IF((r2.NE.curvec).AND.(r2.NE.r1))EXIT
      END DO

!     GET THIRD RANDOM INDIVIDUAL
      IF(archive_yn.EQ.1)THEN
        DO 
          r3=intrand(1,nbp+narchive)
          IF((r3.NE.curvec).AND.(r3.NE.r2).AND.(r2.NE.r1))EXIT
        END DO
      ELSE
        DO 
          r3=intrand(1,nbp)
          IF((r3.NE.curvec).AND.(r3.NE.r2).AND.(r2.NE.r1))EXIT
        END DO
      END IF

!     DE/CURRENT-TO-PBEST/1/BIN
!     INCLUDE ARCHIVE
      IF(r3.GT.nbp)THEN
          r3=r3-nbp
          mutvector(:)=xneighbour(curvec,:)+&
     &       mutation*(xneighbour(r1,:)-xneighbour(curvec,:))+&
     &       mutation*(xneighbour(r2,:)-xarchive(r3,:))
      ELSE
          mutvector(:)=xneighbour(curvec,:)+&
     &       mutation*(xneighbour(r1,:)-xneighbour(curvec,:))+&
     &       mutation*(xneighbour(r2,:)-xneighbour(r3,:))
      END IF



      CONTAINS
      INTEGER FUNCTION intrand(a,b)
      INTEGER, INTENT(IN):: a,b
      DOUBLE PRECISION:: rand
      CALL RANDOM_NUMBER(rand)
      intrand=INT(rand*DBLE(b+1-a))+a
      END FUNCTION

      END SUBROUTINE
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      INTEGER FUNCTION BTS(f1,c1,f2,c2)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: f1,f2,c1,c2

!     IF BOTH FEASIBLE THEN ONE WITH BETS FUNCTION VALUE WINS
      IF((c1.LE.0.0d0).AND.(c2.LE.0.0d0))THEN
          IF(f1.LT.f2)THEN
            BTS=1
          ELSEIF(f2.LT.f1)THEN
            BTS=2
          ELSE
            BTS=0
          END IF
!     IF ONE NOT FEASIBLE THEN THE ONLY FEASIBLE ONE WINS
      ELSEIF(c1.LE.0.0d0)THEN
          BTS=1
      ELSEIF(c2.LE.0.0d0)THEN
          BTS=2
!     IF BOTH NOT FEASIBLE THEN THE MINIMUM CONSTRAINT VIOLATION WINS
      ELSEIF(c1.LT.c2)THEN
          BTS=1
      ELSEIF(c2.LT.c1)THEN
          BTS=2
      ELSE
          BTS=0
      END IF

      END FUNCTION
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE trialvec(dims,xcurrent,cross_pool,mutvector,&
     &   trialvector)
      IMPLICIT NONE


      INTEGER, INTENT(IN):: dims
      DOUBLE PRECISION, INTENT(IN):: xcurrent(dims),&
     &   mutvector(dims),cross_pool
      DOUBLE PRECISION, INTENT(OUT):: trialvector(dims)

      INTEGER:: testr, j
      DOUBLE PRECISION:: rand1,crossover

!     MUTATION CONSTANT FROM POOL
      crossover=cross_pool

!     BINOMIAL CROSSOVER
      testr=intrand(1,dims)
      DO j=1,dims
           CALL RANDOM_NUMBER(rand1)
           IF((rand1 .LE. crossover) .OR. (j .EQ. testr))THEN
             trialvector(j)=mutvector(j)
           ELSE
             trialvector(j)=xcurrent(j)
           END IF
      END DO

      CONTAINS
      INTEGER FUNCTION intrand(a,b)
      INTEGER, INTENT(IN):: a,b
      DOUBLE PRECISION:: rand
      CALL RANDOM_NUMBER(rand)
      intrand=INT(rand*DBLE(b+1-a))+a
      END FUNCTION


      END SUBROUTINE
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE boundaries(dims,xupper,xlower,trialvector)
      IMPLICIT NONE


      INTEGER, INTENT(IN):: dims
      DOUBLE PRECISION, INTENT(IN):: xupper(dims), xlower(dims)
      DOUBLE PRECISION, INTENT(INOUT):: trialvector(dims)

      INTEGER:: j
      DOUBLE PRECISION:: rand1

!      RANDOM REINITIALISATION
       DO j=1,dims
           IF((trialvector(j) .GT. xupper(j)) .OR. &
     &     (trialvector(j) .LT. xlower(j)))THEN
             CALL RANDOM_NUMBER(rand1)
             trialvector(j)=xlower(j)+rand1*(xupper(j)-xlower(j))
           END IF
       END DO

      END SUBROUTINE
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE selection(i,dims,cons,np,xtrial,&
     &      ftrial,ctrial,xparticle,fparticle,cparticle,&
     &      crossover,mutation,nsuc_cross,nsuc_mut,success_cross,&
     &      success_mut,narchive,xarchive)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: dims,cons,np,i
      INTEGER, INTENT(INOUT):: nsuc_cross,nsuc_mut,narchive
      DOUBLE PRECISION, INTENT(IN):: xtrial(dims), ftrial(2), &
     &   ctrial(cons),crossover,mutation
      DOUBLE PRECISION, INTENT(INOUT):: xparticle(np,dims),&
     &   fparticle(np,2),cparticle(np,cons),success_cross(np),&
     &   success_mut(np), xarchive(np,dims)
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: xold
      INTEGER:: r1
      LOGICAL:: selected

      ALLOCATE(xold(dims))
      xold(:)=xparticle(i,:)

!     SELECTION
!     NEIGHBOURHOOD SPECIES DE
      CALL deselect(np,i,dims,cons,xtrial,ftrial,ctrial,&
     &     xparticle,fparticle,cparticle,selected)

      IF(selected)THEN
!       UPDATE HISTORICAL MEMORY
        nsuc_cross=nsuc_cross+1
        nsuc_mut=nsuc_mut+1
        success_cross(nsuc_cross)=crossover
        success_mut(nsuc_mut)=mutation

!       UPDATE ARCHIVE - MUST NOT BE >NP SO PUT IN RANDOM PLACE
        IF(narchive.GE.np)THEN
          r1=intrand(1,np)
          xarchive(r1,:)=xold(:)
        ELSE
          narchive=narchive+1
          xarchive(narchive,:)=xold(:)
        END IF
      END IF

      CONTAINS
      INTEGER FUNCTION intrand(a,b)
      INTEGER, INTENT(IN):: a,b
      DOUBLE PRECISION:: rand
      CALL RANDOM_NUMBER(rand)
      intrand=INT(rand*DBLE(b+1-a))+a
      END FUNCTION

      END SUBROUTINE
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE deselect(np,i,dims,cons,xtrial,ftrial,ctrial,xparticle,&
     &      fparticle,cparticle,selected)

      IMPLICIT NONE
      INTEGER, INTENT(IN):: dims,cons,np,i
      DOUBLE PRECISION, INTENT(IN):: xtrial(dims), ftrial(2), &
     &   ctrial(cons)
      DOUBLE PRECISION, INTENT(INOUT):: xparticle(np,dims), &
     &   fparticle(np,2), cparticle(np,cons)
      LOGICAL, INTENT(OUT):: selected

      INTEGER:: BTS, iwin

!     NORMAL DE SELECTION PROCEDURE:
!       IF TRIAL IS BETTER THAN CURRENT, REPLACE CURRENT

      iwin=BTS(ftrial(1),ftrial(2),fparticle(i,1),fparticle(i,2))
      IF(iwin.EQ.1)THEN
            xparticle(i,:)=xtrial(:)
            fparticle(i,:)=ftrial(:)
            cparticle(i,:)=ctrial(:)
            selected=.TRUE.
      ELSE
            xparticle(i,:)=xparticle(i,:)
            fparticle(i,:)=fparticle(i,:)
            cparticle(i,:)=cparticle(i,:)
            selected=.FALSE.
      END IF

      END SUBROUTINE
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE lshade_F_CR(np,mem_cross,mem_mut,crossover,mutation)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: np
      DOUBLE PRECISION, INTENT(IN):: mem_cross(np), mem_mut(np)
      DOUBLE PRECISION, INTENT(OUT):: crossover,mutation

      INTEGER:: ri

      ri=intrand(1,np)
      IF(mem_cross(ri).GE. 1.0E9)THEN
        crossover=0.0d0
      ELSE
        crossover=rand_normal(mem_cross(ri),0.1d0)
        IF(crossover.GT.1.0d0)THEN
          crossover=1.0d0
        ELSEIF(crossover.LT.0.0d0)THEN
          crossover=0.0d0
        END IF
      END IF

      mutation=-1.0d0
      DO
        mutation=rand_cauchy(mem_mut(ri),0.1d0)
        IF(mutation.GT.1.0d0)THEN
          mutation=1.0d0
          EXIT
        ELSEIF(mutation.GT.0.0d0)THEN
          EXIT
        END IF
      END DO

      CONTAINS
      INTEGER FUNCTION intrand(a,b)
      INTEGER, INTENT(IN):: a,b
      DOUBLE PRECISION:: rand
      CALL RANDOM_NUMBER(rand)
      intrand=INT(rand*DBLE(b+1-a))+a
      END FUNCTION

      FUNCTION rand_normal(mean,stdev) RESULT(c)
      DOUBLE PRECISION :: mean,stdev,c,temp(2),theta,r,pi
      pi=4.0d0*ATAN(1.0d0)
      IF(stdev <= 0.0d0) THEN
        WRITE(*,*) "Standard Deviation must be +ve"
      ELSE
        CALL RANDOM_NUMBER(temp)
        r=(-2.0d0*log(temp(1)))**0.5
        theta = 2.0d0*PI*temp(2)
        c= mean+stdev*r*sin(theta)
      END IF
      END FUNCTION
      FUNCTION rand_cauchy(median, SCALE1) RESULT(ans)
      DOUBLE PRECISION ans,median,scale1,p,pi
      pi=4.0d0*ATAN(1.0d0)
      IF (scale1 <= 0.0d0) THEN
        WRITE(*,*) "Scale PARAMETER must be positive"
      END IF
      CALL RANDOM_NUMBER(p)
      ans = median + SCALE1*tan(PI*(p - 0.5))
      END FUNCTION

      END SUBROUTINE
!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!
      SUBROUTINE lshade_memupdate(np,fparticle,ftrial,fold,&
     &   nsuc_cross,nsuc_mut,success_cross,success_mut,ksuccess,&
     &      mem_cross,mem_mut)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: np, nsuc_cross, nsuc_mut
      DOUBLE PRECISION, INTENT(IN):: fparticle(np,2), ftrial(np,2), &
     &   success_cross(np), success_mut(np), fold(np,2)
      INTEGER, INTENT(INOUT):: ksuccess 
      DOUBLE PRECISION, INTENT(INOUT):: mem_cross,mem_mut

      INTEGER:: ii,icount
      DOUBLE PRECISION:: meanwa_cr, weight, denom, meanwl_f, upp, low
      DOUBLE PRECISION:: upp2, low2

      IF((nsuc_cross.NE.0).AND.(nsuc_mut.NE.0))THEN

        denom=0.0d0
        DO ii=1,np
          IF(ABS(fparticle(ii,1)-ftrial(ii,1)).LE.1.0E-15)THEN
            denom=denom+ABS(fold(ii,1)-ftrial(ii,1))
          END IF
        END DO
        denom=denom+1E-15

        icount=0
        upp=0.0d0
        low=0.0d0
        upp2=0.0d0
        low2=0.0d0
        DO ii=1,np
          IF(ABS(fparticle(ii,1)-ftrial(ii,1)).LE.1.0E-15)THEN
            weight=ABS(fold(ii,1)-ftrial(ii,1))/denom
            icount=icount+1
            upp=upp+weight*success_mut(icount)*success_mut(icount)
            low=low+weight*success_mut(icount)
            upp2=upp2+weight*success_cross(icount)*success_cross(icount)
            low2=low2+weight*success_cross(icount)
          END IF
        END DO
        meanwl_f=upp/(low+1.0E-25)
        meanwa_cr=upp2/(low2+1.0E-25)
        mem_mut=meanwl_f

        IF((mem_cross.GE.1.0E9).OR.(MAXVAL(success_cross(:)).LE.1E-10))THEN
          mem_cross=1.0E10
        ELSE
          mem_cross=meanwa_cr
        END IF

        ksuccess=ksuccess+1
        IF(ksuccess.GT.np)ksuccess=1

      END IF

      END SUBROUTINE

!
!
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!
!



