      SUBROUTINE de_config(particles, nb_size,&
     &     seedrad, cross_pool, mut_pool)

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

   INTEGER, INTENT(OUT):: particles, nb_size
   DOUBLE PRECISION, INTENT(OUT):: seedrad, cross_pool, mut_pool

   INTEGER:: i, n, j, params
   CHARACTER(LEN=200):: linestring,linestring1
   CHARACTER(LEN=50),DIMENSION(:),ALLOCATABLE:: label
   LOGICAL,DIMENSION(:),ALLOCATABLE:: READIN

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   params=5

   ALLOCATE(label(params))
   ALLOCATE(READIN(params))

   DO n=1,params
      READIN(n)=.FALSE.
   END DO

   label(1)='PARTICLES'
   label(2)='CROSSOVER'
   label(3)='MUTATION'
   label(4)='NEIGHBOURHOOD_SIZE'
   label(5)='SEEDING_RADIUS'

   !Comment lines start with '#' so make these blank
   OPEN(10,FILE='DE.conf')
88   linestring(:)=' '
   linestring1(:)=' '
   READ(10,888,ERR=2000,END=2000) linestring
   linestring1=linestring
   CALL TOUC(linestring)
   DO i=1,199
      IF(linestring(i:i).eq.'#') THEN
         DO j=i,200
            linestring(j:j)=' '
         END DO
      END IF
   END DO

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   DO i=1,100

    IF(linestring(i:i+8).eq.'PARTICLES') THEN
        DO j=1,i+8
           linestring(j:j)=' '
        END DO
        DO j=1,100
           IF(linestring(j:j).eq.'=') then
              linestring(j:j)=' '
           END IF
        END DO
        READ(linestring,*) particles
        READIN(1)=.TRUE.

    ELSEIF(linestring(i:i+8).eq.'CROSSOVER') THEN
        DO j=1,i+8
           linestring(j:j)=' '
        END DO
        DO j=1,100
           IF(linestring(j:j).eq.'=') then
              linestring(j:j)=' '
           END IF
        END DO
        READ(linestring,*) cross_pool
        READIN(2)=.TRUE.

    ELSEIF(linestring(i:i+7).eq.'MUTATION') THEN
        DO j=1,i+7
           linestring(j:j)=' '
        END DO
        DO j=1,100
           IF(linestring(j:j).eq.'=') then
              linestring(j:j)=' '
           END IF
        END DO
        READ(linestring,*) mut_pool
        READIN(3)=.TRUE.

    ELSEIF(linestring(i:i+17).eq.'NEIGHBOURHOOD_SIZE') THEN
        DO j=1,i+17
           linestring(j:j)=' '
        END DO
        DO j=1,100
           IF(linestring(j:j).eq.'=') then
              linestring(j:j)=' '
           END IF
        END DO
        READ(linestring,*) nb_size
        READIN(4)=.TRUE.

    ELSEIF(linestring(i:i+13).eq.'SEEDING_RADIUS') THEN
        DO j=1,i+13
           linestring(j:j)=' '
        END DO
        DO j=1,100
           IF(linestring(j:j).eq.'=') then
              linestring(j:j)=' '
           END IF
        END DO
        READ(linestring,*) seedrad
        READIN(5)=.TRUE.

    END IF
   END DO

   GOTO 88
2000    CONTINUE
   CLOSE(10)

   ! CHECK SOLUTION TYPE
   DO n=1,params
      IF(READIN(n).eqv..FALSE.) THEN
         PRINT*,'ERROR: MASTER VARIABLE:' 
         PRINT*,label(n)
         PRINT*,'NOT READ IN CORRECTLY'
         STOP
      END IF
   END DO 

   RETURN

888   FORMAT((A))

CONTAINS

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE touc(AIN)
   IMPLICIT NONE
   CHARACTER(LEN=200):: ain
   INTEGER:: l,lenain,ascbot,asctop,asclow,iain

   DATA ascbot/97/,asctop/122/,asclow/32/

   lenain = LEN(ain)
   DO l=1,lenain
      iain = ichar(ain(l:l))
      if (iain.ge.ascbot.and.iain.le.asctop) then
         ain(l:l) = char(iain-asclow)
      endif
   enddo
   return
   end subroutine touc 
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END SUBROUTINE
