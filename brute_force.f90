                           PROGRAM bruteForce
!**********************************************************************!
!   This program computes icospheres at any order.                     !
!   More time consuming, but more efficient about memory,              !
!   the comparation between burte force and actualizing the faces is   !
!   in the pdf in the repository.                                      !
!**********************************************************************!
!                                          ✨ Victoria Castor, 2023 ✨ !
!**********************************************************************!
IMPLICIT NONE
!----------------------------------------------------------------------!
INTEGER                   :: i, j, k, l, order, nver, nverb, nvertb
REAL(KIND=8)              :: golden, radius, edgeanal, norma
REAL(KIND=8)              :: disatorder, dis, timestart, timeend
REAL(KIND=8)              :: diff(3), midPoint(3)
REAL(KIND=8), ALLOCATABLE :: vertex(:,:)
!----------------------------------------------------------------------!

READ(*,*)  order
golden   = 0.5d0*(1.d0 + DSQRT(5.d0))
radius   = DSQRT(2.d0 + golden)
norma    = 2.d0*golden
nver     = 10*4**order + 2
edgeanal = 2.d0
ALLOCATE(vertex(3,nver)); vertex = 0.d0
!  At order zero, i. e., icosahedron;
! (0, ±1, ±φ); (±1, ±φ, 0); (±φ, 0, ±1)
vertex(:, 1) = (/0.d0,  1.d0,  golden/)
vertex(:, 2) = (/0.d0, -1.d0,  golden/)
vertex(:, 3) = (/0.d0,  1.d0, -golden/)
vertex(:, 4) = (/0.d0, -1.d0, -golden/)
vertex(:, 5) = (/ golden, 0.d0,  1.d0/)
vertex(:, 6) = (/-golden, 0.d0,  1.d0/)
vertex(:, 7) = (/ golden, 0.d0, -1.d0/)
vertex(:, 8) = (/-golden, 0.d0, -1.d0/)
vertex(:, 9) = (/ 1.d0,  golden, 0.d0/)
vertex(:,10) = (/-1.d0,  golden, 0.d0/)
vertex(:,11) = (/ 1.d0, -golden, 0.d0/)
vertex(:,12) = (/-1.d0, -golden, 0.d0/)

CALL CPU_TIME(timestart)
DO i = 1, order
  nverb      = 10*4**(i-1) +2
  disatorder = 1.2d0*edgeanal
  l          = nverb
  DO j = 1, nverb-1
    DO k = j+1, nverb
      diff = vertex(:,j) - vertex(:,k)
      dis  = DSQRT(DOT_PRODUCT(diff,diff))
      IF (dis .LT. disatorder) THEN
        l           = l + 1
        midPoint(:) = vertex(:,j) + vertex(:,k)
        vertex(:,l) = radius*midPoint(:)/norma
      ENDIF
    ENDDO
  ENDDO
  diff     = vertex(:,1) - vertex(:,nverb+1)
  edgeanal = DSQRT(DOT_PRODUCT(diff,diff))
ENDDO
CALL CPU_TIME(timeend)
WRITE(*,*) timeend-timestart

                           END PROGRAM
