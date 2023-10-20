                           PROGRAM icosphere
!**********************************************************************!
!   This program computes icospheres at any order.                     !
!   This is not the most efficient way to know in which points the     !
!   new vertices should be added, we do it by brute force.             !
!                                                                      !
!   Future versions will have a better way.                            !
!                                                                      !
!**********************************************************************!
!                                          ✨ Victoria Castor, 2023 ✨ !
!**********************************************************************!
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                   :: i, j, k, l, m, order, nver, nverb
REAL(KIND=8)              :: golden, radius, edgeanal, norma, disatorder, dis
REAL(KIND=8)              :: centred(3), diff(3), midPoint(3)
REAL(KIND=8), ALLOCATABLE :: vertex(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------

order    = 7! The order we want
golden   = 0.5d0*(1.d0 + DSQRT(5.d0))
radius   = DSQRT(2.d0 + golden)
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

DO i = 1, order
  nverb      = 10*4**(i-1) +2
  disatorder = 1.2d0*edgeanal
  l          = nverb
  DO j = 1, nverb-1
    m = 0
    DO k = j+1, nverb
      diff = vertex(:,j) - vertex(:,k)
      dis  = DSQRT(DOT_PRODUCT(diff,diff))
      IF (dis .LT. disatorder) THEN
        m           = m + 1
        l           = l + 1
        midPoint(:) = vertex(:,j) + vertex(:,k)
        norma       = DSQRT(DOT_PRODUCT(midPoint,midPoint))
        vertex(:,l) = radius*midPoint(:)/norma
      ENDIF
      IF (m .EQ. 5) EXIT
    ENDDO
  ENDDO
  diff     = vertex(:,1) - vertex(:,nverb+1)
  edgeanal = DSQRT(DOT_PRODUCT(diff,diff))
ENDDO

                           END PROGRAM
