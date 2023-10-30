                           PROGRAM icosphere
!**********************************************************************!
!   This program computes icospheres at any order.                     !
!   Less time consuming, but we need more memory, the comparision      !
!   between this algorithm and the brute force is in the pdf file      !
!   in the repository.                                                 !
!**********************************************************************!
!                                          ✨ Victoria Castor, 2023 ✨ !
!**********************************************************************!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
INTEGER                   :: i, j, k, l, order, nver, nfaces, nvertb, nfacesb
INTEGER                   :: vera, verb, verc, verd, vere, verf
INTEGER,      ALLOCATABLE :: atvertex(:,:), faces(:,:)
REAL(KIND=8)              :: golden, radius, norma, timestart, timeend
REAL(KIND=8)              :: midPoint(3)
REAL(KIND=8), ALLOCATABLE :: vertex(:,:)
LOGICAL,      ALLOCATABLE :: edgecomputed(:,:)
!----------------------------------------------------------------------------------------------------------------------------------!

READ(*,*)  order
golden   = 0.5d0*(1.d0 + DSQRT(5.d0))
radius   = DSQRT(2.d0 + golden)
nver     = 10*4**order + 2
nfaces   = 20*4**order
ALLOCATE(vertex(3,nver), faces(3,nfaces)); vertex = 0.d0
ALLOCATE(atvertex(nver,nver), edgecomputed(nver,nver)); atvertex = 0; edgecomputed = .FALSE.
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
faces(:, 1) = [1, 2,  5]; faces(:, 2) = [1, 2,  6]; faces(:, 3) = [1,  5,  9]; faces(:, 4) = [1, 6, 10]; faces(:, 5) = [1, 9, 10]
faces(:, 6) = [2, 5, 11]; faces(:, 7) = [2, 6, 12]; faces(:, 8) = [2, 11, 12]
faces(:, 9) = [3, 4,  7]; faces(:,10) = [3, 4,  8]; faces(:,11) = [3,  7,  9]; faces(:,12) = [3, 8, 10]; faces(:,13) = [3, 9, 10]
faces(:,14) = [4, 7, 11]; faces(:,15) = [4, 8, 12]; faces(:,16) = [4, 11, 12]
faces(:,17) = [5, 7,  9]; faces(:,18) = [5, 7, 11]
faces(:,19) = [6, 8, 10]; faces(:,20) = [6, 8, 12]

norma = DSQRT(DOT_PRODUCT(vertex(:,1)+vertex(:,2),vertex(:,1)+vertex(:,2)))
CALL CPU_TIME(timestart)
DO i = 1, order
  nvertb  = 10*4**(i-1) + 2
  nfacesb = 20*4**(i-1)
  k       = nvertb  + 1
  l       = nfacesb + 1
  DO j = 1, nfacesb
    vera = faces(1,j); verb = faces(2,j); verc = faces(3,j)
    ! Get the vertice number
    CALL get_vertexn(vera, verb, verd, k, edgecomputed, atvertex, vertex, norma, radius)
    CALL get_vertexn(verb, verc, vere, k, edgecomputed, atvertex, vertex, norma, radius)
    CALL get_vertexn(vera, verc, verf, k, edgecomputed, atvertex, vertex, norma, radius)
    ! Actualise the triangles
    faces(2,j) = verd; faces(3,j) = verf
    faces(1,l) = verd; faces(2,l) = verb; faces(3,l) = vere; l = l + 1
    faces(1,l) = verf; faces(2,l) = vere; faces(3,l) = verc; l = l + 1
    faces(1,l) = verd; faces(2,l) = vere; faces(3,l) = verf; l = l + 1
  ENDDO
  norma   = DSQRT(DOT_PRODUCT(vertex(:,verf)+vertex(:,verd),vertex(:,verf)+vertex(:,verd)))
ENDDO
CALL CPU_TIME(timeend)
WRITE(*,*) timeend-timestart

                            CONTAINS
SUBROUTINE get_vertexn(verx, very, verz, k, edgecomputed, atvertex, vertex, norm, radi)
  INTEGER, INTENT(IN)         :: verx, very
  INTEGER, INTENT(OUT)        :: verz
  INTEGER, INTENT(INOUT)      :: k
  INTEGER, INTENT(INOUT)      :: atvertex(:,:)
  LOGICAL, INTENT(INOUT)      :: edgecomputed(:,:)
  REAL(KIND=8), INTENT(IN)    :: radi
  REAL(KIND=8), INTENT(IN)   :: norm
  !REAL(KIND=8), INTENT(OUT)   :: norm
  REAL(KIND=8), INTENT(INOUT) :: vertex(:,:)

  IF (edgecomputed(verx,very)) THEN
    verz = atvertex(verx,very)
  ELSE
    midPoint(:) = vertex(:,verx) + vertex(:,very)
    !norm        = DSQRT(DOT_PRODUCT(midPoint,midPoint))
    verz = k; k = k + 1
    vertex(:,verz) = radi*midPoint(:)/norm
    edgecomputed(verx,very) = .TRUE.
    edgecomputed(very,verx) = .TRUE.
    atvertex(verx,very) = verz
    atvertex(very,verx) = verz
  ENDIF

ENDSUBROUTINE
                           END PROGRAM
