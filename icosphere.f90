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
INTEGER                   :: i, j, k, l, m, n, order, nver, nverb, nedges, nedgesb, nfaces, nvertb, nfacesb
INTEGER                   :: vera, verb, verc, verd, vere, verf
INTEGER, ALLOCATABLE      :: edges(:,:), atvertex(:,:), auxfaces(:,:), faces(:,:)
REAL(KIND=8)              :: golden, radius, edgeanal, norma, disatorder, dis, timestart, timeend
REAL(KIND=8)              :: centred(3), diff(3), midPoint(3)
REAL(KIND=8), ALLOCATABLE :: vertex(:,:), expvert(:,:)
LOGICAL                   :: found, matrices_are_equal
LOGICAL, ALLOCATABLE      :: already_checked(:), edgecomputed(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------

order    = 7! The order we want
golden   = 0.5d0*(1.d0 + DSQRT(5.d0))
radius   = DSQRT(2.d0 + golden)
nver     = 10*4**order + 2
nedges   = 30*4**order
nfaces   = 20*4**order
edgeanal = 2.d0
WRITE(*,*) nver, nedges, nfaces
ALLOCATE(vertex(3,nver), expvert(3,nver), faces(3,nfaces), edges(2,nedges)); vertex = 0.d0
ALLOCATE(atvertex(nver,nver), edgecomputed(nver,nver)); atvertex = 0
ALLOCATE(auxfaces(3,nfaces))
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

expvert = vertex; edgecomputed = .FALSE.
CALL CPU_TIME(timestart)
DO i = 1, order
  nvertb  = 10*4**(i-1) + 2
  nfacesb = 20*4**(i-1)
  nedgesb = 30*4**(i-1)
  k       = nvertb  + 1
  l       = nfacesb + 1
  DO j = 1, nfacesb
    vera = faces(1,j); verb = faces(2,j); verc = faces(3,j)
    IF (edgecomputed(vera,verb)) THEN
      verd = atvertex(vera,verb)
    ELSE
      midPoint(:) = expvert(:,vera) + expvert(:,verb)
      norma       = DSQRT(DOT_PRODUCT(midPoint,midPoint))
      verd = k; k = k + 1
      expvert(:,verd) = radius*midPoint(:)/norma
      edgecomputed(vera,verb) = .TRUE.
      edgecomputed(verb,vera) = .TRUE.
      atvertex(vera,verb) = verd
      atvertex(verb,vera) = verd
    ENDIF
    IF (edgecomputed(verb,verc)) THEN
      vere = atvertex(verb,verc)
    ELSE
      midPoint(:) = expvert(:,verb) + expvert(:,verc)
      norma       = DSQRT(DOT_PRODUCT(midPoint,midPoint))
      vere = k; k = k + 1
      expvert(:,vere) = radius*midPoint(:)/norma
      edgecomputed(verb,verc) = .TRUE.
      edgecomputed(verc,verb) = .TRUE.
      atvertex(verb,verc) = vere
      atvertex(verc,verb) = vere
    ENDIF
    IF (edgecomputed(vera,verc)) THEN
      verf = atvertex(vera,verc)
    ELSE
      midPoint(:) = expvert(:,vera) + expvert(:,verc)
      norma       = DSQRT(DOT_PRODUCT(midPoint,midPoint))
      verf = k; k = k + 1
      expvert(:,verf) = radius*midPoint(:)/norma
      edgecomputed(vera,verc) = .TRUE.
      edgecomputed(verc,vera) = .TRUE.
      atvertex(vera,verc) = verf
      atvertex(verc,vera) = verf
    ENDIF
    ! Actualise the triangles
    faces(2,j) = verd; faces(3,j) = verf
    faces(1,l) = verd; faces(2,l) = verb; faces(3,l) = vere; l = l + 1
    faces(1,l) = verf; faces(2,l) = vere; faces(3,l) = verc; l = l + 1
    faces(1,l) = verd; faces(2,l) = vere; faces(3,l) = verf; l = l + 1
  ENDDO
ENDDO
CALL CPU_TIME(timeend)
WRITE(*,*) timeend-timestart

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
        norma       = DSQRT(DOT_PRODUCT(midPoint,midPoint))
        vertex(:,l) = radius*midPoint(:)/norma
      ENDIF
    ENDDO
  ENDDO
  diff     = vertex(:,1) - vertex(:,nverb+1)
  edgeanal = DSQRT(DOT_PRODUCT(diff,diff))
ENDDO
CALL CPU_TIME(timeend)
WRITE(*,*) timeend-timestart

allocate(already_checked(SIZE(vertex,2)))
matrices_are_equal = .TRUE.; already_checked = .FALSE.

l = 0
DO i = 1, SIZE(vertex,2)
    found = .FALSE.
    DO j = 1, SIZE(vertex,2)
        IF (already_checked(j)) CYCLE
        IF (ALL(vertex(:,i) == expvert(:,j))) THEN
            found = .TRUE.
            already_checked(j) = .TRUE.
            !EXIT
        ENDIF
    ENDDO

    IF (.NOT. found) THEN
        l = l +1
        matrices_are_equal = .FALSE.
        !RETURN
    ENDIF
ENDDO
write(*,*) 'ended', l

VERTEX=3.0*VERTEX; EXPVERT=3.0*EXPVERT
OPEN(13, FILE='vertex.xyz')
    WRITE(UNIT=13,FMT=*) SIZE(VERTEX,2)
    WRITE(UNIT=13,FMT=*) 
    DO I = 1, SIZE(VERTEX,2)
        WRITE(UNIT=13,FMT='(A,3(F9.4))') 'H', VERTEX(:,I)
    ENDDO
CLOSE(13)
OPEN(14, FILE='expvert.xyz')
    WRITE(UNIT=14,FMT=*) SIZE(EXPVERT,2)
    WRITE(UNIT=14,FMT=*) 
    DO I = 1, 12
        WRITE(UNIT=14,FMT='(A,3(F9.4))') 'He', EXPVERT(:,I)
    ENDDO
    DO I = 13, SIZE(EXPVERT,2)
        WRITE(UNIT=14,FMT='(A,3(F9.4))') 'H', EXPVERT(:,I)
    ENDDO
CLOSE(14)

                           END PROGRAM
