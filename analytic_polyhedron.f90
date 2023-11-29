SUBROUTINE icosphereAt(centred, vertex, step, icotype)
    !------------------------------------------------------------------!
    ! Dummy SUBROUTINE, it just call the correct icosphere             !
    !------------------------------------------------------------------!
    INTEGER, INTENT(IN)                                   :: icotype
    REAL(KREAL), INTENT(IN)                               :: step
    REAL(KREAL), DIMENSION(3), INTENT(IN)                 :: centred
    REAL(KREAL), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: vertex
    SELECTCASE(icotype)
       CASE(12)
          ALLOCATE(vertex(12,3))
          CALL icosahedron_at(centred, vertex, step)
       CASE(42)
          ALLOCATE(vertex(24,3))
          CALL rhombicuboctahedron_at(centred, vertex, step)
       CASE(162)
          ALLOCATE(vertex(60,3))
          CALL truncatedicosahedron_at(centred, vertex, step)
    ENDSELECT
ENDSUBROUTINE
!----------------------------------------------------------------------!
SUBROUTINE icosahedron_at(centred, vertex, step)
    !------------------------------------------------------------------!
    ! Dummy SUBROUTINE, it just computes the vertices of a icosahedron !
    !                  (0, ±1, ±φ); (±1, ±φ, 0); (±φ, ±1, 0)           !
    !------------------------------------------------------------------!
    REAL(KREAL)                               :: goldenstep
    REAL(KREAL), INTENT(IN)                   :: step
    REAL(KREAL), DIMENSION(3), INTENT(IN)     :: centred
    REAL(KREAL), DIMENSION(12,3), INTENT(OUT) :: vertex
    goldenstep   = step*0.5_KREAL*(1.0_KREAL + DSQRT(5.0_KREAL))
    vertex(1 ,:) = centred(:) + (/0.0_KREAL,  step,  goldenstep/)
    vertex(2 ,:) = centred(:) + (/0.0_KREAL, -step,  goldenstep/)
    vertex(3 ,:) = centred(:) + (/0.0_KREAL,  step, -goldenstep/)
    vertex(4 ,:) = centred(:) + (/0.0_KREAL, -step, -goldenstep/)
    vertex(5 ,:) = centred(:) + (/ goldenstep, 0.0_KREAL,  step/)
    vertex(6 ,:) = centred(:) + (/-goldenstep, 0.0_KREAL,  step/)
    vertex(7 ,:) = centred(:) + (/ goldenstep, 0.0_KREAL, -step/)
    vertex(8 ,:) = centred(:) + (/-goldenstep, 0.0_KREAL, -step/)
    vertex(9 ,:) = centred(:) + (/ step,  goldenstep, 0.0_KREAL/)
    vertex(10,:) = centred(:) + (/-step,  goldenstep, 0.0_KREAL/)
    vertex(11,:) = centred(:) + (/ step, -goldenstep, 0.0_KREAL/)
    vertex(12,:) = centred(:) + (/-step, -goldenstep, 0.0_KREAL/)
ENDSUBROUTINE
!----------------------------------------------------------------------!
SUBROUTINE rhombicuboctahedron_at(centred, vertex, step)
    !------------------------------------------------------------------!
    ! Dummy SUBROUTINE, it just computes the vertices of a             !
    !                  rhombicuboctahedron                             !
    !                  all even permutations of (±1, ±1, ±(1+√2))      !
    !------------------------------------------------------------------!
    REAL(KREAL)                               :: oneplusstep
    REAL(KREAL), INTENT(IN)                   :: step
    REAL(KREAL), DIMENSION(3), INTENT(IN)     :: centred
    REAL(KREAL), DIMENSION(24,3), INTENT(OUT) :: vertex
    oneplusstep  = step*(1.0_KREAL + DSQRT(2.0_KREAL))
    vertex(1 ,:) = centred(:) + (/ step,  step,  oneplusstep/)
    vertex(2 ,:) = centred(:) + (/ step, -step,  oneplusstep/)
    vertex(3 ,:) = centred(:) + (/ step,  step, -oneplusstep/)
    vertex(4 ,:) = centred(:) + (/ step, -step, -oneplusstep/)
    vertex(5 ,:) = centred(:) + (/-step,  step,  oneplusstep/)
    vertex(6 ,:) = centred(:) + (/-step, -step,  oneplusstep/)
    vertex(7 ,:) = centred(:) + (/-step,  step, -oneplusstep/)
    vertex(8 ,:) = centred(:) + (/-step, -step, -oneplusstep/)
    vertex(9 ,:) = centred(:) + (/ oneplusstep,  step,  step/)
    vertex(10,:) = centred(:) + (/ oneplusstep, -step,  step/)
    vertex(11,:) = centred(:) + (/ oneplusstep,  step, -step/)
    vertex(12,:) = centred(:) + (/ oneplusstep, -step, -step/)
    vertex(13,:) = centred(:) + (/-oneplusstep,  step,  step/)
    vertex(14,:) = centred(:) + (/-oneplusstep, -step,  step/)
    vertex(15,:) = centred(:) + (/-oneplusstep,  step, -step/)
    vertex(16,:) = centred(:) + (/-oneplusstep, -step, -step/)
    vertex(17,:) = centred(:) + (/ step,  oneplusstep,  step/)
    vertex(18,:) = centred(:) + (/ step, -oneplusstep,  step/)
    vertex(19,:) = centred(:) + (/ step,  oneplusstep, -step/)
    vertex(20,:) = centred(:) + (/ step, -oneplusstep, -step/)
    vertex(21,:) = centred(:) + (/-step,  oneplusstep,  step/)
    vertex(22,:) = centred(:) + (/-step, -oneplusstep,  step/)
    vertex(23,:) = centred(:) + (/-step,  oneplusstep, -step/)
    vertex(24,:) = centred(:) + (/-step, -oneplusstep, -step/)
ENDSUBROUTINE
!----------------------------------------------------------------------!
SUBROUTINE truncatedicosahedron_at(centred, vertex, step)
    !------------------------------------------------------------------!
    ! Dummy SUBROUTINE,                                                !
    !       it just computes the vertices of a truncated icosahedron   !
    !       all even permutations of:                                  !
    !                                     (  0,    ±1,      ±3φ  )     !
    !                                     ( ±1,  ±(2+φ),    ±2φ  )     !
    !                                     ( ±φ,    ±2,    ±(2φ+1))     !
    !------------------------------------------------------------------!
    REAL(KREAL)                               :: gr, gs, ts
    REAL(KREAL)                               :: dgs, tgs, dpg, dgp1
    REAL(KREAL), INTENT(IN)                   :: step
    REAL(KREAL), DIMENSION(3), INTENT(IN)     :: centred
    REAL(KREAL), DIMENSION(60,3), INTENT(OUT) :: vertex
    gr  = 0.5_KREAL*(1.0_KREAL + DSQRT(5.0_KREAL))
    gs  = step*gr;          ts   = 2.0_KREAL*step
    dgs = 2.0_KREAL*gs;     dpg  = step*(2.0_KREAL + gr)
    tgs = 3.0_KREAL*gs;     dgp1 = step*(2.0_KREAL*gr + 1.0_KREAL)
    vertex(1 ,:) = centred(:) + (/ 0.0_KREAL,  step,  tgs/)
    vertex(2 ,:) = centred(:) + (/ 0.0_KREAL,  step, -tgs/)
    vertex(3 ,:) = centred(:) + (/ 0.0_KREAL, -step,  tgs/)
    vertex(4 ,:) = centred(:) + (/ 0.0_KREAL, -step, -tgs/)
    vertex(5 ,:) = centred(:) + (/ step,  0.0_KREAL,  tgs/)
    vertex(6 ,:) = centred(:) + (/ step,  0.0_KREAL, -tgs/)
    vertex(7 ,:) = centred(:) + (/-step,  0.0_KREAL,  tgs/)
    vertex(8 ,:) = centred(:) + (/-step,  0.0_KREAL, -tgs/)
    vertex(9 ,:) = centred(:) + (/ step,  tgs,  0.0_KREAL/)
    vertex(10,:) = centred(:) + (/ step, -tgs,  0.0_KREAL/)
    vertex(11,:) = centred(:) + (/-step,  tgs,  0.0_KREAL/)
    vertex(12,:) = centred(:) + (/-step, -tgs,  0.0_KREAL/)
    vertex(13,:) = centred(:) + (/ step,  dpg,  dgs/)
    vertex(14,:) = centred(:) + (/ step,  dpg, -dgs/)
    vertex(15,:) = centred(:) + (/ step, -dpg,  dgs/)
    vertex(16,:) = centred(:) + (/ step, -dpg, -dgs/)
    vertex(17,:) = centred(:) + (/-step,  dpg,  dgs/)
    vertex(18,:) = centred(:) + (/-step,  dpg, -dgs/)
    vertex(19,:) = centred(:) + (/-step, -dpg,  dgs/)
    vertex(20,:) = centred(:) + (/-step, -dpg, -dgs/)
    vertex(21,:) = centred(:) + (/ dpg,  step,  dgs/)
    vertex(22,:) = centred(:) + (/ dpg,  step, -dgs/)
    vertex(23,:) = centred(:) + (/ dpg, -step,  dgs/)
    vertex(24,:) = centred(:) + (/ dpg, -step, -dgs/)
    vertex(25,:) = centred(:) + (/-dpg,  step,  dgs/)
    vertex(26,:) = centred(:) + (/-dpg,  step, -dgs/)
    vertex(27,:) = centred(:) + (/-dpg, -step,  dgs/)
    vertex(28,:) = centred(:) + (/-dpg, -step, -dgs/)
    vertex(29,:) = centred(:) + (/ dpg,  dgs,  step/)
    vertex(30,:) = centred(:) + (/ dpg,  dgs, -step/)
    vertex(31,:) = centred(:) + (/ dpg, -dgs,  step/)
    vertex(32,:) = centred(:) + (/ dpg, -dgs, -step/)
    vertex(33,:) = centred(:) + (/-dpg,  dgs,  step/)
    vertex(34,:) = centred(:) + (/-dpg,  dgs, -step/)
    vertex(35,:) = centred(:) + (/-dpg, -dgs,  step/)
    vertex(36,:) = centred(:) + (/-dpg, -dgs, -step/)
    vertex(37,:) = centred(:) + (/ gs,  ts,  dgp1/)
    vertex(38,:) = centred(:) + (/ gs,  ts, -dgp1/)
    vertex(39,:) = centred(:) + (/ gs, -ts,  dgp1/)
    vertex(40,:) = centred(:) + (/ gs, -ts, -dgp1/)
    vertex(41,:) = centred(:) + (/-gs,  ts,  dgp1/)
    vertex(42,:) = centred(:) + (/-gs,  ts, -dgp1/)
    vertex(43,:) = centred(:) + (/-gs, -ts,  dgp1/)
    vertex(44,:) = centred(:) + (/-gs, -ts, -dgp1/)
    vertex(45,:) = centred(:) + (/ ts,  gs,  dgp1/)
    vertex(46,:) = centred(:) + (/ ts, -gs,  dgp1/)
    vertex(47,:) = centred(:) + (/ ts,  gs, -dgp1/)
    vertex(48,:) = centred(:) + (/ ts, -gs, -dgp1/)
    vertex(49,:) = centred(:) + (/-ts,  gs,  dgp1/)
    vertex(50,:) = centred(:) + (/-ts, -gs,  dgp1/)
    vertex(51,:) = centred(:) + (/-ts,  gs, -dgp1/)
    vertex(52,:) = centred(:) + (/-ts, -gs, -dgp1/)
    vertex(53,:) = centred(:) + (/ ts,  dgp1,  gs/)
    vertex(54,:) = centred(:) + (/ ts,  dgp1, -gs/)
    vertex(55,:) = centred(:) + (/ ts, -dgp1,  gs/)
    vertex(56,:) = centred(:) + (/ ts, -dgp1, -gs/)
    vertex(57,:) = centred(:) + (/-ts,  dgp1,  gs/)
    vertex(58,:) = centred(:) + (/-ts,  dgp1, -gs/)
    vertex(59,:) = centred(:) + (/-ts, -dgp1,  gs/)
    vertex(60,:) = centred(:) + (/-ts, -dgp1, -gs/)
ENDSUBROUTINE
!----------------------------------------------------------------------!
SUBROUTINE truncatedicosidodecahedron_at(centred, vertex, step)
    !------------------------------------------------------------------!
    ! Dummy SUBROUTINE,                                                !
    !       it just computes the vertices of a truncated               !
    !       icosidodecahedron, all even permutations of:               !
    !                                   (  ±1/φ,    ±1/φ,  ±(3+φ)  )   !
    !                                   (  ±2/φ,     ±φ,   ±(1+2φ) )   !
    !                                   (  ±1/φ,    ±φ^2,  ±(3φ-1) )   !
    !                                   ( ±(2φ−1),   ±2,   ±(2+φ)  )   !
    !                                   (   ±φ,      ±3,    ±2φ    )   !
    !------------------------------------------------------------------!
    REAL(KREAL)                                :: gr, igr, digr, tpgr
    REAL(KREAL)                                :: dpgr, op2gr
    REAL(KREAL)                                :: tgrm1, dgrm1, gr2
    REAL(KREAL)                                :: dstep, tstep, dgr
    REAL(KREAL), INTENT(IN)                    :: step
    REAL(KREAL), DIMENSION(3), INTENT(IN)      :: centred
    REAL(KREAL), DIMENSION(120,3), INTENT(OUT) :: vertex
    gr    = 0.5_KREAL*(1.0_KREAL + DSQRT(5.0_KREAL))
    igr   = step/gr
    digr  = 2.0_KREAL*igr
    tpgr  = step*(3.0_KREAL + gr)
    dpgr  = step*(2.0_KREAL + gr)
    op2gr = step*(1.0_KREAL + 2.0_KREAL*gr)
    tgrm1 = step*(3.0_KREAL*gr - 1.0_KREAL)
    dgrm1 = step*(2.0_KREAL*gr - 1.0_KREAL)
    gr2   = step*gr*gr
    dstep = 2.0_KREAL*step
    tstep = 3.0_KREAL*step
    dgr   = 2.0_KREAL*step*gr
    gr    = step*gr
    vertex(1  ,:) = centred(:) + (/ igr,  igr,  tpgr /)
    vertex(2  ,:) = centred(:) + (/ igr,  igr, -tpgr /)
    vertex(3  ,:) = centred(:) + (/ igr, -igr,  tpgr /)
    vertex(4  ,:) = centred(:) + (/ igr, -igr, -tpgr /)
    vertex(5  ,:) = centred(:) + (/-igr,  igr,  tpgr /)
    vertex(6  ,:) = centred(:) + (/-igr,  igr, -tpgr /)
    vertex(7  ,:) = centred(:) + (/-igr, -igr,  tpgr /)
    vertex(8  ,:) = centred(:) + (/-igr, -igr, -tpgr /)
    vertex(9  ,:) = centred(:) + (/ igr,  tpgr,  igr /)
    vertex(10 ,:) = centred(:) + (/ igr,  tpgr, -igr /)
    vertex(11 ,:) = centred(:) + (/ igr, -tpgr,  igr /)
    vertex(12 ,:) = centred(:) + (/ igr, -tpgr, -igr /)
    vertex(13 ,:) = centred(:) + (/-igr,  tpgr,  igr /)
    vertex(14 ,:) = centred(:) + (/-igr,  tpgr, -igr /)
    vertex(15 ,:) = centred(:) + (/-igr, -tpgr,  igr /)
    vertex(16 ,:) = centred(:) + (/-igr, -tpgr, -igr /)
    vertex(17 ,:) = centred(:) + (/ tpgr,  igr,  igr /)
    vertex(18 ,:) = centred(:) + (/ tpgr,  igr, -igr /)
    vertex(19 ,:) = centred(:) + (/ tpgr, -igr,  igr /)
    vertex(20 ,:) = centred(:) + (/ tpgr, -igr, -igr /)
    vertex(21 ,:) = centred(:) + (/-tpgr,  igr,  igr /)
    vertex(22 ,:) = centred(:) + (/-tpgr,  igr, -igr /)
    vertex(23 ,:) = centred(:) + (/-tpgr, -igr,  igr /)
    vertex(24 ,:) = centred(:) + (/-tpgr, -igr, -igr /)
    vertex(25 ,:) = centred(:) + (/ digr,  gr,  op2gr/)
    vertex(26 ,:) = centred(:) + (/ digr,  gr, -op2gr/)
    vertex(27 ,:) = centred(:) + (/ digr, -gr,  op2gr/)
    vertex(28 ,:) = centred(:) + (/ digr, -gr, -op2gr/)
    vertex(29 ,:) = centred(:) + (/-digr,  gr,  op2gr/)
    vertex(30 ,:) = centred(:) + (/-digr,  gr, -op2gr/)
    vertex(31 ,:) = centred(:) + (/-digr, -gr,  op2gr/)
    vertex(32 ,:) = centred(:) + (/-digr, -gr, -op2gr/)
    vertex(33 ,:) = centred(:) + (/ op2gr,  digr,  gr/)
    vertex(34 ,:) = centred(:) + (/ op2gr,  digr, -gr/)
    vertex(35 ,:) = centred(:) + (/ op2gr, -digr,  gr/)
    vertex(36 ,:) = centred(:) + (/ op2gr, -digr, -gr/)
    vertex(37 ,:) = centred(:) + (/-op2gr,  digr,  gr/)
    vertex(38 ,:) = centred(:) + (/-op2gr,  digr, -gr/)
    vertex(39 ,:) = centred(:) + (/-op2gr, -digr,  gr/)
    vertex(40 ,:) = centred(:) + (/-op2gr, -digr, -gr/)
    vertex(41 ,:) = centred(:) + (/ gr,  op2gr,  digr/)
    vertex(42 ,:) = centred(:) + (/ gr,  op2gr, -digr/)
    vertex(43 ,:) = centred(:) + (/ gr, -op2gr,  digr/)
    vertex(44 ,:) = centred(:) + (/ gr, -op2gr, -digr/)
    vertex(45 ,:) = centred(:) + (/-gr,  op2gr,  digr/)
    vertex(46 ,:) = centred(:) + (/-gr,  op2gr, -digr/)
    vertex(47 ,:) = centred(:) + (/-gr, -op2gr,  digr/)
    vertex(48 ,:) = centred(:) + (/-gr, -op2gr, -digr/)
    vertex(49 ,:) = centred(:) + (/ igr,  gr2,  tgrm1/)
    vertex(50 ,:) = centred(:) + (/ igr,  gr2, -tgrm1/)
    vertex(51 ,:) = centred(:) + (/ igr, -gr2,  tgrm1/)
    vertex(52 ,:) = centred(:) + (/ igr, -gr2, -tgrm1/)
    vertex(53 ,:) = centred(:) + (/-igr,  gr2,  tgrm1/)
    vertex(54 ,:) = centred(:) + (/-igr,  gr2, -tgrm1/)
    vertex(55 ,:) = centred(:) + (/-igr, -gr2,  tgrm1/)
    vertex(56 ,:) = centred(:) + (/-igr, -gr2, -tgrm1/)
    vertex(57 ,:) = centred(:) + (/ tgrm1,  igr,  gr2/)
    vertex(58 ,:) = centred(:) + (/ tgrm1,  igr, -gr2/)
    vertex(59 ,:) = centred(:) + (/ tgrm1, -igr,  gr2/)
    vertex(60 ,:) = centred(:) + (/ tgrm1, -igr, -gr2/)
    vertex(61 ,:) = centred(:) + (/-tgrm1,  igr,  gr2/)
    vertex(62 ,:) = centred(:) + (/-tgrm1,  igr, -gr2/)
    vertex(63 ,:) = centred(:) + (/-tgrm1, -igr,  gr2/)
    vertex(64 ,:) = centred(:) + (/-tgrm1, -igr, -gr2/)
    vertex(65 ,:) = centred(:) + (/ gr2,  tgrm1,  igr/)
    vertex(66 ,:) = centred(:) + (/ gr2,  tgrm1, -igr/)
    vertex(67 ,:) = centred(:) + (/ gr2, -tgrm1,  igr/)
    vertex(68 ,:) = centred(:) + (/ gr2, -tgrm1, -igr/)
    vertex(69 ,:) = centred(:) + (/-gr2,  tgrm1,  igr/)
    vertex(70 ,:) = centred(:) + (/-gr2,  tgrm1, -igr/)
    vertex(71 ,:) = centred(:) + (/-gr2, -tgrm1,  igr/)
    vertex(72 ,:) = centred(:) + (/-gr2, -tgrm1, -igr/)
    vertex(73 ,:) = centred(:) + (/ dgrm1,  dstep,  dpgr/)
    vertex(74 ,:) = centred(:) + (/ dgrm1,  dstep, -dpgr/)
    vertex(75 ,:) = centred(:) + (/ dgrm1, -dstep,  dpgr/)
    vertex(76 ,:) = centred(:) + (/ dgrm1, -dstep, -dpgr/)
    vertex(77 ,:) = centred(:) + (/-dgrm1,  dstep,  dpgr/)
    vertex(78 ,:) = centred(:) + (/-dgrm1,  dstep, -dpgr/)
    vertex(79 ,:) = centred(:) + (/-dgrm1, -dstep,  dpgr/)
    vertex(80 ,:) = centred(:) + (/-dgrm1, -dstep, -dpgr/)
    vertex(81 ,:) = centred(:) + (/ dpgr,  dgrm1,  dstep/)
    vertex(82 ,:) = centred(:) + (/ dpgr, -dgrm1,  dstep/)
    vertex(83 ,:) = centred(:) + (/ dpgr,  dgrm1, -dstep/)
    vertex(84 ,:) = centred(:) + (/ dpgr, -dgrm1, -dstep/)
    vertex(85 ,:) = centred(:) + (/-dpgr,  dgrm1,  dstep/)
    vertex(86 ,:) = centred(:) + (/-dpgr, -dgrm1,  dstep/)
    vertex(87 ,:) = centred(:) + (/-dpgr,  dgrm1, -dstep/)
    vertex(88 ,:) = centred(:) + (/-dpgr, -dgrm1, -dstep/)
    vertex(89 ,:) = centred(:) + (/ dstep,  dpgr,  dgrm1/)
    vertex(90 ,:) = centred(:) + (/ dstep,  dpgr, -dgrm1/)
    vertex(91 ,:) = centred(:) + (/ dstep, -dpgr,  dgrm1/)
    vertex(92 ,:) = centred(:) + (/ dstep, -dpgr, -dgrm1/)
    vertex(93 ,:) = centred(:) + (/-dstep,  dpgr,  dgrm1/)
    vertex(94 ,:) = centred(:) + (/-dstep,  dpgr, -dgrm1/)
    vertex(95 ,:) = centred(:) + (/-dstep, -dpgr,  dgrm1/)
    vertex(96 ,:) = centred(:) + (/-dstep, -dpgr, -dgrm1/)
    vertex(97 ,:) = centred(:) + (/ gr,  tstep,  dgr/)
    vertex(98 ,:) = centred(:) + (/ gr,  tstep, -dgr/)
    vertex(99 ,:) = centred(:) + (/ gr, -tstep,  dgr/)
    vertex(100,:) = centred(:) + (/ gr, -tstep, -dgr/)
    vertex(101,:) = centred(:) + (/-gr,  tstep,  dgr/)
    vertex(102,:) = centred(:) + (/-gr,  tstep, -dgr/)
    vertex(103,:) = centred(:) + (/-gr, -tstep,  dgr/)
    vertex(104,:) = centred(:) + (/-gr, -tstep, -dgr/)
    vertex(105,:) = centred(:) + (/ dgr,  gr,  tstep/)
    vertex(106,:) = centred(:) + (/ dgr,  gr, -tstep/)
    vertex(107,:) = centred(:) + (/ dgr, -gr,  tstep/)
    vertex(108,:) = centred(:) + (/ dgr, -gr, -tstep/)
    vertex(109,:) = centred(:) + (/-dgr,  gr,  tstep/)
    vertex(110,:) = centred(:) + (/-dgr,  gr, -tstep/)
    vertex(111,:) = centred(:) + (/-dgr, -gr,  tstep/)
    vertex(112,:) = centred(:) + (/-dgr, -gr, -tstep/)
    vertex(113,:) = centred(:) + (/ tstep,  dgr,  gr/)
    vertex(114,:) = centred(:) + (/ tstep,  dgr, -gr/)
    vertex(115,:) = centred(:) + (/ tstep, -dgr,  gr/)
    vertex(116,:) = centred(:) + (/ tstep, -dgr, -gr/)
    vertex(117,:) = centred(:) + (/-tstep,  dgr,  gr/)
    vertex(118,:) = centred(:) + (/-tstep,  dgr, -gr/)
    vertex(119,:) = centred(:) + (/-tstep, -dgr,  gr/)
    vertex(120,:) = centred(:) + (/-tstep, -dgr, -gr/)
ENDSUBROUTINE
