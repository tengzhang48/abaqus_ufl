subroutine drive_uel(rhs, amatrx, svars, coords, u, du, props, jprops, &
                     time, dtime, pnewdt, lflags, params, jtype, &
                     kstep, kinc, jelem, period, ndofel, nsvars, &
                     mcrd, nnode, nprops, njprop)
  implicit none
  integer, intent(in) :: ndofel, nsvars, mcrd, nnode, nprops, njprop
  integer, intent(in) :: jtype, kstep, kinc, jelem
  integer, intent(in) :: lflags(6)
  integer, intent(in) :: jprops(njprop)
  double precision, intent(out) :: rhs(ndofel)
  double precision, intent(out) :: amatrx(ndofel, ndofel)
  double precision, intent(inout) :: svars(nsvars)
  double precision, intent(in) :: coords(mcrd, nnode)
  double precision, intent(in) :: u(ndofel), du(ndofel)
  double precision, intent(in) :: props(nprops)
  double precision, intent(in) :: time(2), dtime, params(3), period
  double precision, intent(inout) :: pnewdt

  integer, parameter :: nrhs = 1
  integer, parameter :: mdload = 1
  integer, parameter :: ndload = 0
  integer, parameter :: npredf = 1
  integer :: i
  integer :: jdltYP(mdload, 1)
  double precision :: rhs2(ndofel, nrhs)
  double precision :: du2(ndofel, nrhs)
  double precision :: energy(8)
  double precision :: vel(ndofel), accel(ndofel)
  double precision :: adlmag(mdload, 1), ddlmag(mdload, 1)
  double precision :: predef(2, npredf, nnode)

!f2py intent(out) :: rhs, amatrx
!f2py intent(in,out) :: svars, pnewdt
!f2py intent(in) :: coords, u, du, props, jprops, time, dtime, lflags, params
!f2py intent(in) :: jtype, kstep, kinc, jelem, period
!f2py integer intent(hide), depend(u) :: ndofel = shape(u,0)
!f2py integer intent(hide), depend(svars) :: nsvars = shape(svars,0)
!f2py integer intent(hide), depend(coords) :: mcrd = shape(coords,0)
!f2py integer intent(hide), depend(coords) :: nnode = shape(coords,1)
!f2py integer intent(hide), depend(props) :: nprops = shape(props,0)
!f2py integer intent(hide), depend(jprops) :: njprop = shape(jprops,0)

  rhs2 = 0.0d0
  amatrx = 0.0d0
  energy = 0.0d0
  vel = 0.0d0
  accel = 0.0d0
  adlmag = 0.0d0
  ddlmag = 0.0d0
  predef = 0.0d0
  jdltYP = 0

  do i = 1, ndofel
    du2(i, 1) = du(i)
  end do

  call uel(rhs2, amatrx, svars, energy, ndofel, nrhs, nsvars, &
           props, nprops, coords, mcrd, nnode, u, du2, vel, accel, &
           jtype, time, dtime, kstep, kinc, jelem, params, ndload, &
           jdltYP, adlmag, predef, npredf, lflags, ndofel, ddlmag, &
           mdload, pnewdt, jprops, njprop, period)

  do i = 1, ndofel
    rhs(i) = rhs2(i, 1)
  end do
end subroutine drive_uel


subroutine drive_uel_batch(rhs, amatrx, svars, coords, u, du, props, jprops, &
                           time, dtime, pnewdt, lflags, params, jtype, &
                           kstep, kinc, period, nelem, ndofel, nsvars, &
                           mcrd, nnode, nprops, njprop)
  ! Batched driver: one f2py call evaluates ALL nelem elements (the per-element
  ! marshalling across the f2py boundary -- ~70% of assembly time -- is the cost
  ! this removes).  Same UEL, looped in Fortran.  jelem = element index (1-based).
  implicit none
  integer, intent(in) :: nelem, ndofel, nsvars, mcrd, nnode, nprops, njprop
  integer, intent(in) :: jtype, kstep, kinc
  integer, intent(in) :: lflags(6), jprops(njprop)
  double precision, intent(out) :: rhs(ndofel, nelem)
  double precision, intent(out) :: amatrx(ndofel, ndofel, nelem)
  double precision, intent(inout) :: svars(nsvars, nelem)
  double precision, intent(in) :: coords(mcrd, nnode, nelem)
  double precision, intent(in) :: u(ndofel, nelem), du(ndofel, nelem)
  double precision, intent(in) :: props(nprops)
  double precision, intent(in) :: time(2), dtime, params(3), period
  double precision, intent(inout) :: pnewdt

  integer, parameter :: nrhs = 1, mdload = 1, ndload = 0, npredf = 1
  integer :: ie, i, j
  integer :: jdltYP(mdload, 1)
  double precision :: rhs2(ndofel, nrhs), du2(ndofel, nrhs), amx(ndofel, ndofel)
  double precision :: svloc(nsvars), uloc(ndofel), crd(mcrd, nnode)
  double precision :: energy(8), vel(ndofel), accel(ndofel)
  double precision :: adlmag(mdload, 1), ddlmag(mdload, 1)
  double precision :: predef(2, npredf, nnode), pn, pn_in

!f2py intent(out) :: rhs, amatrx
!f2py intent(in,out) :: svars, pnewdt
!f2py intent(in) :: coords, u, du, props, jprops, time, dtime, lflags, params
!f2py intent(in) :: jtype, kstep, kinc, period
!f2py integer intent(hide), depend(u) :: ndofel = shape(u,0)
!f2py integer intent(hide), depend(u) :: nelem = shape(u,1)
!f2py integer intent(hide), depend(svars) :: nsvars = shape(svars,0)
!f2py integer intent(hide), depend(coords) :: mcrd = shape(coords,0)
!f2py integer intent(hide), depend(coords) :: nnode = shape(coords,1)
!f2py integer intent(hide), depend(props) :: nprops = shape(props,0)
!f2py integer intent(hide), depend(jprops) :: njprop = shape(jprops,0)

  pn_in = pnewdt
  do ie = 1, nelem
    rhs2 = 0.0d0; amx = 0.0d0; energy = 0.0d0; vel = 0.0d0; accel = 0.0d0
    adlmag = 0.0d0; ddlmag = 0.0d0; predef = 0.0d0; jdltYP = 0
    do i = 1, ndofel
      du2(i, 1) = du(i, ie)
      uloc(i) = u(i, ie)
    end do
    do i = 1, nsvars
      svloc(i) = svars(i, ie)
    end do
    do j = 1, nnode
      do i = 1, mcrd
        crd(i, j) = coords(i, j, ie)
      end do
    end do
    pn = pn_in
    call uel(rhs2, amx, svloc, energy, ndofel, nrhs, nsvars, &
             props, nprops, crd, mcrd, nnode, uloc, du2, vel, accel, &
             jtype, time, dtime, kstep, kinc, ie, params, ndload, &
             jdltYP, adlmag, predef, npredf, lflags, ndofel, ddlmag, &
             mdload, pn, jprops, njprop, period)
    do i = 1, ndofel
      rhs(i, ie) = rhs2(i, 1)
    end do
    do i = 1, nsvars
      svars(i, ie) = svloc(i)
    end do
    do j = 1, ndofel
      do i = 1, ndofel
        amatrx(i, j, ie) = amx(i, j)
      end do
    end do
    if (pn < pnewdt) pnewdt = pn
  end do
end subroutine drive_uel_batch
