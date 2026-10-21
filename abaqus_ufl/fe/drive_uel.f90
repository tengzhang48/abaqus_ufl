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

! One f2py boundary for a whole batch. Element IDs are global and explicit.
subroutine drive_uel_batch(rhs, amatrx, svars, coords, u, du, props, jprops, &
                           time, dtime, pnewdt, lflags, params, jtype, &
                           kstep, kinc, period, elem_ids, ndofel, nsvars, &
                           mcrd, nnode, nprops, njprop, nelem)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, intent(in) :: ndofel, nsvars, mcrd, nnode, nprops, njprop, nelem
  integer, intent(in) :: jtype, kstep, kinc, lflags(6), jprops(njprop)
  integer, intent(in) :: elem_ids(nelem)
  double precision, intent(out) :: rhs(ndofel, nelem)
  double precision, intent(out) :: amatrx(ndofel, ndofel, nelem)
  double precision, intent(inout) :: svars(nsvars, nelem)
  double precision, intent(in) :: coords(mcrd, nnode, nelem)
  double precision, intent(in) :: u(ndofel, nelem), du(ndofel, nelem)
  double precision, intent(in) :: props(nprops), time(2), dtime, params(3), period
  double precision, intent(inout) :: pnewdt
  double precision :: pn
  integer :: e
!f2py intent(out) :: rhs, amatrx
!f2py intent(in,out) :: svars, pnewdt
!f2py intent(in) :: coords, u, du, props, jprops, time, dtime, lflags, params
!f2py intent(in) :: jtype, kstep, kinc, period, elem_ids
!f2py integer intent(hide), depend(u) :: ndofel = shape(u,0)
!f2py integer intent(hide), depend(u) :: nelem = shape(u,1)
!f2py integer intent(hide), depend(svars) :: nsvars = shape(svars,0)
!f2py integer intent(hide), depend(coords) :: mcrd = shape(coords,0)
!f2py integer intent(hide), depend(coords) :: nnode = shape(coords,1)
!f2py integer intent(hide), depend(props) :: nprops = shape(props,0)
!f2py integer intent(hide), depend(jprops) :: njprop = shape(jprops,0)
  do e = 1, nelem
    pn = pnewdt
    call drive_uel(rhs(:,e), amatrx(:,:,e), svars(:,e), coords(:,:,e), &
                   u(:,e), du(:,e), props, jprops, time, dtime, pn, &
                   lflags, params, jtype, kstep, kinc, elem_ids(e), period, &
                   ndofel, nsvars, mcrd, nnode, nprops, njprop)
    if (.not. ieee_is_finite(pn)) then
      pnewdt = 0.0d0
    else
      pnewdt = min(pnewdt, pn)
    end if
  end do
end subroutine drive_uel_batch
