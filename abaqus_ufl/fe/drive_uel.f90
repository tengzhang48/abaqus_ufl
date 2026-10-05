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
