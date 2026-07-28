subroutine drive_umat(stress, statev, ddsdde, dfgrd1, dfgrd0, dstran, &
                      stran, props, time, dtime, temp, dtemp, &
                      predef, dpred, pnewdt, coords, drot, &
                      ntens, nstatv, nprops, ndi, nshr, cmname)
  implicit none
  integer, intent(in) :: ntens, nstatv, nprops, ndi, nshr
  double precision, intent(inout) :: stress(ntens)
  double precision, intent(inout) :: statev(nstatv)
  double precision, intent(out) :: ddsdde(ntens, ntens)
  double precision, intent(in) :: dfgrd1(3,3), dfgrd0(3,3)
  double precision, intent(in) :: dstran(ntens), stran(ntens)
  double precision, intent(in) :: props(nprops)
  double precision, intent(in) :: time(2), dtime, temp, dtemp
  double precision, intent(in) :: predef(1), dpred(1)
  double precision, intent(inout) :: pnewdt
  double precision, intent(in) :: coords(3), drot(3,3)
  character(len=80), intent(in) :: cmname

  double precision :: sse, spd, scd, rpl, drpldt
  double precision :: ddsddt(ntens), drplde(ntens)
  integer :: noel, npt, layer, kspt, kinc
  integer :: jstep(4)
  double precision :: celent

!f2py intent(in,out) :: stress, statev, pnewdt
!f2py intent(out) :: ddsdde
!f2py intent(in) :: dfgrd1, dfgrd0, dstran, stran, props
!f2py intent(in) :: time, dtime, temp, dtemp, predef, dpred
!f2py intent(in) :: coords, drot, cmname
!f2py integer intent(hide), depend(dstran) :: ntens = shape(dstran,0)
!f2py integer intent(hide), depend(statev) :: nstatv = shape(statev,0)
!f2py integer intent(hide), depend(props) :: nprops = shape(props,0)
!f2py integer intent(in) :: ndi, nshr

  sse = 0.0d0
  spd = 0.0d0
  scd = 0.0d0
  rpl = 0.0d0
  ddsddt = 0.0d0
  drplde = 0.0d0
  drpldt = 0.0d0
  noel = 1
  npt = 1
  layer = 1
  kspt = 1
  kinc = 1
  jstep = (/1, 1, 0, 0/)
  celent = 1.0d0

  call umat(stress, statev, ddsdde, sse, spd, scd, &
            rpl, ddsddt, drplde, drpldt, &
            stran, dstran, time, dtime, temp, dtemp, &
            predef, dpred, cmname, ndi, nshr, ntens, nstatv, &
            props, nprops, coords, drot, pnewdt, celent, &
            dfgrd0, dfgrd1, noel, npt, layer, kspt, jstep, kinc)

end subroutine drive_umat
