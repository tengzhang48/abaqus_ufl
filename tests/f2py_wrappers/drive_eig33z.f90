subroutine drive_sqrtm33z_cs(a0, ii, jj, h, ds)
  implicit none
  double precision, intent(in) :: a0(3,3), h
  integer, intent(in) :: ii, jj
  double precision, intent(out) :: ds(3,3)
  double complex :: az(3,3), sz(3,3)
  integer :: i, j

  do i = 1, 3
    do j = 1, 3
      az(i,j) = dcmplx(a0(i,j), 0.0d0)
    end do
  end do
  az(ii,jj) = az(ii,jj) + dcmplx(0.0d0, h)
  if (ii .ne. jj) then
    az(jj,ii) = az(jj,ii) + dcmplx(0.0d0, h)
  end if

  call sqrtm33z(az, sz)
  do i = 1, 3
    do j = 1, 3
      ds(i,j) = aimag(sz(i,j)) / h
    end do
  end do
end subroutine drive_sqrtm33z_cs

subroutine drive_sqrtm33z_direction_cs(a0, direction, h, ds)
  implicit none
  double precision, intent(in) :: a0(3,3), direction(3,3), h
  double precision, intent(out) :: ds(3,3)
  double complex :: az(3,3), sz(3,3)
  integer :: i, j

  do i = 1, 3
    do j = 1, 3
      az(i,j) = dcmplx(a0(i,j), h * direction(i,j))
    end do
  end do

  call sqrtm33z(az, sz)
  do i = 1, 3
    do j = 1, 3
      ds(i,j) = aimag(sz(i,j)) / h
    end do
  end do
end subroutine drive_sqrtm33z_direction_cs

subroutine drive_sqrtm33z_value(a0, s)
  implicit none
  double precision, intent(in) :: a0(3,3)
  double precision, intent(out) :: s(3,3)
  double complex :: az(3,3), sz(3,3)
  integer :: i, j

  do i = 1, 3
    do j = 1, 3
      az(i,j) = dcmplx(a0(i,j), 0.0d0)
    end do
  end do

  call sqrtm33z(az, sz)
  do i = 1, 3
    do j = 1, 3
      s(i,j) = dble(sz(i,j))
    end do
  end do
end subroutine drive_sqrtm33z_value
