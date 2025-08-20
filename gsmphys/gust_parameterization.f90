subroutine compute_gust(im, u_10m, v_10m, ustar, zol, z1, gust_parameter, gust)
  use machine, only : kind_phys
  implicit none

  integer, intent(in) :: im
  real(kind=kind_phys), intent(in) :: u_10m(im), v_10m(im), ustar(im), &
                                      zol(im), z1(im), gust_parameter(im)
  real(kind=kind_phys), intent(out) :: gust(im)

  real(kind=kind_phys), parameter :: PBL_HEIGHT = 1000.0_kind_phys ! in meters
  real(kind=kind_phys) :: H_o_z(im)
  integer :: i

  ! Follows the IFS gust implementation (see Eq. 3.109 in physics documentation for cycle Cy49r1)

  do i = 1,im
    gust(i) = gust_parameter * ustar(i)

    H_o_z(i) = PBL_HEIGHT / z1(i)

    ! where zol > 0, multiply by similarity function f(H/L)
    if (zol(i) > 0.0_kind_phys) then
      gust(i) = gust(i) * (max(0.0_kind_phys, 1.0_kind_phys - 1.0_kind_phys / 24.0_kind_phys * H_o_z(i) * zol(i)))**(1.0_kind_phys/3.0_kind_phys)
    end if
    gust(i) = gust(i) + sqrt(u_10m(i)**2 + v_10m(i)**2)
  end do

  

end subroutine compute_gust