subroutine compute_gust(im, u_10m, v_10m, ustar, zol, z1, gust_parameter, gust)
  use machine, only : kind_phys
  implicit none

  integer, intent(in) :: im
  real(kind=kind_phys), intent(in) :: u_10m(im), v_10m(im), ustar(im), &
                                      zol(im), z1(im), gust_parameter(im)
  real(kind=kind_phys), intent(out) :: gust(im)

  real(kind=kind_phys), parameter :: PBL_HEIGHT = 1000.0_kind_phys ! in meters
  real(kind=kind_phys) :: H_o_z(im)

  ! Follows the IFS gust implementation (see Eq. 3.109 in physics documentation for cycle Cy49r1)

  gust = gust_parameter * ustar

  H_o_z = PBL_HEIGHT / z1

  ! where zol > 0, multiply by similarity function f(H/L)
  where (zol > 0.0_kind_phys)
    gust = gust * &
           (max(0.0_kind_phys, 1.0_kind_phys - 1.0_kind_phys / 24.0_kind_phys * H_o_z * zol))**(1.0_kind_phys/3.0_kind_phys)
  end where


  gust = gust + sqrt(u_10m**2 + v_10m**2)

end subroutine compute_gust