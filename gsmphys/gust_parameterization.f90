subroutine compute_gust(im, u_10m, v_10m, ustar, monin_obukhov_length, gust_parameter, gust)
  use machine, only : kind_phys
  implicit none

  integer, intent(in) :: im
  real(kind=kind_phys), intent(in) :: u_10m(im), v_10m(im), ustar(im), &
                                      monin_obukhov_length(im), gust_parameter(im)
  real(kind=kind_phys), intent(out) :: gust(im)

  real(kind=kind_phys), parameter :: PBL_HEIGHT = 1000.0_kind_phys ! in meters

  ! Follows the IFS gust implementation (see Eq. 3.109 in physics documentation for cycle Cy49r1)

  gust = sqrt(u_10m**2 + v_10m**2) + gust_parameter * ustar

  ! where Monin-Obukov lenght > 0, multiply by similarity function
  where (monin_obukhov_length > 0.0_kind_phys)
    gust = gust * &
           (1.0_kind_phys - 1.0_kind_phys / 24.0_kind_phys * PBL_HEIGHT / monin_obukhov_length)**(1.0_kind_phys/3.0_kind_phys)
  end where

end subroutine compute_gust