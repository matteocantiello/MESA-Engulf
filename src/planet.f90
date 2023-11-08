module planet_module 

use star_lib
use star_def
use const_def
use math_lib

implicit none 

public
contains




subroutine calculate_intercepted_area (id, Orbital_separation, R_influence, f_disruption, area)
     implicit none
     integer, intent(in) :: id
     real(dp), intent(in) :: Orbital_separation, R_influence, f_disruption
     real(dp), intent(out):: area
     real(dp) :: penetration_depth
     integer :: ierr

     type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


        ! Calculate area used for drag calculation 
        penetration_depth = 0d0
        area = 0d0 ! Initialize cross section of companion (physical or Bondi) for calculating aerodynamic or gravitational drag

      ! Do the calculation only if this is a grazing collision and if the planet has not been destroyed yet
        if (Orbital_separation > s% r(1) + R_influence) then
            penetration_depth = 0.d0
        else
            penetration_depth = calculate_penetration_depth(R_influence,s% r(1), Orbital_separation)
        endif

        if (penetration_depth >= 0.0 .and. (Orbital_separation >= (s% r(1) - R_influence)) .and. (f_disruption <= 1d0)) then
            ! Calculate intersected area. Rstar-rr is x in sketch
              area = intercepted_area (penetration_depth, R_influence)
            !  write(*,*) 'Grazing Collision. Engulfed area fraction: ', s% model_number, area/(pi * pow(R_influence, 2.0))
        else
            ! Full engulfment. Cross section area = Planet area
              area = pi * pow(R_influence, 2d0)
            !  write(*,*) 'Full engulfment. R_influence, area',s% model_number,R_influence/Rsun,area
        end if

end subroutine calculate_intercepted_area



function intercepted_area(x, radius) result(area)

  ! Calculate 2D Intercepted area of planet grazing host star noting that the radius is not
  ! necessarily the radius of the planet, it could be the Bondi radius if it is larger.

  implicit none
  
  real(dp), intent(in) :: x, radius    
  real(dp) :: area
  real(dp) :: alpha, y
  
  if (x < radius) then ! Case when less than half of the planet is engulfed
  
    y = radius - x
    alpha = acos(y/radius)
    area = radius * (radius*alpha - y*sin(alpha))
    
  else                  ! Case when more than half of the planet is engulfed

    y = x - radius
    alpha = acos(y/radius)
    area = pi*radius**2 - radius*(radius*alpha - y*sin(alpha))
    
  end if

end function intercepted_area



real(dp) function calculate_penetration_depth(R_influence, R_star, Orbital_separation) result(penetration_depth)
     real(dp) :: R_influence, R_star, Orbital_separation
     penetration_depth = R_influence + R_star - Orbital_separation
     if (penetration_depth < 0d0) penetration_depth = 0d0
    ! write(*,*)'From penetration_depth, R_influence,Orbital_separation,penetration_depth' &
    ! ,R_influence,Orbital_separation,penetration_depth
end function calculate_penetration_depth

real(dp) function check_disruption(M_companion,R_companion,v_planet,rho_ambient) result(f)
   ! f > 1 means disruption. This is expected when the ram pressure integrated
   ! over the planet cross section approaches the planet binding energy
     real(dp), intent(in) :: M_companion,R_companion,v_planet,rho_ambient
     real(dp) :: v_esc_planet_square, rho_planet
     rho_planet = 3d0*M_companion/(4d0*pi*pow(R_companion, 3d0))
     v_esc_planet_square = standard_cgrav*M_companion/R_companion
   ! Eq.5 in Jia & Spruit 2018  https://arxiv.org/abs/1808.00467
     f = (rho_ambient*pow(v_planet, 2d0)) / (rho_planet*v_esc_planet_square)
    ! write(*,*)'From Check_Disruption rho_ambient, rho_planet, v_planet, f', &
    !           rho_ambient, rho_planet, v_planet, f
end function check_disruption

      
end module 