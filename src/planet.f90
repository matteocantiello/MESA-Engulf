module planet_module 

use star_lib
use star_def
use const_def
use math_lib

implicit none 

public
contains



! subroutine calculate_intercepted_area(id, separation, radius, disrupt, area)
  
!   ! Calculate intercepted area
  
!   implicit none

!   integer, intent(in) :: id
!   real(dp), intent(in) :: separation, radius, disrupt 
!   real(dp), intent(out) :: area
!   real(dp) :: depth
!   integer :: ierr

!   type(star_info), pointer :: s
!     include 'formats'
!     ierr = 0
!     call star_ptr(id, s, ierr)
!     if (ierr /= 0) return

!   ! Initialization, pointers, checks

!   depth = calculate_penetration_depth(radius, s%r(1), separation)

!   if (depth > 0 .and. separation > s%r(1)-radius .and. disrupt <= 1) then
!     area = intercepted_area(depth, radius) 
!   else
!     area = pi*radius**2
!   end if

! contains

!   function calculate_penetration_depth(radius, rstar, separation)
!     real(dp) :: calculate_penetration_depth
!     real(dp), intent(in) :: radius, rstar, separation
!     calculate_penetration_depth = radius + rstar - separation 
!   end function

!   function intercepted_area(x, radius)
!     real(dp) :: intercepted_area
!     real(dp), intent(in) :: x, radius  
!     ! Calculation
!   end function

! end subroutine


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

     
function calculate_penetration_depth(r_inf, r_star, separation) result(depth)
  
  real(dp) :: depth
  real(dp), intent(in) :: r_inf, r_star, separation

  depth = r_inf + r_star - separation
  depth = max(depth, 0.0_dp)

end function calculate_penetration_depth


function check_disruption(m_comp, r_comp, v_planet, rho_ambient) result(f_disrupt)

  real(dp) :: f_disrupt
  real(dp), intent(in) :: m_comp, r_comp, v_planet, rho_ambient

  real(dp) :: v_esc2, rho_planet

  rho_planet = 3.0_dp*m_comp/(4.0_dp*pi*powe(r_comp,3.0))
  v_esc2 = standard_cgrav*m_comp/r_comp

  ! See Jia & Spruit 2018, Eq. 5 
  f_disrupt = (rho_ambient*v_planet**2) / (rho_planet*v_esc2)

end function check_disruption

end module 