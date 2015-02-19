! lennard_jones function calculates the forces between particles and total potential energy

! lennard_jones(r_min, eps, L, r_c, p, N, d, force, Ep)

! input:
!   r_min   -> lj parameter
!   eps     -> lj parameter
!   L       -> size of the box
!   r_c     -> cutoff length for interaction
!   p       -> position components matrix
!   N       -> number of particles (automatic)
!   d       -> dimensionality (automatic)

! output:
!   force   -> force components matrix
!   Ep      -> potential energy
!   vir     -> factor for the verial theorem pressure calculation



subroutine lennard_jones(r_min, eps, L, r_c, p, N, d, force, Ep)
    implicit none
    
    integer, intent(inout) :: N, d; real(8), intent(in) :: r_min, eps, r_c, L
    real(8), intent(in), dimension(N, d) :: p

    real(8), intent(out), dimension(N, d) :: force
    real(8), intent(out) :: Ep

    real(8), dimension(d) :: r
    real(8) :: r2
    integer :: i, j
    
    force = 0
    Ep = 0
    vir = 0
    do i=1,N
        do j=1,N
            r = p(i,:) - p(j,:) - L*ANINT((p(i,:) - p(j,:))/L)
            r2 = sum(r**2)

            if ((r2 .NE. 0) .and. (r2 < r_c**2))  then
                ! potential energy
                Ep = Ep + &
                    0.5*eps*((r_min**12/r2**6) - 2*(r_min**6/r2**3))
                
                ! force matrix
                force(i,:) = force(i,:) + &
                    6.0*eps*((r_min**12/r2**7) - (r_min**6/r2**4))*r
            end if
        end do
    end do

end subroutine lennard_jones
