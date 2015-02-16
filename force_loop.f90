subroutine lennard_jones(r_min, eps, p, p_i, force, potential_energy, N, N_i, d)
    implicit none
    
    integer, intent(inout) :: N, N_i, d; real(8), intent(in) :: r_min, eps
    real(8), intent(in), dimension(N, d) :: p
    real(8), intent(in), dimension(N_i, d) :: p_i    

    real(8), intent(out), dimension(N, d) :: force
    real(8), intent(out) :: potential_energy

    real(8), dimension(d) :: comp_distance
    real(8) :: r2
    integer :: i, j
    
    do i=1,N
        
        ! due to particles in box
        do j=1,N
            comp_distance = (p(i,:)-p(j,:))
            r2 = sum(comp_distance**2)
            if (r2 .NE. 0)  then
                potential_energy = potential_energy + &
                    0.5*eps*((r_min**12/r2**6) - (r_min**6/r2**3))

                force(i,:) = force(i,:) + &
                    6.0*eps*((r_min**12/r2**7) - (r_min**6/r2**4))*comp_distance
            end if
        end do

        ! due to imaged particles outside box
        do j=1,N_i
            comp_distance = (p(i,:)-p_i(j,:))
            r2 = sum(comp_distance**2)
            if (r2 .NE. 0)  then
                potential_energy = potential_energy + &
                    0.5*eps*((r_min**12/r2**6) - (r_min**6/r2**3))

                force(i,:) = force(i,:) + &
                    6.0*eps*((r_min**12/r2**7) - (r_min**6/r2**4))*comp_distance
            end if
        end do
    end do

end subroutine lennard_jones
