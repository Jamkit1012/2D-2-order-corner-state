subroutine ham_ribbon(k, Hamk_ribbon)

     use para
     implicit none

! loop index  
     integer :: i1,j1,i2,j2

! wave vector in 2d
     real(Dp) :: k

! Hamiltonian of slab system
     complex(Dp),intent(out) :: Hamk_ribbon(Num_wann*nslab1*nslab2, &
        Num_wann*nslab1*nslab2)

! the factor 2 is induced by spin
     complex(Dp), allocatable :: Hij(:, :, :, :)

     allocate(Hij(-ijmax:ijmax,-ijmax:ijmax,Num_wann,Num_wann))
     call ham_qlayer2qlayerribbon(k,Hij)    
    ! Initialize Hamiltonian
    Hamk_ribbon = 0.0d0
    
    ! Construct the Hamiltonian
    do i1 = 1, nslab1
        do j1 = 1, nslab2
            do i2 = max(1, i1 - ijmax), min(nslab1, i1 + ijmax)
                do j2 = max(1, j1 - ijmax), min(nslab2, j1 + ijmax) 
                    Hamk_ribbon(idx(i1, j1):idx(i1, j1) + Num_wann - 1, &
                                idx(i2, j2):idx(i2, j2) + Num_wann - 1) &
                                  = Hij(i2 - i1, j2 - j1, :, :)
                end do
            end do
        end do
    end do
    
    ! Deallocate
    deallocate(Hij)
    
contains
    ! Function to calculate the flattened index
    pure function idx(i, j) result(ind)
        integer, intent(in) :: i, j
        integer :: ind
        ind = (i - 1) * nslab2 * Num_wann + (j - 1) * Num_wann + 1
    end function idx
    
end subroutine ham_ribbon
