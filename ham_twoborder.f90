subroutine ham_twoborder(k, Hamk_ribbon)

     use para
     implicit none
  
     integer :: i1,j1,i2,j2
     real(Dp) :: k
     
     complex(Dp),intent(out) :: ham_twoborder(Num_wann*nslab1*nslab2, &
        Num_wann*nslab1*nslab2)

     complex(Dp), allocatable :: Hij(:, :, :, :)

    allocate(Hij(-ijmax:ijmax,-ijmax:ijmax,Num_wann,Num_wann))
    call ham_qlayer2qlayerribbon(k,Hij)    
    ham_twoborder = 0.0d0
    
    ! Construct the Hamiltonian
    do i1 = 1, nslab1
        do j1 = 1, nslab2
            do i2 = max(1, i1 - ijmax), min(nslab1, i1 + ijmax)
                do j2 = max(1, j1 - ijmax), min(nslab2, j1 + ijmax) 
                    ham_twoborder(idx(i1, j1):idx(i1, j1) + Num_wann - 1, &
                                idx(i2, j2):idx(i2, j2) + Num_wann - 1) &
                                  = Hij(i2 - i1, j2 - j1, :, :)
                end do
            end do
        end do
    end do
    
    deallocate(Hij)
    
contains
    ! Function to calculate the flattened index
    pure function idx(i, j) result(ind)
        integer, intent(in) :: i, j
        integer :: ind
        ind = (i - 1) * nslab2 * Num_wann + (j - 1) * Num_wann + 1
    end function idx
    
end subroutine ham_twoborder
