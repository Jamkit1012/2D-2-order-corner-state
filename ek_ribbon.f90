! This subroutine is used to caculate 2D second order corner state.

subroutine ek_ribbon
    use wmpi
    use para
    implicit none

    ! Define variables
    integer :: ndim1, i, j, m, l,mdim,k,b
    real(dp) :: kmax, n, s
    integer :: lwork, il, iu, info, ierr
    integer, allocatable :: iwork(:)
    real(dp) :: vl, vu, abstol
!    real(dp), allocatable :: ekribbon(:,:), ekribbon_mpi(:,:)
    real(dp), allocatable :: rwork(:), eigenvec(:,:)
    real(dp), allocatable :: eigenvalue(:), wf(:), eigenvalue_mpi(:)
    complex(dp), allocatable :: z(:,:), CHamk(:,:), CHamk_t(:,:)
    complex(dp), allocatable :: work(:), psi(:,:)
    real(dp), allocatable :: wavefunction(:,:), wavefunction_mpi(:,:)
    integer :: x, y
!    character(len=30) :: filename

    ! Set parameters
    ndim1 = Num_wann * nslab1 * nslab2
    lwork = 128 * ndim1
    kmax = 0.5d0
    abstol = 1e-10
    vl = omegamin
    vu = omegamax
!    il = (NumOccupied - 2) * nslab1 * nslab2
!    iu = (NumOccupied + 2) * nslab1 * nslab2
!    mdim=iu-il+1
!    il = max(1, (NumOccupied - 2) * Nslab1 * Nslab2)
!    iu = min(Ndim1, (NumOccupied + 2) * Nslab1 * Nslab2)
    il = 1
    iu = ndim1/2
    m=0

    if (il < 1 .or. iu > Ndim1 .or. il > iu) then
        print *, 'Error: Invalid il or iu values.'
        stop
    endif

    ! Allocate arrays
!    allocate(ekribbon(ndim1, Nk1))
!    allocate(ekribbon_mpi(ndim1, Nk1))
    allocate(eigenvalue_mpi(ndim1))
    allocate(z(ndim1, ndim1))
    allocate(CHamk(ndim1, ndim1), CHamk_t(ndim1, ndim1))
    allocate(rwork(7 * ndim1))
    allocate(work(lwork))
    allocate(eigenvalue(ndim1))
    allocate(eigenvec(ndim1, ndim1))
    allocate(wf(ndim1))
    allocate(wavefunction(nslab1, nslab2), wavefunction_mpi(nslab1,nslab2))
    allocate(iwork(3 * ndim1))  ! Allocate iwork array

    ierr = 0
!    ekribbon = 0.0d0
!    ekribbon_mpi = 0.0d0
    eigenvalue = 0.0d0
    eigenvalue_mpi = 0.0d0
    ! Initialize wavefunction arrays
    wf = 0.0d0
    CHamk_t = 0.0d0
    CHamk = 0.0d0
    eigenvec = 0.0d0
    s = 0.0d0
    wavefunction=0.0d0
    wavefunction_mpi=0.0d0
    ! Loop over k points
    do i=1+cpuid, Nk1, num_cpu
        if (cpuid==0) write(stdout, *) "Ribbonek the i'th kpoint", i, Nk1
        k=kmax*real(i-1)/(Nk1-1)
    ! Construct Hamiltonian at k-point
        call ham_ribbon(k, CHamk)
        CHamk_t = CHamk
    
        ! Ensure the matrix is Hermitian
        do j = 1, ndim1
             do l = 1, j
                 if (CHamk_t(j, l) /= conjg(CHamk_t(l, j))) then
                      print *, 'Matrix is not Hermitian at:', j, l
                      stop 'CHamk_t matrix is not Hermitian'
                 endif
             end do
        end do

    ! Solve eigenvalue problem
!        call zheevx('V', 'I', 'U', ndim1, CHamk_t, ndim1, vl, vu, il, iu, abstol, ndim1/2, &
!                eigenvalue, eigenvec, ndim1, work, lwork, rwork, iwork, info)

!        call zheevx_pack('V', 'U', ndim1, il, iu, CHamk_t, eigenvalue, eigenvec)

        call eigensystem_c('V', 'U', Ndim1, CHamk_t, eigenvalue)
!        ekribbon(:,i)=eigenvalue
     ! Loop over each eigenvalue to calculate the energy and wavefunction
        do b = 1, ndim1
            s = eigenvalue(b)
    !        call append(energy, s)
    !        print *, s
            ! Check the condition for the energy range and update wf
            if (-0.0005d0 < s .and. s < 0.0005d0) then
!            if (3168 < b .and. b < 3177) then
                wf = wf + abs(CHamk_t(:, b))**2
            end if
        end do
    
        ! Construct the wavefunction from wf
        do y = 1, nslab2
            do x = 1, nslab1
                n = 0.0d0
                do m = 1, Num_wann
                    n = n + wf((y - 1) * nslab1 * Num_wann + (x - 1) * Num_wann + m)
                end do
                print *, n
                wavefunction(x, y) = n
            end do
        end do

        if (cpuid.eq.0)write(stdout,'(a2,i4,f12.5,f10.2,a2)')'k',i,eigenvalue(i)
    enddo


#if defined (MPI)
     call mpi_allreduce(eigenvalue,eigenvalue_mpi,size(eigenvalue),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(wavefunction,wavefunction_mpi,size(wavefunction),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     eigenvalue_mpi= eigenvalue
     wavefunction_mpi=wavefunction
#endif

    if (cpuid.eq.0) then
        open(unit=200, file='corner.dat', status='unknown')
        do y = 1, nslab2
            do x = 1, nslab1
                write(200, '(I5, I5, F20.10)') x, y, wavefunction_mpi(x, y)
            end do
            write(200, *)
        end do
    
        close(200)
    
        open(unit=200, file='eigenvalues.dat', status='unknown')
        do j = 1, Ndim1
            write(200, *) j,  eigenvalue_mpi(j)
        end do
        close(200)
        write(stdout,*) 'calculate ribbon_eigenval  done'
    endif

     if (cpuid==0) then
        open(unit=118, file='DOS_corner.gnu')
        write(118, '(a)')'set encoding iso_8859_1'
        write(118, '(a)')'set terminal postscript enhanced color'
        write(118, '(3a)')"set output 'corner.eps'"
        write(118, '(3a)')'set palette cubehelix start 3 cycles 0.5 saturation 1 negative '
        write(118, '(a)')'set size square'
        write(118,'(a, I0, a)')'set xrange [1 : ', nslab1, ']'
        write(118, '(a, I0, a)')'set yrange [1 :', nslab2, ']'
        write(118, '(a)')'set pm3d'
        write(118, '(a)')'set view equal xyz'
        write(118, '(a)')'set view map'
        write(118, '(a)')'set border lw 1'
        write(118, '(a)')'unset key'
        write(118, '(a)')'unset colorbox'
        write(118, '(a)')'set pm3d interpolate 10,10'
        write(118, '(a)')'splot "corner.dat" u 1:2:($3) w pm3d'

        close(118)
     endif

    ! Combine results from all CPUs using MPI_Reduce
    ! Deallocate arrays
    deallocate(eigenvalue_mpi, z, CHamk, CHamk_t, rwork, work, &
               eigenvalue, eigenvec, wavefunction, iwork)

    return
end subroutine ek_ribbon
