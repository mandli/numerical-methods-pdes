real function matrix_multiply_test(N,method)

    use mod_rand
    implicit none
    
    external DGEMM,DDOT
    
    double precision :: DDOT
    integer, intent(in) :: N,method
    integer :: start,finish,count_rate
    double precision, dimension(:,:), allocatable :: A,B,C
    
    ! Local
    integer :: i,j,k
    
    ! Create the random arrays
    call init_random_seed()
    allocate(A(N,N),B(N,N),C(N,N))
    call random_number(A)
    call random_number(B)
    
    ! Start the timer and start multiplying
    call system_clock(start,count_rate)
    select case(method)
        case(1) ! Default method provided as an intrinsic method
            C = matmul(A,B)
        case(2) ! Simple three loop multiplication
            !$OMP parallel do private(j,k)
            do i=1,N
                do j=1,N
                    do k=1,N
                        C(i,j) = C(i,j) + A(i,k) * B(k,j)
                    enddo
                enddo
            enddo
        case(3) ! OpenMP parallelized double loop
            !$OMP parallel do private(j)
            do i=1,N
                do j=1,N
                    C(i,j) = DDOT(N, A(i,:), 1, B(:,j), 1)
                enddo
            enddo
        case(4) ! BLAS Routine call
            ! call DGEMM(transa,transb,l,n,m,alpha,a,lda,b,ldb,beta,c,ldc)
            call DGEMM('N', 'N', N, N, N, 1.d0, A, N, B, N, 0.d0, C, N)
        case default
            print *, "***ERROR*** Invalid multiplication method!"
            matrix_multiply_test = -1
            return
    end select
    call system_clock(finish,count_rate)
    
    matrix_multiply_test = float(finish - start) / float(count_rate)
    
end function matrix_multiply_test
    
program matrix_multiply
    
    use omp_lib

    implicit none
    
    integer :: N, method, threads
    character(len=10) :: input_N, input_method, input_threads
    real :: matrix_multiply_test, time
    
    select case(iargc())
        case(0)
            N = 1000
            method = 1
            threads = 1
        case(1)
            call getarg(1,input_N)
            read(input_N,'(I10)') N
            method = 1
        case(2)
            call getarg(1,input_N)
            call getarg(2,input_method)
            read(input_N,'(I10)') N
            read(input_method,'(I10)') method
        case(3)
            call getarg(1,input_N)
            call getarg(2,input_method)
            call getarg(3,input_threads)
            read(input_N,'(I10)') N
            read(input_method,'(I10)') method
            read(input_threads,'(I10)') threads
        case default
            print *, "***ERROR*** Too many arguments!"
            stop
    end select
    
    !$ call omp_set_num_threads(threads)

    time = matrix_multiply_test(N, method)
    
    print '("Timing for ",i5,"x",i5," matrices: ",f10.5," s")',N,N,time
    
end program matrix_multiply