program fine_grain
   
    use omp_lib
    
    implicit none

    integer :: i, thread_num
    integer(kind=8), parameter :: n = 2**10
 
    real(kind=8), dimension(n) :: x, y
    real(kind=8) :: norm,ynorm
 
    integer :: nthreads
    
    ! Specify number of threads to use:
    nthreads = 1       ! need this value in serial mode
    !$ nthreads = 1    
    !$ call omp_set_num_threads(nthreads)
    !$ print "('Using OpenMP with ',i3,' threads')", nthreads

    ! Specify number of threads to use:
    !$ call omp_set_num_threads(4)
 
    ! initialize x:
    !$omp parallel do 
    do i=1,n
        x(i) = real(i, kind=8)  ! convert to double float
    enddo

    norm = 0.d0
    ynorm = 0.d0

    !$omp parallel private(i)

    !$omp do reduction(+ : norm)
    do i=1,n
        norm = norm + abs(x(i))
    enddo

     !$omp barrier   ! not needed (implicit)

    !$omp do reduction(+ : ynorm)
    do i=1,n
        y(i) = x(i) / norm
        ynorm = ynorm + abs(y(i))
    enddo
    
    !$omp end parallel

    print *, "norm of x = ",norm, "  n(n+1)/2 = ",n*(n+1)/2
    print *, 'ynorm should be 1.0:   ynorm = ', ynorm

end program fine_grain
