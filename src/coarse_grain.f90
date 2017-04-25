program coarse_grain
    
    use omp_lib
    implicit none
    integer(kind=8), parameter :: n = 2**10
    real(kind=8), dimension(n) :: x,y
    real(kind=8) :: norm,norm_thread,ynorm,ynorm_thread
    integer :: nthreads, points_per_thread,thread_num
    integer :: i,istart,iend

    ! Specify number of threads to use:
    nthreads = 1       ! need this value in serial mode
    !$ nthreads = 4    
    !$ call omp_set_num_threads(nthreads)
    !$ print "('Using OpenMP with ',i3,' threads')", nthreads

    ! Determine how many points to handle with each thread.
    ! Note that dividing two integers and assigning to an integer will
    ! round down if the result is not an integer.  
    ! This, together with the min(...) in the definition of iend below,
    ! insures that all points will get distributed to some thread.
    points_per_thread = (n + nthreads - 1) / nthreads
    print *, "points_per_thread = ",points_per_thread

    ! initialize x:
    do i=1,n
        x(i) = dble(i)  ! convert to double float
        enddo

    norm = 0.d0
    ynorm = 0.d0

    !$omp parallel private(i,norm_thread, &
    !$omp                  istart,iend,thread_num,ynorm_thread) 

    thread_num = 0     ! needed in serial mode
    !$ thread_num = omp_get_thread_num()    ! unique for each thread

    ! Determine start and end index for the set of points to be 
    ! handled by this thread:
    istart = thread_num * points_per_thread + 1
    iend = min((thread_num+1) * points_per_thread, n)

    !$omp critical
    print '("Thread ",i2," will take i = ",i6," through i = ",i6)', thread_num, istart, iend
    !$omp end critical

    norm_thread = 0.d0
    do i=istart,iend
        norm_thread = norm_thread + abs(x(i))
        enddo

    ! update global norm with value from each thread:
    !$omp critical
      norm = norm + norm_thread
      print *, "norm updated to: ",norm
    !$omp end critical

    ! make sure all have updated norm before proceeding:
    !$omp barrier

    ynorm_thread = 0.d0
    do i=istart,iend
        y(i) = x(i) / norm
        ynorm_thread = ynorm_thread + abs(y(i))
        enddo

    ! update global ynorm with value from each thread:
    !$omp critical
      ynorm = ynorm + ynorm_thread
      print *, "ynorm updated to: ",ynorm
    !$omp end critical
    !$omp barrier

    !$omp end parallel 

    print *, "norm of x = ",norm, "  n(n+1)/2 = ",n*(n+1)/2
    print *, 'ynorm should be 1.0:   ynorm = ', ynorm

end program coarse_grain