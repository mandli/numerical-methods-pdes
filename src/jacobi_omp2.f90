! Solve Poisson equation
!  u_{xx} = f(x)    x \in [a, b]
! with 
!  u(a) = alpha, u(b) = beta
! using Jacobi iteration and a coarse-grain parallelism approach using OpenMP.
program jacobi_omp2

    use omp_lib
    
    implicit none
    
    ! Problem specification and storage
    real(kind=8) :: a, b, alpha, beta
    real(kind=8), dimension(:), allocatable :: x, u, u_new, f
    real(kind=8) :: dx, tolerance, du_max, du_max_thread

    integer(kind=4) :: n, num_threads, points_per_thread, thread_id, start, end
    integer(kind=8) :: i, num_iterations
    integer(kind=8), parameter :: MAX_ITERATIONS = 100000
    real(kind=8) :: time(2)

    ! Boundaries
    a = 0.d0
    b = 1.d0
    alpha = 0.d0
    beta = 3.d0

    ! Specify number of threads to use:
    num_threads = 2
    !$ call omp_set_num_threads(num_threads)
    !$ print "('Using OpenMP with ',i3,' threads')", num_threads

    N = 100

    ! Allocate storage for boundary points too
    allocate(x(0:N + 1), u(0:N + 1), u_new(0:N + 1), f(0:N + 1))

    call cpu_time(time(1))

    ! grid spacing:
    dx = (b - a) / (N + 1.d0)

    ! Tolerance
    tolerance = 0.1d0 * dx**2

    ! Determine how many points to handle with each thread.
    ! Note that dividing two integers and assigning to an integer will
    ! round down if the result is not an integer.  
    ! This, together with the min(...) in the definition of iend below,
    ! insures that all points will get distributed to some thread.
    points_per_thread = (n + num_threads - 1) / num_threads
    print *, "points_per_thread = ", points_per_thread

    ! Start of the parallel block... 
    ! ------------------------------

    ! This is the only time threads are forked in this program:
    !$omp parallel private(thread_id, num_iterations, start, end, i,   &
    !$OMP                  du_max_thread)

    ! Set thread is, default to 0 if in serial
    thread_id = 0
    !$ thread_id = omp_get_thread_num()

    ! Determine start and end index
    start = thread_id * points_per_thread + 1
    end = min((thread_id + 1) * points_per_thread, N)

    ! Output some thread information and indexing
    !$omp critical
    print '("Thread ",i2," will take i = ",i6," through i = ",i6)', &
          thread_id, start, end
    !$omp end critical

    ! Set iniital guess and construct the grid
    do i=start, end
        ! grid points:
        x(i) = i * dx
        ! source term:
        f(i) = exp(x(i))
        ! initial guess (satisfies boundary conditions and sets them)
        u(i) = alpha + x(i) * (beta - alpha)
    enddo
    ! Note that the above does not set the boundaries, do this in a single thread
    !$omp single
    x(0) = a
    x(N + 1) = b
    u(0) = alpha
    u(N + 1) = beta
    !$omp end single nowait

    ! Main Jacobi iteration loop
    do num_iterations=1, MAX_ITERATIONS
        
        ! Make one thread reset the global du_max
        !$omp single
        du_max = 0.d0
        !$omp end single

        ! Private to each thread
        du_max_thread = 0.d0
        do i=start, end
            u_new(i) = 0.5d0 * (u(i-1) + u(i+1) - dx**2 * f(i))
            du_max_thread = max(du_max_thread, abs(u_new(i) - u(i)))
        end do

        ! Compute global du_max
        !$omp critical
        du_max = max(du_max, du_max_thread)
        !$omp end critical

        ! Make sure all threads are done contributing to du_max
        !$omp barrier

        ! Have one thread print out the convergence info
        !$omp single
        if (mod(num_iterations, 1000) == 0) then
            print '("After ",i8," iterations, dumax = ",d16.6,/)', num_iterations, du_max
        end if
        !$omp end single nowait

        ! Copy new data into u
        do i=start, end
            u(i) = u_new(i)
        end do

        ! Check exit criteria
        if (du_max < tolerance) then
            exit
        end if

        ! Make sure all threads are caught up to this point before starting
        ! another iteration
        !$omp barrier
    end do

    !$omp end parallel

    call cpu_time(time(2))
    if (num_iterations > MAX_ITERATIONS) then
        print *, "Iteration failed!"
        stop
    end if
    print '("CPU time = ",f12.8, " seconds")', time(2) - time(1)
    print *, "Total number of iterations: ", num_iterations

    ! Write out solution for checking
    open(unit=20, file="jacobi_omp2.txt", status="unknown")
    do i=0,n+1
        write(20,'(2e20.10)') x(i), u(i)
    enddo
    close(20)

end program jacobi_omp2