! Solve Poisson equation
!  u_{xx} = f(x)    x \in [a, b]
! with 
!  u(a) = alpha, u(b) = beta
! using Jacobi iteration and a fine-grain parallelism approach using OpenMP.
program jacobi_omp1

    use omp_lib
    
    implicit none
    
    ! Problem specification and storage
    real(kind=8) :: a, b, alpha, beta
    real(kind=8), dimension(:), allocatable :: x, u, u_new, f
    real(kind=8) :: dx, tolerance, du_max

    integer(kind=4) :: n, num_threads
    integer(kind=8) :: i, num_iterations
    integer(kind=8), parameter :: MAX_ITERATIONS = 100000
    real(kind=8) :: time(2)

    ! Boundaries
    a = 0.d0
    b = 1.d0
    alpha = 0.d0
    beta = 3.d0

    ! Specify number of threads to use:
    num_threads = 4
    !$ call omp_set_num_threads(num_threads)
    !$ print "('Using OpenMP with ',i3,' threads')", num_threads

    N = 100

    ! Allocate storage for boundary points too
    allocate(x(0:N + 1), u(0:N + 1), u_new(0:N + 1), f(0:N + 1))

    call cpu_time(time(1))

    ! grid spacing:
    dx = (b - a) / (N + 1.d0)

    ! Set iniital guess and construct the grid
    ! Note that here we are breaking up the problem already to the threads
    !$omp parallel do
    do i=0, N + 1
        ! grid points:
        x(i) = i * dx
        ! source term:
        f(i) = exp(x(i))
        ! initial guess (satisfies boundary conditions and sets them)
        u(i) = alpha + x(i) * (beta - alpha)
    enddo

    ! Tolerance
    tolerance = 0.1d0 * dx**2

    ! Main Jacobi iteration loop
    ! Copy old values to new
    num_iterations = 0
    du_max = 1d99
    do while (du_max >= tolerance .and. num_iterations <= MAX_ITERATIONS)
        du_max = 0.d0
        !$omp parallel do reduction(max : du_max)
        do i=1, N
            u_new(i) = 0.5d0 * (u(i-1) + u(i+1) - dx**2 * f(i))
            ! print *, abs(u(i) - u_old(i))
            du_max = max(du_max, abs(u_new(i) - u(i)))
        end do

        if (mod(num_iterations, 1000) == 0) then
            print *, "du_max, iteration = ", du_max, num_iterations
        end if

        ! Copy old data into new
        !$omp parallel do 
        do i=1, N
            u(i) = u_new(i)
        end do
        num_iterations = num_iterations + 1
    end do

    call cpu_time(time(2))
    if (num_iterations > MAX_ITERATIONS) then
        print *, "Iteration failed!"
        stop
    end if
    print '("CPU time = ",f12.8, " seconds")', time(2) - time(1)
    print *, "Total number of iterations: ", num_iterations

    ! Write out solution for checking
    open(unit=20, file="jacobi_omp1.txt", status="unknown")
    do i=0,n+1
        write(20,'(2e20.10)') x(i), u(i)
    enddo
    close(20)

end program jacobi_omp1