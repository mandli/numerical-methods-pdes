! Example program that computes pi in parallel
program compute_pi

    use mpi
    
    implicit none

    integer :: error, num_procs, proc_id, points_per_proc, n, i, start, end
    real(kind=8) :: x, dx, pi_sum, pi_sum_proc

    real(kind=8), parameter :: pi = 3.1415926535897932384626433832795d0

    call mpi_init(error)
    call mpi_comm_size(MPI_COMM_WORLD, num_procs, error)
    call mpi_comm_rank(MPI_COMM_WORLD, proc_id, error)

    ! Proc 0 will ask the user for the number of points
    if (proc_id == 0) then
        print *, "Using ",num_procs," processors"
        print *, "Input n ... "
        read *, n
    end if

    ! Broadcast to all procs; everybody gets the value of n from proc 0
    call mpi_bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error)

    dx = 1.d0 / n

    ! Determine how many points to handle with each proc
    points_per_proc = (n + num_procs - 1) / num_procs
    ! Only print out the number of points per proc by process 0
    if (proc_id == 0) then   
        print *, "points_per_proc = ", points_per_proc
    end if

    ! Determine start and end index for this proc's points
    start = proc_id * points_per_proc + 1
    end = min((proc_id + 1) * points_per_proc, n)

    ! Diagnostic: tell the user which points will be handled by which proc
    print '("Process ",i2," will take i = ",i8," through i = ",i8)', &
          proc_id, start, end

    pi_sum_proc = 0.d0
    do i=start,end
        x = (i - 0.5d0) * dx
        pi_sum_proc = pi_sum_proc + 1.d0 / (1.d0 + x**2)
    enddo

    call MPI_REDUCE(pi_sum_proc, pi_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                        MPI_COMM_WORLD, error)

    if (proc_id == 0) then
        print *, "The approximation to pi is ", 4.d0 * dx * pi_sum
        print *, "Difference = ", abs(pi - 4.d0 * dx * pi_sum)
    endif

    call mpi_finalize(error)

end program compute_pi