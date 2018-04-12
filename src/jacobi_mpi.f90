! Solve Poisson equation
!  u_{xx} = f(x)    x \in [a, b]
! with 
!  u(a) = alpha, u(b) = beta
! using Jacobi iteration using MPI.
program jacobi_mpi

    use mpi

    implicit none

    real (kind=8), parameter :: alpha = 0.d0, beta = 3.d0

    integer :: i, start, end, points_per_task, task_id, N
    integer :: error, num_tasks, my_id, req1, req2
    integer, dimension(MPI_STATUS_SIZE) :: mpi_status
    real(kind=8), dimension(:), allocatable :: f, u, u_old
    real(kind=8) :: x, du_max_task, du_max_global, dx, tolerance

    integer(kind=8) :: num_iterations
    integer(kind=8), parameter :: MAX_ITERATIONS = 2**32_8
    integer, parameter :: PRINT_INTERVAL = 5000

    ! Initialize MPI; get total number of tasks and ID of this task
    call mpi_init(error)
    call mpi_comm_size(MPI_COMM_WORLD, num_tasks, error)
    call mpi_comm_rank(MPI_COMM_WORLD, my_id, error)

    ! Ask the user for the number of points
    ! if (me == 0) then
    !     print *, "Input n ... "
    !     read *, n
    ! end if
    N = 100
    ! Broadcast to all tasks; everybody gets the value of n from task 0
    call mpi_bcast(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error)

    dx = 1.d0 / real(N + 1, kind=8)
    tolerance = 0.1d0*dx**2

    ! Determine how many points to handle with each task
    points_per_task = (N + num_tasks - 1) / num_tasks
    if (my_id == 0) then   ! Only one task should print to avoid clutter
        print *, "points_per_task = ", points_per_task
    end if

    ! Determine start and end index for this task's points
    start = my_id * points_per_task + 1
    end = min((my_id + 1) * points_per_task, n)

    ! Diagnostic: tell the user which points will be handled by which task
    print '("Task ",i2," will take i = ",i6," through i = ",i6)', &
        my_id, start, end


    ! Initialize:
    ! -----------

    ! This makes the indices run from istart-1 to end+1
    ! This is more or less cosmetic, but makes things easier to think about
    allocate(f(start - 1:end + 1))
    allocate(u(start - 1:end + 1))
    allocate(u_old(start - 1:end + 1))

    ! Each task sets its own, independent array
    do i = start, end
        ! Each task is a single thread with all its variables private
        ! so re-using the scalar variable x from one loop iteration to
        ! the next does not produce a race condition.
        x = dx * real(i, kind=8)
        f(i) = exp(x)                        ! RHS
        u(i) = alpha + x * (beta - alpha)    ! Initial guess
    end do
    
    ! Set boundary conditions if this task is keeping track of a boundary
    ! point
    if (my_id == 0) then
        u(start - 1) = alpha
    end if
    if (my_id == num_tasks - 1) then
        u(end + 1) = beta
    end if


    ! Jacobi iteratation:
    ! -------------------

    do num_iterations = 1, MAX_ITERATIONS
        u_old = u

        ! Send endpoint values to tasks handling neighboring sections
        ! of the array.  Note that non-blocking sends are used; note
        ! also that this sends from u_old, so the buffer we're sending
        ! from won't be modified while it's being sent.
        !
        ! tag = 1 is used for messages sent to the left
        ! tag = 2 is used for messages sent to the right

        if (my_id > 0) then
            ! Send left endpoint value to process to the "left"
            ! Note that this is a "non-blocking send"
            call mpi_isend(u_old(start), 1, MPI_DOUBLE_PRECISION, my_id - 1, &
                1, MPI_COMM_WORLD, req1, error)
        end if
        if (my_id < num_tasks - 1) then
            ! Send right endpoint value to process on the "right"
            ! Note that this is a "non-blocking send"
            call mpi_isend(u_old(end), 1, MPI_DOUBLE_PRECISION, my_id + 1, &
                2, MPI_COMM_WORLD, req2, error)
        end if

        ! Accept incoming endpoint values from other tasks.  Note that
        ! these are blocking receives, because we can't run the next step
        ! of the Jacobi iteration until we've received all the
        ! incoming data.

        if (my_id < num_tasks-1) then
            ! Receive right endpoint value
            call mpi_recv(u_old(end+1), 1, MPI_DOUBLE_PRECISION, my_id + 1, &
                1, MPI_COMM_WORLD, mpi_status, error)
        end if
        if (my_id > 0) then
            ! Receive left endpoint value
            call mpi_recv(u_old(start-1), 1, MPI_DOUBLE_PRECISION, my_id - 1, &
                2, MPI_COMM_WORLD, mpi_status, error)
        end if

        du_max_task = 0.d0   ! Max seen by this task

        ! Apply Jacobi iteration on this task's section of the array
        do i = start, end
            u(i) = 0.5d0*(u_old(i-1) + u_old(i+1) - dx**2*f(i))
            du_max_task = max(du_max_task, abs(u(i) - u_old(i)))
        end do

        ! Take global maximum of dumax values
        call mpi_allreduce(du_max_task, du_max_global, 1,       &
                           MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, error)
        ! Note that this MPI_ALLREDUCE call acts as an implicit barrier,
        ! since no process can return from it until all processes
        ! have called it.  Because of this, after this call we know
        ! that all the send and receive operations initiated at the
        ! top of the loop have finished -- all the MPI_RECV calls have
        ! finished in order for each process to get here, and if the
        ! MPI_RECV calls have finished, the corresponding MPI_ISEND
        ! calls have also finished.  Thus we can safely modify u_old
        ! again.

        ! Also periodically report progress to the user
        if (my_id == 0) then
            if (mod(num_iterations, PRINT_INTERVAL)==0) then
                print '("After ",i8," iterations, dumax = ",d16.6,/)', &
                    num_iterations, du_max_global
            end if
        end if

        ! All tasks now have du_max_global, and can check for convergence
        if (du_max_global < tolerance) exit
    end do

    print '("Task number ",i2," finished after ",i9," iterations, dumax = ",&
            e16.6)', my_id, num_iterations, du_max_global


    ! Output result:
    ! --------------

    ! Note: this only works if all processes share a file system
    ! and can all open and write to the same file!

    ! Synchronize to keep the next part from being non-deterministic
    call mpi_barrier(MPI_COMM_WORLD, error)

    ! Check to make sure that we were successful
    if (num_iterations >= MAX_ITERATIONS) then
        if (my_id == 0) then
            print *, "Jacobi failed to converge!"
            print *, "Reached du_max = ", du_max_global
            print *, "Tolerance = ", tolerance
        end if
        call mpi_finalize(error)
        stop
    endif

    ! Have each task output to a file in sequence, using messages to
    ! coordinate

    if (my_id == 0) then    ! Task 0 goes first
        ! Open file for writing, replacing any previous version:
        open(unit=20, file="jacobi_mpi.txt", status="replace")
        write(20, '(2e20.10)') 0.d0, u(0)    ! Boundary value at left end

        do i = start, end
            write(20, '(2e20.10)') i * dx, u(i)
        end do

        close(unit=20)
        ! Closing the file should guarantee that all the ouput 
        ! will be written to disk.
        ! If the file isn't closed before the next process starts writing,
        ! output may be jumbled or missing.

        ! Send go-ahead message to next task
        ! Only the fact that the message was sent is important, not its contents
        ! so we send the special address MPI_BOTTOM and length 0.
        ! tag=4 is used for this message.

        if (num_tasks > 1) then
            call mpi_send(MPI_BOTTOM, 0, MPI_INTEGER, 1, 4, &
                          MPI_COMM_WORLD, error)
            endif

    else
        ! Wait for go-ahead message from previous task
        call mpi_recv(MPI_BOTTOM, 0, MPI_INTEGER, my_id - 1, 4, &
                          MPI_COMM_WORLD, mpi_status, error)
        ! Open file for appending; do not destroy previous contents
        open(unit=20, file="jacobi_mpi.txt", status="old", access="append")
        do i = start, end
            write(20, '(2e20.10)') i * dx, u(i)
        end do

        ! Boundary value at right end:
        if (my_id == num_tasks - 1) write(20, '(2e20.10)') 1.d0, u(end + 1)    

        ! Flush all pending writes to disk
        close(unit=20)

        if (my_id < num_tasks - 1) then
            ! Send go-ahead message to next task
            call mpi_send(MPI_BOTTOM, 0, MPI_INTEGER, my_id + 1, 4, &
                          MPI_COMM_WORLD, error)
        end if
    end if

    ! Close out MPI
    call mpi_finalize(error)

end program jacobi_mpi
