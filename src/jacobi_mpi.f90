! Solve Poisson equation
!  u_{xx} = f(x)    x \in [a, b]
! with 
!  u(a) = alpha, u(b) = beta
! using Jacobi iteration using MPI.
program jacobi_mpi
    use mpi

    implicit none

    integer, parameter :: maxiter = 100000, nprint = 5000
    real (kind=8), parameter :: alpha = 0.d0, beta = 3.d0

    integer :: i, iter, istart, iend, points_per_task, itask, n
    integer :: ierr, ntasks, me, req1, req2
    integer, dimension(MPI_STATUS_SIZE) :: mpistatus
    real (kind = 8), dimension(:), allocatable :: f, u, uold
    real (kind = 8) :: x, dumax_task, dumax_global, dx, tol

    ! Initialize MPI; get total number of tasks and ID of this task
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, ntasks, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, me, ierr)

    ! Ask the user for the number of points
    ! if (me == 0) then
    !     print *, "Input n ... "
    !     read *, n
    ! end if
    n = 100
    ! Broadcast to all tasks; everybody gets the value of n from task 0
    call mpi_bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    dx = 1.d0/real(n+1, kind=8)
    tol = 0.1d0*dx**2

    ! Determine how many points to handle with each task
    points_per_task = (n + ntasks - 1)/ntasks
    if (me == 0) then   ! Only one task should print to avoid clutter
        print *, "points_per_task = ", points_per_task
    end if

    ! Determine start and end index for this task's points
    istart = me * points_per_task + 1
    iend = min((me + 1)*points_per_task, n)

    ! Diagnostic: tell the user which points will be handled by which task
    print '("Task ",i2," will take i = ",i6," through i = ",i6)', &
        me, istart, iend


    ! Initialize:
    ! -----------

    ! This makes the indices run from istart-1 to iend+1
    ! This is more or less cosmetic, but makes things easier to think about
    allocate(f(istart-1:iend+1), u(istart-1:iend+1), uold(istart-1:iend+1))

    ! Each task sets its own, independent array
    do i = istart, iend
        ! Each task is a single thread with all its variables private
        ! so re-using the scalar variable x from one loop iteration to
        ! the next does not produce a race condition.
        x = dx*real(i, kind=8)
        f(i) = exp(x)               ! Source term
        u(i) = alpha + x*(beta - alpha)    ! Initial guess
    end do
    
    ! Set boundary conditions if this task is keeping track of a boundary
    ! point
    if (me == 0)        u(istart-1) = alpha
    if (me == ntasks-1) u(iend+1)   = beta


    ! Jacobi iteratation:
    ! -------------------

    do iter = 1, maxiter
        uold = u

        ! Send endpoint values to tasks handling neighboring sections
        ! of the array.  Note that non-blocking sends are used; note
        ! also that this sends from uold, so the buffer we're sending
        ! from won't be modified while it's being sent.
        !
        ! tag=1 is used for messages sent to the left
        ! tag=2 is used for messages sent to the right

        if (me > 0) then
            ! Send left endpoint value to process to the "left"
            call mpi_isend(uold(istart), 1, MPI_DOUBLE_PRECISION, me - 1, &
                1, MPI_COMM_WORLD, req1, ierr)
        end if
        if (me < ntasks-1) then
            ! Send right endpoint value to process on the "right"
            call mpi_isend(uold(iend), 1, MPI_DOUBLE_PRECISION, me + 1, &
                2, MPI_COMM_WORLD, req2, ierr)
        end if

        ! Accept incoming endpoint values from other tasks.  Note that
        ! these are blocking receives, because we can't run the next step
        ! of the Jacobi iteration until we've received all the
        ! incoming data.

        if (me < ntasks-1) then
            ! Receive right endpoint value
            call mpi_recv(uold(iend+1), 1, MPI_DOUBLE_PRECISION, me + 1, &
                1, MPI_COMM_WORLD, mpistatus, ierr)
        end if
        if (me > 0) then
            ! Receive left endpoint value
            call mpi_recv(uold(istart-1), 1, MPI_DOUBLE_PRECISION, me - 1, &
                2, MPI_COMM_WORLD, mpistatus, ierr)
        end if

        dumax_task = 0.d0   ! Max seen by this task

        ! Apply Jacobi iteration on this task's section of the array
        do i = istart, iend
            u(i) = 0.5d0*(uold(i-1) + uold(i+1) - dx**2*f(i))
            dumax_task = max(dumax_task, abs(u(i) - uold(i)))
        end do

        ! Take global maximum of dumax values
        call mpi_allreduce(dumax_task, dumax_global, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, MPI_COMM_WORLD, ierr)
        ! Note that this MPI_ALLREDUCE call acts as an implicit barrier,
        ! since no process can return from it until all processes
        ! have called it.  Because of this, after this call we know
        ! that all the send and receive operations initiated at the
        ! top of the loop have finished -- all the MPI_RECV calls have
        ! finished in order for each process to get here, and if the
        ! MPI_RECV calls have finished, the corresponding MPI_ISEND
        ! calls have also finished.  Thus we can safely modify uold
        ! again.

        ! Also periodically report progress to the user
        if (me == 0) then
            if (mod(iter, nprint)==0) then
                print '("After ",i8," iterations, dumax = ",d16.6,/)', &
                    iter, dumax_global
            end if
        end if

        ! All tasks now have dumax_global, and can check for convergence
        if (dumax_global < tol) exit
    end do

    print '("Task number ",i2," finished after ",i9," iterations, dumax = ",&
            e16.6)', me, iter, dumax_global


    ! Output result:
    ! --------------

    ! Note: this only works if all processes share a file system
    ! and can all open and write to the same file!

    ! Synchronize to keep the next part from being non-deterministic
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    ! Have each task output to a file in sequence, using messages to
    ! coordinate

    if (me == 0) then    ! Task 0 goes first
        ! Open file for writing, replacing any previous version:
        open(unit=20, file="jacobi_mpi.txt", status="replace")
        write(20, '(2e20.10)') 0.d0, u(0)    ! Boundary value at left end

        do i = istart, iend
            write(20, '(2e20.10)') i*dx, u(i)
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

        if (ntasks > 1) then
            call mpi_send(MPI_BOTTOM, 0, MPI_INTEGER, 1, 4, &
                          MPI_COMM_WORLD, ierr)
            endif

    else
        ! Wait for go-ahead message from previous task
        call mpi_recv(MPI_BOTTOM, 0, MPI_INTEGER, me - 1, 4, &
                          MPI_COMM_WORLD, mpistatus, ierr)
        ! Open file for appending; do not destroy previous contents
        open(unit=20, file="jacobi_mpi.txt", status="old", access="append")
        do i = istart, iend
            write(20, '(2e20.10)') i*dx, u(i)
        end do

        ! Boundary value at right end:
        if (me == ntasks - 1) write(20, '(2e20.10)') 1.d0, u(iend+1)    

        ! Flush all pending writes to disk
        close(unit=20)

        if (me < ntasks - 1) then
            ! Send go-ahead message to next task
            call mpi_send(MPI_BOTTOM, 0, MPI_INTEGER, me + 1, 4, &
                          MPI_COMM_WORLD, ierr)
        end if
    end if

    ! Close out MPI
    call mpi_finalize(ierr)

end program jacobi_mpi