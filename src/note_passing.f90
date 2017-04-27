! Demonstration for message passing between MPI processes
program note_passing

    use mpi

    implicit none

    integer :: proc_id, num_procs, error, tag
    real(kind=8) :: important_note
    integer, dimension(MPI_STATUS_SIZE) :: status

    call mpi_init(error)
    call mpi_comm_size(MPI_COMM_WORLD, num_procs, error)
    call mpi_comm_rank(MPI_COMM_WORLD, proc_id, error)

    if (num_procs == 1) then
        print *, "Only one process, cannot do anything."
        call MPI_FINALIZE(error)
        stop
    endif

    ! Not really important in this case but important to sort through messages
    tag = 42

    ! If we are process 0 then set the value to be passed
    if (proc_id == 0) then
        important_note = 2.718281828d0
        print '("Process ",i3," note = ",e18.8)', proc_id, important_note

        call mpi_send(important_note, 1, MPI_DOUBLE_PRECISION, 1, tag, &
                      MPI_COMM_WORLD, error)

    ! If we are one of the processes in between pass it on to the next process
    else if (proc_id < num_procs - 1) then

        call MPI_RECV(important_note, 1, MPI_DOUBLE_PRECISION, proc_id-1, tag, &
                      MPI_COMM_WORLD, status, error)

        print '("Process ",i3," received note = ",e18.8)', proc_id, important_note

        call mpi_send(important_note, 1, MPI_DOUBLE_PRECISION, proc_id + 1, &
                      tag, MPI_COMM_WORLD, error)

        print '("Process ",i3,"    sent note = ",e18.8)', proc_id, important_note

    ! If we are the last process in the class to find out, well...
    else if (proc_id == num_procs - 1) then

        call mpi_recv(important_note, 1, MPI_DOUBLE_PRECISION, proc_id - 1, &
                      tag, MPI_COMM_WORLD, status, error)

        print '("Process ",i3," ends up with note = ",e18.8)', proc_id, important_note
      endif

    call MPI_FINALIZE(error)

end program note_passing
