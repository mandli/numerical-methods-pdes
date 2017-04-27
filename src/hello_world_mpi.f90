program hello_world_mpi
    
    use mpi

    implicit none
    integer :: error, num_procs, proc_num

    call mpi_init(error)
    call mpi_comm_size(MPI_COMM_WORLD, num_procs, error)
    call mpi_comm_rank(MPI_COMM_WORLD, proc_num, error)

    print *, 'Hello from Process number', proc_num + 1, &
             ' of ', num_procs, ' processes'

    call mpi_finalize(error)

end program hello_world_mpi
