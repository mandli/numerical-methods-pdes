program hello_world_omp
    
    use omp_lib

    implicit none
    integer :: num_threads, thread_id

    !$OMP parallel private(num_threads, thread_id)
    !$ num_threads = omp_get_num_threads()
    !$ thread_id = omp_get_thread_num()
    print *, 'Hello from thread number', thread_id + 1, &
             ' of ', num_threads, ' processes'

    !$OMP end parallel

end program hello_world_omp
