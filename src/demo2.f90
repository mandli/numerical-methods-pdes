program demo2
   
   use omp_lib
   implicit none
   integer :: i
   integer(kind=8), parameter :: n = 2**16
   real(kind=8), dimension(n) :: x,y,z
   
   ! Specify number of threads to use:
   !$ call omp_set_num_threads(2)

   !$omp parallel  ! spawn two threads
   !$omp sections  ! split up work between them

     !$omp section
     x = 1.d0   ! one thread initializes x array

     !$omp section
     y = 1.d0   ! another thread initializes y array

   !$omp end sections
   !$omp barrier   ! not needed, implied at end of sections

   !$omp single    ! only want to print once:
   print *, "Done initializing x and y"
   !$omp end single nowait  ! ok for other thread to continue

   !$omp do   ! split work between threads:
   do i=1,n
       z(i) = x(i) + y(i)
       enddo

   !$omp end parallel
   print *, "max value of z is ",maxval(z)
    

end program demo2