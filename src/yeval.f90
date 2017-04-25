program yeval
   
   use omp_lib

   implicit none

   integer(kind=8), parameter :: n = 2**16
   integer(kind=4) :: i, nthreads
   real(kind=8), dimension(n) :: y
   real(kind=8) :: dx, x

   ! Specify number of threads to use:
   !$ print *, "How many threads to use? "
   !$ read *, nthreads
   !$ call omp_set_num_threads(nthreads)
   !$ print "('Using OpenMP with ',i3,' threads')", nthreads

   dx = 1.d0 / (n+1.d0)

   !$omp parallel do private(x) 
   do i=1, n
      x = i * dx
      y(i) = exp(x) * cos(x) * sin(x) * sqrt(5.d0 * x + 6.d0)
   enddo
   !$omp end parallel do

   print *, "Filled vector y of length", n

end program yeval