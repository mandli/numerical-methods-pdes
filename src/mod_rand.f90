! ============================================================================
!  Program:     
!  File:        untitled.f95
!  Created:     2009-03-25
!  Author:      Kyle Mandli
! ============================================================================
!  Description:
!
! ============================================================================

module mod_rand

    implicit none
    
contains

    SUBROUTINE init_random_seed()
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        CALL SYSTEM_CLOCK(COUNT=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)

        DEALLOCATE(seed)
    END SUBROUTINE

end module mod_rand
