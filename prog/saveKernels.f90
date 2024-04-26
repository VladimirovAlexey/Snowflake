!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Code that save the kernels for future use.      A.Vladimirov 01.04.2024                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program SaveKernelsToFile
use SnowFlake
use EvolutionKernels

implicit none

real*8::t1,t2
!$ real*8::omp_get_wtime

!!!! initialize the Snowflake with particular ini file
! call  SnowFlake_Initialize("kernels10x15.ini","prog/")
!
! call cpu_time(t1)
! !$ t1=omp_get_wtime()
! call SaveKernels("kernels/10x15/")
!
!
! call ReadKernels("kernels/10x15/")

!!! initialize the Snowflake with particular ini file
call  SnowFlake_Initialize("kernels20x25.ini","prog/")

call cpu_time(t1)
!$ t1=omp_get_wtime()
call SaveKernels("kernels/20x25/")


call ReadKernels("kernels/20x25/")
call cpu_time(t2)
!$ t2=omp_get_wtime()
write(*,*) "Time for computation of evolution",t2-t1

end program SaveKernelsToFile
