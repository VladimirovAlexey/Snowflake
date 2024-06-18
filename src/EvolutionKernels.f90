!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! This module computes and saves the matrices of evolution kernels for the tw-3 evolution  !!
!!                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module EvolutionKernels
use IO_snowflake
use HexGrid
implicit none

private

INCLUDE 'commonVariables.f90'

!!!! GK_ORDER=1 = G7K15
!!!! GK_ORDER=1 = G20K41
#define GK_ORDER 2

#if GK_ORDER==1
INCLUDE 'Tables/G7K15.f90'
#elif GK_ORDER==2
INCLUDE 'Tables/G20K41.f90'
#endif

!!! QCD constants
!!!! color coefficient CA
real(dp),parameter::CA=3._dp
!!!! color coefficient CF
real(dp),parameter::CF=4._dp/3._dp
!!!! color coefficient CF-CA/2
real(dp),parameter::CFminusCA2=-1/6._dp
!!!! color coefficient ColorQG = - (n^2-4)/n
real(dp),parameter::ColorQG=-5/3._dp

type SparseM2
    integer, allocatable :: i(:)
    real(dp), allocatable :: val(:)
end type SparseM2

!!!! Evolution kernels for the quark. Hqq_J is the mixture contribution
type(SparseM2),allocatable,dimension(:):: Hqq, Hqq_J, Hqq_CO
!!!! Evolution kernels for the quark/gluons. Plus and minus components
type(SparseM2),allocatable,dimension(:):: Hqg_PLUS, Hgq_PLUS, Hgg_PLUS
type(SparseM2),allocatable,dimension(:):: Hqg_MINUS, Hgq_MINUS, Hgg_MINUS

public:: EvolutionKernels_Initialize,EvSingletPLUS,EvSingletMINUS,EvNonSinglet,EvChiralOdd
public:: SaveKernels,ReadKernels

contains

INCLUDE 'ExpressionsForKernels.f90'

subroutine EvolutionKernels_Initialize(path)
real*8::t1,t2,tt1,tt2
!$ real*8::omp_get_wtime
integer::sparsity
character(len=*)::path
logical::readKERNEL
character(len=300)::pathToKernels

!----------------- reading ini-file --------------------------------------
OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
!!! Number of grid-nodes in the sector-angle
call MoveTO(51,'0   :')
read(51,*) showINI

!!! Number of grid-nodes in the sector-angle
call MoveTO(51,'1   :')
read(51,*) NUM_phi
NUM_perimeter=6*NUM_phi-1

!!! Number of grid-nodes in the radius
call MoveTO(51,'2   :')
read(51,*) NUM_R
NUM_TOT=(NUM_perimeter+1)*(NUM_R+1)

!!! Value of minimal radius
call MoveTO(51,'3   :')
read(51,*) xMIN

!!! Value of zero (used for comparison)
call MoveTO(51,'4   :')
read(51,*) zero

!!! Integration tolerance
call MoveTO(51,'5   :')
read(51,*) toleranceINT

!!! Maximum number of processors to use
call MoveTO(51,'6   :')
read(51,*) RGh_max

!!! Maximal size of step in the Runge-Kutta iteraction (for t=log[mu])
call MoveTO(51,'7   :')
read(51,*) allowedNumProcessor

!!! Initialize Chiral-Even kernels
call MoveTO(51,'8   :')
read(51,*) IncludeChiralEvenEvolution
!!! Initialize Chiral-Odd kernels
call MoveTO(51,'9   :')
read(51,*) IncludeChiralOddEvolution
!!! Take into account gluon and flavor mixing (singlet evolution)
call MoveTO(51,'10  :')
read(51,*) useSingletEvolution
!!! Mass of CHARM threshold [GeV]
call MoveTO(51,'11  :')
read(51,*) massCHARM
!!! Mass of BOTTOM threshold [GeV]
call MoveTO(51,'12  :')
read(51,*) massBOTTOM

!!! Read kernels from files
call MoveTO(51,'13  :')
read(51,*) readKERNEL
if(readKERNEL) then
    !!! Path to the kernels
    call MoveTO(51,'14  :')
    read(51,'(A)') pathToKernels
end if

CLOSE (51, STATUS='KEEP')

!$ call OMP_set_num_threads(allowedNumProcessor)

if(showINI) then
    write(*,*) "Estimated zero:",zero
    write(*,*) "Tolerance of integration:",toleranceINT
    write(*,*) "Maximal step of Runge-Kutta:",RGh_max
#if GK_ORDER==1
write(*,*) ">>> Using GK7/15 for integration"
#elif GK_ORDER==2
write(*,*) ">>> Using GK20/41 for integration"
#else
write(*,*) ">>> Using GK7/15 for integration"
#endif
!$ write(*,'("Parallel computation is ",A," with ",I3," processors.")') color("ON",c_yellow),allowedNumProcessor
end if

if(readKERNEL) then
    call ReadKernels(pathToKernels)
else

    call cpu_time(t1)
    !$ t1=omp_get_wtime()


    if(IncludeChiralEvenEvolution) then
            !---------------------------- Hqq
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hqq kernels ..."
        call PreComputeMatrix(kernelHqq, Hqq, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hqq =",real(sparsity,dp)/NUM_TOT/NUM_TOT

        if(useSingletEvolution) then
        !!!!! in thecase of Singlet evolution mixing with gluon is accounted

        !---------------------------- Hqq - J
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hqq' kernels ..."
        call PreComputeMatrix(kernelHqq_J, Hqq_J, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hqq =",real(sparsity,dp)/NUM_TOT/NUM_TOT


        !---------------------------- Hqg plus
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hqg+ kernels ..."
        call PreComputeMatrix(kernelHqg_plus, Hqg_PLUS, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hqg+ =",real(sparsity,dp)/NUM_TOT/NUM_TOT

        !---------------------------- Hgq plus
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hgq+ kernels ..."
        call PreComputeMatrix(kernelHgq_plus, Hgq_PLUS, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hgq+ =",real(sparsity,dp)/NUM_TOT/NUM_TOT

        !---------------------------- Hgg plus
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hgg+ kernels ..."
        call PreComputeMatrix(kernelHgg_plus, Hgg_PLUS, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hgg+ =",real(sparsity,dp)/NUM_TOT/NUM_TOT

        !---------------------------- Hqg plus
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hqg- kernels ..."
        call PreComputeMatrix(kernelHqg_minus, Hqg_MINUS, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hqg- =",real(sparsity,dp)/NUM_TOT/NUM_TOT

        !---------------------------- Hgq plus
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hgq- kernels ..."
        call PreComputeMatrix(kernelHgq_minus, Hgq_MINUS, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hgq- =",real(sparsity,dp)/NUM_TOT/NUM_TOT

        !---------------------------- Hgg plus
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hgg- kernels ..."
        call PreComputeMatrix(kernelHgg_minus, Hgg_MINUS, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hgg- =",real(sparsity,dp)/NUM_TOT/NUM_TOT
        end if
    end if

    if(IncludeChiralOddEvolution) then
        !---------------------------- Hqq - J
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hqq (chiral odd)  kernels ..."
        call PreComputeMatrix(kernelHqq_CO, Hqq_CO, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hqq_CO =",real(sparsity,dp)/NUM_TOT/NUM_TOT
    end if

    call cpu_time(t2)
    !$ t2=omp_get_wtime()

    if(showINI) write(*,*) "--------- All necesary kernels computed. Total time =",t2-t1
end if

end subroutine EvolutionKernels_Initialize

!!!! This subroutine computes the kernel-matrix by the function Mfunction
!!!! and stores it into the COMPACT-matrix M
subroutine PreComputeMatrix(Mfunction,M,sparsity)
type(SparseM2),allocatable, dimension(:)::M
real(dp)::Mfunction
integer::n,k,nn,kk,counter,c1,c2,sparsity
real(dp)::x1,x2
integer::OMP_GET_THREAD_NUM,thread_id

real(dp),dimension(0:NUM_R,0:NUM_perimeter):: Hinter

sparsity=0

allocate(M(0:NUM_TOT-1))

!$OMP PARALLEL DO private(counter,c1,thread_id,x1,x2,k,nn,kk,Hinter)
do n=0,NUM_R
do k=0,NUM_perimeter
    thread_id = OMP_GET_THREAD_NUM()


    !!!! first: compute the matrix for element (n,k)
    counter=0
    do nn=0,NUM_R
    do kk=0,NUM_perimeter

        call NKtoX12(n,k,x1,x2)
        Hinter(nn,kk)=Mfunction(nn,kk,x1,x2)

        if(abs(Hinter(nn,kk))>zero) counter=counter+1
    end do
    end do

    !!!! second: allocate the value of M
    c1=Index1D(n,k)
    allocate(M(c1)%i(1:counter))
    allocate(M(c1)%val(1:counter))

    !!!! third: fill the elements of M
    counter=1
    do nn=0,NUM_R
    do kk=0,NUM_perimeter
        if(abs(Hinter(nn,kk))>zero) then
            M(c1)%i(counter)=Index1D(nn,kk)
            M(c1)%val(counter)=Hinter(nn,kk)
            counter=counter+1
        end if
    end do
    end do

end do
end do
!$OMP END PARALLEL DO

sparsity=0
do n=0,NUM_TOT-1
    sparsity=sparsity+size(M(n)%i)
end do

end subroutine PreComputeMatrix

!!!!!! Save Kernel into a directory
subroutine SaveKernels(path)
character(len=*)::path
write(*,'(A)',advance="no") "Saving Hqq ..."
call SaveKernelMatrix(Hqq,trim(path)//trim("Hqq.ker"))
write(*,'(A)',advance="no") "Saving Hqq_J ..."
call SaveKernelMatrix(Hqq_J,trim(path)//trim("Hqq_J.ker"))
write(*,'(A)',advance="no") "Saving Hqq_CO ..."
call SaveKernelMatrix(Hqq_CO,trim(path)//trim("Hqq_CO.ker"))
write(*,'(A)',advance="no") "Saving Hqg_PLUS ..."
call SaveKernelMatrix(Hqg_PLUS,trim(path)//trim("Hqg_PLUS.ker"))
write(*,'(A)',advance="no") "Saving Hgq_PLUS ..."
call SaveKernelMatrix(Hgq_PLUS,trim(path)//trim("Hgq_PLUS.ker"))
write(*,'(A)',advance="no") "Saving Hgg_PLUS ..."
call SaveKernelMatrix(Hgg_PLUS,trim(path)//trim("Hgg_PLUS.ker"))
write(*,'(A)',advance="no") "Saving Hqg_MINUS ..."
call SaveKernelMatrix(Hqg_MINUS,trim(path)//trim("Hqg_MINUS.ker"))
write(*,'(A)',advance="no") "Saving Hgq_MINUS ..."
call SaveKernelMatrix(Hgq_MINUS,trim(path)//trim("Hgq_MINUS.ker"))
write(*,'(A)',advance="no") "Saving Hgg_MINUS ..."
call SaveKernelMatrix(Hgg_MINUS,trim(path)//trim("Hgg_MINUS.ker"))
write(*,'(A)') "All kernels are stored in "//trim(path)

end subroutine SaveKernels

!!!!!! Save matrix For the Kernel into a file
subroutine SaveKernelMatrix(M,path)
type(SparseM2),dimension(0:NUM_TOT-1),intent(in)::M
character(len=*)::path
integer::i,j

write(*,'(A)') trim(path)

OPEN(UNIT=51, FILE=path, ACTION="write", STATUS="replace")
write(51,*) NUM_phi
write(51,*) NUM_R
do i=0,NUM_TOT-1
write(51,*) i, size(M(i)%i)
end do
do i=0,NUM_TOT-1
do j=1,size(M(i)%i)
write(51,*) M(i)%i(j), M(i)%val(j)
end do
end do

CLOSE(51, STATUS='KEEP')

end subroutine SaveKernelMatrix

!!!!!! Read all Kernels from the directory
subroutine ReadKernels(path)
type(SparseM2),dimension(0:NUM_TOT-1)::M22
character(len=*)::path


if(showINI) write(*,*) "Reading Hqq ..."
if(allocated(Hqq)) deallocate(Hqq)
allocate(Hqq(0:NUM_TOT-1))
call ReadKernelMatrix(Hqq,trim(path)//trim("Hqq.ker"))

if(showINI) write(*,*) "Reading Hqq_J ..."
if(allocated(Hqq_J)) deallocate(Hqq_J)
allocate(Hqq_J(0:NUM_TOT-1))
call ReadKernelMatrix(Hqq_J,trim(path)//trim("Hqq_J.ker"))

if(showINI) write(*,*) "Reading Hqq_CO ..."
if(allocated(Hqq_CO)) deallocate(Hqq_CO)
allocate(Hqq_CO(0:NUM_TOT-1))
call ReadKernelMatrix(Hqq_CO,trim(path)//trim("Hqq_CO.ker"))

if(showINI) write(*,*) "Reading Hqg_PLUS ..."
if(allocated(Hqg_PLUS)) deallocate(Hqg_PLUS)
allocate(Hqg_PLUS(0:NUM_TOT-1))
call ReadKernelMatrix(Hqg_PLUS,trim(path)//trim("Hqg_PLUS.ker"))

if(showINI) write(*,*) "Reading Hgq_PLUS ..."
if(allocated(Hgq_PLUS)) deallocate(Hgq_PLUS)
allocate(Hgq_PLUS(0:NUM_TOT-1))
call ReadKernelMatrix(Hgq_PLUS,trim(path)//trim("Hgq_PLUS.ker"))

if(showINI) write(*,*) "Reading Hgg_PLUS ..."
if(allocated(Hgg_PLUS)) deallocate(Hgg_PLUS)
allocate(Hgg_PLUS(0:NUM_TOT-1))
call ReadKernelMatrix(Hgg_PLUS,trim(path)//trim("Hgg_PLUS.ker"))

if(showINI) write(*,*) "Reading Hqg_MINUS ..."
if(allocated(Hqg_MINUS)) deallocate(Hqg_MINUS)
allocate(Hqg_MINUS(0:NUM_TOT-1))
call ReadKernelMatrix(Hqg_MINUS,trim(path)//trim("Hqg_MINUS.ker"))

if(showINI) write(*,*) "Reading Hgq_MINUS ..."
if(allocated(Hgq_MINUS)) deallocate(Hgq_MINUS)
allocate(Hgq_MINUS(0:NUM_TOT-1))
call ReadKernelMatrix(Hgq_MINUS,trim(path)//trim("Hgq_MINUS.ker"))

if(showINI) write(*,*) "Reading Hgg_MINUS ..."
if(allocated(Hgg_MINUS)) deallocate(Hgg_MINUS)
allocate(Hgg_MINUS(0:NUM_TOT-1))
call ReadKernelMatrix(Hgg_MINUS,trim(path)//trim("Hgg_MINUS.ker"))
end subroutine ReadKernels

!!!!!! Read matrix For the Kernel from a file
subroutine ReadKernelMatrix(M,path)
type(SparseM2),dimension(0:NUM_TOT-1)::M
character(len=*)::path
integer::i,j
integer::int1,int2
real(dp)::real1
logical::fileExist

inquire( file=trim(path), exist=fileExist )
if(.not.fileExist) then
    write(*,*) ErrorString("file for kernel is not found. Check the path."," ")
    write(*,*) trim(path)
    stop
end if

OPEN(UNIT=51, FILE=trim(path), ACTION="read", STATUS="old")
!!! read size
read(51,*) int1
read(51,*) int2
if(NUM_phi/=int1 .or. NUM_R/=int2) then
    write(*,*) ErrorString("Saved kernels use different grid-sizes than decalred in INI-file"," ")
    write(*,'("Files = (",I4 ," x ",I4,") vs. INI-file (",I4 ," x ",I4,")")') int1,int2,NUM_phi,NUM_R
    write(*,*) "Setup INI-file appropriately. Evaluation STOP"
    CLOSE(51, STATUS='KEEP')
    stop
end if

!!! read sizes of individual SparseM2, and allocate them
do i=0,NUM_TOT-1
    if(allocated(M(i)%i)) deallocate(M(i)%i)
    if(allocated(M(i)%val)) deallocate(M(i)%val)
    read(51,*) int1, int2
    if(int1/=i) then
        write(*,*) ErrorString("Kernel file corrupted. ERROR1"," ")
        write(*,*) trim(path)
        write(*,*) "Evaluation STOP"
        CLOSE(51, STATUS='KEEP')
        stop
    end if

    allocate(M(i)%i(1:int2),M(i)%val(1:int2))
end do

!!! Finally read enties
do i=0,NUM_TOT-1
   do j=1,size(M(i)%i)
    read(51,*) int1,real1
    M(i)%i(j)=int1
    M(i)%val(j)=real1
   end do
end do

CLOSE(51, STATUS='KEEP')

end subroutine ReadKernelMatrix


!!!!!!! composition of elementary kernels
!!!----- quark-quark part
function kernelHqq(i,j,x1,x2)
real(dp)::kernelHqq
integer,intent(in)::i,j
real(dp),intent(in)::x1,x2
     kernelHqq=CA*(H12hat(i,j,x1,x2)+H23hat(i,j,x1,x2)-2._dp*H12plus(i,j,x1,x2))&
        -1._dp/CA*(H13hat(i,j,x1,x2)-H13plus(i,j,x1,x2)-H23eP(i,j,x1,x2)+2*H12minus(i,j,x1,x2))&
        -3._dp*CF*Windex(i,j,x1,x2)
end function kernelHqq

!!!----- quark-quark mixture part
function kernelHqq_J(i,j,x1,x2)
real(dp)::kernelHqq_J
integer,intent(in)::i,j
real(dp),intent(in)::x1,x2
    kernelHqq_J=4*H13d(i,j,x1,x2)
end function kernelHqq_J

!!!----- quark-quark part (chiral odd)
function kernelHqq_CO(i,j,x1,x2)
real(dp)::kernelHqq_CO
integer,intent(in)::i,j
real(dp),intent(in)::x1,x2
    kernelHqq_CO=CA*(H12hat(i,j,x1,x2)+H23hat(i,j,x1,x2)-2*H12plus(i,j,x1,x2)-2*H23plus(i,j,x1,x2))&
       -1._dp/CA*(H13hat(i,j,x1,x2)+2*H12minus(i,j,x1,x2)+2*H23minus(i,j,x1,x2))&
       -3._dp*CF*Windex(i,j,x1,x2)
end function kernelHqq_CO

!!!----- gluon-gluon mixture part [WITHOUT BETA0 PART]
function kernelHgg_plus(i,j,x1,x2)
real(dp)::kernelHgg_plus
integer,intent(in)::i,j
real(dp),intent(in)::x1,x2
    kernelHgg_plus=CA*(H12hat_FFF(i,j,x1,x2)+H23hat_FFF(i,j,x1,x2)+H13hat_FFF(i,j,x1,x2) &
       -4*(H12plus_FFF(i,j,x1,x2)+H13plus_FFF(i,j,x1,x2)) &
       -2*(H12tilde_FFF(i,j,x1,x2)+H13tilde_FFF(i,j,x1,x2)) &
       +6*(H12minus_FFF(i,j,x1,x2)+H13minus_FFF(i,j,x1,x2)))
end function kernelHgg_plus

!!!----- gluon-gluon mixture part [WITHOUT BETA0 PART]
function kernelHgg_minus(i,j,x1,x2)
real(dp)::kernelHgg_minus
integer,intent(in)::i,j
real(dp),intent(in)::x1,x2
    kernelHgg_minus=CA*(H12hat_FFF(i,j,x1,x2)+H23hat_FFF(i,j,x1,x2)+H13hat_FFF(i,j,x1,x2) &
        -4*(H12plus_FFF(i,j,x1,x2)+H13plus_FFF(i,j,x1,x2)) &
        -2*(H12tilde_FFF(i,j,x1,x2)+H13tilde_FFF(i,j,x1,x2)) &
        -6*(H12minus_FFF(i,j,x1,x2)+H13minus_FFF(i,j,x1,x2)))
end function kernelHgg_minus

!!!----- quark-gluon mixture part
function kernelHqg_plus(i,j,x1,x2)
real(dp)::kernelHqg_plus
integer,intent(in)::i,j
real(dp),intent(in)::x1,x2
    kernelHqg_plus=V13plus(i,j,x1,x2)-V13minus(i,j,x1,x2)
end function kernelHqg_plus

!!!----- quark-gluon mixture part
function kernelHqg_minus(i,j,x1,x2)
real(dp)::kernelHqg_minus
integer,intent(in)::i,j
real(dp),intent(in)::x1,x2
    kernelHqg_minus=V13plus(i,j,x1,x2)+V13minus(i,j,x1,x2)
end function kernelHqg_minus

!!!----- gluon-quark mixture part
function kernelHgq_plus(i,j,x1,x2)
real(dp)::kernelHgq_plus
integer,intent(in)::i,j
real(dp),intent(in)::x1,x2
    kernelHgq_plus=CA*(Wplus(i,j,x1,x2)+Wminus(i,j,x1,x2)-2*DeltaW(i,j,x1,x2)&
                       - WplusP(i,j,x1,x2)-WminusP(i,j,x1,x2)+2*DeltaWP(i,j,x1,x2))
end function kernelHgq_plus

!!!----- gluon-quark mixture part
function kernelHgq_minus(i,j,x1,x2)
real(dp)::kernelHgq_minus
integer,intent(in)::i,j
real(dp),intent(in)::x1,x2
    kernelHgq_minus=ColorQG*(Wplus(i,j,x1,x2)+Wminus(i,j,x1,x2)+ WplusP(i,j,x1,x2)+WminusP(i,j,x1,x2))
end function kernelHgq_minus

!!!! multiplies the sparse ``H'' by ``vector'' F
function HxF(H,F)
    real(dp),dimension(0:NUM_TOT-1)::HxF
    real(dp),dimension(0:NUM_TOT-1),intent(in)::F
    type(SparseM2),dimension(0:NUM_TOT-1),intent(in)::H

    integer::c,n

    !$OMP PARALLEL DO private(n)
    do c=0,NUM_TOT-1
        HxF(c)=0._dp
        do n=1,size(H(c)%i)
            HxF(c)=HxF(c)+H(c)%val(n)*F(H(c)%i(n))
        end do
    end do
    !$OMP END PARALLEL DO
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Runge Kuta-4 within vector space !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! various combinations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! All functions have the same interface
!!!! F is the initial boundary condition as the vector in S-space
!!!! alpha is the a(t) is alpha_s/4pi(log mu^2)
!!!! t0, t1 are initial and final scales of evolution
!!!! nf is the active number of flavors
!!!! output is placed into F0

subroutine EvNonSinglet(F,alpha,t0,t1)
    real(dp),dimension(0:NUM_TOT-1)::F
    real(dp)::alpha
    real(dp),intent(in)::t0,t1

    real(dp)::tn,RGh
    real(dp),dimension(0:NUM_TOT-1)::k1,k2,k3,k4
    integer::i,numSteps

   !!! number of steps, and the step size
    numSteps=int((t1-t0)/RGh_max)!+1
    RGh=(t1-t0)/numSteps

    tn=t0
    do i=1,numSteps
        !!! the factor -1 is the factor in the evolution equation

        k1=-RGh*alpha(tn)*HxF(Hqq,F)
        k2=-RGh*alpha(tn+RGh/2)*HxF(Hqq,F+k1/2)
        k3=-RGh*alpha(tn+RGh/2)*HxF(Hqq,F+k2/2)
        k4=-RGh*alpha(tn+RGh)*HxF(Hqq,F+k3)

        F=F+(k1+2*k2+2*k3+k4)/6
        tn=tn+RGh
    end do
end subroutine EvNonSinglet

subroutine EvSingletPLUS(F,Fg,alpha,t0,t1,nf)
    real(dp),dimension(0:NUM_TOT-1)::F,Fg
    real(dp)::alpha
    real(dp),intent(in)::t0,t1
    integer,intent(in)::nf

    real(dp)::tn,RGh
    real(dp),dimension(0:NUM_TOT-1)::k1,k2,k3,k4,k1g,k2g,k3g,k4g,F_dum,Fg_dum
    integer::i,numSteps

    !!! number of steps, and the step size
    numSteps=int((t1-t0)/RGh_max)+1
    RGh=(t1-t0)/numSteps
    tn=t0
    do i=1,numSteps
        !!! the factor -1 is the factor in the evolution equation
        k1=-RGh*alpha(tn)*(HxF(Hqq,F)+Nf*(HxF(Hqq_J,F)+HxF(Hqg_PLUS,Fg)))
        k1g=-RGh*alpha(tn)*(HxF(Hgq_PLUS,F)+HxF(Hgg_PLUS,Fg)-(11._dp-2._dp*nf/3._dp)*Fg)

        F_dum=F+k1/2
        Fg_dum=Fg+k1g/2

        k2=-RGh*alpha(tn+RGh/2)*(HxF(Hqq,F_dum)+Nf*(HxF(Hqq_J,F_dum)+HxF(Hqg_PLUS,Fg_dum)))
        k2g=-RGh*alpha(tn+RGh/2)*(HxF(Hgq_PLUS,F_dum)+HxF(Hgg_PLUS,Fg_dum)-(11._dp-2._dp*nf/3._dp)*Fg_dum)

        F_dum=F+k2/2
        Fg_dum=Fg+k2g/2

        k3=-RGh*alpha(tn+RGh/2)*(HxF(Hqq,F_dum)+Nf*(HxF(Hqq_J,F_dum)+HxF(Hqg_PLUS,Fg_dum)))
        k3g=-RGh*alpha(tn+RGh/2)*(HxF(Hgq_PLUS,F_dum)+HxF(Hgg_PLUS,Fg_dum)-(11._dp-2._dp*nf/3._dp)*Fg_dum)

        F_dum=F+k3
        Fg_dum=Fg+k3g

        k4=-RGh*alpha(tn+RGh)*(HxF(Hqq,F_dum)+Nf*(HxF(Hqq_J,F_dum)+HxF(Hqg_PLUS,Fg_dum)))
        k4g=-RGh*alpha(tn+RGh)*(HxF(Hgq_PLUS,F_dum)+HxF(Hgg_PLUS,Fg_dum)-(11._dp-2._dp*nf/3._dp)*Fg_dum)

        F=F+(k1+2*k2+2*k3+k4)/6
        Fg=Fg+(k1g+2*k2g+2*k3g+k4g)/6
        tn=tn+RGh
     end do

end subroutine EvSingletPLUS

subroutine EvSingletMINUS(F,Fg,alpha,t0,t1,nf)
    real(dp),dimension(0:NUM_TOT-1)::F,Fg
    real(dp)::alpha
    real(dp),intent(in)::t0,t1

    real(dp)::tn,RGh
    real(dp),dimension(0:NUM_TOT-1)::k1,k2,k3,k4,k1g,k2g,k3g,k4g,F_dum,Fg_dum
    integer::nf
    integer::i,numSteps

    !!! number of steps, and the step size
    numSteps=int((t1-t0)/RGh_max)+1
    RGh=(t1-t0)/numSteps
    tn=t0

   do i=1,numSteps
        !!! the factor -1 is the factor in the evolution equation

        k1=-RGh*alpha(tn)*(HxF(Hqq,F)+Nf*HxF(Hqg_MINUS,Fg))
        k1g=-RGh*alpha(tn)*(HxF(Hgq_MINUS,F)+HxF(Hgg_MINUS,Fg)-(11._dp-2._dp*nf/3._dp)*Fg)

        F_dum=F+k1/2
        Fg_dum=Fg+k1g/2

        k2=-RGh*alpha(tn+RGh/2)*(HxF(Hqq,F_dum)+Nf*HxF(Hqg_MINUS,Fg_dum))
        k2g=-RGh*alpha(tn+RGh/2)*(HxF(Hgq_MINUS,F_dum)+HxF(Hgg_MINUS,Fg_dum)-(11._dp-2._dp*nf/3._dp)*Fg_dum)

        F_dum=F+k2/2
        Fg_dum=Fg+k2g/2

        k3=-RGh*alpha(tn+RGh/2)*(HxF(Hqq,F_dum)+Nf*HxF(Hqg_MINUS,Fg_dum))
        k3g=-RGh*alpha(tn+RGh/2)*(HxF(Hgq_MINUS,F_dum)+HxF(Hgg_MINUS,Fg_dum)-(11._dp-2._dp*nf/3._dp)*Fg_dum)

        F_dum=F+k3
        Fg_dum=Fg+k3g

        k4=-RGh*alpha(tn+RGh)*(HxF(Hqq,F_dum)+Nf*HxF(Hqg_MINUS,Fg_dum))
        k4g=-RGh*alpha(tn+RGh)*(HxF(Hgq_MINUS,F_dum)+HxF(Hgg_MINUS,Fg_dum)-(11._dp-2._dp*nf/3._dp)*Fg_dum)

        F=F+(k1+2*k2+2*k3+k4)/6
        Fg=Fg+(k1g+2*k2g+2*k3g+k4g)/6

        tn=tn+RGh

    end do

end subroutine EvSingletMINUS

subroutine EvChiralOdd(F,alpha,t0,t1)
    real(dp),dimension(0:NUM_TOT-1)::F
    real(dp)::alpha
    real(dp),intent(in)::t0,t1

    real(dp)::tn,RGh
    real(dp),dimension(0:NUM_TOT-1)::k1,k2,k3,k4
    integer::i,numSteps

    !!! number of steps, and the step size
    numSteps=int((t1-t0)/RGh_max)+1
    RGh=(t1-t0)/numSteps
    tn=t0

    !!! all steps except the last
    do i=1,numSteps
        k1=-RGh*alpha(tn)*HxF(Hqq_CO,F)
        k2=-RGh*alpha(tn+RGh/2)*HxF(Hqq_CO,F+k1/2)
        k3=-RGh*alpha(tn+RGh/2)*HxF(Hqq_CO,F+k2/2)
        k4=-RGh*alpha(tn+RGh)*HxF(Hqq_CO,F+k3)

        F=F+(k1+2*k2+2*k3+k4)/6
        tn=tn+RGh
    end do

end subroutine EvChiralOdd

!!!!-------------------------------------------------------------------------------------------------------------
!!! Gauss-Kronrod 7/15 adaptive
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! tolerance is relative (it is wieghted by approximate value of integral)
function Integrate_GK(f,xMin,xMax)
    real(dp)::f,Integrate_GK
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,g7,k15,eps,fI
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp

    if(delta<zero) then
        Integrate_GK=0._dp
        return
    end if

    g7=0._dp
    k15=0._dp
#if GK_ORDER==1
    do i=1,15
        fI=f(Xi_k15(i)*delta+av)
        g7=g7+Wi_g7(i)*fI
        k15=k15+Wi_k15(i)*fI
    end do
#elif  GK_ORDER==2
    do i=1,41
        fI=f(Xi_k41(i)*delta+av)
        g7=g7+Wi_g20(i)*fI
        k15=k15+Wi_k41(i)*fI
    end do
#endif

    eps=delta*abs(k15)*toleranceINT

!    write(*,*) "-->",delta*abs(k15-g7),eps,delta*k15
    if(abs(delta*k15)<zero) then
        Integrate_GK=delta*k15
    else if(delta*abs(k15-g7)>eps) then
        Integrate_GK=GK_Rec(f,xMin,av,eps)+GK_Rec(f,av,xMax,eps)
    else
        Integrate_GK=delta*k15
    end if

end function Integrate_GK

recursive function GK_Rec(f,xMin,xMax,eps) result(res)
    real(dp)::f,res
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,g7,k15,eps,fI
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp

    if(delta<zero) then
        res=0._dp
        return
    end if

    g7=0._dp
    k15=0._dp
#if GK_ORDER==1
    do i=1,15
        fI=f(Xi_k15(i)*delta+av)
        g7=g7+Wi_g7(i)*fI
        k15=k15+Wi_k15(i)*fI
    end do
#elif  GK_ORDER==2
    do i=1,41
        fI=f(Xi_k41(i)*delta+av)
        g7=g7+Wi_g20(i)*fI
        k15=k15+Wi_k41(i)*fI
    end do
#endif

    if(abs(delta*k15)<zero) then
        res=delta*k15
    else if(delta*abs(k15-g7)>eps) then
        res=GK_Rec(f,xMin,av,eps)+GK_Rec(f,av,xMax,eps)
    else
        res=delta*k15
    end if

end function GK_Rec
!!!!-------------------------------------------------------------------------------------------------------------

end module EvolutionKernels
