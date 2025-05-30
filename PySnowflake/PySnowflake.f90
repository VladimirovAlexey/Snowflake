!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! This module is the python interface to the snoflake, designed for fits                   !!
!!                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module PySnowFlake
use SnowFlake
use Snowflake_Model

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GENERAL
subroutine Initialize(file)
character(len=*)::file

call  SnowFlake_Initialize(file)
call SnowFlake_Model_Initialize()

end subroutine Initialize

!!!!! update the evolution tables
subroutine UpdateEvolutionTable(mu0,mu1)
real*8,intent(in)::mu0,mu1

!call ComputeEvolution(mu0,mu1,alpha,U1=Tu,D1=Td,S1=Ts,U2=dTu,D2=dTd,S2=dTs,G1=Tp,G2=Tm,inputQ="T",inputG="T")
call ComputeEvolution(mu0,mu1,alpha,U1=SplusU,D1=SplusD,S1=SplusR,U2=SminusU,D2=SminusD,S2=SminusR,G1=Tp,G2=Tm,&
    inputQ="C",inputG="T")

end subroutine UpdateEvolutionTable

subroutine UpdateNPparameters(lambdaIN)
real*8,intent(in)::lambdaIN(:)

call SetNPparameters(lambdaIN)

end subroutine UpdateNPparameters


!!!!!!!!!!!!!!!!!!!!!!!!!!!!! call routines
function SnowFlake_G2_List(X,Q,f,ListLength)
integer,intent(in)::ListLength
real*8,intent(in),dimension(:)::X               !x-Bjorken
real*8,intent(in),dimension(:)::Q               !Q
integer,intent(in),dimension(:)::f              !flavor/process
real*8,dimension(1:ListLength)::SnowFlake_G2_List

call G2_List(SnowFlake_G2_List,X,Q,f)

end function SnowFlake_G2_List

function SnowFlake_D2_List(Q,f,ListLength)
integer,intent(in)::ListLength
real*8,intent(in),dimension(:)::Q               !Q
integer,intent(in),dimension(:)::f              !flavor/process
real*8,dimension(1:ListLength)::SnowFlake_D2_List

call D2_List(SnowFlake_D2_List,Q,f)

end function SnowFlake_D2_List

end module PySnowFlake
