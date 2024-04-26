!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! This module operates with the grid and grid variables on the hexagon                     !!
!!                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module HexGrid
use IO_snowflake
implicit none

private

interface NKtoX12
    module procedure NKtoX12_int, NKtoX12_real
end interface

INCLUDE 'commonVariables.f90'

real(dp),parameter::pi=3.141592653589793238462643383279502884197169399375105820974944592_dp        !!pi

!!!!!! Pre-processor parameter that specifies the grid type
!!!!!! GRIDTYPE_R=1;  log-exp grid
!!!!!!      r_j=xMin^{(n-j)/n}, logarithm grid from xMin to 1.
!!!!!!
!!!!!! GRIDTYPE_R=2; root-power
!!!!!!      r_j=((j/n+c)/(1+c))^a, power grid from xMin to 1.
!!!!!!
!!!!!! GRIDTYPE_R=3; log-log
!!!!!!      r_j=2 exp(j-n/n/c)/(1+exp(..)), logarith+logarithm grid from xMin to 1.
!!!!!!
!!!!!! GRIDTYPE_R=4; cosh
!!!!!!      r_j=1/cosh^p(r), inverse cosh grid
!!!!!!
!!!!!! GRIDTYPE_PHI=1; linear
!!!!!!      phi_j=n*phi
!!!!!!
!!!!!! GRIDTYPE_PHI=2; cos-like
!!!!!!      phi_j=n/pi*arcos(1-2phi)
!!!!!!
#define GRIDTYPE_R 4
#define GRIDTYPE_PHI 1

real(dp)::grid_paramC
real(dp)::grid_param0

public::Initialize_HexGrid
public::Index1D,Index2D,F2Dto1D,F1Dto2D,FXYto1D,X12toNK,NKtoX12,NtoX12
public::FatXfrom2D,Windex,LimitsX1,LimitsX2,LimitsX3



contains

subroutine Initialize_HexGrid(path)
character(len=*)::path

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

CLOSE (51, STATUS='KEEP')


!!!! parameters of the grid
#if GRIDTYPE_R==1
grid_paramC=1/log(xMIN)
grid_param0=1._dp

#elif GRIDTYPE_R==2
grid_param0=3._dp
grid_paramC=1._dp/(xMIN**(-1/grid_param0)-1)

#elif GRIDTYPE_R==3
grid_param0=1._dp
grid_paramC=1/log(xMIN/(2-xMIN))

#elif GRIDTYPE_R==4
grid_param0=3._dp
grid_paramC=1/acosh(xMIN**(-1/grid_param0))

#else
grid_paramC=1/log(xMIN)
grid_param0=1._dp

#endif

if(showINI) then
#if GRIDTYPE_R==1
write(*,*) ">>> Using log-exp radius-grid."
#elif GRIDTYPE_R==2
write(*,*) ">>> Using root-power radius-grid."
#elif GRIDTYPE_R==3
write(*,*) ">>> Using log-log radius-grid."
#elif GRIDTYPE_R==4
write(*,*) ">>> Using cosh-cosh radius-grid."
#else
write(*,*) ">>> Using log-exp radius-grid."
#endif
#if GRIDTYPE_PHI==1
write(*,*) ">>> Using linear angular-grid."
#elif GRIDTYPE_PHI==2
write(*,*) ">>> Using cos-like angular-grid."
#else
write(*,*) ">>> Using linear angular-grid."
#endif
end if

end subroutine Initialize_HexGrid

!!!!!!-----------------------------------------------------------------------------------------------------
!!!! --------------------------- OPERATIONS WITH INDICES AND COORDINATES-----------------------------------
!!!!!!-----------------------------------------------------------------------------------------------------

!!!! Radius of the point on the hexagon
pure function PointRadius(x1,x2)
    real(dp),intent(in)::x1,x2
    real(dp)::PointRadius
    PointRadius=max(abs(x1),abs(x2),abs(x1+x2))
end function PointRadius

!!!! Angle of the point on the hexagon
pure function PointAngle(x1,x2)
    real(dp),intent(in)::x1,x2
    real(dp)::PointAngle
    real(dp)::r,x3
    x3=-x1-x2
    r=max(abs(x1),abs(x2),abs(x3))

    if(x1>0 .and. x2>=0 .and. x3<0) then
        PointAngle=x2/r
    else if(x1<=0 .and. x2>0 .and. x3<0) then
        PointAngle=(1-x1/r)
    else if(x1<0 .and. x2>0 .and. x3>=0) then
        PointAngle=(3-x2/r)
    else if(x1<0 .and. x2<=0 .and. x3>0) then
        PointAngle=(3-x2/r)
    else if(x1>=0 .and. x2<0 .and. x3>0) then
        PointAngle=(4+x1/r)
    else if(x1>0 .and. x2<0 .and. x3<=0) then
        PointAngle=(6+x2/r)
    else
        PointAngle=0._dp
    end if
end function PointAngle

!!!! the function that make the grid for radious variable
!!!! at r=1,    n->NUM_R
!!!! at r=xMIN, n->0
pure function RadiusFunction(r)
    real(dp),intent(in)::r
    real(dp)::RadiusFunction
#if GRIDTYPE_R==1
    RadiusFunction=NUM_R*(1-Log(r)*grid_paramC)

#elif GRIDTYPE_R==2
    RadiusFunction=NUM_R*((grid_paramC+1)*r**(1._dp/grid_param0)-grid_paramC)

#elif GRIDTYPE_R==3
    RadiusFunction=NUM_R*(1-Log(r/(2-r))*grid_paramC)

#elif GRIDTYPE_R==4
    RadiusFunction=NUM_R*(1-acosh(r**(-1._dp/grid_param0))*grid_paramC)

#else
    RadiusFunction=NUM_R*(1-Log(r)*grid_paramC)

#endif

end function RadiusFunction

!!!! the function that make the radious from the grid variable
!!!! at n=0,     r->xMIN
!!!! at n=NUM_R, r->1
pure function InvRadiusFunction(n)
    real(dp),intent(in)::n
    real(dp)::InvRadiusFunction
#if GRIDTYPE_R==1
    InvRadiusFunction=exp(real(NUM_R-n,dp)/NUM_R/grid_paramC)

#elif GRIDTYPE_R==2
    InvRadiusFunction=((real(n,dp)/NUM_R+grid_paramC)/(grid_paramC+1))**grid_param0

#elif GRIDTYPE_R==3
    InvRadiusFunction=2._dp*exp(real(NUM_R-n,dp)/NUM_R/grid_paramC)/(1+exp(real(NUM_R-n,dp)/NUM_R/grid_paramC))

#elif GRIDTYPE_R==4
    InvRadiusFunction=cosh(real(NUM_R-n,dp)/NUM_R/grid_paramC)**(-grid_param0)

#else
    InvRadiusFunction=exp(real(NUM_R-n,dp)/NUM_R/grid_paramC)

#endif

end function InvRadiusFunction

!!!! the function that make the grid for angle variable
!!!! at phi=0,    n->0
!!!! at phi=1,    n->NUM_phi
!!!! at phi=2,    n->2NUM_phi
!!!!  ...
!!!! at phi=5,    n->5NUM_phi
!!!! at phi=6,    n->6NUM_phi->0
pure function AngleFunction(phi)
    real(dp),intent(in)::phi
    real(dp)::ANGLEFunction
#if GRIDTYPE_PHI==1
    AngleFunction=NUM_phi*phi

#elif GRIDTYPE_PHI==2
    AngleFunction=int(mod(phi,6._dp))*NUM_phi+NUM_phi/pi*acos(1-2*mod(phi,1._dp))

#else
    AngleFunction=NUM_phi*phi

#endif

end function AngleFunction

!!!! the function that make the angle from the grid variable within one sector(!)
!!!! at k=0,            f->0
!!!! at k=NUM_phi,      f->1
pure function InvAngleFunction(k)
    real(dp),intent(in)::k
    real(dp)::k_pos
    real(dp)::InvAngleFunction

    !!! first project to [0,NUM_perimeter]
    !!! fortran mod does not account circles. Thus for negative numbers we need to acount it by hands
    if(k>=0) then
        k_pos=mod(k,real(NUM_perimeter))
    else
        k_pos=NUM_perimeter+mod(k,real(NUM_perimeter))
    end if

#if GRIDTYPE_PHI==1
    InvAngleFunction=mod(k,real(NUM_phi,dp))/NUM_phi

#elif GRIDTYPE_PHI==2
    InvAngleFunction=0.5_dp*(1-cos(pi*mod(k,real(NUM_phi,dp))/NUM_phi))

#else
    InvAngleFunction=mod(k,real(NUM_phi,dp))/NUM_phi
#endif


end function InvAngleFunction


!!! the function which takes (x1,x2), and return the value of (n=radious grid, k=perimeter grid)
pure subroutine X12toNK(x1,x2,n,k)
    real(dp),intent(in)::x1,x2
    real(dp),intent(out)::n,k

    n=RadiusFunction(PointRadius(x1,x2))
    k=AngleFunction(PointAngle(x1,x2))

end subroutine X12toNK

!!! the function which takes (n,k), and return the value of (x1,x2)
subroutine NKtoX12_real(n,k,x1,x2)
    real(dp),intent(out)::x1,x2
    real(dp),intent(in)::n,k
    real(dp)::r,phi,kk
    r=InvRadiusFunction(n)

    if(k>=0) then
        kk=mod(k,real(NUM_perimeter+1,dp))
    else
        kk=mod(NUM_perimeter+1-Abs(k),real(NUM_perimeter+1,dp))
    end if

    !!! this 0<phi<1
    phi=InvAngleFunction(kk)
    !!! the sector is selected here:
    if(kk<NUM_phi) then
        x1=r*(1-phi)
        x2=r*phi
    else if(kk<2*NUM_phi) then
        x1=-r*phi
        x2=r
    else if(kk<3*NUM_phi) then
        x1=-r
        x2=r*(1-phi)
    else if(kk<4*NUM_phi) then
        x1=-r*(1-phi)
        x2=-r*phi
    else if(kk<5*NUM_phi) then
        x1=r*phi
        x2=-r
    else
        x1=r
        x2=-r*(1-phi)
    end if

end subroutine NKtoX12_real
!!!!! same as above but for integer parameter
subroutine NKtoX12_int(n,k,x1,x2)
    real(dp),intent(out)::x1,x2
    integer,intent(in)::n,k

    call NKtoX12_real(real(n,dp),real(k,dp),x1,x2)

end subroutine NKtoX12_int


!!! Maps (n,k) to 1D array
pure function Index1D(n,k)
    integer::Index1D
    integer,intent(in)::n,k

    Index1D=(NUM_perimeter+1)*n+mod(k,6*NUM_phi)

end function Index1D

!!! Maps 1D variable to (n,k) coordinates
pure function Index2D(c)
    integer,dimension(1:2)::index2D
    integer,intent(in)::c

    Index2D(2)=mod(c,NUM_perimeter+1)
    Index2D(1)=(c-Index2D(2))/(NUM_perimeter+1)
end function Index2D

!!!!! Translation of 1D index directly to the (x1,x2)
subroutine NtoX12(n,x1,x2)
    real(dp),intent(out)::x1,x2
    integer,intent(in)::n
    integer,dimension(1:2)::ind

    ind=Index2D(n)

    call NKtoX12_real(real(ind(1),dp),real(ind(2),dp),x1,x2)

end subroutine NtoX12

!!!!!!-----------------------------------------------------------------------------------------------------
!!!!------------------------------------------------- INTERPOLATION FUNCTIONS------------------------------
!!!!!!-----------------------------------------------------------------------------------------------------

!!!!! the interpolating function
!!! i,j is the number of node (int)
!!! x, y is the requested point (in the position space)
pure function Windex_plus(n,k,x1,x2)
real(dp)::Windex_plus
integer,intent(in)::n,k
real(dp),intent(in)::x1,x2
real(dp)::dr,dphi,dd,nn,kk
integer::k_red

if(k>=0) then
    k_red=mod(k,NUM_perimeter+1)
else
    k_red=mod(NUM_perimeter+1-Abs(k),NUM_perimeter+1)
end if


call X12toNK(x1,x2,nn,kk)
if(k==0 .and. kk>NUM_perimeter) kk=kk-NUM_perimeter-1
dr=nn-n
dphi=kk-k_red
dd=dphi+dr

Windex_plus=max(1._dp-max(abs(dr),abs(dphi),abs(dd)),0._dp)

end function Windex_plus

!!!!! the interpolating function
!!! i,j is the number of node (int)
!!! x, y is the requested point (in the position space)
pure function Windex_minus(n,k,x1,x2)
real(dp)::Windex_minus
integer,intent(in)::n,k
real(dp),intent(in)::x1,x2
real(dp)::dr,dphi,dd,nn,kk
integer::k_red

if(k>=0) then
    k_red=mod(k,NUM_perimeter+1)
else
    k_red=mod(NUM_perimeter+1-Abs(k),NUM_perimeter+1)
end if

call X12toNK(x1,x2,nn,kk)
if(k==0 .and. kk>NUM_perimeter) kk=kk-NUM_perimeter-1

dr=nn-n
dphi=kk-k_red
dd=dphi-dr

Windex_minus=max(1._dp-max(abs(dr),abs(dphi),abs(dd)),0._dp)

end function Windex_minus

!!!!! the interpolating function, which is (W_+ + W_-)/2
!!! i,j is the number of node (int)
!!! x, y is the requested point (in the position space)
pure function Windex_av(n,k,x1,x2)
real(dp)::Windex_av
integer,intent(in)::n,k
real(dp),intent(in)::x1,x2
real(dp)::dr,dphi,d1,d2,nn,kk,w1,w2
integer::k_red

if(k>=0) then
    k_red=mod(k,NUM_perimeter+1)
else
    k_red=mod(NUM_perimeter+1-Abs(k),NUM_perimeter+1)
end if


call X12toNK(x1,x2,nn,kk)
if(k==0 .and. kk>NUM_perimeter) kk=kk-NUM_perimeter-1
dr=nn-n
dphi=kk-k_red
d1=dphi+dr
d2=dphi-dr

w1=max(1._dp-max(abs(dr),abs(dphi),abs(d1)),0._dp)
w2=max(1._dp-max(abs(dr),abs(dphi),abs(d2)),0._dp)

Windex_av=(w1+w2)/2

end function Windex_av

!!! gridTYPE is the T or F
!!! if gridTYPE=T, the grid is with \-section i.e. (+) grid
!!! if gridTYPE=F, the grid is with /-section i.e. (-) grid
pure function Windex(n,k,x1,x2)
real(dp)::Windex
integer,intent(in)::n,k
real(dp),intent(in)::x1,x2
logical::t

Windex=Windex_av(n,k,x1,x2)

end function Windex

!!! Project Function over 2D grid to 1D grid
pure function F2Dto1D(F)
    real(dp),dimension(0:NUM_R,0:NUM_perimeter),intent(in)::F
    real(dp),dimension(0:NUM_TOT-1)::F2Dto1D
    integer::n,k,c1

    do n=0,NUM_R
    do k=0,NUM_perimeter

    c1=Index1D(n,k)
    F2Dto1D(c1)=F(n,k)
    end do
    end do

end function F2Dto1D

!!! Project Function over X,Y to 1D grid
function FXYto1D(F)
    real(dp)::F
    real(dp),dimension(0:NUM_TOT-1)::FXYto1D
    integer::n,k,c1
    real(dp)::x1,x2

    do n=0,NUM_R
    do k=0,NUM_perimeter
    c1=Index1D(n,k)
    call NKtoX12_int(n,k,x1,x2)
    FXYto1D(c1)=F(x1,x2)
    end do
    end do

end function FXYto1D

!!! Project Function over 1D grid to 2D grid
pure function F1Dto2D(F)
    real(dp),dimension(0:NUM_R,0:NUM_perimeter)::F1Dto2D
    real(dp),dimension(0:NUM_TOT-1),intent(in)::F
    integer::i(1:2),c1

    F1Dto2D=0._dp
    do c1=0,NUM_TOT-1
        i=Index2D(c1)
        F1Dto2D(i(1),i(2))=F(c1)
    end do

end function F1Dto2D

!!! actually intepolation procedure
function FatXfrom2D(F,x1,x2)
real(dp)::FatXfrom2D
real(dp),intent(in),dimension(0:NUM_R,0:NUM_perimeter)::F
real(dp),intent(in)::x1,x2
real(dp)::n_r,k_r
integer::n,k

if(abs(x1)>=1 .or. abs(x2)>=1 .or. abs(x1+x2)>=1) then
        FatXfrom2D=0._dp
else

    call X12toNK(x1,x2,n_r,k_r)
    n=int(n_r)
    k=int(k_r)

    FatXfrom2D=F(n,k)*Windex(n,k,x1,x2)&
            +F(n+1,k)*Windex(n+1,k,x1,x2)&
            +F(n,k+1)*Windex(n,k+1,x1,x2)&
            +F(n+1,k+1)*Windex(n+1,k+1,x1,x2)
end if

end function FatXfrom2D

!!!!!!-----------------------------------------------------------------------------------------------------
!!!!!!-----------------------------------------------------OPERATIONS WITH BOUNDARIES ---------------------
!!!!!!-----------------------------------------------------------------------------------------------------

!!!! This function returns the min/max values of x1,x2,x3 for the hexagon specified by n,k
!!!! It is needed to estimate the integration limits for kernels
!!!! this function is very straightforward, it is made intentionally, to be sure that the boundary is correct.
!!!! the output is (x1Min,x1Max,x2Min,x2Max,x3Min,x3Max)
function BoundaryValues(n,k)
real(dp),dimension(1:6)::BoundaryValues
integer,intent(in)::n,k
logical::t
real(dp)::x11,x12,x13,x14,x15,x16,x21,x22,x23,x24,x25,x26,x31,x32,x33,x34,x35,x36
real(dp),parameter::delta=1.3_dp


call NKtoX12_real(real(n,dp),real(k+delta,dp),x12,x22)
call NKtoX12_real(real(n,dp),real(k-delta,dp),x15,x25)

if(n>0) then
    call NKtoX12_real(real(n-1,dp),real(k+delta,dp),x13,x23)
    call NKtoX12_real(real(n-1,dp),real(k-delta,dp),x14,x24)
else
    !!!!! I must put a little margin in this case, because the line passes exactly on the border of hexagon
    !!!!! and sometimes it does not cross it by 10^{-...} and produces incorrect result
    call NKtoX12_real(real(n-0.001,dp),real(k+delta,dp),x13,x23)
    call NKtoX12_real(real(n-0.001,dp),real(k-delta,dp),x14,x24)
    !x13=x12
    !x23=x22
    !x14=x12
    !x24=x22
end if

if(n<NUM_R) then
    call NKtoX12_real(real(n+1,dp),real(k-delta,dp),x16,x26)
    call NKtoX12_real(real(n+1,dp),real(k+delta,dp),x11,x21)
else
    !!!! same here, line passes just by the border of hexagon
    call NKtoX12_real(real(n+0.001,dp),real(k-delta,dp),x16,x26)
    call NKtoX12_real(real(n+0.001,dp),real(k+delta,dp),x11,x21)
    !x16=x12
    !x26=x22
    !x11=x12
    !x21=x22
end if

!!! also the x3's fro each corner
x31=-x11-x21
x32=-x12-x22
x33=-x13-x23
x34=-x14-x24
x35=-x15-x25
x36=-x16-x26

BoundaryValues(1)=min(x11,x12,x13,x14,x15,x16)
BoundaryValues(2)=max(x11,x12,x13,x14,x15,x16)
BoundaryValues(3)=min(x21,x22,x23,x24,x25,x26)
BoundaryValues(4)=max(x21,x22,x23,x24,x25,x26)
BoundaryValues(5)=min(x31,x32,x33,x34,x35,x36)
BoundaryValues(6)=max(x31,x32,x33,x34,x35,x36)

end function BoundaryValues

!!!!!! Subroutine checks if the line (x1-v,x2+v) intersects the approximate hexagon n,k
!!!!!! if negativeV, then (x1+v,x2-v)
!!!!!! it returns intersect=T, if it does; set vMin and vMax
!!!!!! if there is no intersection then intersect=F, and vMin vMax are undefined.
subroutine LimitsX3(n,k,x1,x2,intersect,vMin,vMax,negativeV)
integer,intent(in)::n,k
real(dp),intent(in)::x1,x2
real(dp),intent(out)::vMin,vMax
logical,intent(out)::intersect
logical,optional,intent(in)::negativeV
logical::nV

real(dp),dimension(1:6)::limits
real(dp)::x3

limits=BoundaryValues(n,k)

if(present(negativeV)) then
    nV=negativeV
else
    nV=.false.
end if

x3=-x1-x2

!!! the intersection is possible if the line is in the strip.
if(limits(5)<=x3 .and. x3<=limits(6)) then
    intersect=.true.
else
    intersect=.false.
    return
end if

if(nV) then
    vMin=max(x2-limits(4),limits(1)-x1)
    vMax=min(x2-limits(3),limits(2)-x1)
else
    vMin=max(x1-limits(2),limits(3)-x2)
    vMax=min(x1-limits(1),limits(4)-x2)
end if

if(vMax<vMin) then
    !write(*,*) "ERROR 1: in LimitsX3"
    intersect=.false.
end if
end subroutine LimitsX3

!!!!!! Subroutine checks if the line (x2-v,x3+v) intersects the approximate hexagon n,k
!!!!!! it returns intersect=T, if it does; set vMin and vMax
!!!!!! if there is no intersection then intersect=F, and vMin vMax are undefined.
subroutine LimitsX1(n,k,x1,x2,intersect,vMin,vMax,negativeV)
integer,intent(in)::n,k
real(dp),intent(in)::x1,x2
real(dp),intent(out)::vMin,vMax
logical,intent(out)::intersect
logical,optional,intent(in)::negativeV
logical::nV

real(dp),dimension(1:6)::limits
real(dp)::x3


limits=BoundaryValues(n,k)

if(present(negativeV)) then
    nV=negativeV
else
    nV=.false.
end if

x3=-x1-x2

!!! the intersection is possible if the line is in the strip.
if(limits(1)<=x1 .and. x1<=limits(2)) then
    intersect=.true.
else
    intersect=.false.
    return
end if

if(nV) then
    vMin=max(x3-limits(6),limits(3)-x2)
    vMax=min(x3-limits(5),limits(4)-x2)
else
    vMin=max(x2-limits(4),limits(5)-x3)
    vMax=min(x2-limits(3),limits(6)-x3)
end if
if(vMax<vMin) then
    !write(*,*) "ERROR 1: in LimitsX1"
    intersect=.false.
end if
end subroutine LimitsX1


!!!!!! Subroutine checks if the line (x1-v,x3+v) intersects the approximate hexagon n,k
!!!!!! it returns intersect=T, if it does; set vMin and vMax
!!!!!! if there is no intersection then intersect=F, and vMin vMax are undefined.
subroutine LimitsX2(n,k,x1,x2,intersect,vMin,vMax,negativeV)
integer,intent(in)::n,k
real(dp),intent(in)::x1,x2
real(dp),intent(out)::vMin,vMax
logical,intent(out)::intersect
logical,optional,intent(in)::negativeV
logical::nV

real(dp),dimension(1:6)::limits
real(dp)::x3


limits=BoundaryValues(n,k)

if(present(negativeV)) then
    nV=negativeV
else
    nV=.false.
end if

x3=-x1-x2

!!! the intersection is possible if the line is in the strip.
if(limits(3)<=x2 .and. x2<=limits(4)) then
    intersect=.true.
else
    intersect=.false.
    return
end if

if(nV) then
    vMin=max(x3-limits(6),limits(1)-x1)
    vMax=min(x3-limits(5),limits(2)-x1)
else
    vMin=max(x1-limits(2),limits(5)-x3)
    vMax=min(x1-limits(1),limits(6)-x3)
end if
if(vMax<vMin) then
    !write(*,*) "ERROR 1: in LimitsX2", n,k,x1,x2
    intersect=.false.
end if
end subroutine LimitsX2

end module HexGrid
