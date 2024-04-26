!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! This module is the central module for SnowFlake library.                                 !!
!! It decompose the input and call correponsing evolution routines                          !!
!!                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module SnowFlake
use IO_snowflake
use HexGrid
use EvolutionKernels
implicit none

private

INCLUDE 'commonVariables.f90'

!!!!! These global variables save the result of current evolution
!!!!! These are C-odd functions
real(dp),allocatable:: Uplus(:,:),Dplus(:,:),Splus(:,:),Cplus(:,:),Bplus(:,:),Gplus(:,:)
real(dp),allocatable:: Uminus(:,:),Dminus(:,:),Sminus(:,:),Cminus(:,:),Bminus(:,:),Gminus(:,:)
!!!!! These are C-even functions
real(dp),allocatable:: Uodd(:,:),Dodd(:,:),Sodd(:,:),Codd(:,:),Bodd(:,:)

public:: SnowFlake_Initialize,ComputeEvolution, ComputeEvolutionChiralOdd, GetPDF,GetPDFChiralOdd

contains


subroutine SnowFlake_Initialize(file,prefix)
character(len=*)::file
character(len=*),optional::prefix
character(len=300)::path
!$ real*8::omp_get_wtime
real*8::t1,t2

if(present(prefix)) then
    path=trim(adjustl(prefix))//trim(adjustr(file))
else
    path=trim(adjustr(file))
end if

call cpu_time(t1)
!$ t1=omp_get_wtime()

write(*,*) color("----------------------- SNOWFLAKE V1.0 -----------------------",c_cyan)
write(*,*) color("                           /\‾‾‾/\                            ",c_cyan)
write(*,*) color("                          /  \ /  \                           ",c_cyan)
write(*,*) color("                          ----o----                           ",c_cyan)
write(*,*) color("                          \  / \  /                           ",c_cyan)
write(*,*) color("                           \/___\/                            ",c_cyan)
write(*,*) color("- please, cite [S.Rodini, L.Rossi, A.Vladimirov, [2404.????] -",c_cyan)

write(*,*) "SnowFlake initlization with INI-file:"//trim(path)

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

!!! Value of zero
call MoveTO(51,'4   :')
read(51,*) zero

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

CLOSE (51, STATUS='KEEP')

if(showINI) then
    write(*,*) "SnowFlake parameters:"
    write(*,'("Grid size = (",I4," x",I4,") = ",I8, " nodes.")') NUM_perimeter,NUM_R,NUM_TOT
    write(*,'(A)',advance="no") "Include Chiral-odd evolution:"
    if(IncludeChiralEvenEvolution) then
        write(*,*) color(" YES",c_green)
    else
        write(*,*) color(" NO",c_red)
    end if
    write(*,'(A)',advance="no") "Include Chiral-odd evolution:"
    if(IncludeChiralOddEvolution) then
        write(*,*) color(" YES",c_green)
    else
        write(*,*) color(" NO",c_red)
    end if
    write(*,'(A)',advance="no") "Include gluon/flavor mixing :"
    if(useSingletEvolution) then
        write(*,*) color(" YES",c_green)
    else
        write(*,*) color(" NO",c_red)
    end if
    write(*,'(A,F8.3)') "Mass of CHARM threshold [GeV] :",massCHARM
    write(*,'(A,F8.3)') "Mass of BOTTOM threshold [GeV]:",massBOTTOM
end if

if(.not.(IncludeChiralEvenEvolution .or. IncludeChiralOddEvolution)) then
    write(*,*) ErrorString("Non of evolution types are included. EVALUTION TERMINATED.","main")
    stop
end if

allocate(Uplus(0:NUM_R,0:NUM_perimeter))
allocate(Dplus(0:NUM_R,0:NUM_perimeter))
allocate(Splus(0:NUM_R,0:NUM_perimeter))
allocate(Cplus(0:NUM_R,0:NUM_perimeter))
allocate(Bplus(0:NUM_R,0:NUM_perimeter))
allocate(Gplus(0:NUM_R,0:NUM_perimeter))

allocate(Uminus(0:NUM_R,0:NUM_perimeter))
allocate(Dminus(0:NUM_R,0:NUM_perimeter))
allocate(Sminus(0:NUM_R,0:NUM_perimeter))
allocate(Cminus(0:NUM_R,0:NUM_perimeter))
allocate(Bminus(0:NUM_R,0:NUM_perimeter))
allocate(Gminus(0:NUM_R,0:NUM_perimeter))

allocate(Uodd(0:NUM_R,0:NUM_perimeter))
allocate(Dodd(0:NUM_R,0:NUM_perimeter))
allocate(Sodd(0:NUM_R,0:NUM_perimeter))
allocate(Codd(0:NUM_R,0:NUM_perimeter))
allocate(Bodd(0:NUM_R,0:NUM_perimeter))

call Initialize_HexGrid(path)
call EvolutionKernels_Initialize(path)

call cpu_time(t2)
!$ t2=omp_get_wtime()
write(*,*) "--------- SnowFlake initialization complete. Initialization time =",t2-t1

end subroutine SnowFlake_Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! INTERFACES GLOBAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compute the evolution of chiral-even distributions
!!! mu0, mu1 initial and final scales of evolution (GEV)
!!! alpha = alpha_s(mu)
!!! inputQ,inputG defines the type of input functions
!!! --------
!!! inputQ = 'C'  C-definite input for quarks, i.e. Q1=\mathcal{S}+, Q2=\mathcal{S}-
!!! inputQ = 'S'  S-function input for quarks, i.e. Q1=S+, Q2=S-
!!! inputQ = 'T'  T-function input for quarks, i.e. Q1=T, Q2=DeltaT
!!! --------
!!! inputQ = 'C'  C-definite input for gluons, i.e. G1=\mathcal{G}+, G2=\mathcal{G}-
!!! inputQ = 'T'  T-function input for gluons, i.e. G1=T, G2=DeltaT
!!! --------
!!! the result of evolution is stored in the global variable and can be called separately.
subroutine ComputeEvolution(mu0,mu1,alpha,G1,U1,D1,S1,C1,B1,G2,U2,D2,S2,C2,B2,inputQ,inputG)
real(dp),external,optional::G1,U1,D1,S1,C1,B1,G2,U2,D2,S2,C2,B2
real(dp),external::alpha
real(dp),dimension(0:NUM_TOT-1)::G_p,U_p,D_p,S_p,C_p,B_p
real(dp),dimension(0:NUM_TOT-1)::G_m,U_m,D_m,S_m,C_m,B_m
real(dp),intent(in)::mu0,mu1
character(len=1),optional,intent(in)::inputQ,inputG
character(len=1)::inQ,inG

integer::n,k,cc
real(dp)::x1,x2,x3

if(present(inputQ)) then
    inQ=inputQ
else
    inQ='C'
end if
if(present(inputG)) then
    inG=inputG
else
    inG='C'
end if

!!!! if some function is absent it is interpreted as zero


!!! deshifrating input
SELECT CASE (inQ)
    CASE ('T')
        !!!! \mathcal{S}+=T(123)+T(321)-DeltaT(123)+DeltaT(321)
        !!!! \mathcal{S}-=T(123)-T(321)-DeltaT(123)-DeltaT(321)
        do n=0,NUM_R
        do k=0,NUM_perimeter
        cc=Index1D(n,k)
        call NKtoX12(n,k,x1,x2)
        x3=-x1-x2
            !!!! U-quark
            if(present(U1)) then !!!! T
                U_p(cc)=U1(x1,x2)+U1(x3,x2)
                U_m(cc)=U1(x1,x2)-U1(x3,x2)
            else
                U_p(cc)=0._dp
                U_m(cc)=0._dp
            end if
            if(present(U2)) then !!!! DeltaT
                U_p(cc)=U_p(cc)-U2(x1,x2)+U2(x3,x2)
                U_m(cc)=U_m(cc)-U2(x1,x2)-U2(x3,x2)
            end if

            !!!! D-quark
            if(present(D1)) then !!!! T
                D_p(cc)=D1(x1,x2)+D1(x3,x2)
                D_m(cc)=D1(x1,x2)-D1(x3,x2)
            else
                D_p(cc)=0._dp
                D_m(cc)=0._dp
            end if
            if(present(D2)) then !!!! DeltaT
                D_p(cc)=D_p(cc)-D2(x1,x2)+D2(x3,x2)
                D_m(cc)=D_m(cc)-D2(x1,x2)-D2(x3,x2)
            end if

            !!!! S-quark
            if(present(S1)) then !!!! T
                S_p(cc)=S1(x1,x2)+S1(x3,x2)
                S_m(cc)=S1(x1,x2)-S1(x3,x2)
            else
                S_p(cc)=0._dp
                S_m(cc)=0._dp
            end if
            if(present(S2)) then !!!! DeltaT
                S_p(cc)=S_p(cc)-S2(x1,x2)+S2(x3,x2)
                S_m(cc)=S_m(cc)-S2(x1,x2)-S2(x3,x2)
            end if

            !!!! C-quark
            if(present(C1)) then !!!! T
                C_p(cc)=C1(x1,x2)+C1(x3,x2)
                C_m(cc)=C1(x1,x2)-C1(x3,x2)
            else
                C_p(cc)=0._dp
                C_m(cc)=0._dp
            end if
            if(present(C2)) then !!!! DeltaT
                C_p(cc)=C_p(cc)-C2(x1,x2)+C2(x3,x2)
                C_m(cc)=C_m(cc)-C2(x1,x2)-C2(x3,x2)
            end if

            !!!! B-quark
            if(present(B1)) then !!!! T
                B_p(cc)=B1(x1,x2)+B1(x3,x2)
                B_m(cc)=B1(x1,x2)-B1(x3,x2)
            else
                B_p(cc)=0._dp
                B_m(cc)=0._dp
            end if
            if(present(B2)) then !!!! DeltaT
                B_p(cc)=B_p(cc)-B2(x1,x2)+B2(x3,x2)
                B_m(cc)=B_m(cc)-B2(x1,x2)-B2(x3,x2)
            end if

        end do
        end do

    CASE ('S')
        !!!! \mathcal{S}+=-2(S+(123)+S-(321))
        !!!! \mathcal{S}-=-2(S+(123)-S-(321))
        do n=0,NUM_R
        do k=0,NUM_perimeter
        cc=Index1D(n,k)
        call NKtoX12(n,k,x1,x2)
        x3=-x1-x2
            !!!! U-quark
            if(present(U1)) then !!!! S+
                U_p(cc)=-2*U1(x1,x2)
                U_m(cc)=-2*U1(x1,x2)
            else
                U_p(cc)=0._dp
                U_m(cc)=0._dp
            end if
            if(present(U2)) then !!!! S-
                U_p(cc)=U_p(cc)-2*U2(x3,x2)
                U_m(cc)=U_m(cc)+2*U2(x3,x2)
            end if

            !!!! D-quark
            if(present(D1)) then !!!! S+
                D_p(cc)=-2*D1(x1,x2)
                D_m(cc)=-2*D1(x1,x2)
            else
                D_p(cc)=0._dp
                D_m(cc)=0._dp
            end if
            if(present(D2)) then !!!! S-
                D_p(cc)=D_p(cc)-2*D2(x3,x2)
                D_m(cc)=D_m(cc)+2*D2(x3,x2)
            end if

            !!!! S-quark
            if(present(S1)) then !!!! S+
                S_p(cc)=-2*S1(x1,x2)
                S_m(cc)=-2*S1(x1,x2)
            else
                S_p(cc)=0._dp
                S_m(cc)=0._dp
            end if
            if(present(S2)) then !!!! S-
                S_p(cc)=S_p(cc)-2*S2(x3,x2)
                S_m(cc)=S_m(cc)+2*S2(x3,x2)
            end if

            !!!! C-quark
            if(present(C1)) then !!!! S+
                C_p(cc)=-2*C1(x1,x2)
                C_m(cc)=-2*C1(x1,x2)
            else
                C_p(cc)=0._dp
                C_m(cc)=0._dp
            end if
            if(present(C2)) then !!!! S-
                C_p(cc)=C_p(cc)-2*C2(x3,x2)
                C_m(cc)=C_m(cc)+2*C2(x3,x2)
            end if

            !!!! B-quark
            if(present(B1)) then !!!! S+
                B_p(cc)=-2*B1(x1,x2)
                B_m(cc)=-2*B1(x1,x2)
            else
                B_p(cc)=0._dp
                B_m(cc)=0._dp
            end if
            if(present(B2)) then !!!! S-
                B_p(cc)=B_p(cc)-2*B2(x3,x2)
                B_m(cc)=B_m(cc)+2*B2(x3,x2)
            end if

        end do
        end do

    CASE ('C')  !!! the definition is C-definite
        if(present(U1)) then
            U_p=FXYto1D(U1)
        else
            U_p=0._dp
        end if
        if(present(D1)) then
            D_p=FXYto1D(D1)
        else
            D_p=0._dp
        end if
        if(present(S1)) then
            S_p=FXYto1D(S1)
        else
            S_p=0._dp
        end if
        if(present(C1)) then
            C_p=FXYto1D(C1)
        else
            C_p=0._dp
        end if
        if(present(B1)) then
            B_p=FXYto1D(B1)
        else
            B_p=0._dp
        end if
        !--------------
        if(present(U2)) then
            U_m=FXYto1D(U2)
        else
            U_m=0._dp
        end if
        if(present(D2)) then
            D_m=FXYto1D(D2)
        else
            D_m=0._dp
        end if
        if(present(S2)) then
            S_m=FXYto1D(S2)
        else
            S_m=0._dp
        end if
        if(present(C2)) then
            C_m=FXYto1D(C2)
        else
            C_m=0._dp
        end if
        if(present(B2)) then
            B_m=FXYto1D(B2)
        else
            B_m=0._dp
        end if
    CASE DEFAULT
        write(*,*) ErrorString("unknown specification for inputQ. Evaluation STOP"," ")
        stop
END SELECT

!!!! if below thresholds nullify baoundary values of C and B
if(mu0<massCHARM) then
    C_p=0._dp
    C_m=0._dp
end if
if(mu0<massBOTTOM) then
    B_p=0._dp
    B_m=0._dp
end if

!!! deshifrating gluon input
SELECT CASE(inG)
    CASE ('T')
        !!! \mathcal{G}+=T+(123)-T+(132)+T+(213)
        !!! \mathcal{G}-=T-(123)+T-(132)-T-(213)
        do n=0,NUM_R
        do k=0,NUM_perimeter
        cc=Index1D(n,k)
        call NKtoX12(n,k,x1,x2)
        x3=-x1-x2
            !!!! G+
            if(present(G1)) then !!!! T
                G_p(cc)=G1(x1,x2)-G1(x1,x3)+G1(x2,x1)
            else
                G_p(cc)=0._dp
            end if
            !!!! G-
            if(present(G2)) then !!!! T
                G_m(cc)=G2(x1,x2)+G2(x1,x3)-G2(x2,x1)
            else
                G_m(cc)=0._dp
            end if
        end do
        end do
    CASE ('C') !!! the definition is C-definite
        if(present(G1)) then
            G_p=FXYto1D(G1)
        else
            G_p=0._dp
        end if
        if(present(G2)) then
            G_m=FXYto1D(G2)
        else
            G_m=0._dp
        end if
    CASE DEFAULT
        write(*,*) ErrorString("unknown specification for inputG. Evaluation STOP"," ")
        stop
END SELECT


!!! finally call evolution
call EvolvePLUS(alpha,mu0,mu1,G_p,U_p,D_p,S_p,C_p,B_p)
call EvolveMINUS(alpha,mu0,mu1,G_m,U_m,D_m,S_m,C_m,B_m)

end subroutine ComputeEvolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compute the evolution of chiral-odd distributions
!!! mu0, mu1 initial and final scales of evolution (GEV)
!!! alpha = alpha_s(mu)
!!! input are functions H, or E
!!! --------------------------
!!! the result of evolution is stored in the global variable and can be called separately.
subroutine ComputeEvolutionChiralOdd(mu0,mu1,alpha,U1,D1,S1,C1,B1)
real(dp),external,optional::U1,D1,S1,C1,B1
real(dp),external::alpha
real(dp),dimension(0:NUM_TOT-1)::U_p,D_p,S_p,C_p,B_p
real(dp),intent(in)::mu0,mu1

integer::n,k,cc
real(dp)::x1,x2,x3,t0,t1

    ! no evolution option
    if(abs(mu0-mu1)<zero) then
        if(present(U1)) then
            U_p=FXYto1D(U1)
            Uodd=F1Dto2D(U_p)
        else
            Uodd=0._dp
        end if

        if(present(D1)) then
            D_p=FXYto1D(D1)
            Dodd=F1Dto2D(D_p)
        else
            Dodd=0._dp
        end if

        if(present(S1)) then
            S_p=FXYto1D(S1)
            Sodd=F1Dto2D(S_p)
        else
            Sodd=0._dp
        end if

        if(present(C1)) then
            C_p=FXYto1D(C1)
            Codd=F1Dto2D(C_p)
        else
            Codd=0._dp
        end if

        if(present(B1)) then
            B_p=FXYto1D(B1)
            Bodd=F1Dto2D(B_p)
        else
            Bodd=0._dp
        end if
    return
    end if

    !!!! t=log[mu^2]
    t0=2*log(mu0)
    t1=2*log(mu1)

    !!! the chiral odd functions do not mix in flavor.
    !!! Thus if something is zero, it remains zero. And there is no need to compute it

    if(present(U1)) then
        U_p=FXYto1D(U1)
        call EvChiralOdd(U_p,as,t0,t1)
        Uodd=F1Dto2D(U_p)
    else

        Uodd=0._dp
    end if

    if(present(D1)) then
        D_p=FXYto1D(D1)
        call EvChiralOdd(D_p,as,t0,t1)
        Dodd=F1Dto2D(D_p)
    else
        Dodd=0._dp
    end if

    if(present(S1)) then
        S_p=FXYto1D(S1)
        call EvChiralOdd(S_p,as,t0,t1)
        Sodd=F1Dto2D(S_p)
    else
        Sodd=0._dp
    end if

    if(present(C1)) then
        C_p=FXYto1D(C1)
        call EvChiralOdd(C_p,as,t0,t1)
        Codd=F1Dto2D(C_p)
    else
        Codd=0._dp
    end if

    if(present(B1)) then
        B_p=FXYto1D(B1)
        call EvChiralOdd(B_p,as,t0,t1)
        Bodd=F1Dto2D(B_p)
    else
        Bodd=0._dp
    end if

contains

    !!!!! as=alpha_s(mu)/(4pi)
    function as(t)
        real(dp)::as,t
        !!! 1/(4 pi)
        real(dp),parameter::pi4minus1=real(0.0795774715459476678844418816862571810,dp)
        as=alpha(Exp(t/2))*pi4minus1
    end function as

end subroutine ComputeEvolutionChiralOdd

!!!!! Computes interpolation flavor f
!!!!! f=0=10=gluon
!!!!! f=1=d
!!!!! f=2=u
!!!!! f=3=s
!!!!! f=4=c
!!!!! f=5=b
!!!!!
function GetPDF(x1,x2,f,outputT)
    real(dp)::GetPDF
    real(dp),intent(in)::x1,x2
    integer,intent(in)::f
    character(len=1),optional,intent(in)::outputT
    character(len=1)::outT
    real(dp)::x3,r
    x3=-x1-x2

    r=max(abs(x1),abs(x2),abs(x1+x2))
    if(r<xMIN) then
        write(*,*) ErrorString("requested value with |x|<xMIN"," ")
        write(*,*) "Requested (x1,x2,x3)=(",x1,x2,x3,") with |x|=",r
        write(*,*) "Current xMIN =",xMIN
        write(*,*) "Increase lower limit in INI-file. Evaluation STOP"
        stop
    end if
    if(abs(x1)>1.or.abs(x2)>1.or.abs(x1+x2)>1) then
        write(*,*) WarningString("requested values at |x|>1. Zero returned."," ")
    GetPDF=0._dp
    end if


    if(present(outputT)) then
    outT=outputT
    else
    outT='C'
    end if

    SELECT CASE (outT)
    CASE ('C')
        if(f==10.or.f==0) then
            GetPDF=FatXfrom2D(Gplus,x1,x2)
        else if(f==1) then
            GetPDF=FatXfrom2D(Dplus,x1,x2)
        else if(f==2) then
            GetPDF=FatXfrom2D(Uplus,x1,x2)
        else if(f==3) then
            GetPDF=FatXfrom2D(Splus,x1,x2)
        else if(f==4) then
            GetPDF=FatXfrom2D(Cplus,x1,x2)
        else if(f==5) then
            GetPDF=FatXfrom2D(Bplus,x1,x2)
        else if(f==-10) then
            GetPDF=FatXfrom2D(Gminus,x1,x2)
        else if(f==-1) then
            GetPDF=FatXfrom2D(Dminus,x1,x2)
        else if(f==-2) then
            GetPDF=FatXfrom2D(Uminus,x1,x2)
        else if(f==-3) then
            GetPDF=FatXfrom2D(Sminus,x1,x2)
        else if(f==-4) then
            GetPDF=FatXfrom2D(Cminus,x1,x2)
        else if(f==-5) then
            GetPDF=FatXfrom2D(Bminus,x1,x2)
        else
            write(*,*) WarningString("the flavor can be only -5,...,5, 10, -10. Zero returned."," ")
            GetPDF=0._dp
        end if
    CASE ('T')
        if(f==10.or.f==0) then
            GetPDF=(FatXfrom2D(Gplus,x1,x2)-FatXfrom2D(Gplus,x3,x2))/2
        else if(f==1) then
            GetPDF=(FatXfrom2D(Dplus,x1,x2)+FatXfrom2D(Dplus,x3,x2)+FatXfrom2D(Dminus,x1,x2)-FatXfrom2D(Dminus,x3,x2))/4
        else if(f==2) then
            GetPDF=(FatXfrom2D(Uplus,x1,x2)+FatXfrom2D(Uplus,x3,x2)+FatXfrom2D(Uminus,x1,x2)-FatXfrom2D(Uminus,x3,x2))/4
        else if(f==3) then
            GetPDF=(FatXfrom2D(Splus,x1,x2)+FatXfrom2D(Splus,x3,x2)+FatXfrom2D(Sminus,x1,x2)-FatXfrom2D(Sminus,x3,x2))/4
        else if(f==4) then
            GetPDF=(FatXfrom2D(Cplus,x1,x2)+FatXfrom2D(Cplus,x3,x2)+FatXfrom2D(Cminus,x1,x2)-FatXfrom2D(Cminus,x3,x2))/4
        else if(f==5) then
            GetPDF=(FatXfrom2D(Bplus,x1,x2)+FatXfrom2D(Bplus,x3,x2)+FatXfrom2D(Bminus,x1,x2)-FatXfrom2D(Bminus,x3,x2))/4
        else if(f==-10) then
            GetPDF=(FatXfrom2D(Gminus,x1,x2)+FatXfrom2D(Gminus,x3,x2))/2
        else if(f==-1) then
            GetPDF=-(FatXfrom2D(Dplus,x1,x2)-FatXfrom2D(Dplus,x3,x2)+FatXfrom2D(Dminus,x1,x2)+FatXfrom2D(Dminus,x3,x2))/4
        else if(f==-2) then
            GetPDF=-(FatXfrom2D(Uplus,x1,x2)-FatXfrom2D(Uplus,x3,x2)+FatXfrom2D(Uminus,x1,x2)+FatXfrom2D(Uminus,x3,x2))/4
        else if(f==-3) then
            GetPDF=-(FatXfrom2D(Splus,x1,x2)-FatXfrom2D(Splus,x3,x2)+FatXfrom2D(Sminus,x1,x2)+FatXfrom2D(Sminus,x3,x2))/4
        else if(f==-4) then
            GetPDF=-(FatXfrom2D(Cplus,x1,x2)-FatXfrom2D(Cplus,x3,x2)+FatXfrom2D(Cminus,x1,x2)+FatXfrom2D(Cminus,x3,x2))/4
        else if(f==-5) then
            GetPDF=-(FatXfrom2D(Bplus,x1,x2)-FatXfrom2D(Bplus,x3,x2)+FatXfrom2D(Bminus,x1,x2)+FatXfrom2D(Bminus,x3,x2))/4
        else
            write(*,*) WarningString("the flavor can be only -5,...,5, 10, -10. Zero returned."," ")
            GetPDF=0._dp
        end if
     CASE ('S')
        if(f==10.or.f==0) then
            GetPDF=(FatXfrom2D(Gplus,x1,x2)-FatXfrom2D(Gplus,x3,x2))/2
        else if(f==1) then
            GetPDF=-(FatXfrom2D(Dplus,x1,x2)+FatXfrom2D(Dminus,x1,x2))/4
        else if(f==2) then
            GetPDF=-(FatXfrom2D(Uplus,x1,x2)+FatXfrom2D(Uminus,x1,x2))/4
        else if(f==3) then
            GetPDF=-(FatXfrom2D(Splus,x1,x2)+FatXfrom2D(Sminus,x1,x2))/4
        else if(f==4) then
            GetPDF=-(FatXfrom2D(Cplus,x1,x2)+FatXfrom2D(Cminus,x1,x2))/4
        else if(f==5) then
            GetPDF=-(FatXfrom2D(Bplus,x1,x2)+FatXfrom2D(Bminus,x1,x2))/4
        else if(f==-10) then
            GetPDF=(FatXfrom2D(Gminus,x1,x2)+FatXfrom2D(Gminus,x3,x2))/2
        else if(f==-1) then
            GetPDF=-(FatXfrom2D(Dplus,x3,x2)-FatXfrom2D(Dminus,x3,x2))/4
        else if(f==-2) then
            GetPDF=-(FatXfrom2D(Uplus,x3,x2)-FatXfrom2D(Uminus,x3,x2))/4
        else if(f==-3) then
            GetPDF=-(FatXfrom2D(Splus,x3,x2)-FatXfrom2D(Sminus,x3,x2))/4
        else if(f==-4) then
            GetPDF=-(FatXfrom2D(Cplus,x3,x2)-FatXfrom2D(Cminus,x3,x2))/4
        else if(f==-5) then
            GetPDF=-(FatXfrom2D(Bplus,x3,x2)-FatXfrom2D(Bminus,x3,x2))/4
        else
            write(*,*) WarningString("the flavor can be only -5,...,5, 10, -10. Zero returned."," ")
            GetPDF=0._dp
        end if
    CASE DEFAULT
        write(*,*) ErrorString("unknown outputT. Evaluation STOP."," ")
        stop
    END SELECT

end function GetPDF

!!!!! Computes interpolation flavor f for chiral odd-case
!!!!! f=1=d
!!!!! f=2=u
!!!!! f=3=s
!!!!! f=4=c
!!!!! f=5=b
!!!!!
function GetPDFChiralOdd(x1,x2,f)
    real(dp)::GetPDFChiralOdd
    real(dp),intent(in)::x1,x2
    integer,intent(in)::f
    real(dp)::x3,r
    x3=-x1-x2

    r=max(abs(x1),abs(x2),abs(x1+x2))
    if(r<xMIN) then
        write(*,*) ErrorString("requested value with |x|<xMIN"," ")
        write(*,*) "Requested (x1,x2,x3)=(",x1,x2,x3,") with |x|=",r
        write(*,*) "Current xMIN =",xMIN
        write(*,*) "Increase lower limit in INI-file. Evaluation STOP"
        stop
    end if
    if(abs(x1)>1.or.abs(x2)>1.or.abs(x1+x2)>1) then
        write(*,*) WarningString("requested values at |x|>1. Zero returned."," ")
    GetPDFChiralOdd=0._dp
    end if

    if(f==1) then
        GetPDFChiralOdd=FatXfrom2D(Dodd,x1,x2)
    else if(f==2) then
        GetPDFChiralOdd=FatXfrom2D(Uodd,x1,x2)
    else if(f==3) then
        GetPDFChiralOdd=FatXfrom2D(Sodd,x1,x2)
    else if(f==4) then
        GetPDFChiralOdd=FatXfrom2D(Codd,x1,x2)
    else if(f==5) then
        GetPDFChiralOdd=FatXfrom2D(Bodd,x1,x2)
    else
        write(*,*) WarningString("the flavor can be only -5,...,5, 10, -10. Zero returned."," ")
        GetPDFChiralOdd=0._dp
    end if

end function GetPDFChiralOdd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! INTERFACES FOR EVOLUTION OF DITRIBUTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Evolution of PLUS components
!!! G,U,... = 1D flavor arrays to evolve
!!! alpha = alpha_s(mu)
!!! mu0, mu1 = initial/final scales in GeV
subroutine EvolvePLUS(alpha,mu0,mu1,G,U,D,S,C,B)
    real(dp),external::alpha
    real(dp),intent(in)::mu0,mu1
    real(dp)::t0,t1
    real(dp),dimension(0:NUM_TOT-1)::G,U,D,S,C,B
    real(dp),allocatable,dimension(:)::Singlet,N1,N2,N3,N4!!! singlet, and 3 non-singlet combinations
    integer::i

     if(useSingletEvolution) then

        if(mu0>mu1) then
            write(*,*) ErrorString("evolution from high-to-low scale is not implemented yet"," ")
            write(*,*) "       computation terminated."
            stop
        end if

        !!!! depending on the range of mu one need to allocate different set of variables
        !!!! in anycase singlet Singlet,N1,N2
        allocate(Singlet(0:NUM_TOT-1))
        allocate(N1(0:NUM_TOT-1))
        allocate(N2(0:NUM_TOT-1))
        if(mu1>massCHARM) allocate(N3(0:NUM_TOT-1))
        if(mu1>massBOTTOM) allocate(N4(0:NUM_TOT-1))

!         N1=U
!         N2=U
!         call EvNonSinglet(N1,as,0.d0,1.d0)
!         call EvNonSingletOLD(N2,as,0.d0,1.d0)
!         write(*,*) N1-N2
!         stop



        !!!! The boundary condition of the singlet combination depend on the range of mu
        if(mu0<massCHARM) then
            Singlet=U+D+S
        else if(mu0<massBOTTOM) then
            Singlet=U+D+S+C
        else
            Singlet=U+D+S+C+B
        end if
        !!! boundary conditions for non-singlet
        N1=U-D
        N2=U+D-2*S

        !!!!! evolving the singlet+gluon part
        !!!!! depending on range of mu, there are several thresholds
        !!!!! Evolution through the threshold defines boundary for N3 and N4 (if necesary)
        if(mu1<massCHARM) then
            !!!! CASE :  mu0<mu1<massCHARM  Nf=3

            !!!! t=log[mu^2]
            t0=2*log(mu0)
            t1=2*log(mu1)
            call EvSingletPLUS(Singlet,G,as,t0,t1,3)

        else if(mu1<massBOTTOM) then

            if(mu0<massCHARM) then
                !!!! CASE :  mu0<massCharm<mu1<massBOTTOM  Nf=3, then Nf=4
                !!!! t=log[mu^2]
                t0=2*log(mu0)
                t1=2*log(massCHARM)
                call EvSingletPLUS(Singlet,G,as,t0,t1,3)

                N3=Singlet !!! boundary conidition for N3-combination

                t0=t1
                t1=2*log(mu1)
                call EvSingletPLUS(Singlet,G,as,t0,t1,4)
            else
                N3=U+D+S-3*C !!! boundary condition for N3-combination

                !!!! CASE :  massCharm<mu0<mu1<massBOTTOM  Nf=4
                !!!! t=log[mu^2]
                t0=2*log(mu0)
                t1=2*log(mu1)
                call EvSingletPLUS(Singlet,G,as,t0,t1,4)
            end if
        else
            if(mu0<massCHARM) then

                !!!! CASE :  mu0<massCharm<massBOTTOM<mu1  Nf=3, then Nf=4, then Nf=5
                !!!! t=log[mu^2]
                t0=2*log(mu0)
                t1=2*log(massCHARM)
                call EvSingletPLUS(Singlet,G,as,t0,t1,3)

                N3=Singlet !!! boundary conidition for N3-combination

                t0=t1
                t1=2*log(massBOTTOM)
                call EvSingletPLUS(Singlet,G,as,t0,t1,4)

                N4=Singlet !!! boundary conidition for N4-combination

                t0=t1
                t1=2*log(mu1)
                call EvSingletPLUS(Singlet,G,as,t0,t1,5)

            else if(mu0<massBOTTOM) then

                !!!! CASE :  massCharm<mu0<massBOTTOM<mu1  Nf=4 then Nf=5
                !!!! t=log[mu^2]
                t0=2*log(mu0)
                t1=2*log(massBOTTOM)
                call EvSingletPLUS(Singlet,G,as,t0,t1,4)

                N3=U+D+S-3*C !!! boundary conidition for N3-combination
                N4=Singlet !!! boundary conidition for N4-combination

                t0=t1
                t1=2*log(mu1)
                call EvSingletPLUS(Singlet,G,as,t0,t1,5)

            else
                N3=U+D+S-3*C !!! boundary conidition for N3-combination
                N4=U+D+S+C-4*B !!! boundary conidition for N4-combination

                !!!! CASE :  massCharm<massBOTTOM<mu0<mu1  Nf=5
                !!!! t=log[mu^2]
                t0=2*log(mu0)
                t1=2*log(mu1)
                call EvSingletPLUS(Singlet,G,as,t0,t1,5)
            end if
        end if

        !!! non-singletes evolve as they are

        !!!! t=log[mu^2]
        t0=2*log(mu0)
        t1=2*log(mu1)

        call EvNonSinglet(N1,as,t0,t1)
        call EvNonSinglet(N2,as,t0,t1)

        !!!! the evolution for N3, and N4 is used only above thresholds.
        if(mu1>massCHARM) then
            if(mu0<massCHARM) then
                t0=2*log(massCHARM)
            end if
            call EvNonSinglet(N3,as,t0,t1)
        end if

        if(mu1>massBOTTOM) then
            if(mu0<massBOTTOM) then
                t0=2*log(massBOTTOM)
            else
                t0=2*log(mu0)
            end if
            call EvNonSinglet(N4,as,t0,t1)
        end if

        !!! FINALY we uncover the quark contributions
        if(mu1<massCHARM) then
            U=(3*N1+N2+2*Singlet)/6
            D=(-3*N1+N2+2*Singlet)/6
            S=(-N2+Singlet)/3
            C=0._dp
            B=0._dp
        else if(mu1<massBOTTOM) then
            U=(6*N1+2*N2+N3+3*Singlet)/12
            D=(-6*N1+2*N2+N3+3*Singlet)/12
            S=(-4*N2+N3+3*Singlet)/12
            C=(-N3+Singlet)/4
            B=0._dp
        else
            U=(30*N1+10*N2+5*N3+3*n4+12*Singlet)/60
            D=(-30*N1+10*N2+5*N3+3*n4+12*Singlet)/60
            S=(-20*N2+5*N3+3*n4+12*Singlet)/60
            C=(-5*N3+n4+4*Singlet)/20
            B=(-N4+Singlet)/5
        end if

    else

        !!! non-singlet evolution is Nf-independent
        !!! thus boundary as they are

        !!!! t=log[mu^2]
        t0=2*log(mu0)
        t1=2*log(mu1)

        call EvNonSinglet(U,as,t0,t1)
        call EvNonSinglet(D,as,t0,t1)
        call EvNonSinglet(S,as,t0,t1)
        call EvNonSinglet(C,as,t0,t1)
        call EvNonSinglet(B,as,t0,t1)

    end if

    !!!! store result globally as 2D
    Uplus=F1Dto2D(U)
    Dplus=F1Dto2D(D)
    Splus=F1Dto2D(S)
    Cplus=F1Dto2D(C)
    Bplus=F1Dto2D(B)
    Gplus=F1Dto2D(G)

    contains

    !!!!! as=alpha_s(mu)/(4pi)
    function as(t)
        real(dp)::as,t
        !!! 1/(4 pi)
        real(dp),parameter::pi4minus1=real(0.0795774715459476678844418816862571810,dp)
        as=alpha(Exp(t/2))*pi4minus1
    end function as

end subroutine EvolvePLUS

!!! Evolution of MINUS components
!!! G,U,... = 1D flavor arrays to evolve
!!! alpha = alpha_s(mu)
!!! mu0, mu1 = initial/final scales in GeV
subroutine EvolveMINUS(alpha,mu0,mu1,G,U,D,S,C,B)
    real(dp),external::alpha
    real(dp),intent(in)::mu0,mu1
    real(dp)::t0,t1
    real(dp),dimension(0:NUM_TOT-1)::G,U,D,S,C,B
    real(dp),allocatable,dimension(:)::Singlet,N1,N2,N3,N4!!! singlet, and 3 non-singlet combinations
    integer::i

    if(useSingletEvolution) then
        !!! singlet evolution

        if(mu0>mu1) then
            write(*,*) ErrorString("evolution from high-to-low scale is not implemented yet"," ")
            write(*,*) "       computation terminated."
            stop
        end if

        !!!! depending on the range of mu one need to allocate different set of variables
        !!!! in anycase singlet Singlet,N1,N2
        allocate(Singlet(0:NUM_TOT-1))
        allocate(N1(0:NUM_TOT-1))
        allocate(N2(0:NUM_TOT-1))
        if(mu1>massCHARM) allocate(N3(0:NUM_TOT-1))
        if(mu1>massBOTTOM) allocate(N4(0:NUM_TOT-1))


        !!!! The boundary condition of the singlet combination depend on the range of mu
        if(mu0<massCHARM) then
            Singlet=U+D+S
        else if(mu0<massBOTTOM) then
            Singlet=U+D+S+C
        else
            Singlet=U+D+S+C+B
        end if
        !!! boundary conditions for non-singlet
        N1=U-D
        N2=U+D-2*S

        !!!!! evolving the singlet+gluon part
        !!!!! depending on range of mu, there are several thresholds
        !!!!! Evolution through the threshold defines boundary for N3 and N4 (if necesary)
        if(mu1<massCHARM) then
            !!!! CASE :  mu0<mu1<massCHARM  Nf=3

            !!!! t=log[mu^2]
            t0=2*log(mu0)
            t1=2*log(mu1)
            call EvSingletMINUS(Singlet,G,as,t0,t1,3)

        else if(mu1<massBOTTOM) then

            if(mu0<massCHARM) then
                !!!! CASE :  mu0<massCharm<mu1<massBOTTOM  Nf=3, then Nf=4
                !!!! t=log[mu^2]
                t0=2*log(mu0)
                t1=2*log(massCHARM)
                call EvSingletMINUS(Singlet,G,as,t0,t1,3)

                N3=Singlet !!! boundary conidition for N3-combination

                t0=t1
                t1=2*log(mu1)
                call EvSingletMINUS(Singlet,G,as,t0,t1,4)
            else
                N3=U+D+S-3*C !!! boundary condition for N3-combination

                !!!! CASE :  massCharm<mu0<mu1<massBOTTOM  Nf=4
                !!!! t=log[mu^2]
                t0=2*log(mu0)
                t1=2*log(mu1)
                call EvSingletMINUS(Singlet,G,as,t0,t1,4)
            end if
        else
            if(mu0<massCHARM) then
                !!!! CASE :  mu0<massCharm<massBOTTOM<mu1  Nf=3, then Nf=4, then Nf=5
                !!!! t=log[mu^2]
                t0=2*log(mu0)
                t1=2*log(massCHARM)
                call EvSingletMINUS(Singlet,G,as,t0,t1,3)

                N3=Singlet !!! boundary conidition for N3-combination

                t0=t1
                t1=2*log(massBOTTOM)
                call EvSingletMINUS(Singlet,G,as,t0,t1,4)

                N4=Singlet !!! boundary conidition for N4-combination

                t0=t1
                t1=2*log(mu1)
                call EvSingletMINUS(Singlet,G,as,t0,t1,5)

            else if(mu0<massBOTTOM) then

                !!!! CASE :  massCharm<mu0<massBOTTOM<mu1  Nf=4 then Nf=5
                !!!! t=log[mu^2]
                t0=2*log(mu0)
                t1=2*log(massBOTTOM)
                call EvSingletMINUS(Singlet,G,as,t0,t1,4)

                N3=U+D+S-3*C !!! boundary conidition for N3-combination
                N4=Singlet !!! boundary conidition for N4-combination

                t0=t1
                t1=2*log(mu1)
                call EvSingletMINUS(Singlet,G,as,t0,t1,5)

            else
                N3=U+D+S-3*C !!! boundary conidition for N3-combination
                N4=U+D+S+C-4*B !!! boundary conidition for N4-combination

                !!!! CASE :  massCharm<massBOTTOM<mu0<mu1  Nf=5
                !!!! t=log[mu^2]
                t0=2*log(mu0)
                t1=2*log(mu1)
                call EvSingletMINUS(Singlet,G,as,t0,t1,5)
            end if
        end if

        !!! non-singletes evolve as they are

        !!!! t=log[mu^2]
        t0=2*log(mu0)
        t1=2*log(mu1)
        call EvNonSinglet(N1,as,t0,t1)
        call EvNonSinglet(N2,as,t0,t1)

        !!!! the evolution for N3, and N4 is used only above thresholds.
        if(mu1>massCHARM) then
            if(mu0<massCHARM) then
                t0=2*log(massCHARM)
            end if

            call EvNonSinglet(N3,as,t0,t1)
        end if

        if(mu1>massBOTTOM) then
            if(mu0<massBOTTOM) then
                t0=2*log(massBOTTOM)
            else
                t0=2*log(mu0)
            end if
            call EvNonSinglet(N4,as,t0,t1)
        end if

        !!! FINALY we uncover the quark contributions
        if(mu1<massCHARM) then
            U=(3*N1+N2+2*Singlet)/6
            D=(-3*N1+N2+2*Singlet)/6
            S=(-N2+Singlet)/3
            C=0._dp
            B=0._dp
        else if(mu1<massBOTTOM) then
            U=(6*N1+2*N2+N3+3*Singlet)/12
            D=(-6*N1+2*N2+N3+3*Singlet)/12
            S=(-4*N2+N3+3*Singlet)/12
            C=(-N3+Singlet)/4
            B=0._dp
        else
            U=(30*N1+10*N2+5*N3+3*N4+12*Singlet)/60
            D=(-30*N1+10*N2+5*N3+3*N4+12*Singlet)/60
            S=(-20*N2+5*N3+3*N4+12*Singlet)/60
            C=(-5*N3+N4+4*Singlet)/20
            B=(-N4+Singlet)/5
        end if

    else

        !!! non-singlet evolution is Nf-independent
        !!! thus boundary as they are

        !!!! t=log[mu^2]
        t0=2*log(mu0)
        t1=2*log(mu1)

        call EvNonSinglet(U,as,t0,t1)
        call EvNonSinglet(D,as,t0,t1)
        call EvNonSinglet(S,as,t0,t1)
        call EvNonSinglet(C,as,t0,t1)
        call EvNonSinglet(B,as,t0,t1)

    end if

    !!!! store result globally as 2D
    Uminus=F1Dto2D(U)
    Dminus=F1Dto2D(D)
    Sminus=F1Dto2D(S)
    Cminus=F1Dto2D(C)
    Bminus=F1Dto2D(B)
    Gminus=F1Dto2D(G)

    contains

    !!!!! as=alpha_s(mu)/(4pi)
    function as(t)
        real(dp)::as,t
        !!! 1/(4 pi)
        real(dp),parameter::pi4minus1=real(0.0795774715459476678844418816862571810,dp)
        as=alpha(Exp(t/2))*pi4minus1
    end function as

end subroutine EvolveMINUS


end module SnowFlake
