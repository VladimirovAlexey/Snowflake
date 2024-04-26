!!!! Working precision of variables
integer, parameter :: dp = selected_real_kind(15, 307)

!!!! Print initialization information
logical::showINI=.false.

!!!! GRID type is defined in the top of the HexGrid.f90

!!!! size of the grid-points along perimeter [k=0,...,6 NUM_phi-1]
integer::NUM_phi=10
!!!! total size of the perimeter [k=0,...,6 NUM_phi-1]
integer::NUM_perimeter=59!6*NUM_phi-1
!!!! size of the grid-points along radius [n=0,...,NUM_R]
integer::NUM_R=15
!!!! total size of the grid
integer::NUM_TOT=960!(NUM_perimeter+1)*(NUM_R+1)
!!!! minimal accesable x (r=1)
real(dp)::xMIN=0.01_dp



!!!! Value of zero to compare with
real(dp)::zero=1d-12
!!!! Integration tolerance
real(dp)::toleranceINT=1d-8
!!!! Runge-Kutta maximal step
real(dp)::RGh_max=0.001_dp
!!!! Number of processeor allows for parralel computation
integer::allowedNumProcessor=10

!!!! Setup chiral even evolution
logical:: IncludeChiralEvenEvolution=.true.
!!!! Setup chiral odd evolution
logical:: IncludeChiralOddEvolution=.true.
!!!! Include mixing with gluon, other flavors etc.
logical:: useSingletEvolution=.true.

!!!! Mass of the CHARM threshold [GeV]
real(dp):: massCHARM=1.27_dp
!!!! Mass of the BOTTOM threshold [GeV]
real(dp):: massBOTTOM=4.18_dp
