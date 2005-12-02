module grid_paw_variables
  !
  !   WARNINGS:
  !
  ! NO spin-orbit
  ! NO EXX
  ! NO Parallelism
  ! NO rinner > 0
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : lqmax, nbrx, npsx, nqfx, ndmx
  !
  implicit none
  public!              <===
  save

  ! Analogous to okvan in  "uspp_param" (Modules/uspp.f90)
  LOGICAL :: &
       okpaw              ! if .TRUE. at least one pseudo is PAW

  ! Analogous to tvanp in "uspp_param" (Modules/uspp.f90)
  LOGICAL :: &
       tpawp(npsx)            ! if .TRUE. the atom is of PAW type

  ! Analogous to qfunc in "uspp_param" (Modules/uspp.f90)
  REAL(DP), TARGET :: &
       pfunc(ndmx,nbrx,nbrx,npsx), &! AE: \phi_{mu}(|r|)-\phi_{nu}(|r|)
       ptfunc(ndmx,nbrx,nbrx,npsx)  ! PS: \tilde{\phi}_{mu}(|r|)-\tilde{\phi}_{nu}(|r|)

  ! Analogous to qq in "uspp_param" (Modules/uspp.f90)
  REAL(DP), ALLOCATABLE, TARGET :: &
       pp(:,:,:),             &! the integrals of p functions in the solid
       ppt(:,:,:)              ! the integrals of pt functions in the solid

  ! Analogous to qrad in "us" (PW/pwcom.f90)
  REAL(DP), ALLOCATABLE, TARGET :: &
       prad(:,:,:,:),         &! radial FT of P functions
       ptrad(:,:,:,:)          ! radial FT of \tilde{P} functions

  ! Analogous to rho in "scf" (PW/pwcom.f90) + index scanning atoms
  REAL(DP), ALLOCATABLE, TARGET :: &
       rho1(:,:,:),             &! 1center AE charge density in real space
       rho1t(:,:,:),            &! 1center PS charge density in real space
       rho1h(:,:,:)              ! 1center compensation charge in real space

  ! Analogous to vr in "scf" (PW/pwcom.f90) + index scanning atoms
  REAL(DP), ALLOCATABLE, TARGET :: &
       vr1(:,:,:),        &! the Hartree+XC potential in real space of rho1
       vr1t(:,:,:)         ! the Hartree+XC potential in real space of rho1t

  ! Analogous to qq in "uspp_param" (Modules/uspp.f90)
  REAL(DP), ALLOCATABLE, TARGET :: &
       int_r2pfunc(:,:,:),   &! Integrals of r^2 * pfunc(r) (AE)
       int_r2ptfunc(:,:,:)    ! Integrals of r^2 * pfunc(r) (PS)

  ! Analogous to vloc_at in "uspp_param" (Modules/uspp.f90)
  REAL(DP), TARGET :: &
      aevloc_at(ndmx,npsx),               &! AE descreened potential
      psvloc_at(ndmx,npsx)                 ! PS descreened potential

end module grid_paw_variables
