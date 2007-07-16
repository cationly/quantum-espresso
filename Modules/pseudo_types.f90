!
! Copyright (C) 2002-2003 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  extracted from module "pseudo_types" of FPMD
!

      MODULE pseudo_types

!  this module contains the definitions of several TYPE structures,
!  together with their allocation/deallocation routines

        USE kinds, ONLY: DP
        USE parameters, ONLY: cp_lmax, lmaxx
        use radial_grids, ONLY: ndmx, radial_grid_type
        !  USE ld1_parameters, ONLY: ndm, nwfsx
        USE parameters, ONLY: nwfsx=>nchix

        IMPLICIT NONE
        SAVE


TYPE :: paw_t
   !
   ! Type describing a PAW dataset (temporary).
   ! Functions are defined on a logarithmic radial mesh.
   !
   CHARACTER(LEN=2) :: symbol
   REAL (DP) :: zval
   REAL (DP) :: z
   CHARACTER(LEN=80) :: dft
   TYPE(radial_grid_type) :: grid
!   INTEGER        :: mesh      ! the size of the mesh
!   REAL (DP), POINTER :: r(:) !r (ndmx)     ! the mesh
!   REAL (DP), POINTER :: r2(:) !r2 (ndmx)    ! r^2
!   REAL (DP), POINTER :: sqrtr(:) !sqrtr (ndmx) ! sqrt(r)
!   REAL (DP) :: dx          ! log(r(i+1))-log(r(i))
!   REAL (DP) :: zmesh
   REAL (DP) :: rmatch_augfun  ! the matching radius for augmentation charges
   LOGICAL :: nlcc ! nonlinear core correction
   INTEGER :: nwfc ! number of wavefunctions/projectors
   INTEGER :: lmax ! maximum angular momentum of projectors
   INTEGER, POINTER :: l(:) !l(nwfsx) ! angular momentum of projectors
   INTEGER, POINTER :: ikk(:) !ikk(nwfsx) ! cutoff radius for the projectors
   INTEGER :: irc ! r(irc) = radius of the augmentation sphere
   REAL (DP), POINTER :: &
        oc(:), &!oc (nwfsx), & ! the occupations
        enl(:), &!enl (nwfsx), & ! the energy of the wavefunctions
        aewfc(:,:), &!aewfc (ndmx,nwfsx), &  ! all-electron wavefunctions
        pswfc(:,:), &!pswfc (ndmx,nwfsx),        & ! pseudo wavefunctions
        proj(:,:), &!proj (ndmx,nwfsx),     & ! projectors
        augfun(:,:,:,:), &!augfun(ndmx,nwfsx,nwfsx,0:2*lmaxx+1),      & ! augmentation functions, augfun(:,:,:,0) are AE augm. functions
        augmom(:,:,:), &!augmom(nwfsx,nwfsx,0:2*lmaxx) , & ! moments of the augmentation functions
        aeccharge(:), &!aeccharge (ndmx),  & ! AE core charge * 4PI r^2
        psccharge(:), &!psccharge (ndmx),  & ! PS core charge * 4PI r^2
        pscharge(:), &!pscharge (ndmx),  & ! PS charge * 4PI r^2
        aeloc(:), &!aeloc (ndmx),     & ! descreened AE potential: v_AE-v_H[n1]-v_XC[n1+nc]
        psloc(:), &!psloc (ndmx),     & ! descreened local PS potential: v_PS-v_H[n~+n^]-v_XC[n~+n^+n~c]
        kdiff(:,:), &!kdiff (nwfsx,nwfsx) ,&      ! kinetic energy differences
        dion(:,:) !dion(nwfsx,nwfsx)
!!!  Notes about screening:
!!!       Without nlcc, the local PSpotential is descreened with n~+n^ only.
!!!       The local AEpotential is descreened ALWAYS with n1+nc. This improves
!!!       the accuracy, and will not cost in the plane wave code (atomic
!!!       contribution only).
END TYPE paw_t

!
!============================================================================

!  BEGIN manual
!  TYPE DEFINITIONS

        TYPE pseudo_upf
          CHARACTER(LEN=80):: generated   ! 
          CHARACTER(LEN=80):: date_author ! Misc info
          CHARACTER(LEN=80):: comment     !
          CHARACTER(LEN=2) :: psd       ! Element label
          CHARACTER(LEN=20) :: typ      ! Pseudo type ( NC or US )
          LOGICAL  :: tvanp             ! .true. if Ultrasoft
          LOGICAL :: nlcc               ! Non linear core corrections
          CHARACTER(LEN=20) :: dft      ! Exch-Corr type
          REAL(DP) :: zp               ! z valence
          REAL(DP) :: etotps           ! total energy
          REAL(DP) :: ecutwfc          ! suggested cut-off for wfc
          REAL(DP) :: ecutrho          ! suggested cut-off for rho

          LOGICAL :: has_so             ! if .true. includes spin-orbit
          REAL(DP) :: xmin             ! the minimum x of the linear mesh
          REAL(DP) :: rmax             ! the maximum radius of the mesh
          REAL(DP) :: zmesh            ! the nuclear charge used for mesh
          REAL(DP) :: dx               ! the deltax of the linear mesh
          INTEGER, POINTER :: nn(:)      ! nn(nwfc)
          REAL(DP), POINTER :: rcut(:)  ! cut-off radius(nbeta)
          REAL(DP), POINTER :: rcutus(:)! cut-off ultrasoft radius (nbeta)
          REAL(DP), POINTER :: epseu(:) ! energy (nwfc)
          REAL(DP), POINTER :: jchi(:)  ! jchi(nwfc)
          REAL(DP), POINTER :: jjj(:)   ! jjj(nbeta)

          INTEGER :: nv                 ! UPF file version number
          INTEGER :: lmax               ! maximum angular momentum component
          INTEGER :: mesh               ! number of point in the radial mesh
          INTEGER :: nwfc               ! number of wavefunctions
          INTEGER :: nbeta              ! number of projectors
          CHARACTER(LEN=2), POINTER :: els(:)  ! els(nwfc)
          CHARACTER(LEN=2), POINTER :: els_beta(:)  ! els(nbeta)
          INTEGER, POINTER :: lchi(:)   ! lchi(nwfc)
          REAL(DP), POINTER :: oc(:)   ! oc(nwfc)
          REAL(DP), POINTER :: r(:)    ! r(mesh)
          REAL(DP), POINTER :: rab(:)  ! rab(mesh)
          REAL(DP), POINTER :: rho_atc(:) ! rho_atc(mesh)
          REAL(DP), POINTER :: vloc(:)    ! vloc(mesh)
          INTEGER, POINTER :: lll(:)       ! lll(nbeta)
          INTEGER, POINTER :: kkbeta(:)    ! kkbeta(nbeta)
          REAL(DP), POINTER :: beta(:,:)  ! beta(mesh,nbeta)
          INTEGER :: nd
          REAL(DP), POINTER :: dion(:,:)  ! dion(nbeta,nbeta)
          INTEGER :: nqf
          INTEGER :: nqlc
          REAL(DP), POINTER :: rinner(:)  ! rinner(0:2*lmax)
          REAL(DP), POINTER :: qqq(:,:)   ! qqq(nbeta,nbeta)
          REAL(DP), POINTER :: qfunc(:,:,:) ! qfunc(mesh,nbeta,nbeta)
          REAL(DP), POINTER :: qfcoef(:,:,:,:) ! qfcoef(nqf,0:2*lmax,nbeta,nbeta)
          REAL(DP), POINTER :: chi(:,:) !  chi(mesh,nwfc)
          REAL(DP), POINTER :: rho_at(:) !  rho_at(mesh)
        END TYPE


        TYPE pseudo_ncpp
          CHARACTER(LEN=4) :: psd         ! Element label
          CHARACTER(LEN=20) :: pottyp     ! Potential type
          LOGICAL :: tmix
          LOGICAL :: tnlcc
          INTEGER :: igau
          INTEGER :: lloc
          INTEGER :: nbeta 
          INTEGER :: lll(cp_lmax)
          INTEGER :: nchan
          INTEGER :: mesh
          REAL(DP) ::  zv
          REAL(DP) ::  dx            ! r(i) = cost * EXP( xmin + dx * (i-1) )
          REAL(DP) ::  rab(ndmx)
          REAL(DP) ::  rw(ndmx)
          REAL(DP) ::  vnl(ndmx, cp_lmax)
          REAL(DP) ::  vloc(ndmx)
          REAL(DP) ::  vrps(ndmx, cp_lmax)
          REAL(DP) ::  wgv(cp_lmax)
          REAL(DP) ::  rc(2)
          REAL(DP) ::  wrc(2)
          REAL(DP) ::  rcl(3,3)
          REAL(DP) ::  al(3,3)
          REAL(DP) ::  bl(3,3)
          INTEGER :: nrps                     ! number of atomic wave function
          INTEGER :: lrps(cp_lmax)            ! angular momentum
          REAL(DP) :: oc(cp_lmax)            ! occupation for each rps
          REAL(DP) :: rps(ndmx, cp_lmax)  ! atomic pseudo wave function
          REAL(DP) :: rhoc(ndmx)          ! core charge
        END TYPE pseudo_ncpp

!  ----------------------------------------------
!  END manual

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

SUBROUTINE nullify_pseudo_paw( paw )
  TYPE( paw_t ), INTENT(INOUT) :: paw
  NULLIFY( paw%l, paw%ikk )
  NULLIFY( paw%oc, paw%enl, paw%aewfc, paw%pswfc, paw%proj )
  NULLIFY( paw%augfun, paw%augmom, paw%aeccharge, paw%psccharge, paw%pscharge )
  NULLIFY( paw%aeloc, paw%psloc, paw%kdiff, paw%dion )
  RETURN
END SUBROUTINE nullify_pseudo_paw

SUBROUTINE allocate_pseudo_paw( paw, size_mesh, size_nwfc, size_lmax )
  TYPE( paw_t ), INTENT(INOUT) :: paw
  INTEGER, INTENT(IN) :: size_mesh, size_nwfc, size_lmax
  ALLOCATE ( paw%l(size_nwfc) )
  ALLOCATE ( paw%ikk(size_nwfc) )
  ALLOCATE ( paw%oc(size_nwfc) )
  ALLOCATE ( paw%enl(size_nwfc) )
  ALLOCATE ( paw%aewfc(size_mesh,size_nwfc) )
  ALLOCATE ( paw%pswfc(size_mesh,size_nwfc) )
  ALLOCATE ( paw%proj (size_mesh,size_nwfc) )
  ALLOCATE ( paw%augfun(size_mesh,size_nwfc,size_nwfc,0:2*size_lmax+2) )
  ALLOCATE ( paw%augmom(size_nwfc,size_nwfc,0:2*size_lmax+1) )
  ALLOCATE ( paw%aeccharge(size_mesh) )
  ALLOCATE ( paw%psccharge(size_mesh) )
  ALLOCATE ( paw%pscharge(size_mesh) )
  ALLOCATE ( paw%aeloc(size_mesh) )
  ALLOCATE ( paw%psloc(size_mesh) )
  ALLOCATE ( paw%kdiff(size_nwfc,size_nwfc) )
  ALLOCATE ( paw%dion (size_nwfc,size_nwfc) )
END SUBROUTINE allocate_pseudo_paw

SUBROUTINE deallocate_pseudo_paw( paw )
  TYPE( paw_t ), INTENT(INOUT) :: paw
  IF( ASSOCIATED( paw%l ) ) DEALLOCATE( paw%l )
  IF( ASSOCIATED( paw%ikk ) ) DEALLOCATE( paw%ikk )
  IF( ASSOCIATED( paw%oc ) ) DEALLOCATE( paw%oc )
  IF( ASSOCIATED( paw%enl ) ) DEALLOCATE( paw%enl )
  IF( ASSOCIATED( paw%aewfc ) ) DEALLOCATE( paw%aewfc )
  IF( ASSOCIATED( paw%pswfc ) ) DEALLOCATE( paw%pswfc )
  IF( ASSOCIATED( paw%proj ) ) DEALLOCATE( paw%proj )
  IF( ASSOCIATED( paw%augfun ) ) DEALLOCATE( paw%augfun )
  IF( ASSOCIATED( paw%augmom ) ) DEALLOCATE( paw%augmom )
  IF( ASSOCIATED( paw%aeccharge ) ) DEALLOCATE( paw%aeccharge )
  IF( ASSOCIATED( paw%psccharge ) ) DEALLOCATE( paw%psccharge )
  IF( ASSOCIATED( paw%pscharge ) ) DEALLOCATE( paw%pscharge )
  IF( ASSOCIATED( paw%aeloc ) ) DEALLOCATE( paw%aeloc )
  IF( ASSOCIATED( paw%psloc ) ) DEALLOCATE( paw%psloc )
  IF( ASSOCIATED( paw%kdiff ) ) DEALLOCATE( paw%kdiff )
  IF( ASSOCIATED( paw%dion ) ) DEALLOCATE( paw%dion )
  RETURN
END SUBROUTINE deallocate_pseudo_paw

!  subroutines

!  ----------------------------------------------

        SUBROUTINE nullify_pseudo_upf( upf )
          TYPE( pseudo_upf ), INTENT(INOUT) :: upf
          NULLIFY( upf%els, upf%lchi, upf%jchi, upf%oc )
          NULLIFY( upf%r, upf%rab )  
          NULLIFY( upf%rho_atc, upf%vloc )  
          NULLIFY( upf%nn, upf%rcut)
          NULLIFY( upf%els_beta)
          NULLIFY( upf%rcutus, upf%epseu)
          NULLIFY( upf%lll, upf%jjj, upf%kkbeta, upf%beta, upf%dion )  
          NULLIFY( upf%rinner, upf%qqq, upf%qfunc, upf%qfcoef )  
          NULLIFY( upf%chi )  
          NULLIFY( upf%rho_at )  
          RETURN
        END SUBROUTINE nullify_pseudo_upf

        SUBROUTINE deallocate_pseudo_upf( upf )
          TYPE( pseudo_upf ), INTENT(INOUT) :: upf
          IF( ASSOCIATED( upf%els ) ) DEALLOCATE( upf%els )
          IF( ASSOCIATED( upf%lchi ) ) DEALLOCATE( upf%lchi )
          IF( ASSOCIATED( upf%jchi ) ) DEALLOCATE( upf%jchi )
          IF( ASSOCIATED( upf%oc ) ) DEALLOCATE( upf%oc )
          IF( ASSOCIATED( upf%r ) ) DEALLOCATE( upf%r )
          IF( ASSOCIATED( upf%rab ) ) DEALLOCATE( upf%rab )
          IF( ASSOCIATED( upf%nn ) ) DEALLOCATE( upf%nn )
          IF( ASSOCIATED( upf%els_beta ) ) DEALLOCATE( upf%els_beta )
          IF( ASSOCIATED( upf%rcut ) ) DEALLOCATE( upf%rcut )
          IF( ASSOCIATED( upf%rcutus ) ) DEALLOCATE( upf%rcutus )
          IF( ASSOCIATED( upf%epseu ) ) DEALLOCATE( upf%epseu )
          IF( ASSOCIATED( upf%rho_atc ) ) DEALLOCATE( upf%rho_atc )
          IF( ASSOCIATED( upf%vloc ) ) DEALLOCATE( upf%vloc )
          IF( ASSOCIATED( upf%lll ) ) DEALLOCATE( upf%lll )
          IF( ASSOCIATED( upf%jjj ) ) DEALLOCATE( upf%jjj )
          IF( ASSOCIATED( upf%kkbeta ) ) DEALLOCATE( upf%kkbeta )
          IF( ASSOCIATED( upf%beta ) ) DEALLOCATE( upf%beta )
          IF( ASSOCIATED( upf%dion ) ) DEALLOCATE( upf%dion )
          IF( ASSOCIATED( upf%rinner ) ) DEALLOCATE( upf%rinner )
          IF( ASSOCIATED( upf%qqq ) ) DEALLOCATE( upf%qqq )
          IF( ASSOCIATED( upf%qfunc ) ) DEALLOCATE( upf%qfunc )
          IF( ASSOCIATED( upf%qfcoef ) ) DEALLOCATE( upf%qfcoef )
          IF( ASSOCIATED( upf%chi ) ) DEALLOCATE( upf%chi )
          IF( ASSOCIATED( upf%rho_at ) ) DEALLOCATE( upf%rho_at )
          RETURN
        END SUBROUTINE deallocate_pseudo_upf

      END MODULE pseudo_types

