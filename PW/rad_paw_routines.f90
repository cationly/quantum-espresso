!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"
! FIXME: modularize code and remove the following includes
#include "../atomic/hartree.f90"
#include "../atomic/series.f90"
!
MODULE rad_paw_routines
  !
  IMPLICIT NONE
  PUBLIC
  SAVE
CONTAINS

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! This function computes energy from potential and density on a radial grid
!!! if e_na is provided it's filled with per atom energies
!!

SUBROUTINE PAW_energy(becsum)
    USE kinds,                  ONLY : DP
    USE lsda_mod,               ONLY : nspin
    USE parameters,             ONLY : ndmx
    USE ions_base,              ONLY : nat, ityp

    USE grid_paw_variables,     ONLY : pfunc, ptfunc, tpawp
    USE uspp_param,             ONLY : augfun, nhm, lmaxq

    REAL(DP), INTENT(IN)         :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations

    INTEGER, PARAMETER           :: AE = 1, PS = 2   ! All-Electron and Pseudo
    INTEGER                      :: i_what           ! counter on AE and PS
    INTEGER                      :: na,nt,first_nat,last_nat ! atoms counters and indexes

    REAL(DP), ALLOCATABLE   :: rho_lm(:,:,:)      ! radial density expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: v_h_lm(:,:,:)      ! hartree potential
    REAL(DP)                :: e_h(lmaxq**2,nspin)! hartree energy components
    REAL(DP)                :: e                  ! placeholder

    CALL start_clock ('PAW_energy')
    whattodo: DO i_what = AE, PS
        ! The following operations will be done, first on AE, then on PS part
        ! furthermore they will be done for one atom at a time (to reduce memory usage)
        ! in the future code will be parallelized on atoms:
        !
        ! 1. build rho_lm (PAW_rho_lm)
        ! 2. compute Hartree energy
        !   a. use rho_lm to compute hartree potential (PAW_v_h)
        !   b. use v_h & rho_lm to build hartree energy (PAW_h_energy)
        ! 3. compute XC energy
        !   a. build rho1_rad(theta, phi) from rho_lm
        !   b. build v_xc(theta, phi)
        !   c. compute E_xc(theta, phi)
        !   d. repeat a,b,c integrating om theta & phi
        ! 4. free memory (dealloc rho_lm) and return
        !
        first_nat = 1
        last_nat  = nat
        ! Operations from here on are (will be) parallelized on atoms
        atoms: DO na = first_nat, last_nat
        ifpaw: IF (tpawp(nt)) THEN
            nt = ityp(na) ! the type of atom "na"
            ! STEP: 1 [ build rho_lm (PAW_rho_lm) ]
            ALLOCATE(rho_lm(ndmx,lmaxq**2,nspin))
            IF (i_what == AE) THEN
                ! passing "na" as an argument is dirtyer but faster and
                ! uses less memory than passing only a hyperslice of the array
                CALL PAW_rho_lm(na, becsum, pfunc, rho_lm)
            ELSE
                CALL PAW_rho_lm(na, becsum, ptfunc, rho_lm, augfun)
                !      optional argument for pseudo part --> ^^^^^^
            ENDIF
            ! STEP: 2 [ compute Hartree energy ]
            ALLOCATE(v_h_lm(ndmx,lmaxq**2,nspin))
            !   2a. use rho_lm to compute hartree potential (PAW_v_h)
            CALL PAW_v_h(na, rho_lm, v_h_lm)
            !   2b. use v_h & rho_lm to build hartree energy (PAW_h_energy)
            e = PAW_h_energy(na, rho_lm, v_h_lm, e_h)
            WRITE(6,*) "******************************"
            WRITE(6,*) e
            WRITE(6,*) e_h
            !
            DEALLOCATE(v_h_lm)
            ! STEP: 3 [ compute XC energy ]
            DEALLOCATE(rho_lm)

        ENDIF ifpaw
        ENDDO atoms
    ENDDO whattodo
    CALL stop_clock ('PAW_energy')


END SUBROUTINE PAW_energy

FUNCTION PAW_h_energy(na, pot_lm, rho_lm, e_lm)
    USE kinds,                  ONLY : DP
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq
    USE atom,                   ONLY : r, rab, mesh, msh
    USE parameters,             ONLY : ndmx
    USE ions_base,              ONLY : nat, ityp

    INTEGER,  INTENT(IN)             :: na ! atom index
    REAL(DP), INTENT(IN)             :: rho_lm(ndmx,lmaxq**2,nspin)! in:  rho
    REAL(DP), INTENT(IN)             :: pot_lm(ndmx,lmaxq**2,nspin)! in:  potential 
    REAL(DP),OPTIONAL,INTENT(OUT)    :: e_lm(lmaxq**2,nspin)    ! out: energy components 
    REAL(DP)                         :: PAW_h_energy            ! total hartree energy

                                        ! optional output: energy per atom

    REAL(DP)                :: aux(ndmx)   ! workspace
    REAL(DP)                :: par_energy  ! workspace
    INTEGER                 :: ispin,lm !counters on atoms, spins, angular momentum
    CALL start_clock ('PAW_h_energy')

    PAW_h_energy = 0._dp
    IF (present(e_lm)) e_lm(:,:) = 0._dp

    spins: DO ispin = 1,nspin
        DO lm = 1,lmaxq**2
            aux(:) = rho_lm(:,lm,ispin) * pot_lm(:,lm,ispin)
            ! V*rho is integrated using simpson algorithm
            CALL simpson (msh(ityp(na)),aux,rab(1,ityp(na)),par_energy)
            PAW_h_energy = PAW_h_energy + par_energy
            IF (present(e_lm)) e_lm(lm,ispin) = par_energy
        ENDDO
    ENDDO spins

    CALL stop_clock ('PAW_h_energy')

END FUNCTION PAW_h_energy

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! use the density produced by sum_rad_rho to compute potential (only hartree, at this moment)
!!
SUBROUTINE PAW_v_h(na, rho_lm, pot_lm)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi
    USE parameters,             ONLY : ndmx, npsx
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq
    USE ions_base,              ONLY : nat, ityp, ntyp => nsp
    USE atom,                   ONLY : r, msh, mesh

    INTEGER,  INTENT(IN)  :: na
    REAL(DP), INTENT(IN)  :: rho_lm(ndmx,lmaxq**2,nspin)! charge density as lm components
    REAL(DP), INTENT(OUT) :: pot_lm(ndmx,lmaxq**2,nspin)! potential as lm components
    REAL(DP)              :: auxrho(ndmx)               ! workspace

    REAL(DP)              :: r2(ndmx,npsx)  ! r**2
    REAL(DP)              :: rs(ndmx,npsx)  ! r**.5
    REAL(DP)              :: dx(npsx)       ! integration step used in atomic code
    INTEGER               :: nt,&           ! atom type
                             ispin, &       ! counter on spins
                             lm,l           ! counter on composite angmom lm = l**2 +m
!     INTEGER :: k !DEBUG

!     REAL(DP)                     :: e(nat,2) !DEBUG
!     REAL(DP)                     :: ecomps(lmaxq**2,nspin,nat,2) !DEBUG
    CALL start_clock ('PAW_v_h')

    ! get type of atom
    nt = ityp(na)

    ! prepare r**2 and r**.5 arrays                     !FIXME: move to init
    r2(:,:) = r(:,:)**2
    rs(:,:) = sqrt(r(:,:))

    !  I have to initialize the integration step dx:    !FIXME: move to init
    dx(:) = 0._dp
!     DO nt = 1,ntyp
    IF (r(1,nt) > 0.0_dp) THEN
        ! r(i+1) = exp(xmin)/zmesh * exp(i*dx)
        dx(nt)=log(r(2,nt)/r(1,nt))
    ELSE
        ! r(i+1) = exp(xmin)/zmesh * ( exp(i*dx) - 1 )
        dx(nt)=log( (r(3,nt)-r(2,nt)) / r(2,nt) )
    ENDIF
!     ENDDO


    ! this loop computes the hartree potential using the following formula:
    !               l is the first argument in hartree subroutine
    !               r1 = min(r,r'); r2 = MAX(r,r')
    ! V_h(r) = \sum{lm} Y_{lm}(\hat{r})/(2L+1) \int dr' 4\pi r'^2 \rho^{lm}(r') (r1^l/r2^{l+1})
    !     done here --> ^^^^^^^^^^^^^^^^^^^^^           ^^^^^^^^^^^^^^^^^^^^^^ <-- input to the hartree subroutine
    !                 output from the h.s. --> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !pot(:,:,:) = 0._dp
    spins: DO ispin = 1,nspin
        DO lm = 1, lmaxq**2
            l = INT(sqrt(DBLE(lm-1)))     ! l has to start from *zero*
                auxrho(:) = fpi/(2*l+1)*rho_lm(:,lm,ispin)
                CALL hartree(l, 2*l+2, mesh(nt), r(:,nt),r2(:,nt),rs(:,nt),&
                                dx(nt),auxrho(:), pot_lm(:,lm,ispin))
        ENDDO ! lm
    ENDDO spins

    CALL stop_clock ('PAW_v_h')

END SUBROUTINE PAW_v_h

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! sum up pfuncs x occupation in order to build radial density's angular momentum components
!!
! CALL PAW_rho_lm(becsum(:,na,:), ptfunc(:,:,:,nt), rho_lm, augfun(:,:,:,:,nt))
SUBROUTINE PAW_rho_lm(na, becsum, pfunc, rho_lm, augfun)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : eps8
    USE atom,                   ONLY : r, rab, mesh, msh
    USE ions_base,              ONLY : ityp, ntyp => nsp, nat 
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq!, augfun
    USE uspp,                   ONLY : indv, ap, nhtolm,lpl,lpx
    USE parameters,             ONLY : ndmx, nbrx, lqmax
    USE grid_paw_variables,     ONLY : okpaw, tpawp!, pfunc, ptfunc

    INTEGER,  INTENT(IN)         :: na     ! index of atom to use
    REAL(DP), INTENT(IN)         :: becsum(nhm*(nhm+1)/2,nat,nspin)    ! cross band occupation
    REAL(DP), INTENT(IN)         :: pfunc(ndmx,nbrx,nbrx,ntyp)         ! psi_i * psi_j
    REAL(DP), INTENT(OUT)        :: rho_lm(ndmx,lmaxq**2,nspin)       ! AE charge density on rad. grid
    REAL(DP), OPTIONAL,INTENT(IN):: augfun(ndmx,nbrx,nbrx,0:lqmax,ntyp)! augmentation functions (only for PS part)

    REAL(DP)                :: pref ! workspace (ap*becsum)

    INTEGER                 :: nb,mb, &     ! counters for pfunc nb,mb = 1, nh
                               nmb, &       ! composite "triangular" index for pfunc nmb = 1,nh*(nh+1)/2
                               lm,lp,l,m, & ! counters for angular momentum lm = l**2+m
                               ispin,&      ! counter for spin (FIXME: may be unnecessary)
                               nt           ! type of atom na
    INTEGER :: k !DEBUG

    ! This functions computes the angular momentum components of rho
    ! using the following formula:
    !   rho(\vec{r}) = \sum_{LM} Y_{LM} \sum_{i,j} (\hat{r}) a_{LM}^{(lm)_i(lm)_j} becsum_ij pfunc_ij(r)
    !
    ! actually different angular momentum components are stored separately:
    !   rho^{LM}(\vec{r}) = \sum_{i,j} (\hat{r}) a_{LM}^{(lm)_i(lm)_j} becsum_ij pfunc_ij(r)
    !
    ! notice that pfunc's are already multiplied by r^2 and they are indexed on the atom
    ! (they only depends on l, not on m), the augmentation charge depend only on l
    ! but the becsum depend on both l and m

    CALL start_clock ('PAW_rho_lm')

    ! get type of atom
    nt = ityp(na)

    ! initialize density
    rho_lm(:,:,:) = 0._dp

    !ifpaw: IF (tpawp(nt)) THEN !FIXME:: do this check outside
    spins: DO ispin = 1, nspin
    nmb = 0
        ! loop on all pfunc for this kind of pseudo
        DO nb = 1, nh(nt)
        DO mb = nb, nh(nt)
            nmb = nmb+1 ! nmb = 1, nh*(nh+1)/2
!           DO lm = 1, lmaxq**2  ! loop on all lm has been replaced by lpl+lpx trick
            DO lp = 1, lpx (nhtolm(mb,nt), nhtolm(nb,nt)) !lmaxq**2
                lm = lpl (nhtolm(mb,nt), nhtolm(nb,nt), lp)
                ! becsum already contains a factor 2 for off-diagonal pfuncs
                pref = becsum(nmb,na,ispin) * ap(lm, nhtolm(nb,nt), nhtolm(mb,nt))
                !
                rho_lm(1:msh(nt),lm,ispin) = rho_lm(1:msh(nt),lm,ispin) +&
                                pref * pfunc(1:msh(nt), indv(nb,nt), indv(mb,nt), nt)
                IF (present(augfun)) THEN
                    ! if I'm doing the pseudo part I have to add the augmentation charge
                    l = INT(sqrt(DBLE(lm-1))) ! l has to start from zero
                    rho_lm(1:msh(nt),lm,ispin) = rho_lm(1:msh(nt),lm,ispin) +&
                                pref * augfun(1:msh(na), indv(nb,nt), indv(mb,nt), l, nt)
                ENDIF ! augfun
            ENDDO ! lm
!                     DO k = 1, msh(nt)
!                         WRITE(90000+1000*i_what+100*na+10*nb+mb,"(2f15.8)")&
!                                r(k,nt), pfunc_(k,indv(nb,nt),indv(mb,nt),nt)
!                     ENDDO
        ENDDO !mb
        ENDDO !nb
    ENDDO spins
!   ENDIF ifpaw

    CALL stop_clock ('PAW_rho_lm')

END SUBROUTINE PAW_rho_lm

! REMOVE ME:
#include "rad_paw_trash.f90"

END MODULE rad_paw_routines
