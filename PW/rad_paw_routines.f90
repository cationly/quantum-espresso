!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE rad_paw_routines
  !
  IMPLICIT NONE
  PUBLIC
  SAVE
CONTAINS

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! This is the main driver of PAW routines, it uses directly or indirectly
!!! all the other routines of the module
!!
SUBROUTINE PAW_energy(becsum)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : pi
    USE lsda_mod,               ONLY : nspin
    USE radial_grids,           ONLY : ndmx
    USE ions_base,              ONLY : nat, ityp

    USE grid_paw_variables,     ONLY : pfunc, ptfunc, tpawp
    USE uspp_param,             ONLY : augfun, nhm, lmaxq

    REAL(DP), INTENT(IN)    :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations

    INTEGER, PARAMETER      :: AE = 1, PS = 2   ! All-Electron and Pseudo
    INTEGER                 :: i_what           ! counter on AE and PS
    INTEGER                 :: na,nt,first_nat,last_nat ! atoms counters and indexes

    ! hartree energy scalar fields expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: rho_lm(:,:,:)      ! radial density expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: v_h_lm(:,:,:)      ! hartree potential
    REAL(DP)                :: e_h(lmaxq**2,nspin)! hartree energy components
    REAL(DP)                :: e                  ! placeholder
    INTEGER                 :: n!debug
    ! BEWARE THAT HARTREE ONLY DEPENDS ON THE TOTAL RHO NOT ON RHOUP AND RHODW SEPARATELY...
    ! TREATMENT OF NSPIN>1 MUST BE CHECKED AND CORRECTED
    
    ! XC energy scalar fields sampled in radial coordinates r,th(eta),ph(i)
    REAL(DP)                :: th,ph ! placeholder
    REAL(DP), ALLOCATABLE   :: rho_rad(:,:)
    INTEGER :: lm

    CALL start_clock ('PAW_energy')

    WRITE(6,*) "***************************************************************"

    ! CHECK: it may be better to switch "whattodo" and "atoms" loops: it 
    !        saves one IF and may have additional pros with parallel
    ! CHECk: maybe we don't need to alloc/dealloc rho_lm every time
    whattodo: DO i_what = AE, PS
        ! The following operations will be done, first on AE, then on PS part
        ! furthermore they will be done for one atom at a time (to reduce memory usage)
        ! in the future code will be parallelized on atoms:
        !
        ! 1. build rho_lm (PAW_rho_lm)
        ! 2. compute Hartree energy
        !   a. use rho_lm to compute hartree potential (PAW_v_h)
        !   b. use v_h_lm & rho_lm to build hartree energy (PAW_h_energy)
        ! 3. compute XC energy
        !   a. build rho_rad(theta, phi) from rho_lm
        !   b. build v_xc(theta, phi)
        !   c. compute E_xc(theta, phi)
        !   d. repeat a,b,c integrating on theta & phi
        ! 4. free memory (dealloc rho_lm) and return
        !
        first_nat = 1
        last_nat  = nat
        ! Operations from here on are (will be) parallelized on atoms
        atoms: DO na = first_nat, last_nat
        nt = ityp(na) ! the type of atom na
        ifpaw: IF (tpawp(nt)) THEN
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
            !   2b. use v_h_lm & rho_lm to build hartree energy (PAW_h_energy)
            e = PAW_h_energy(na, rho_lm, v_h_lm, e_h)
            WRITE(6,*) "******************************"
!            WRITE(6,*) e
            WRITE(6,*) "==PAW RADIAL ENERGY: ", e
            WRITE(6,'(a,i1,a,f15.7)') ("==RADIAL PAW ENERGY (LM=",lm,"):",e_h(lm,1),lm=1,lmaxq**2)
            !
            ! STEP: 3 [ compute XC energy ]
            th = pi/3._dp
            ph   = pi/4._dp
            ALLOCATE(rho_rad(ndmx,nspin))
!            CALL PAW_rho_rad(th,ph, rho_lm, rho_rad)
!            DO n = 10,10000,10
               CALL start_clock('sph10')
                n = 10
                e = PAW_sph_integral(n, rho_lm, v_h_lm)
               CALL stop_clock('sph10')
                write(6,*) n,e
!            ENDDO
                !
                n = 100
                CALL start_clock('sph100')
                e = PAW_sph_integral(n, rho_lm, v_h_lm)
                CALL stop_clock('sph100')
                write(6,*) n,e
                !
                n = 1000
                CALL start_clock('sph1000')
                e = PAW_sph_integral(n, rho_lm, v_h_lm)
                CALL stop_clock('sph1000')
                write(6,*) n,e
                !
                n = 10000
                CALL start_clock('sph10000')
                e = PAW_sph_integral(n, rho_lm, v_h_lm)
                CALL stop_clock('sph10000')
                write(6,*) n,e


            DEALLOCATE(rho_rad)
            DEALLOCATE(v_h_lm)
            !
            DEALLOCATE(rho_lm)

        ENDIF ifpaw
        ENDDO atoms
    ENDDO whattodo

    CALL print_clock('sph10')
    CALL print_clock('sph100')
    CALL print_clock('sph1000')
    CALL print_clock('sph10000')

    WRITE(6,*) "***************************************************************"

    CALL stop_clock ('PAW_energy')


END SUBROUTINE PAW_energy

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! This function computes energy from potential and density on a radial grid
!!! if e_lm is provided it's filled with energy of single lm components
!!
FUNCTION PAW_h_energy(na, pot_lm, rho_lm, e_lm)
    USE kinds,                  ONLY : DP
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq
    USE atom,                   ONLY : r, rab, mesh, msh
    USE radial_grids,           ONLY : ndmx
    USE ions_base,              ONLY : ityp

    INTEGER,  INTENT(IN)           :: na ! atom index
    REAL(DP), INTENT(IN)           :: rho_lm(ndmx,lmaxq**2,nspin)! in:  rho
    REAL(DP), INTENT(IN)           :: pot_lm(ndmx,lmaxq**2,nspin)! in:  potential 
    REAL(DP), OPTIONAL,INTENT(OUT) :: e_lm(lmaxq**2,nspin)    ! out: energy components 
    REAL(DP)                       :: PAW_h_energy            ! total hartree energy

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
    USE parameters,             ONLY : npsx
    USE radial_grids,           ONLY : ndmx, hartree
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq
    USE ions_base,              ONLY : ityp
    USE atom,                   ONLY : rgrid

    INTEGER,  INTENT(IN)  :: na
    REAL(DP), INTENT(IN)  :: rho_lm(ndmx,lmaxq**2,nspin)! charge density as lm components
    REAL(DP), INTENT(OUT) :: pot_lm(ndmx,lmaxq**2,nspin)! potential as lm components
    REAL(DP)              :: auxrho(ndmx)               ! workspace

    INTEGER               :: nt,&           ! atom type
                             ispin, &       ! counter on spins
                             lm,l           ! counter on composite angmom lm = l**2 +m
!     INTEGER :: k !DEBUG

!     REAL(DP)                     :: e(nat,2) !DEBUG
!     REAL(DP)                     :: ecomps(lmaxq**2,nspin,nat,2) !DEBUG
    CALL start_clock ('PAW_v_h')

    ! get type of atom
    nt = ityp(na)
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
                CALL hartree(l, 2*l+2, rgrid(nt)%mesh, rgrid(nt), auxrho(:), pot_lm(:,lm,ispin))
        ENDDO ! lm
    ENDDO spins

    CALL stop_clock ('PAW_v_h')

END SUBROUTINE PAW_v_h

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! sum up pfuncs x occupation in order to build radial density's angular momentum components
!!
SUBROUTINE PAW_rho_lm(na, becsum, pfunc, rho_lm, augfun)
    USE kinds,                  ONLY : DP
    USE atom,                   ONLY : r, rab, mesh, msh
    USE ions_base,              ONLY : ityp, ntyp => nsp, nat 
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq!, augfun
    USE uspp,                   ONLY : indv, ap, nhtolm,lpl,lpx
    USE parameters,             ONLY : nbrx, lqmax
    USE radial_grids,           ONLY : ndmx
!    USE grid_paw_variables,     ONLY : okpaw, tpawp!, pfunc, ptfunc

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
        ENDDO !mb
        ENDDO !nb
    ENDDO spins
!   ENDIF ifpaw

    CALL stop_clock ('PAW_rho_lm')

END SUBROUTINE PAW_rho_lm

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! build radial charge distribution from its spherical harmonics expansion
!!
SUBROUTINE PAW_rho_rad(th,ph, rho_lm, rho_rad)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : eps8, pi
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE parameters,             ONLY : nbrx, lqmax
    USE radial_grids,           ONLY : ndmx
!     USE atom,                   ONLY : r

    REAL(DP), INTENT(IN)        :: rho_lm(ndmx,lmaxq**2,nspin)! Y_lm expansion of rho
    REAL(DP), INTENT(IN)        :: th, ph                     ! direction on which rho_rad is summed
    REAL(DP), INTENT(OUT)       :: rho_rad(ndmx,nspin)        ! charge density on rad. grid

    REAL(DP)                    :: v(3,1)           ! the versor pointed by angles th,ph
    REAL(DP)                    :: sin_ph           ! aux
    REAL(DP),PARAMETER          :: nv(1) = (/1._dp/)! it has to be a vector
    REAL(DP)                    :: ylm(1,lmaxq**2)  ! the spherical harmonics

    INTEGER                     :: ispin, lm ! counters on angmom and spin
!     INTEGER                     :: k !debug

    CALL start_clock ('PAW_rho_rad')
    rho_rad(:,:) = 0._dp

    ! prepare the versor
    sin_ph = sin(ph)
    v(1,1) = cos(th) * sin_ph
    v(2,1) = sin(th) * sin_ph
    v(3,1) = cos(ph)
!     WRITE(100,*) 0._dp,0._dp,0._dp
!     WRITE(100,*) v(:,1)
!     WRITE(101,*) v(:,1)
!     WRITE(100,*) 0._dp,0._dp,0._dp

    CALL ylmr2(lmaxq**2, 1, v, nv, ylm)
    rho_rad(:,:) = 0._dp
    spins: DO ispin = 1,nspin
        DO lm = 1, lmaxq**2
            !IF (ABS(ylm(1,lm)) < eps8 ) CONTINUE
            rho_rad(:,ispin) = rho_rad(:,ispin) +&
                    ylm(1,lm)*rho_lm(:,lm,ispin)
        ENDDO ! lm
    ENDDO spins

    CALL stop_clock ('PAW_rho_rad')

END SUBROUTINE PAW_rho_rad
!------------------------------------------------------------
SUBROUTINE PAW_rho_rad2(x, rho_lm, rho_rad)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : eps8, pi
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE parameters,             ONLY : nbrx, lqmax
    USE radial_grids,           ONLY : ndmx
!     USE atom,                   ONLY : r

    REAL(DP), INTENT(IN)        :: rho_lm(ndmx,lmaxq**2,nspin)! Y_lm expansion of rho
    REAL(DP), INTENT(IN)        :: x(3)
    REAL(DP), INTENT(OUT)       :: rho_rad(ndmx,nspin)        ! charge density on rad. grid

    REAL(DP)                    :: v(3,1)           ! the versor pointed by angles th,ph
    REAL(DP)                    :: sin_ph           ! aux
    REAL(DP),PARAMETER          :: nv(1) = (/1._dp/)! it has to be a vector
    REAL(DP)                    :: ylm(1,lmaxq**2)  ! the spherical harmonics

    INTEGER                     :: ispin, lm ! counters on angmom and spin
!     INTEGER                     :: k !debug

    CALL start_clock ('PAW_rho_rad')
    rho_rad(:,:) = 0._dp

    ! prepare the versor
    v(1,1) = x(1)
    v(2,1) = x(2)
    v(3,1) = x(3)
!     WRITE(100,*) 0._dp,0._dp,0._dp
!     WRITE(100,*) v(:,1)
!     WRITE(101,*) v(:,1)
!     WRITE(100,*) 0._dp,0._dp,0._dp

    CALL ylmr2(lmaxq**2, 1, v, nv, ylm)
    rho_rad(:,:) = 0._dp
    spins: DO ispin = 1,nspin
        DO lm = 1, lmaxq**2
            !IF (ABS(ylm(1,lm)) < eps8 ) CONTINUE
            rho_rad(:,ispin) = rho_rad(:,ispin) +&
                    ylm(1,lm)*rho_lm(:,lm,ispin)
        ENDDO ! lm
    ENDDO spins

    CALL stop_clock ('PAW_rho_rad')

END SUBROUTINE PAW_rho_rad2
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! integrate
!!
FUNCTION PAW_sph_integral(n, f1_lm, f2_lm)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : pi,eps8,fpi
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE parameters,             ONLY : nbrx, lqmax
    USE radial_grids,           ONLY : ndmx
    USE atom,                   ONLY : r, rab
    USE uspp,                   ONLY : gen_rndm_r


    REAL(DP), INTENT(IN)        :: f1_lm(ndmx,lmaxq**2,nspin)! Y_lm expansion of rho
    REAL(DP), INTENT(IN)        :: f2_lm(ndmx,lmaxq**2,nspin)! Y_lm expansion of rho
    INTEGER,  INTENT(IN)        :: n ! 4n**2 ~= number of integration directions
    REAL(DP)                    :: PAW_sph_integral

    REAL(DP)                    :: f1_rad(ndmx,nspin)
    REAL(DP)                    :: f2_rad(ndmx,nspin)
    REAL(DP)                    :: aux(ndmx)

    REAL(DP)                    :: th, ph    ! direction on which rho_rad is summed
    REAL(DP)                    :: d         ! surface element
    REAL(DP)                    :: pref      ! 
    INTEGER                     :: i         ! counters on angmom and spin
    INTEGER                     :: ispin, lm ! counters on angmom and spin
    REAL(DP)                    :: integral  ! aux

    INTEGER                     :: nx
    INTEGER                     :: dum
    REAL(DP)                    :: x(3,(lmaxq+1)**2)  !
    REAL(DP)                    :: xx((lmaxq+1)**2)  ! 
    REAL(DP)                    :: ylm((lmaxq+1)**2, (lmaxq+1)**2)  ! 
    REAL(DP)                    :: mly((lmaxq+1)**2, (lmaxq+1)**2)  ! 
    REAL(DP)                    :: w((lmaxq+1)**2)  ! 
    REAL(DP)                    :: o((lmaxq+1)**2)  ! 

    CALL start_clock ('PAW_sph_int')


    ! Angular integration is done on a spiral path on the unitary sphere
    ! (a) use a uniform step for azimuth angle ph
    ! (b) use a variable step so that the surface element is constant
    ! 
    ! it looks like (a) is 10 times more efficient, I don't know why,
    ! tests on different atoms and configurations are necessary
    !
    ! pref is the number of complete loops (times pi) that theta
    ! does while phi goes from 0 to pi. Precision is very sensitive
    ! to it's  choice. While the dependence on sqrt(n) is quite 
    ! clear the extra factor sqrt(pi) is choosen empirically
    ! in order to provide:
    ! 1. as little periodicity as possible 
    !    (n may be a perfect square)
    ! 2. a grid as homogeneous as possible
    !
    ! integration on radial grid has been tryed but proved worst
    ! 
    ! TODO: it may be smart to choose n as a function of max(lm)

    o(:) = 0._dp
    o(1) = sqrt(fpi)
    nx = (lmaxq+1)**2
    ylm(:,:) = 0._dp
    mly(:,:) = 0._dp
    w(:) = 0._dp

!     write(6,*) "lmaxq ::", lmaxq
    CALL gen_rndm_r(nx, x, xx)
    WRITE(6,*) "nx ",1,nx
    CALL ylmr2(nx, nx, x, xx, ylm)
    WRITE(6,*) "nx ",2,nx
!     write(6,*) "ylm ::::::"
!     write(6,"(6f12.6)") ylm
!     write(6,*)
    CALL invmat(nx, ylm, mly, dum)
    WRITE(6,*) "nx ",3,nx
!     write(6,*) "mly ::::::"
!     write(6,"(6f12.6)") mly
!     write(6,*)

!     write(6,*) "1 ::::::"
!     write(6,"(6f12.6)") MATMUL(mly,ylm)
!     write(6,*)


    w = MATMUL(mly,o)
    WRITE(6,*) "nx ",4,nx
!      write(6,*) "w :::::::"
!      WRITE(6,"(16f11.6)") w(:)
!      write(6,*)
!    o = MATMUL(ylm,w)
!      write(6,*) "o :::::::"
!      WRITE(6,"(16f11.6)") o(:)
!      write(6,*)


    PAW_sph_integral = 0._dp
    WRITE(6,*) "nx ",5,nx
    DO i = 1,16 !nx
        !
        CALL PAW_rho_rad2(x(:,i), f1_lm, f1_rad)
        CALL PAW_rho_rad2(x(:,i), f2_lm, f2_rad)
        aux(:) = f1_rad(:,1) * f2_rad(:,1)
        CALL simpson (ndmx,aux,rab(:,1),integral)
        WRITE(6,*) "-->", i,integral,w(i)
        PAW_sph_integral = PAW_sph_integral + integral*w(i)
    ENDDO

    CALL stop_clock ('PAW_sph_int')

END FUNCTION PAW_sph_integral

!***********************
   FUNCTION int2char(i)
   !***********************
      INTEGER, INTENT(in) :: i
      CHARACTER(len=15)   :: int2char
      WRITE(int2char,"(i15)") i
      int2char = ADJUSTL(int2char)
   END FUNCTION int2char


! REMOVE ME:
#include "rad_paw_trash.f90"

END MODULE rad_paw_routines
