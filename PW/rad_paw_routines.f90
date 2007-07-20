!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


MODULE rad_paw_routines
    !
    USE kinds,      ONLY : DP
    !
    IMPLICIT NONE
    PUBLIC
    SAVE
    !
    LOGICAL              :: is_init = .false.
    ! the following variables are used to convert spherical harmonics expansion
    ! to radial sampling, they are initialized for an angular momentum up to
    ! l = max_l and (l+1)**2 = max_lm = nx
    ! see function PAW_rad_init for details
    INTEGER              :: l_max  = 0
    INTEGER              :: lm_max = 0
    INTEGER              :: nx     = 0
    REAL(DP),ALLOCATABLE :: w(:)
    REAL(DP),ALLOCATABLE :: ylm(:,:)
    
CONTAINS
! these has to be modularized too:
#include "../atomic/vxc_t.f90"
#include "../atomic/exc_t.f90"

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

    USE grid_paw_variables,     ONLY : pfunc, ptfunc, tpawp, aerho_atc, psrho_atc
    USE uspp_param,             ONLY : augfun, nhm, lmaxq

    REAL(DP), INTENT(IN)    :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations

    INTEGER, PARAMETER      :: AE = 1, PS = 2   ! All-Electron and Pseudo
    INTEGER                 :: i_what           ! counter on AE and PS
    INTEGER                 :: na,nt,first_nat,last_nat ! atoms counters and indexes

    ! hartree energy scalar fields expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: rho_lm(:,:,:)      ! radial density expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: v_h_lm(:,:)      ! hartree potential, summed on spins
    REAL(DP)                :: e_h(lmaxq**2,nspin)! hartree energy components
    REAL(DP)                :: e                  ! placeholder
    INTEGER                 :: n,lm!debug
    ! xc variables:
    REAL(DP)                :: e_xc(lmaxq**2,nspin)! hartree energy components
    REAL(DP), POINTER       :: rho_core(:,:)      ! pointer to AE/PS core charge density 

    ! BEWARE THAT HARTREE ONLY DEPENDS ON THE TOTAL RHO NOT ON RHOUP AND RHODW SEPARATELY...
    ! TREATMENT OF NSPIN>1 MUST BE CHECKED AND CORRECTED

    ! initialize for integration on angular momentum up to 2*lmaxq (max angmom in atom +1)
    WRITE(6,*) "2 ***************************************************************"
    CALL PAW_rad_init(6*lmaxq) !4*lmaxq) 

    CALL start_clock ('PAW_energy')

    WRITE(6,*) "***************************************************************"

    ! The following operations will be done, first on AE, then on PS part
    ! furthermore they will be done for one atom at a time (to reduce memory usage)
    ! in the future code will be parallelized on atoms:
    !
    ! 1. build rho_lm (PAW_rho_lm)
    ! 2. compute Hartree energy
    !   a. use rho_lm to compute hartree potential (PAW_v_h)
    !   b. use v_h_lm & rho_lm to build hartree energy (PAW_h_energy)
    !   NOTE: a. & b. have been unified in PAW_h_energy
    ! 3. compute XC energy
    !   a. build rho_rad(theta, phi) from rho_lm
    !   b. build v_xc(theta, phi)
    !   c. compute E_xc(theta, phi)
    !   d. repeat a,b,c integrating on theta & phi
    ! 4. free memory (dealloc rho_lm) and return
    !

    ! CHECK: maybe we don't need to alloc/dealloc rho_lm every time
    WRITE(30,"(f20.10)") aerho_atc(:,1)
    WRITE(31,"(f20.10)") psrho_atc(:,1)

    first_nat = 1
    last_nat  = nat
    ! Operations from here on are (will be) parallelized on atoms
    atoms: DO na = first_nat, last_nat
        !
        nt = ityp(na) ! the type of atom na
        ifpaw: IF (tpawp(nt)) THEN
        whattodo: DO i_what = AE, PS
            ! STEP: 1 [ build rho_lm (PAW_rho_lm) ]
            ALLOCATE(rho_lm(ndmx,lmaxq**2,nspin))
            NULLIFY(rho_core)
            IF (i_what == AE) THEN
                ! passing "na" as an argument is dirtyer but faster and
                ! uses less memory than passing only a hyperslice of the array
                CALL PAW_rho_lm(na, becsum, pfunc, rho_lm)
                ! used later for xc energy:
                rho_core => aerho_atc
            ELSE
                CALL PAW_rho_lm(na, becsum, ptfunc, rho_lm, augfun)
                ! used later for xc energy:
                rho_core => psrho_atc
                !     optional argument for pseudo part --> ^^^^^^
            ENDIF
            ! STEP: 2 [ compute Hartree energy ]
            ALLOCATE(v_h_lm(ndmx,lmaxq**2))
            !   2a. use rho_lm to compute hartree potential (PAW_v_h)
            e = PAW_h_energy(na, rho_lm, v_h_lm, e_h)
            WRITE(6,*) "******************************"
            WRITE(6,*) "==PAW RADIAL ENERGY: ", e
            !WRITE(6,'(a,i1,a,f15.7)') ("==RADIAL PAW ENERGY (LM=",lm,"):",e_h(lm,1),lm=1,lmaxq**2)
            !
            ! STEP: 3 [ compute XC energy ]
            
            WRITE(6,*) "== rho core:",MAXVAL(ABS(rho_core(:,:)))
            !
            CALL PAW_rad_init(2*lmaxq) !4*lmaxq) 
            e = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc)
            WRITE(6,*) "******************************"
            WRITE(6,"(a,2i3,f25.15)") "==PAW RADIAL **XC** ENERGY: ", i_what, na, e
            CALL PAW_rad_init(4*lmaxq) !4*lmaxq) 
            e = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc)
            WRITE(6,*) "******************************"
            WRITE(6,"(a,2i3,f25.15)") "==PAW RADIAL **XC** ENERGY: ", i_what, na, e
            CALL PAW_rad_init(6*lmaxq) !4*lmaxq) 
            e = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc)
            WRITE(6,*) "******************************"
            WRITE(6,"(a,2i3,f25.15)") "==PAW RADIAL **XC** ENERGY: ", i_what, na, e



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             CALL start_clock('sph')
            e = PAW_sph_integral(rho_lm, v_h_lm)
!            CALL stop_clock('sph')
            write(6,*) "==radial hartree integral --> ",n,e
            !
            DEALLOCATE(v_h_lm)
            !
            DEALLOCATE(rho_lm)

        ENDDO whattodo
        ENDIF ifpaw
    ENDDO atoms

!     CALL print_clock('sph')

    WRITE(6,*) "***************************************************************"

    CALL stop_clock ('PAW_energy')


END SUBROUTINE PAW_energy

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! initialize several quantities related to radial integration
!!
SUBROUTINE PAW_rad_init(l)
    USE constants,              ONLY : pi, fpi
!    USE uspp,                   ONLY : gen_rndm_r
    INTEGER,INTENT(IN)          :: l     ! max angular momentum

    REAL(DP),ALLOCATABLE        :: x(:,:)   ! nx versors in smart directions
    REAL(DP),ALLOCATABLE        :: xx(:)    ! the norms of the versors (full of ones)

    INTEGER                     :: i,dum    ! counters
    REAL(DP)                    :: th,ph    ! angles
    REAL(DP)                    :: pref     ! workspace
    REAL(DP),ALLOCATABLE        :: mly(:,:) ! inverse of ylm(direction, lm)

    IF( is_init ) THEN
        IF ( l /= l_max ) THEN
!             CALL errore('PAW_rad_init',&
!               'PAW radial integration already initialized but for a different l',&
!               l+100*l_max)
            CALL infomsg('PAW_rad_init',&
              'PAW radial integration already initialized but for a different l: reinitializing.',&
              -l-100*l_max)
              DEALLOCATE(w, ylm)
        ELSE
            ! if already initialized correctly nothing to be done
            RETURN
        ENDIF
    ENDIF
        

    CALL start_clock ('PAW_rad_init')

    ! Angular integration is done on a spiral path on the unitary sphere
    ! using a uniform step for angles ph and th. 
    ! A set of weight is choosen so that:
    ! \sum_{x} w_{x} y_{lm}(\hat{x}) = \sqrt{4\pi} \delta_{lm,00}
    ! than the integral is computed as
    ! \sum_{x} w_{x} \sum_{lm} y_{lm}(\hat{x}) \int_0^{\inf} f(r) dr
    ! notice that the number of x directions must be equal to max{lm}
    ! so nx will be equal to lm_max, I keep to separate variables as 
    ! it is clearer to use one or the other in different contests
    nx = (l+1)**2
    ! the MAGIC factor is empirical: factors <= 1. or integers > 1
    ! may not work, there are no big differences for anything between 1 and 2
#define MAGIC_FACTOR 1.77245385090551602730_dp
    pref = MAGIC_FACTOR*sqrt(DBLE(nx))*pi
    ALLOCATE(w(nx))
    ALLOCATE(ylm(nx,nx), mly(nx,nx))
    ALLOCATE(x(3,nx), xx(nx))

    DO i = 1, nx
        ph = pi * (DBLE(i)-.75_dp) / DBLE(nx)
        !ph = acos(  2._dp*( (DBLE(i)-.25_dp)-DBLE(nx)/2._dp )/DBLE(nx)  )
        th = MOD(pref * (DBLE(i)-1.0_dp) / DBLE(nx), 2._dp*pi )
        x(1,i) = cos(th) * sin(ph)
        x(2,i) = sin(th) * sin(ph)
        x(3,i) = cos(ph)
        xx(i)  = 1._dp
    ENDDO
!    CALL gen_rndm_r(nx,x,xx)

    CALL ylmr2(nx, nx, x, xx, ylm)
    CALL invmat(nx, ylm, mly, dum)
    w(:) = sqrt(fpi)*mly(1,:) 

    ! DEBUG
!     IF (is_init .eq. .false.) THEN
!     DO i = 1,nx
!         WRITE(26,"(10f20.10)") x(1,i), x(2,i), x(3,i), w(i)
!     ENDDO
!     ENDIF

    DEALLOCATE(mly, x, xx)

    ! global variables
    l_max   = l
    lm_max  = nx
    is_init = .true.

    CALL stop_clock ('PAW_rad_init')

END SUBROUTINE PAW_rad_init 

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! use the density produced by sum_rad_rho to compute xc potential and energy, as
!!! xc functional is not diagonal on angular momentum numerical integartion is performed
!!
FUNCTION PAW_xc_energy(na, rho_lm, rho_core, pot_lm, e_lm)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi, e2
    USE parameters,             ONLY : npsx
    USE radial_grids,           ONLY : ndmx
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq
    USE ions_base,              ONLY : ityp
    USE atom,                   ONLY : rgrid

    REAL(DP)                       :: PAW_xc_energy      ! total xc energy
    !
    INTEGER,  INTENT(IN)  :: na                         ! the number of the atom
    REAL(DP), INTENT(IN)  :: rho_lm(ndmx,lmaxq**2,nspin)! charge density as lm components
    REAL(DP), INTENT(IN)  :: rho_core(ndmx,npsx)        ! core charge, radial and spherical
    REAL(DP), INTENT(OUT) :: pot_lm(ndmx,lmaxq**2)      ! out: potential as lm components
    REAL(DP), OPTIONAL,INTENT(OUT) :: e_lm(lmaxq**2)    ! out: energy components 
    !
    ! PLACEHOLDERS yet to be implemented properly:
!    REAL(DP),PARAMETER    :: rho_core = 0._dp ! core density for correction
    INTEGER               :: lsd          ! switch to control local spin density
    !
    REAL(DP)              :: rho_loc(2) = (/0._dp, 0._dp/) 
    REAL(DP)              :: rho_core_loc
                             ! local density (workspace), up and down
    REAL(DP)              :: e            ! workspace
    REAL(DP)              :: e_rad(ndmx)  ! radial energy (to be integrated)
    REAL(DP)              :: rho_rad(ndmx,nspin)
    INTEGER               :: nt, & ! ityp(na)
                             i,k   ! counters on directions and grid 

    CALL start_clock ('PAW_xc_nrg')
    lsd = nspin-1
    nt = ityp(na)

    PAW_xc_energy = 0._dp
    DO i = 1, nx
        !
        CALL PAW_lm2rad(i, rho_lm, rho_rad)
        DO k = 1,ndmx
            ! rho_loc(2) should remain zero if nspin is 1
            rho_loc(1:nspin) = rho_rad(k,1:nspin)/rgrid(nt)%r2(k)
            rho_core_loc = rho_core(k,nt)/rgrid(nt)%r2(k)
            e_rad(k) = exc_t(rho_loc, rho_core_loc, lsd) * SUM(ABS(rho_rad(k,1:nspin)))
        ENDDO
        !
        CALL simpson (rgrid(nt)%mesh,e_rad,rgrid(nt)%rab,e)
        PAW_xc_energy = PAW_xc_energy + e * w(i)
       ! WRITE(6,*) "xc-->",e ,w(i)
    ENDDO

    CALL stop_clock ('PAW_xc_nrg')

END FUNCTION PAW_xc_energy


!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! use the density produced by sum_rad_rho to compute hartree potential 
!!! the potential is then directly integrated to compute hartree energy
!!
FUNCTION PAW_h_energy(na, rho_lm, pot_lm, e_lm)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi
    USE parameters,             ONLY : npsx
    USE radial_grids,           ONLY : ndmx, hartree
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq
    USE ions_base,              ONLY : ityp
    USE atom,                   ONLY : rgrid

    REAL(DP)                       :: PAW_h_energy      ! total hartree energy
    !
    INTEGER,  INTENT(IN)  :: na                         ! the number of the atom
    REAL(DP), INTENT(IN)  :: rho_lm(ndmx,lmaxq**2,nspin)! charge density as lm components
    REAL(DP), INTENT(OUT) :: pot_lm(ndmx,lmaxq**2)      ! out: potential as lm components
    REAL(DP), OPTIONAL,INTENT(OUT) :: e_lm(lmaxq**2)    ! out: energy components 
    !
    REAL(DP)              :: aux(ndmx)   ! workspace
    REAL(DP)              :: par_energy  ! workspace
    REAL(DP)              :: pref        ! workspace

    INTEGER               :: nt,&        ! atom type
                             ispin, &    ! counter on spins
                             lm,l        ! counter on composite angmom lm = l**2 +m
    CALL start_clock ('PAW_h_energy')

    ! get type of atom
    nt = ityp(na)

    ! init total energy and its lm components
    PAW_h_energy = 0._dp
    IF (present(e_lm)) e_lm(:) = 0._dp

    ! this loop computes the hartree potential using the following formula:
    !               l is the first argument in hartree subroutine
    !               r1 = min(r,r'); r2 = MAX(r,r')
    ! V_h(r) = \sum{lm} Y_{lm}(\hat{r})/(2L+1) \int dr' 4\pi r'^2 \rho^{lm}(r') (r1^l/r2^{l+1})
    !     done here --> ^^^^^^^^^^^^^^^^^^^^^           ^^^^^^^^^^^^^^^^^^^^^^ <-- input to the hartree subroutine
    !                 output from the h.s. --> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    DO lm = 1, lmaxq**2
        l = INT(sqrt(DBLE(lm-1)))     ! l has to start from *zero*
            ! this should be the more efficient way to threat nspin>1 cases
            ! with minimal overhead for nspin=1 and no duplicated code
            pref = fpi/(2*l+1)
            aux(:) = pref * rho_lm(:,lm,1)
            DO ispin = 2,nspin ! if nspin < 2 it jumps to ENDDO
                aux(:) = aux(:)+fpi/(2*l+1)*rho_lm(:,lm,ispin)
            ENDDO
            CALL hartree(l, 2*l+2, rgrid(nt)%mesh, rgrid(nt), aux(:), pot_lm(:,lm))
            !
            ! now energy is computed integrating v_h^{lm} * \sum_{spin} rho^{lm}
            ! and summing on lm, aux already contains \sum_{spin} rho^{lm}
            ! but I have to divide (all the array, or just the energy) by 4pi/(2l+1)
            aux(:) = aux(:) * pot_lm(:,lm) 
            CALL simpson (rgrid(nt)%mesh,aux,rgrid(nt)%rab,par_energy)
            !
            PAW_h_energy = PAW_h_energy + par_energy / pref
            IF (present(e_lm)) e_lm(lm) = par_energy / pref
    ENDDO ! lm

    CALL stop_clock ('PAW_h_energy')

END FUNCTION PAW_h_energy

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

    CALL start_clock ('PAW_rho_rad')
    rho_rad(:,:) = 0._dp

    ! prepare the versor
    sin_ph = sin(ph)
    v(1,1) = cos(th) * sin_ph
    v(2,1) = sin(th) * sin_ph
    v(3,1) = cos(ph)

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
! same as PAW_rho_rad but take versor instead of theta and phi as input
SUBROUTINE PAW_rho_rad2(x, rho_lm, rho_rad)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : eps8, pi
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE parameters,             ONLY : nbrx, lqmax
    USE radial_grids,           ONLY : ndmx

    REAL(DP), INTENT(IN)        :: rho_lm(ndmx,lmaxq**2,nspin)! Y_lm expansion of rho
    REAL(DP), INTENT(IN)        :: x(3)
    REAL(DP), INTENT(OUT)       :: rho_rad(ndmx,nspin)        ! charge density on rad. grid

    REAL(DP)                    :: v(3,1)           ! the versor pointed by angles th,ph
    REAL(DP)                    :: sin_ph           ! aux
    REAL(DP),PARAMETER          :: nv(1) = (/1._dp/)! it has to be a vector
    REAL(DP)                    :: ylm(1,lmaxq**2)  ! the spherical harmonics

    INTEGER                     :: ispin, lm ! counters on angmom and spin

    CALL start_clock ('PAW_rho_rad')
    rho_rad(:,:) = 0._dp

    ! prepare the versor (it has to be a 2D matrix)
    v(1,1) = x(1)
    v(2,1) = x(2)
    v(3,1) = x(3)

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
!------------------------------------------------------------
! same as PAW_rho_rad but take directly the spherical harmonics from
! global module's variables
SUBROUTINE PAW_lm2rad(ix, rho_lm, rho_rad)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : eps8, pi
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE radial_grids,           ONLY : ndmx

    REAL(DP), INTENT(IN)        :: rho_lm(ndmx,lmaxq**2,nspin)! Y_lm expansion of rho
    INTEGER                     :: ix ! line of the ylm matrix to use
                                      ! actually it is one of the lm_max directions
    REAL(DP), INTENT(OUT)       :: rho_rad(ndmx,nspin)        ! charge density on rad. grid

    INTEGER                     :: ispin, lm ! counters on angmom and spin

    CALL start_clock ('PAW_lm2rad')
    rho_rad(:,:) = 0._dp

    rho_rad(:,:) = 0._dp
    ! cycling on spin is a bit less general...
    spins: DO ispin = 1,nspin
        DO lm = 1, lmaxq**2 ! 
            !IF (ABS(ylm(1,lm)) < eps8 ) CONTINUE
            rho_rad(:,ispin) = rho_rad(:,ispin) +&
                    ylm(ix,lm)*rho_lm(:,lm,ispin)
        ENDDO ! lm
    ENDDO spins

    CALL stop_clock ('PAW_lm2rad')

END SUBROUTINE PAW_lm2rad

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! integrate
!!
FUNCTION PAW_sph_integral(f1_lm, f2_lm)
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
    REAL(DP)                    :: PAW_sph_integral

    REAL(DP)                    :: f1_rad(ndmx,nspin)
    REAL(DP)                    :: f2_rad(ndmx,nspin)
    REAL(DP)                    :: aux(ndmx)

    INTEGER                     :: i         ! counters on angmom and spin
    REAL(DP)                    :: integral  ! aux

    INTEGER                     :: dum

    CALL start_clock ('PAW_sph_int')
    !

    PAW_sph_integral = 0._dp
    DO i = 1, nx
        !
        CALL PAW_lm2rad(i, f1_lm, f1_rad)
        CALL PAW_lm2rad(i, f2_lm, f2_rad)
        aux(:) = f1_rad(:,1) * f2_rad(:,1)
        !
        CALL simpson (ndmx,aux,rab(:,1),integral)
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
