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
    ! if set to true (in init) we have to do gradient correction
    LOGICAL              :: do_gcxc = .false.

    ! the following variables are used to convert spherical harmonics expansion
    ! to radial sampling, they are initialized for an angular momentum up to
    ! l = max_l and (l+1)**2 = max_lm = nx
    ! see function PAW_rad_init for details
    INTEGER              :: l_max  = 0
    INTEGER              :: lm_max = 0
    INTEGER              :: nx     = 0
    REAL(DP),ALLOCATABLE :: ww(:)
    ! additional variables for gradient correction
    REAL(DP),PARAMETER   :: delta = 1.e-5_dp ! for numerical derivatives
    REAL(DP),ALLOCATABLE :: ylm(:,:),&
                            dylm2(:,:),&!  |grad(ylm)|**2
                            dylmt(:,:),&! |d(ylm)/dtheta|**2
                            dylmp(:,:),&! |d(ylm)/dphi|**2
                            ylm_dth(:,:,:),&! spherical harmonics r + \delta \hat{\theta}
                            ylm_dph(:,:,:)  ! as above, but s/theta/phi/

CONTAINS
! these has to be modularized too:
#include "../atomic/vxc_t.f90"
#include "../atomic/exc_t.f90"
#include "../atomic/vxcgc.f90"

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
    REAL(DP)                :: e,e1,e2,e3         ! placeholder
    INTEGER                 :: n,lm!debug
    ! xc variables:
    REAL(DP)                :: e_xc(lmaxq**2,nspin)! hartree energy components
    REAL(DP), POINTER       :: rho_core(:,:)      ! pointer to AE/PS core charge density 

    ! BEWARE THAT HARTREE ONLY DEPENDS ON THE TOTAL RHO NOT ON RHOUP AND RHODW SEPARATELY...
    ! TREATMENT OF NSPIN>1 MUST BE CHECKED AND CORRECTED

    CALL start_clock ('PAW_energy')

    ! initialize for integration on angular momentum and gradient
    CALL PAW_rad_init(5*lmaxq)

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
                !     optional argument for pseudo part --> ^^^^^^
                ! used later for xc energy:
                rho_core => psrho_atc
            ENDIF
            ! STEP: 2 [ compute Hartree energy ]
            ALLOCATE(v_h_lm(ndmx,lmaxq**2))
            !   2a. use rho_lm to compute hartree potential (PAW_v_h)
            e = PAW_h_energy(na, rho_lm, v_h_lm, e_h)
            WRITE(6,*) "******************************"
            WRITE(6,*) "==PAW Hartree E: ", e
            !WRITE(6,'(a,i1,a,f15.7)') ("==RADIAL PAW ENERGY (LM=",lm,"):",e_h(lm,1),lm=1,lmaxq**2)

            ! STEP: 3 [ compute XC energy ]
            WRITE(6,*) "******************************"
            e1 = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc,1)
            e2 = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc,2)
            e3 = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc,3)
            WRITE(6,"(a,2i3,3f25.15)") "==XCa: ", i_what, na, e1,e2,e3
            e1 = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc,-1)
            e2 = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc,-2)
            e3 = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc,-3)
            WRITE(6,"(a,2i3,3f25.15)") "==XCb: ", i_what, na, e1,e2,e3
            e1 = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc,8)
            e2 = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc,9)
            e3 = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc,10)
            WRITE(6,"(a,2i3,3f25.15)") "==vxcgc", i_what, na, e1,e2,e3


            ! check if integration is working
            e = PAW_sph_integral(rho_lm, v_h_lm)
            write(6,*) "==radial hartree integral --> ",e
            !
            DEALLOCATE(v_h_lm)
            !
            DEALLOCATE(rho_lm)

        ENDDO whattodo
        ENDIF ifpaw
    ENDDO atoms

    CALL stop_clock ('PAW_energy')


END SUBROUTINE PAW_energy

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! initialize several quantities related to radial integration: spherical harmonics and their 
!!! gradients along several (depending on lmaxq) directions, weights for spherical integration
!!
SUBROUTINE PAW_rad_init(l)
    USE constants,              ONLY : pi, fpi, eps8
    USE funct,                  ONLY : igcx, igcc
    INTEGER,INTENT(IN)          :: l            ! max angular momentum

    REAL(DP),ALLOCATABLE        :: x(:),&       ! nx versors in smart directions
                                   w(:),&       ! temporary integration weights
                                   r(:,:),&     ! integration directions
                                   r2(:),&      ! square modulus of r
                                   ath(:),aph(:)! angles in sph coords for r

    INTEGER                     :: i,ii,n,dum   ! counters
    INTEGER                     :: lm,m         ! indexes for ang.mom
    REAL(DP)                    :: phi,dphi,rho ! spherical coordinates
    REAL(DP)                    :: z            ! cartesian coordinates
    REAL(DP)                    :: int_lm       ! workspace
    ! for gradient corrections:
    INTEGER                     :: ipol
    REAL(DP),ALLOCATABLE        :: aux(:,:),&   ! workspace
                                   s(:,:),&     ! integration directions + delta
                                   s2(:)        ! square modulus of s
    REAL(DP)                    :: vth(3), vph(3) !versors for theta and phi

    ! reinit if necessary
    IF( is_init ) THEN
        IF ( l /= l_max ) THEN
            CALL infomsg('PAW_rad_init',&
              'PAW radial integration already initialized but for a different l: reinitializing.',&
              -l-100*l_max)
              DEALLOCATE(ww, ylm)
              IF (allocated(dylm2)) DEALLOCATE(dylm2)
        ELSE
            ! if already initialized correctly nothing to be done
            RETURN
        ENDIF
    ENDIF

    CALL start_clock ('PAW_rad_init')

    ! maximum value of l correctly integrated
    l_max = l
    ! volume element for angle phi
    dphi = 2.d0*pi/(l_max+1)
    ! number of samples for theta angle
    n = (l_max+2)/2
    ALLOCATE (x(n),w(n))
    ! compute weights for theta integration
    CALL weights(x,w,n)

    ! number of integration directions
    nx = n*(l_max+1)
    ALLOCATE (r(3,nx),r2(nx), ww(nx), ath(nx), aph(nx))

    ! compute real weights multiplying theta and phi weights
    ii = 0
    do i=1,n
        z = x(i)
        rho=sqrt(1.d0-z**2)
        do m=0,l_max
            ii= ii+1
            phi = dphi*m
            r(1,ii) = rho*cos(phi)
            r(2,ii) = rho*sin(phi)
            r(3,ii) = z
            ww(ii) = w(i)*2._dp*pi/(l_max+1)
            r2(ii) = r(1,ii)**2+r(2,ii)**2+r(3,ii)**2
            ! these will be used later:
            ath(ii) = acos(z/sqrt(r2(ii)))
            aph(ii) = phi
        end do
    end do
    ! cleanup
    DEALLOCATE (x,w)

    ! initialize spherical harmonics that will be used
    ! to convert rho_lm to radial grid
    lm_max = (l_max+1)**2
    ALLOCATE(ylm(nx,lm_max))
    CALL ylmr2(lm_max, nx, r,r2,ylm)

    ! if gradient corrections will be used than we need
    ! to initialize the gradient of ylm, as we are working in spherical
    ! coordinates the formula involves \hat{theta} and \hat{phi} 
    gradient: IF (igcx>0 .or. igcc>0) THEN
        do_gcxc = .true.
        ALLOCATE (s(3,nx),s2(nx))
        ALLOCATE(dylm2(nx,lm_max),dylmt(nx,lm_max),dylmp(nx,lm_max),aux(nx,lm_max))
        ! last dimension is +/- delta on theta/phi direction
        ALLOCATE(ylm_dth(nx,lm_max,2), ylm_dph(nx,lm_max,2))
        dylm2(:,:)  = 0._dp
        dylmt(:,:) = 0._dp
        dylmp(:,:) = 0._dp
        ylm_dth(:,:,:) = 0._dp
        ylm_dph(:,:,:) = 0._dp

        ! compute ylm at r +/- delta * vth/vph, code needs massive cleanup
        ! + theta
        DO i = 1,nx
            vph = (/-sin(aph(i)), cos(aph(i)), 0._dp/)
            ! this is the explicit form, but the cross product trick (below) is much faster:
            ! vth = (/cos(aph(i))*cos(ath(i)), sin(aph(i))*cos(ath(i)), -sin(ath(i))/)
            vth = (/vph(2)*r(3,i)-vph(3)*r(2,i), vph(3)*r(1,i)-vph(1)*r(3,i), vph(1)*r(2,i)-vph(2)*r(1,i)/)
            s(:,i) = r(:,i) + delta * vth(:)
            s2(i) = s(1,i)**2 + s(2,i)**2 + s(3,i)**2
        ENDDO
        CALL ylmr2(lm_max, nx, s,s2,ylm_dth(:,:,2))
        ! - theta
        DO i = 1,nx
            vph = (/-sin(aph(i)), cos(aph(i)), 0._dp/)
            vth = (/vph(2)*r(3,i)-vph(3)*r(2,i), vph(3)*r(1,i)-vph(1)*r(3,i), vph(1)*r(2,i)-vph(2)*r(1,i)/)
            s(:,i) = r(:,i) - delta * vth(:)
            s2(i) = s(1,i)**2 + s(2,i)**2 + s(3,i)**2
        ENDDO
        CALL ylmr2(lm_max, nx, s,s2,ylm_dth(:,:,1)) 
        ! + phi
        DO i = 1,nx
            vph = (/-sin(aph(i)), cos(aph(i)), 0._dp/)
            vth = (/vph(2)*r(3,i)-vph(3)*r(2,i), vph(3)*r(1,i)-vph(1)*r(3,i), vph(1)*r(2,i)-vph(2)*r(1,i)/)
            s(:,i) = r(:,i) + delta * vph(:)
            s2(i) = s(1,i)**2 + s(2,i)**2 + s(3,i)**2
        ENDDO
        CALL ylmr2(lm_max, nx, s,s2,ylm_dph(:,:,2)) 
            ! - phi
        DO i = 1,nx
            vph = (/-sin(aph(i)), cos(aph(i)), 0._dp/)
            vth = (/vph(2)*r(3,i)-vph(3)*r(2,i), vph(3)*r(1,i)-vph(1)*r(3,i), vph(1)*r(2,i)-vph(2)*r(1,i)/)
            s(:,i) = r(:,i) - delta * vph(:)
            s2(i) = s(1,i)**2 + s(2,i)**2 + s(3,i)**2
        ENDDO
        CALL ylmr2(lm_max, nx, s,s2,ylm_dph(:,:,1)) 

        ! compute derivative along x, y and z => gradient, then compute the
        ! scalar products with \hat{theta} and \hat{phi} and store them in
        ! dylmt and dylmp respectively
        DO ipol = 1,3 !x,y,z
            CALL dylmr2(lm_max, nx, r,r2, aux, ipol)
            DO lm = 1, lm_max
            DO i = 1,nx
                vph = (/-sin(aph(i)), cos(aph(i)), 0._dp/)
                ! this is the explicit form, but the cross product trick (below) is much faster:
                ! vth = (/cos(aph(i))*cos(ath(i)), sin(aph(i))*cos(ath(i)), -sin(ath(i))/)
                vth = (/vph(2)*r(3,i)-vph(3)*r(2,i), vph(3)*r(1,i)-vph(1)*r(3,i), vph(1)*r(2,i)-vph(2)*r(1,i)/)
                dylm2(i,lm) = dylm2(i,lm) + ABS(aux(i,lm)) ** 2
                dylmt(i,lm) = dylmt(i,lm) + aux(i,lm)*vth(ipol)
                dylmp(i,lm) = dylmp(i,lm) + aux(i,lm)*vph(ipol)/cos(ath(i))
            ENDDO
            ENDDO
        ENDDO
    DEALLOCATE(aux)
    ENDIF gradient

    ! cleanup
    DEALLOCATE (r,r2)

    ! success
    is_init = .true.

    CALL stop_clock ('PAW_rad_init')

CONTAINS
    ! Computes weights for gaussian integrals,
    ! from numerical recipes
    SUBROUTINE weights(x,w,n)
    implicit none
    integer :: n, i,j,m
    real(8), parameter :: eps=1.d-14
    real(8) :: x(n),w(n), z,z1, p1,p2,p3,pp,pi
    
    pi = 4.d0*atan(1.d0)
    m=(n+1)/2
    do i=1,m
        z1 = 2.d0
        z=cos(pi*(i-0.25d0)/(n+0.5d0))
        do while (abs(z-z1).gt.eps)
        p1=1.d0
        p2=0.d0
        do j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
        end do
        pp = n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp
        end do
        x(i) = -z
        x(n+1-i) = z
        w(i) = 2.d0/((1.d0-z*z)*pp*pp)
        w(n+1-i) = w(i)
    end do

    END SUBROUTINE weights
END SUBROUTINE PAW_rad_init 
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! use the density produced by sum_rad_rho to compute xc potential and energy, as
!!! xc functional is not diagonal on angular momentum numerical integartion is performed
FUNCTION PAW_xc_energy(na, rho_lm, rho_core, pot_lm, e_lm, task)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi
    USE parameters,             ONLY : npsx
    USE radial_grids,           ONLY : ndmx
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : lmaxq
    USE ions_base,              ONLY : ityp
    USE atom,                   ONLY : rgrid
    ! for gradient correction:
    USE funct,                  ONLY : igcx, igcc

    REAL(DP)              :: PAW_xc_energy      ! total xc energy
    !
    INTEGER,  INTENT(IN)  :: na                         ! the number of the atom
    INTEGER,  INTENT(IN)  :: task!remove me
    REAL(DP), INTENT(IN)  :: rho_lm(ndmx,lmaxq**2,nspin)! charge density as lm components
    REAL(DP), INTENT(IN)  :: rho_core(ndmx,npsx)        ! core charge, radial and spherical
    ! TODO:
    REAL(DP), INTENT(OUT) :: pot_lm(ndmx,lmaxq**2)      ! out: potential as lm components
    REAL(DP), OPTIONAL,INTENT(OUT) :: e_lm(lmaxq**2)    ! out: energy components 
    !
    INTEGER               :: lsd          ! switch to control local spin density
    !
    REAL(DP)              :: rho_loc(2) = (/0._dp, 0._dp/) 
                             ! local density (workspace), up and down
    REAL(DP)              :: e, e_aux     ! workspace
    REAL(DP)              :: e_rad(ndmx)  ! radial energy (to be integrated)
    REAL(DP)              :: rho_rad(ndmx,nspin) ! workspace (radial slice of rho)
    INTEGER               :: nt, &        ! ityp(na)
                             ix,k          ! counters on directions and radial grid
    ! for gradient correction
    REAL(DP),ALLOCATABLE  :: grho_rad(:,:)! workspace (radial slice of grad(rho))
    REAL(DP)              :: grho_loc(2) = (/0._dp, 0._dp/) !I can afford to waste 16 bytes
    ! for STUPID gradient correction
    REAL(DP)              :: frr(ndmx,2), frc(ndmx), vgc(ndmx), egc(ndmx) !fake rho and rho_core
    REAL(DP)              :: fr(ndmx), fr2(ndmx) !fake rho and rho_core
    INTEGER               :: fns !fake spin number

    CALL start_clock ('PAW_xc_nrg')
    lsd = nspin-1
    nt = ityp(na)

    ! init for gradient correction
    IF (do_gcxc) ALLOCATE(grho_rad(ndmx,nspin))
#define RADIAL_GRADIENT

    PAW_xc_energy = 0._dp
    DO ix = 1, nx
        ! LDA (and LSDA) part (no gradient correction):
        CALL PAW_lm2rad(ix, rho_lm, rho_rad)
        !
        DO k = 1,rgrid(nt)%mesh
            rho_loc(1:nspin) = rho_rad(k,1:nspin)/rgrid(nt)%r2(k)
            !
            e_rad(k) = exc_t(rho_loc, rho_core(k,nt), lsd)&
                     * (SUM(rho_rad(k,1:nspin))+rho_core(k,nt)*rgrid(nt)%r2(k))
        ENDDO
        gradient_correction:& ! add it!
        IF (do_gcxc) THEN
            IF ( task == 1) THEN
                CALL PAW_ngrad(na, ix, rho_lm, rho_rad, rho_core, grho_rad)
                CALL PAW_gcxc(na, rho_rad,rho_core,grho_rad, e_rad,task)
            ELSE IF (task == 2) THEN
                e_rad(:) = 0._dp
                CALL PAW_ngrad(na, ix, rho_lm, rho_rad, rho_core, grho_rad)
                CALL PAW_gcxc(na, rho_rad,rho_core,grho_rad, e_rad,task)
            ELSE IF (task == 3) THEN
                ! nothing to do
            ELSE IF (task == -1) THEN
                CALL PAW_grad(na, ix, rho_lm, rho_rad, rho_core, grho_rad)
                CALL PAW_gcxc(na, rho_rad,rho_core,grho_rad, e_rad,task)
            ELSE IF (task == -2) THEN
                e_rad(:) = 0._dp
                CALL PAW_grad(na, ix, rho_lm, rho_rad, rho_core, grho_rad)
                CALL PAW_gcxc(na, rho_rad,rho_core,grho_rad, e_rad,task)
            ELSE IF (task == -3) THEN
                ! nothing to do
            ELSE IF (task > 7) THEN
                fns = nspin
                DO k = 1,rgrid(nt)%mesh
                    frr(k,1) = fpi*rho_rad(k,1)
                    frr(k,2) = 0._dp
                    frc(k)   = fpi*rgrid(nt)%r2(k)*rho_core(k,nt)
                    fr(k)    = rgrid(nt)%r(k)
                    fr2(k)   = rgrid(nt)%r2(k)
                ENDDO
                write(6,*) "+entering vxcgc"
!              real(DP) :: r(mesh), r2(mesh), rho(ndm,2), rhoc(ndm), &
!                    vxcgc(ndm,mesh, nspin,r,r2,  rho,rhoc, vgc,egc)
                CALL vxcgc(ndmx, ndmx, fns,fr,fr2,frr, frc, vgc,egc)
                write(6,*) "+exiting vxcgc"
                IF (task == 9) e_rad(:) = e_rad(:) + egc(:)
                IF (task == 8) e_rad(:) = egc(:)
                IF (task == 10) THEN
                    e_rad(:) = 0._dp
                    CALL PAW_ngrad(na, ix, rho_lm, rho_rad, rho_core, grho_rad)
                    CALL PAW_gcxc(na, rho_rad,rho_core,grho_rad, e_rad,task)
                    DO k = 1,rgrid(nt)%mesh
                        write(510,"(90f20.12)") rgrid(nt)%r(k), egc(k), e_rad(k), egc(k)-e_rad(k)
                    ENDDO
                    e_rad(:) = egc(:) - e_rad(:)
                ENDIF
                

            ENDIF
        ENDIF gradient_correction
        !
        ! integrate radial slice of xc energy:
        CALL simpson (rgrid(nt)%mesh,e_rad,rgrid(nt)%rab,e)
        ! integrate on sph. surface     v----------------^
        PAW_xc_energy = PAW_xc_energy + e * ww(ix)
    ENDDO

    IF (allocated(grho_rad)) DEALLOCATE(grho_rad)

    CALL stop_clock ('PAW_xc_nrg')

END FUNCTION PAW_xc_energy

! add gradient correction, code adapted from ../atomic/vxcgc.f90
SUBROUTINE PAW_gcxc(na, rho,core,grho,e,task)

USE kinds,                  ONLY : DP
USE ions_base,              ONLY : ityp
USE radial_grids,           ONLY : ndmx
USE lsda_mod,               ONLY : nspin
USE atom,                   ONLY : g => rgrid
USE parameters,             ONLY : npsx
USE constants,              ONLY : fpi,pi,e2

INTEGER, INTENT(IN)    :: na              ! atom index
INTEGER, INTENT(IN)    :: task !DEBUG
REAL(DP), INTENT(IN)   :: rho(ndmx,nspin) ! radial density,
REAL(DP), INTENT(IN)   :: grho(ndmx,nspin)! square gradient of rho
REAL(DP), INTENT(IN)   :: core(ndmx,npsx) ! spherical core density
REAL(DP), INTENT(INOUT):: e(ndmx)         ! radial local xc energy
!
INTEGER                :: i               ! counter for mesh
INTEGER                :: nt              ! atom type
! workspaces:
REAL(DP)               :: arho, sgn
REAL(DP)               :: sx,sc,v1x,v2x,v1c,v2c
REAL(DP),PARAMETER     :: eps = 1.e-12_dp

nt = ityp(na)

IF (nspin.eq.1) THEN
    DO i=1,g(nt)%mesh
        arho = rho(i,1)/g(nt)%r2(i) + core(i,ityp(na))
        sgn  = sign(1.0_dp,arho)
        arho = abs(arho)
        IF (arho.gt.eps.and.abs(grho(i,1)).gt.eps) THEN
            CALL gcxc(arho,grho(i,1),sx,sc,v1x,v2x,v1c,v2c)
            e(i) = e(i) + sgn *e2 *(sx+sc)*g(nt)%r2(i) !&
                    !*(SUM(rho(i,1:nspin))+core(i,nt)*g(nt)%r2(i))
        ELSE IF (i.gt.g(nt)%mesh/2) THEN
            ! asymptotic formula
            e(i) = e(i)-0.0_dp/(2.0_dp*g(nt)%r(i))!&
                !*(SUM(rho(i,1:nspin))+core(i,nt)*g(nt)%r2(i))
        ENDIF
    write(500+task/2,"(90f20.12)") g(nt)%r(i), grho(i,1), e(i), arho,sx,sc,v1x,v2x,v1c,v2c
    ENDDO
ENDIF

END SUBROUTINE PAW_gcxc

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! build *gradient* in a (not so much, after all) ignorant and brutal numerical way
SUBROUTINE PAW_ngrad(na, ix, rho_lm, rho_rad, rho_core, grho_rad)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE radial_grids,           ONLY : ndmx
    USE ions_base,              ONLY : ityp
    USE parameters,             ONLY : npsx
    USE atom,                   ONLY : g => rgrid

    REAL(DP), INTENT(IN)        :: rho_lm(ndmx,lmaxq**2,nspin)! Y_lm expansion of rho
    REAL(DP), INTENT(IN)        :: rho_rad(ndmx,nspin)        ! radial density along direction ix
    REAL(DP), INTENT(IN)        :: rho_core(ndmx,npsx)             ! core density
    INTEGER, INTENT(IN)         :: ix ! line of the dylm2 matrix to use
                                      ! actually it is one of the nx directions
    INTEGER, INTENT(IN)         :: na ! atom index
    REAL(DP), INTENT(OUT)       :: grho_rad(ndmx,nspin)   ! grad of charge density on rad. grid
    !
    REAL(DP)                    :: aux(ndmx,nspin),aux2(ndmx),aux3(ndmx,2) ! workspace
    REAL(DP)                    :: r(ndmx)         ! placeholder (for lazyer programming)
    INTEGER                     :: i, is, lm, nt,k,dir ! counters on angmom and spin

    CALL start_clock ('PAW_grad')
    nt = ityp(na)
    r(:) = g(nt)%r(:) !FIXME: remove this when code works

    ! The problem in computing this gradient is that the radial functions
    ! are sampled on a logarithmic grid, so it is not possible to use simple 
    ! finite difference formulas, furthermore deriving along x,y and z is
    ! almost impossible. I will have to derive along r, theta and phi instead;
    ! it has the additional advantage that I can (actually I have to) consider
    ! lm components of rho constants for (small) displacement along theta and phi,
    ! spherical harmonics are constant for displacements along r

    ! 1. rho_lm has to be divided by r**2
    DO is = 1,nspin
    ! rho/r*2 + rho_core/nspin = real density
    aux(1:g(nt)%mesh,is) = rho_rad(1:g(nt)%mesh,is)/g(nt)%r2(1:g(nt)%mesh) &
                         + rho_core(1:g(nt)%mesh,nt)/nspin
    ENDDO

    ! 2. compute the partial derivative of rho_lm
    ! 3. \sum \partial rho(r) \over \partial r Y_{lm}(th,ph)
    ! 4. compute the square *after* summing the components
    grho_rad(:,:) = 0._dp
    DO is = 1,nspin
        ! numerical derivative by ADC ../atomic/vxcgc.f90
        DO k  = 2,g(nt)%mesh-1
            aux2(k) = (  (r(k+1)-r(k))**2 * (aux(k-1,is)-aux(k,is))  &
                       -(r(k-1)-r(k))**2 * (aux(k+1,is)-aux(k,is))  )&
                    / ( (r(k+1)-r(k)) * (r(k-1)-r(k)) * (r(k+1)-r(k-1)) )
        ENDDO
        ! extremes:
        aux2(g(nt)%mesh)=0.0_dp
        aux2(1)=aux2(2)+(aux2(3)-aux2(2))*(r(1)-r(2))/(r(3)-r(2))
        ! compute the square
        DO k  = 1,g(nt)%mesh
            grho_rad(k,is) = aux2(k)**2
        ENDDO
    ENDDO

#ifdef RADIAL_GRADIENT
    ! now I have to compute numerical derivatives along theta and phi
    ! directions, I will use ylm_dth and ylm_dph that should be already 
    ! initialized
    DO is = 1,nspin
        aux3(:,:) = 0._dp
        DO dir = 1,2 ! positive or negative difference along theta
            ! derivation along theta
            ! core charge has to be added to the spherical component
            aux3(:,dir) = aux3(:,dir) +ylm_dth(ix,1,dir)&
                        *(rho_lm(:,1,is)/g(nt)%r2(:)+rho_core(:,nt)/nspin)
            DO lm = 2, lmaxq**2 ! 
                aux3(:,dir) = aux3(:,dir) + ylm_dth(ix,lm,dir)&
                            * rho_lm(:,lm,is)/g(nt)%r2(:)
            ENDDO ! lm
        ENDDO ! dir
        ! compute the differences and sum them to grho
        DO k  = 1,g(nt)%mesh
            grho_rad(k,is) = grho_rad(k,is)&
                           + ( (aux3(k,1) - aux3(k,2))/ 2/delta/g(nt)%r(k) )**2
        ENDDO
    ENDDO ! spin

    ! the same for phi (yes, this code is duplicated)...
    DO is = 1,nspin
        aux3(:,:) = 0._dp
        DO dir = 1,2 ! positive or negative difference along theta
            aux3(:,dir) = aux3(:,dir) +ylm_dph(ix,1,dir)&
                        *(rho_lm(:,1,is)/g(nt)%r2(:)+rho_core(:,nt)/nspin)
            DO lm = 2, lmaxq**2 ! 
                aux3(:,dir) = aux3(:,dir) + ylm_dph(ix,lm,dir)&
                            * rho_lm(:,lm,is)/g(nt)%r2(:)
            ENDDO ! lm
        ENDDO ! dir
        DO k  = 1,g(nt)%mesh
            grho_rad(k,is) = grho_rad(k,is)&
                           + ( (aux3(k,1) - aux3(k,2))/ 2/delta/g(nt)%r(k) )**2
        ENDDO
    ENDDO ! spin
#endif

    CALL stop_clock ('PAW_grad')

END SUBROUTINE PAW_ngrad
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! build *gradient* of radial charge distribution from its spherical harmonics expansion
SUBROUTINE PAW_grad(na, ix, rho_lm, rho_rad, rho_core, grho_rad)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE radial_grids,           ONLY : ndmx
    USE ions_base,              ONLY : ityp
    USE parameters,             ONLY : npsx
    USE atom,                   ONLY : g => rgrid

    REAL(DP), INTENT(IN)        :: rho_lm(ndmx,lmaxq**2,nspin)! Y_lm expansion of rho
    REAL(DP), INTENT(IN)        :: rho_rad(ndmx,nspin)        ! radial density along direction ix
    REAL(DP), INTENT(IN)        :: rho_core(ndmx,npsx)             ! core density
    INTEGER, INTENT(IN)         :: ix ! line of the dylm2 matrix to use
                                      ! actually it is one of the nx directions
    INTEGER, INTENT(IN)         :: na ! atom index
    REAL(DP), INTENT(OUT)       :: grho_rad(ndmx,nspin)   ! grad of charge density on rad. grid
    !
    REAL(DP)                    :: aux(ndmx),aux2(ndmx),aux_lm(ndmx,lmaxq**2,nspin) ! workspace
    REAL(DP)                    :: r(ndmx)         ! placeholder (for easyer programming)
    INTEGER                     :: i, is, lm, nt,k   ! counters on angmom and spin

    CALL start_clock ('PAW_grad')
    nt = ityp(na)
    r(:) = g(nt)%r(:) !FIXME: remove this when code works
    ! TODO: include core correction charge


    ! from here on \sum => \sum_{l=0,l_max}\sum_{ m=-l,l}
    !                   => \sum_{(lm) = 1,l_max**2}

    ! 1. build real charge density = rho/r**2 + rho_core
    ! 2. compute the partial derivative of rho_lm
    ! 3. \sum \partial rho(r) \over \partial r Y_{lm}(th,ph)
    ! 4. compute the square *after* summing the components
    !
    grho_rad(:,:) = 0._dp
    DO is = 1,nspin
        ! build real charge density
        aux(1:g(nt)%mesh) = rho_rad(1:g(nt)%mesh,is)/g(nt)%r2(1:g(nt)%mesh) &
                          + rho_core(1:g(nt)%mesh,nt)/nspin

        ! numerical derivative by ADC ../atomic/vxcgc.f90
        DO k  = 2,g(nt)%mesh-1
            aux2(k) = (  (r(k+1)-r(k))**2 * (aux(k-1)-aux(k))  &
                        -(r(k-1)-r(k))**2 * (aux(k+1)-aux(k))  )&
                    / ( (r(k+1)-r(k)) * (r(k-1)-r(k)) * (r(k+1)-r(k-1)) )
        ENDDO
        ! extremes:
        aux2(g(nt)%mesh)=0.0_dp
        aux2(1)=aux2(2)+(aux2(3)-aux2(2))*(r(1)-r(2))/(r(3)-r(2))
        ! compute the square
        DO k  = 1,g(nt)%mesh
            grho_rad(k,is) = aux2(k)**2
        ENDDO
    ENDDO

#ifdef RADIAL_GRADIENT
    ! 1. rho_lm has to be divided by r**2
    DO is = 1,nspin
        ! core density has to be added to the spherical component (lm=1)
        ! FIXME: probabily rho_core has to be multiplied/divided by fpi or sqrt(fpi)
        !        depending on Y_00 normalization
        aux_lm(1:g(nt)%mesh,1,is) = rho_lm(1:g(nt)%mesh,1,is) / g(nt)%r2(1:g(nt)%mesh) &
                                + rho_core(1:g(nt)%mesh,nt)/nspin
        DO lm = 2, lmaxq**2
            aux_lm(1:g(nt)%mesh,lm,is) = rho_lm(1:g(nt)%mesh,lm,is) / g(nt)%r2(1:g(nt)%mesh)
        ENDDO
    ENDDO

    ! 5. [ \sum rho(r) (dY_{lm}/dphi /cos(theta))  ]**2
    ! 6. [ \sum rho(r) (dY_{lm}/dtheta)  ]**2
    aux(:)  = 0._dp
    aux2(:) = 0._dp
    DO is = 1,nspin
    DO lm = 1,lmaxq**2
        ! 5:
        aux(:) = aux(1:g(nt)%mesh) + dylmp(ix,lm)* (aux_lm(1:g(nt)%mesh,lm,is))
        ! 6:
        aux2(:) = aux2(1:g(nt)%mesh) + dylmt(ix,lm)* (aux_lm(1:g(nt)%mesh,lm,is))
    ENDDO
    !
    grho_rad(1:g(nt)%mesh,is) = grho_rad(1:g(nt)%mesh,is)&
                              + (aux(1:g(nt)%mesh)**2 + aux2(1:g(nt)%mesh)**2 )&
                                / g(nt)%r2(1:g(nt)%mesh)
    ENDDO
#endif

    CALL stop_clock ('PAW_grad')

END SUBROUTINE PAW_grad

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! use the density produced by sum_rad_rho to compute hartree potential 
!!! the potential is then directly integrated to compute hartree energy
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
            ! but I have to redivide by 4pi/(2l+1)
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
!!! sum up pfuncs x occupation to build radial density's angular momentum components
SUBROUTINE PAW_rho_lm(na, becsum, pfunc, rho_lm, augfun)
    USE kinds,                  ONLY : DP
    USE atom,                   ONLY : r, rab, mesh, msh
    USE ions_base,              ONLY : ityp, ntyp => nsp, nat 
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq!, augfun
    USE uspp,                   ONLY : indv, ap, nhtolm,lpl,lpx
    USE parameters,             ONLY : nbrx, lqmax
    USE radial_grids,           ONLY : ndmx

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
SUBROUTINE PAW_lm2rad(ix, rho_lm, rho_rad)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : eps8, pi
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE radial_grids,           ONLY : ndmx

    REAL(DP), INTENT(IN)        :: rho_lm(ndmx,lmaxq**2,nspin)! Y_lm expansion of rho
    INTEGER                     :: ix ! line of the ylm matrix to use
                                      ! actually it is one of the nx directions
    REAL(DP), INTENT(OUT)       :: rho_rad(ndmx,nspin) ! charge density on rad. grid

    INTEGER                     :: ispin, lm ! counters on angmom and spin

    CALL start_clock ('PAW_lm2rad')
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
!!! example of rdial integration
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
        PAW_sph_integral = PAW_sph_integral + integral*ww(i)
    ENDDO

    CALL stop_clock ('PAW_sph_int')

END FUNCTION PAW_sph_integral

! REMOVE ME:
#include "rad_paw_trash.f90"

END MODULE rad_paw_routines
