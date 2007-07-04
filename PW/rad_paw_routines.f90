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
FUNCTION rad_energy(pot, rho, e_na)
    USE kinds,                  ONLY : DP
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq
    USE atom,                   ONLY : r, rab, mesh, msh
    USE parameters,             ONLY : ndmx
    USE ions_base,              ONLY : nat, ityp

    REAL(DP), INTENT(IN)    :: rho(ndmx,lmaxq**2,nspin,nat) ! input rho
    REAL(DP), INTENT(IN)    :: pot(ndmx,lmaxq**2,nspin,nat)     ! input potential
    REAL(DP)                :: rad_energy ! output
    REAL(DP),OPTIONAL,INTENT(OUT)    :: e_na(nat)  ! optional output: energy per atom

    REAL(DP)   :: aux(ndmx)  ! workspace
    REAL(DP)                :: par_energy  ! workspace
    INTEGER                 :: na,ispin,lm !counters on atoms, spins, angular momentum
    CALL start_clock ('rad_energy')

    rad_energy = 0._dp
    IF (present(e_na)) e_na(:) = 0._dp

    atoms: DO na = 1, nat
    spins: DO ispin = 1,nspin
        DO lm = 1,lmaxq**2
            aux(:) = rho(:,lm,ispin,na) * pot(:,lm,ispin,na)
            ! V*rho is integrated using simpson algorithm
            CALL simpson (msh(ityp(na)),aux,rab(1,ityp(na)),par_energy)
            rad_energy = rad_energy + par_energy
            IF (present(e_na)) e_na(na) = e_na(na) +par_energy
        ENDDO
    ENDDO spins
    ENDDO atoms

    CALL stop_clock ('rad_energy')

END FUNCTION rad_energy

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! use the density produced by sum_rad_rho to compute potential (only hartree, at this moment)
!!
SUBROUTINE rad_potential(rho1rad, rho1trad, pot)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi, eps8
    USE parameters,             ONLY : ndmx, npsx
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq
    USE uspp,                   ONLY : ap
    USE ions_base,              ONLY : nat, ityp, ntyp => nsp
    USE atom,                   ONLY : r, mesh, msh

    REAL(DP), TARGET, INTENT(IN)    :: rho1rad(ndmx,lmaxq**2,nspin,nat) ! AE charge density on radial grid
    REAL(DP), TARGET, INTENT(IN)    :: rho1trad(ndmx,lmaxq**2,nspin,nat)! the same, but pseudo
    REAL(DP), POINTER               :: rho1rad_(:,:,:,:)                ! pointer to both
    REAL(DP), TARGET, INTENT(OUT)   :: pot(ndmx,lmaxq**2,nspin,nat,2)   ! AE charge density on radial grid
    REAL(DP)                        :: auxrho(ndmx)    ! workspace

    REAL(DP)                        :: r2(ndmx,npsx)    ! r**2
    REAL(DP)                        :: rs(ndmx,npsx)    ! r**.5
    REAL(DP)                        :: dx(npsx)         ! integration step used in atomic code
    REAL(DP)                        :: bogus !DEBUG

    INTEGER, PARAMETER              :: AE = 1, PS = 2   ! All-Electron and Pseudo
    INTEGER                         :: i_what, &
                                       na,nt, &     ! counter on atoms and atom types
                                       ispin, &     ! counter on spins
                                       lm,l,&       ! counter on composite angmom lm = l**2 +m
                                       k ! DEBUG

    REAL(DP)                        :: e(2,nat)    ! tmp
    CALL start_clock ('rad_pot')

    WRITE(6,*) "Compunting radial potentials..."
    
    ! prepare r**2 and r**.5 arrays                     !FIXME: move to init
    r2(:,:) = r(:,:)**2
    rs(:,:) = sqrt(r(:,:))

    !  I have to initialize the integration step dx:    !FIXME: move to init
    dx(:) = 0._dp
    DO nt = 1,ntyp
        IF (r(1,nt) > 0.0_dp) THEN
        ! r(i+1) = exp(xmin)/zmesh * exp(i*dx)
        dx(nt)=log(r(2,nt)/r(1,nt))
        ELSE
        ! r(i+1) = exp(xmin)/zmesh * ( exp(i*dx) - 1 )
        dx(nt)=log( (r(3,nt)-r(2,nt)) / r(2,nt) )
        ENDIF
    ENDDO


    ! this loop computes the hartree potential using the following formula:
    !               l is the first argument in hartree subroutine
    !               r1 = min(r,r'); r2 = MAX(r,r')
    ! V_h(r) = \sum{lm} Y_{lm}(\hat{r}) \int dr' 4\pi r'^2 \rho^{lm}(r') (r1^l/r2^{l+1})/(2L+1)
    !                                            ^^^^^^^^^^^^^^^^^^^^^^ <-- input to the hartree subroutine
    !          output from the h.s. --> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    whattodo: DO i_what = AE, PS
    NULLIFY(rho1rad_)
    IF (i_what == AE) THEN
        rho1rad_ => rho1rad
    ELSE IF (i_what == PS) THEN
        rho1rad_ => rho1trad
    ENDIF
        pot(:,:,:,:,:) = 0._dp
        atoms: DO na = 1, nat
        nt = ityp(na)
        spins: DO ispin = 1,nspin
            DO lm = 1, lmaxq**2
                l = INT(sqrt(DBLE(lm-1)))     ! l has to start from *zero*
                    auxrho(:) = fpi*rho1rad_(:,lm,ispin,na)
                    CALL hartree(l, 2*l+2, mesh(nt), r(:,nt),r2(:,nt),rs(:,nt),&
                                 dx(nt),auxrho(:), pot(:,lm,ispin,na,i_what))
            ENDDO ! lm
            ![[-- debug
            DO lm = 1, lmaxq**2
                DO k = 1, mesh(nt)
                    WRITE(10000*i_what+100*na+lm,"(2f15.8)") r(k,nt), pot(k,lm,ispin,na,i_what)
                    WRITE(10000*(i_what+2)+100*na+lm,"(2f15.8)") r(k,nt), rho1rad_(k,lm,ispin,na)
                    WRITE(10000*(i_what+4)+100*na+lm,"(2f15.8)") r(k,nt), rho1rad_(k,lm,ispin,na)*pot(k,lm,ispin,na,i_what)
                ENDDO
                WRITE(10000*i_what+100*na+lm,"()")
                WRITE(10000*(i_what+2)+100*na+lm,"()")
                WRITE(10000*(i_what+4)+100*na+lm,"()")
            ENDDO
            !debug --]]
        ENDDO spins
        ENDDO atoms
        ![[-- debug
        DO na = 1,nat
             bogus = rad_energy(pot(:,:,:,:,i_what), rho1rad_,e(:,i_what))
        ENDDO
        !debug --]]
    ENDDO whattodo

    WRITE(6,*) "PAW RADIAL ENERGIES: "
    WRITE(6,*) "AE 1     :", e(1,1)
    WRITE(6,*) "AE 2     :", e(2,1)
!    WRITE(6,*) "AE 3     :", e(1,3)
    WRITE(6,*) "PS 1     :", e(1,2)
    WRITE(6,*) "PS 2     :", e(2,2)
!    WRITE(6,*) "PS 3     :", e(2,3)
    WRITE(6,*) "AE tot   :", SUM(e(:,1))
    WRITE(6,*) "PS tot   :", SUM(e(:,2))
    WRITE(6,*) "AE-PS 1  :", e(1,1)-e(1,2)
    WRITE(6,*) "AE-PS 2  :", e(2,1)-e(2,2)
!    WRITE(6,*) "AE-PS 3  :", e(1,3)-e(2,3)
    WRITE(6,*) "AE-PS tot:", SUM(e(:,1))-SUM(e(:,2))
        

    CALL stop_clock ('rad_pot')

END SUBROUTINE rad_potential

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! sum up pfuncs x occupation in order to build radial density's angular momentum components
!!
SUBROUTINE sum_rad_rho(bec, rho1rad, rho1trad)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : eps8
    USE atom,                   ONLY : r, rab, mesh, msh
    USE ions_base,              ONLY : nat, ityp, ntyp => nsp
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq, augfun
    USE uspp,                   ONLY : indv, ap, nhtolm
    USE parameters,             ONLY : ndmx
    USE grid_paw_variables,     ONLY : okpaw, tpawp, pfunc, ptfunc


    INTEGER, PARAMETER              :: AE = 1, PS = 2   ! All-Electron and Pseudo

    REAL(DP), INTENT(IN)            :: bec(nhm*(nhm+1)/2,nat,nspin)     ! cross band occupation
    REAL(DP), TARGET, INTENT(OUT)   :: rho1rad(ndmx,lmaxq**2,nspin,nat) ! AE charge density on radial grid
    REAL(DP), TARGET, INTENT(OUT)   :: rho1trad(ndmx,lmaxq**2,nspin,nat)! the same, but pseudo
    REAL(DP), POINTER               :: rho1rad_(:,:,:,:)                ! pointer to rho1rad & rho1trad
    REAL(DP), POINTER               :: pfunc_(:,:,:,:)                  ! pointer to pfunc and ptfunc
    REAL(DP)                        :: pref                             ! workspace

    INTEGER                 :: i_what, &    ! counter on AE/pseudo
                               nb,mb, &     ! counters for pfunc nb,mb = 1, nh
                               nmb, &       ! composite "triangular" index for pfunc nmb = 1,  nh*(nh+1)/2
                               lm,l,m, &    ! counters for angular momentum lm = l**2+m
                               nt,na, &     ! counters for pseudo type and atom nt = indv(na)
                               ispin,&        ! counter for spin (FIXME: may be unnecessary)
                               k !DEBUG

    CALL start_clock ('sum_rad_rho')
    WRITE(6,*) "RADIAL PAW ROUTINES: sum rad rho"

!     WRITE(6,"(a)") "i_what, na, nb,mb,indv(nb,nt),&
!    indv(mb,nt),lm,nhtolm(nb,nt), nhtolm(mb,nt),pref, ap(lm, nhtolm(nb,nt),&
!    nhtolm(mb,nt)), bec(nmb,na,ispin) , bec(nmb, MOD(na,2)+1, ispin)"

    WRITE(21,*) "becsum used in RAD:"
    whattodo: DO i_what = AE, PS
    NULLIFY(rho1rad_, pfunc_)
    IF (i_what == AE) THEN
        rho1rad_ => rho1rad
        pfunc_   => pfunc
    ELSE IF (i_what == PS) THEN
        rho1rad_ => rho1trad
        pfunc_   => ptfunc
    ENDIF
        rho1rad_(:,:,:,:)  = 0._dp
        !
        atoms: DO na = 1, nat
        nt = ityp(na) ! "nt" is the pseudopot of atom "na"
        ifpaw: IF (tpawp(nt)) THEN
            spins: DO ispin = 1, nspin
            nmb = 0
                ! loop on all pfunc for this kind of pseudo
                DO nb = 1, nh(nt)
                DO mb = nb, nh(nt)
                    nmb = nmb+1 ! nmb = 1, nh*(nh+1)/2
                    WRITE(21,"(a,i3,a,4i3,f12.6)") "-->",nmb,":",nb,mb,na,ispin,bec(nmb,na,ispin)
                    DO lm = 1, lmaxq**2
                        IF ( ABS(ap(lm, nhtolm(mb,nt), nhtolm(nb,nt))) > eps8 ) THEN
                            ! becsum already contains a factor 2 for off-diagonal pfuncs
                            pref = bec(nmb,na,ispin) * ap(lm, nhtolm(nb,nt), nhtolm(mb,nt))

!                             WRITE(6,"(a,i2,2i3,i2,i3,i2,i4,i3,i2,4f10.4)") "-->",&
!                                         i_what, na, nb,mb,indv(nb,nt), indv(mb,nt),&
!                                         lm,nhtolm(nb,nt), nhtolm(mb,nt), &
!                                         pref, ap(lm, nhtolm(nb,nt), nhtolm(mb,nt)),&
!                                         bec(nmb,na,ispin) , bec(nmb, MOD(na,2)+1, ispin)

                            rho1rad_(1:msh(nt),lm,ispin,na) = rho1rad_(1:msh(nt),lm,ispin,na) +&
                                            pref * pfunc_(1:msh(nt), indv(nb,nt), indv(mb,nt), nt)
                            IF (i_what == PS) THEN
                                ! if I'm doing the pseudo part I have to add the augmentation charge
                                l = INT(sqrt(DBLE(lm-1))) ! l has to start from zero
                                rho1rad_(1:msh(nt),lm,ispin,na) = rho1rad_(1:msh(nt),lm,ispin,na) +&
                                            pref * augfun(1:msh(na), indv(nb,nt), indv(mb,nt), l, nt)
                            ENDIF ! i_what == PS
                        ENDIF ! |ap| > 1.0e-8
                    ENDDO ! lm
                    DO k = 1, msh(nt)
                        WRITE(90000+1000*i_what+100*na+10*nb+mb,"(2f15.8)")&
                               r(k,nt), pfunc_(k,indv(nb,nt),indv(mb,nt),nt)
                    ENDDO
                ENDDO !mb
                ENDDO !nb
            ENDDO spins
        ENDIF ifpaw
        ENDDO atoms
    ENDDO whattodo

    CALL stop_clock ('sum_rad_rho')

END SUBROUTINE sum_rad_rho




SUBROUTINE integrate_pfunc
    !
    USE kinds,      ONLY : DP
    USE parameters, ONLY : lmaxx, nbrx, lqmax, ndmx
    USE constants,  ONLY : fpi, eps8, eps4
    USE atom,       ONLY : r, rab, mesh, msh
    USE ions_base,  ONLY : ntyp => nsp
    USE cell_base,  ONLY : omega, tpiba
    USE gvect,      ONLY : g, gg
    USE lsda_mod,   ONLY : nspin
    USE us,         ONLY : nqxq, dq, nqx, tab, qrad
    USE uspp 
    USE uspp_param 
    USE spin_orb,   ONLY : lspinorb, rot_ylm, fcoef
    !
    USE grid_paw_variables, ONLY: tpawp, pfunc, ptfunc, pp, ppt, prad, ptrad, okpaw
    ! for FFt method
    USE gvect,         ONLY : gg, gi =>gstart, ngm
    USE grid_paw_routines, ONLY: pvan2
    !
    IMPLICIT NONE
    !, int_pfunc_(:,:,:)
    REAL(DP), POINTER :: pfunc_(:,:,:,:), prad_(:,:,:,:), pp_(:,:,:), int_pfunc_(:,:,:,:,:)
    REAL(DP),TARGET   :: int_pfunc(nbrx,nbrx,nbrx,nbrx,ntyp),&
                         int_ptfunc(nbrx,nbrx,nbrx,nbrx,ntyp)
    REAL(DP)          :: integral, ap2
    !
    INTEGER :: i_what, terms
    REAL(DP) :: aux2(ndmx)
    !
    ! here a few local variables
    !
    INTEGER :: nt, ih, jh, nb, mb, nc, mc, nmb, l, m, lm, ir, iq, is, ndm ! various counters
    REAL(DP), ALLOCATABLE :: aux (:), aux1 (:)
    ! various work space
    INTEGER :: n1, m0, m1, n, li, mi, vi, vj, ijs, is1, is2, &
              lk, mk, vk, kh, lh, sph_ind, nnbrx, ll,j
    COMPLEX(DP) :: coeff, qgm(1)
    REAL(DP) :: ap_tot
    REAL(DP) :: spinor, ji, jk

    ! for FFT method
    REAL(DP), ALLOCATABLE :: qmod (:),  & ! the modulus of G
                             ylmk0 (:,:)  ! the spherical harmonics
    COMPLEX(DP), ALLOCATABLE  :: pft(:,:,:)  ! the AE/PS wfc products
    COMPLEX(DP), TARGET  :: int_pfunc_fft(nbrx,nbrx,nbrx,nbrx, ntyp),&  ! the AE/PS wfc products
                            int_ptfunc_fft(nbrx,nbrx,nbrx,nbrx, ntyp)
    COMPLEX(DP), POINTER :: int_pfunc_fft_(:,:,:,:,:)
    INTEGER :: terms_counter(9,9,9,9,2,1)

WRITE(6,*) "RADIAL PAW ROUTINES: integrate_pfunc (start)"

RETURN
    !
    ! part1: compute P_ij * P_ab on radial, real space, grid
    !--------------------------------------------------------------------------------

whattodo: DO i_what=1, 2
       ! associate a pointer to the AE or PS part
       NULLIFY(pfunc_,int_pfunc_)
       IF (i_what==1) THEN
          pfunc_=> pfunc
          int_pfunc_ => int_pfunc
       ELSE IF (i_what==2) THEN
          pfunc_=> ptfunc
          int_pfunc_ => int_ptfunc
       END IF
       ! Compute the integrals of pfunc
       DO nt = 1, ntyp ! ntype is the # of atomic species (PP's) is .le. than the # of atoms
          IF (tpawp(nt)) THEN
            ! I have to cicle on pfunc TWICE, and each pfunc has 2 indexes => 4 indexes
            ih = 0
            DO nc = 1, nh(nt)
            DO mc = 1, nh(nt)
                DO nb = 1, nh(nt)
                DO mb = 1, nh(nt)
!                    WRITE(6,*) MAXVAL(pfunc_(1:msh(nt), nb, mb, nt))
                    ih = ih+1
                    int_pfunc_(nc,mc,nb,mb,nt) = 0._DP
                    terms = 0
                    DO lm = 1, lmaxq**2 ! FIXME: is this the right upper bound??
                        ap2 = ap(lm,nhtolm(nc,nt),nhtolm(mc,nt))*ap(lm,nhtolm(nb,nt),nhtolm(mb,nt))
                        IF ( ABS(ap2) > eps8 ) THEN
                            terms = terms +1
                            ! if I don't have the augfun the integral have to be computed only once
                            IF ((i_what == 1) .and. (terms == 1)) THEN
                                aux2(1:msh(nt)) = (pfunc_(1:msh(nt), indv(nb,nt), indv(mb,nt), nt)/r(1:msh(nt),nt))*&
                                                  (pfunc_(1:msh(nt), indv(nc,nt), indv(mc,nt), nt)/r(1:msh(nt),nt))
                                CALL simpson (msh(nt),aux2,rab(1,nt),integral)
                            ENDIF
                            ! with the augfun than I have to compute the integral for each value of lm
                            IF ((i_what == 2)) THEN
                                l = INT(sqrt(DBLE(lm-1))) ! the "+1" is not required, as augfun are labelled 0..l
                                aux2(1:msh(nt)) = (pfunc_(1:msh(nt), indv(nb,nt), indv(mb,nt), nt)&
                                                    +augfun(1:msh(nt), indv(nb,nt), indv(mb,nt), l, nt))/r(1:msh(nt),nt)*&
                                                  (pfunc_(1:msh(nt), indv(nc,nt), indv(mc,nt), nt)&
                                                    +augfun(1:msh(nt), indv(nc,nt), indv(mc,nt), l, nt))/r(1:msh(nt),nt)
                                ! the following line is duplicated (better than using two IF..THEN)
                                CALL simpson (msh(nt),aux2,rab(1,nt),integral)
                            ENDIF
                            ! anyway I have to sum
                            int_pfunc_(nc,mc,nb,mb,nt) = int_pfunc_(nc,mc,nb,mb,nt) + integral*ap2
                        ENDIF
                    ENDDO ! l = 1, lmaxq
                    IF (terms > 0.and. i_what==2)&
                       WRITE(1001,"(i4,2i4,3i2,f14.8,i3,2f8.4)"), ih,i_what,nc,mc,nb,mb,  int_pfunc_(nc,mc,nb,mb,nt),&
                                   terms,integral, ap2
                    terms_counter(nc,mc,nb,mb,i_what,nt) = terms
                END DO !mb
                END DO !nb 
            END DO !mc
            END DO !nc 
            !
          END IF ! tpawp
       END DO ! nt
    END DO whattodo

    !
    ! part2: the same in FFT
    !--------------------------------------------------------------------------------
    WRITE(6,*) "done: radial"
    ALLOCATE (qmod(ngm), pft(ngm,nbrx,nbrx), ylmk0(ngm,lmaxq*lmaxq))
    !
    CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
    qmod(:) = SQRT(gg(:))
    !
    whattodo2: DO i_what=1, 2
        !
        NULLIFY(prad_,int_pfunc_fft_)
        IF (i_what==1) THEN
            prad_ => prad
            int_pfunc_fft_ => int_pfunc_fft
        ELSE IF (i_what==2) THEN
            ! ***NOTE: ptrad already has the augmentation charge***
            prad_ => ptrad
            int_pfunc_fft_ => int_ptfunc_fft
        END IF
         pft(:,:,:) = 0._DP !probably unnecessary
        !
        int_pfunc_fft_ (:,:,:,:,:)  = (0.d0, 0.d0)
        !
        DO nt = 1, ntyp
            ih = 0
            pft (:,:,:) = (0.d0, 0.d0)
            IF (tpawp(nt)) THEN
                DO nc = 1, nh (nt)
                DO mc = 1, nh (nt)
                    CALL pvan2 (ngm, mc, nc, nt, qmod, pft(1,nc,mc), ylmk0, prad_, &
                        SIZE(prad_,1),SIZE(prad_,2),SIZE(prad_,3),SIZE(prad_,4))
                ENDDO ! jh
                ENDDO ! ih

                DO nc = 1, nh(nt)
                DO mc = 1, nh(nt)
                    DO nb = 1, nh(nt)
                    DO mb = 1, nh(nt)
                    ih = ih+1
                    int_pfunc_fft_ (mb,nb,mc,nc, nt) = OMEGA *& 
                    SUM( DBLE( CONJG(pft(:,mc,nc))*pft(:,mb,nb) ))
                    !
                    !int_pfunc_fft_ (ijh2, ijh, nt) = CONJG( int_pfunc_fft_ (ijh, ijh2, nt) )
                 !
                    IF (ABS(int_pfunc_fft_(nc,mc,nb,mb,nt))>eps8 .and. i_what==2) &
                        WRITE(1002,"(i4,2i4,3i2,2f14.8)"), ih,i_what,nc,mc,nb,mb,  int_pfunc_fft_(nc,mc,nb,mb,nt)
                    ENDDO ! mb
                    ENDDO ! nb
                ENDDO ! mc
                ENDDO ! nc
            ENDIF ! tpawp
        ENDDO ! nt
     !
    ENDDO whattodo2
    !
    DEALLOCATE (qmod, pft, ylmk0) 
    WRITE(6,*) "done: FFT"

    lm = 1
    whattodo3: DO i_what=1, 2
        !
        NULLIFY(prad_,int_pfunc_fft_,int_pfunc_)
        IF (i_what==1) THEN
            int_pfunc_ => int_pfunc
            int_pfunc_fft_ => int_pfunc_fft
        ELSE IF (i_what==2) THEN
            int_pfunc_ => int_ptfunc
            int_pfunc_fft_ => int_ptfunc_fft
        END IF
        DO nt = 1, ntyp
            ih = 0
            IF (tpawp(nt)) THEN
                DO nc = 1, nh(nt)
                DO mc = 1, nh(nt)
                    DO nb = 1, nh(nt)
                    DO mb = 1, nh(nt)
!                 DO nc = 1, nh(nt)
!                 DO mc = nc, nh(nt)
!                     DO nb = mc, nh(nt)
!                     DO mb = nb, nh(nt)
                        ih = ih+1
                        IF( ABS(int_pfunc_(nc,mc,nb,mb,nt) - int_pfunc_fft_(nc,mc,nb,mb,nt)) > eps4) THEN
                            WRITE(1003,"(3i4,3i2,f14.8,f16.8,f14.8,i3,4i2,a)"), ih,i_what,nc,mc,nb,mb, &
                             int_pfunc_(nc,mc,nb,mb,nt), int_pfunc_fft_(nc,mc,nb,mb,nt),terms_counter(nc,mc,nb,mb,i_what,nt),&
                             nhtolm(nc,nt),nhtolm(mc,nt),nhtolm(nb,nt),nhtolm(mb,nt),"yadda"
                        ELSE IF ( ABS(int_pfunc_(nc,mc,nb,mb,nt)) > eps8 ) THEN
                            WRITE(1004,"(3i4,3i2,f14.8,f16.8,f14.8,i3,4i2,a)"),ih,i_what,nc,mc,nb,mb, &
                             int_pfunc_(nc,mc,nb,mb,nt), int_pfunc_fft_(nc,mc,nb,mb,nt), terms_counter(nc,mc,nb,mb,i_what,nt),&
                             nhtolm(nc,nt),nhtolm(mc,nt),nhtolm(nb,nt),nhtolm(mb,nt),"blah"
                        ENDIF
                    ENDDO ! mb
                    ENDDO ! nb
                ENDDO ! mc
                ENDDO ! nc
            ENDIF ! tpawp
        ENDDO ! nt
    ENDDO whattodo3
WRITE(6,*) "RADIAL PAW ROUTINES: integrate_pfunc (end)"

STOP

END SUBROUTINE integrate_pfunc


END MODULE rad_paw_routines
