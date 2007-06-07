!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
MODULE grid_paw_routines
  !
  !   WARNINGS:
  !
  ! NO spin-orbit
  ! NO EXX
  ! NO Parallelism
  ! NO rinner > 0
  ! NO Gamma (?)
  !
  !USE temp_PAW_variables
  !
  IMPLICIT NONE
  PUBLIC!              <===
  SAVE

!!!=========================================================================
CONTAINS

  ! Analogous to PW/allocate_nlpot.f90 
  SUBROUTINE allocate_paw_internals
    USE gvect,              ONLY : nrxx
    USE lsda_mod,           ONLY : nspin
    USE parameters,         ONLY : nbrx
    USE ions_base,          ONLY : nsp, nat, ntyp => nsp
    USE us,                 ONLY : nqxq
    USE uspp_param,         ONLY : lmaxq, nhm
    USE gvect,              ONLY : ngl
    !
    USE grid_paw_variables, ONLY : pp, ppt, prad, ptrad, rho1, rho1t, &
         vr1, vr1t, int_r2pfunc, int_r2ptfunc, ehart1, etxc1, vtxc1, &
         ehart1t, etxc1t, vtxc1t, aerho_core, psrho_core, &
         radial_distance, &
         dpaw_ae, dpaw_ps, aevloc, psvloc, aevloc_r, psvloc_r, &
         deband_paw, descf_paw, prodp, prodpt, prod0p, prod0pt
    !
    IMPLICIT NONE
    !
    ALLOCATE (pp(   nhm, nhm, nsp))
    ALLOCATE (ppt(  nhm, nhm, nsp))
    !
    ALLOCATE (int_r2pfunc(   nhm, nhm, nsp))
    ALLOCATE (int_r2ptfunc(  nhm, nhm, nsp))
    !
    IF (lmaxq > 0) ALLOCATE (prad( nqxq, nbrx*(nbrx+1)/2, lmaxq, nsp))
    IF (lmaxq > 0) ALLOCATE (ptrad( nqxq, nbrx*(nbrx+1)/2, lmaxq, nsp))
    !
    ALLOCATE (aevloc( ngl, ntyp))
    ALLOCATE (psvloc( ngl, ntyp))  
    !
    ALLOCATE (aevloc_r(nrxx,nat))
    ALLOCATE (psvloc_r(nrxx,nat))
    ALLOCATE (radial_distance(nrxx,nat))
    !
    ALLOCATE(prodp(nhm*(nhm+1)/2,nhm*(nhm+1)/2,ntyp))
    ALLOCATE(prodpt(nhm*(nhm+1)/2,nhm*(nhm+1)/2,ntyp))
    ALLOCATE(prod0p(nhm*(nhm+1)/2,nhm*(nhm+1)/2,ntyp))
    ALLOCATE(prod0pt(nhm*(nhm+1)/2,nhm*(nhm+1)/2,ntyp))
    !
    ALLOCATE(rho1(nrxx, nspin, nat))
    ALLOCATE(rho1t(nrxx, nspin, nat))
    !
    ALLOCATE(vr1(nrxx, nspin, nat))
    ALLOCATE(vr1t(nrxx, nspin, nat))
    !
    ALLOCATE(ehart1 (nat))
    ALLOCATE(etxc1  (nat))
    ALLOCATE(vtxc1  (nat))
    ALLOCATE(ehart1t(nat))
    ALLOCATE(etxc1t (nat))
    ALLOCATE(vtxc1t (nat))
    ALLOCATE (aerho_core(nrxx, nat))
    ALLOCATE (psrho_core(nrxx, nat))
    !
    ALLOCATE(dpaw_ae( nhm, nhm, nat, nspin))
    ALLOCATE(dpaw_ps( nhm, nhm, nat, nspin))
    !
    ALLOCATE(deband_paw( 2, nat))
    ALLOCATE(descf_paw ( 2, nat))
    !
  END SUBROUTINE allocate_paw_internals


!!$=========================================================================
  ! Analogous to part of PW/init_us_1.f90
  !   Notice integration performed not only up to r(kkbeta(:)) but up to r(msh(:))
  !   because pfunc may extend further than qfunc (see comments "!!kk!!")
  ! + Evaluation of int_r2pfunc
!#define __DEBUG_INIT_PRAD
  SUBROUTINE init_prad
    !
    USE kinds,      ONLY : DP
    USE parameters, ONLY : lmaxx, nbrx, lqmax, ndmx
    USE constants,  ONLY : fpi
    USE atom,       ONLY : r, rab, mesh, msh
    USE ions_base,  ONLY : ntyp => nsp
    USE cell_base,  ONLY : omega, tpiba
    USE gvect,      ONLY : g, gg
    USE lsda_mod,   ONLY : nspin
    USE us,         ONLY : nqxq, dq, nqx, tab, qrad
    USE uspp,       ONLY : qq, qq_so
    USE uspp_param, ONLY : lmaxq, betar, qfunc, qfcoef, rinner, nbeta, &
         kkbeta, nqf, nqlc, lll, jjj, lmaxkb, nh, tvanp, nhm, tvanp, augfun
    USE spin_orb,   ONLY : lspinorb, rot_ylm, fcoef
    !
    USE grid_paw_variables, ONLY: tpawp, pfunc, ptfunc, pp, ppt, prad, ptrad, &
         int_r2pfunc, int_r2ptfunc, okpaw
    !
    IMPLICIT NONE
    !
    REAL(DP), POINTER :: pfunc_(:,:,:,:), prad_(:,:,:,:), pp_(:,:,:), int_r2pfunc_(:,:,:)
    !
    INTEGER :: i_what
    REAL(DP) :: aux2(ndmx)
    !
    ! here a few local variables
    !
    INTEGER :: nt, ih, jh, nb, mb, nmb, l, m, ir, iq, is, startq, &
               lastq, ilast, ndm ! various counters
    REAL(DP), ALLOCATABLE :: aux (:), aux1 (:), besr (:), qtot (:,:,:)
    ! various work space
    REAL(DP) :: prefr, & ! the prefactor of the q functions
                pref,  & ! the prefactor of the beta functions
                q,     & ! the modulus of g for each shell
                qi       ! q-point grid for interpolation
    REAL(DP), ALLOCATABLE :: ylmk0 (:) ! the spherical harmonics
    INTEGER :: n1, m0, m1, n, li, mi, vi, vj, ijs, is1, is2, &
              lk, mk, vk, kh, lh, sph_ind, nnbrx, ll
    COMPLEX(DP) :: coeff, qgm(1)
    REAL(DP) :: spinor, ji, jk

!    call start_clock ('init_prad')
    IF (.NOT.okpaw) RETURN
    !
    !    Initialization of the variables
    !
    ndm = MAXVAL (msh(1:ntyp))
    ALLOCATE ( aux(ndm), aux1(ndm), besr(ndm), qtot(ndm,nbrx,nbrx), &
               ylmk0(lmaxq*lmaxq))
    !
    whattodo: DO i_what=1, 2
       ! associate a pointer to the AE or PS part
       NULLIFY(pfunc_,pp_,prad_)
       IF (i_what==1) THEN
          pfunc_=> pfunc
          pp_   => pp
          prad_ => prad
          int_r2pfunc_ => int_r2pfunc
       ELSE IF (i_what==2) THEN
          pfunc_=> ptfunc
          pp_   => ppt
          prad_ => ptrad
          int_r2pfunc_ => int_r2ptfunc
       END IF

       IF (lmaxq > 0) prad_(:,:,:,:)= 0.d0
       !
       prefr = fpi / omega
       pp_ (:,:,:)   = 0.d0
       !
       ! here for the PAW types we compute the Fourier transform of the 
       ! P functions.
       !
       CALL divide (nqxq, startq, lastq)
       DO nt = 1, ntyp
          IF (tpawp (nt) ) THEN
             DO l = 0, nqlc (nt) - 1
                !
                ! for each nb,mb,l we build the total P_l(|r|) function.
                ! l is the total (combined) angular momentum and 
                ! the arrays have dimensions 1..l+1
                !
                DO nb = 1, nbeta (nt)
                   DO mb = nb, nbeta (nt)
                      IF ( (l >= ABS (lll(nb,nt) - lll(mb,nt) ) ) .AND. &
                           (l <= lll(nb,nt) + lll(mb,nt) )        .AND. &
                           (MOD (l + lll(nb,nt) + lll(mb,nt),2) == 0) ) THEN
                         qtot(1:msh(nt),nb,mb) = pfunc_(1:msh(nt),nb,mb,nt)
                      ENDIF
                   ENDDO ! mb
                ENDDO ! nb
                !
                ! here we compute the spherical bessel function for each |g|
                !
                DO iq = startq, lastq
                   q = (iq - 1) * dq * tpiba
                   CALL sph_bes (msh(nt), r(1,nt), q, l, aux)
                   !
                   ! and then we integrate with all the Q functions
                   !
                   DO nb = 1, nbeta (nt)
                      DO mb = nb, nbeta (nt)
                         ! the P are symmetric with respect to indices
                         nmb = mb * (mb - 1) / 2 + nb
                         IF ( (l >= ABS( lll(nb,nt) - lll(mb,nt) ) ) .AND. &
                              (l <= lll(nb,nt) + lll(mb,nt) )        .AND. &
                              (MOD (l+lll(nb,nt) + lll(mb,nt), 2) == 0) ) THEN
                            !!kk!! DO ir = 1, kkbeta (nt)
                            DO ir = 1, msh (nt)
                               aux1 (ir) = aux (ir) * qtot (ir, nb, mb)
                            ENDDO
                            !!kk!! CALL simpson (kkbeta(nt), ...)
                            CALL simpson (msh(nt), aux1, rab(1, nt), &
                                          prad_(iq,nmb,l + 1, nt) ) 
                         ENDIF
                      ENDDO ! mb
                   ENDDO ! nb
                ENDDO ! iq
             ENDDO ! l
             prad_ (:, :, :, nt) = prad_ (:, :, :, nt) * prefr
#ifdef __PARA
             call reduce (nqxq*nbrx*(nbrx+1)/2*lmaxq, prad_(1,1,1,nt) )
#endif
             ! add augmentation charge in the ps-case
             IF (i_what.EQ.2) prad_(:,:,:,nt) = prad_(:,:,:,nt) + qrad(:,:,:,nt)
          ENDIF ! tpawp
       ENDDO ! nt
       !
       !  and finally we compute the pp_ coefficients by integrating the P.
       !  pp are the g=0 components of P.
       !
#ifdef __PARA
       if (gg (1) > 1.0d-8) goto 100
#endif
       CALL ylmr2 (lmaxq * lmaxq, 1, g, gg, ylmk0)
       DO nt = 1, ntyp
          IF (tpawp (nt) ) THEN
          if (lspinorb) call errore ('init_prad','lspinorb not implemented',1)
             DO ih = 1, nh (nt)
                DO jh = ih, nh (nt)
                   !call qvan2 (1, ih, jh, nt, gg, qgm, ylmk0)
                   CALL pvan2 (1, ih, jh, nt, gg, qgm, ylmk0, prad_, &
                        SIZE(prad_,1),SIZE(prad_,2),SIZE(prad_,3),SIZE(prad_,4))
                   pp_ (ih, jh, nt) = omega *  DBLE (qgm (1) )
                   pp_ (jh, ih, nt) = pp_ (ih, jh, nt)
                ENDDO
             ENDDO
          ENDIF
       ENDDO
#ifdef __PARA
100 continue
!    if (lspinorb) then
!       call reduce ( nhm * nhm * ntyp * 8, pp__so )
!    else
       call reduce ( nhm * nhm * ntyp, pp_ )
!    endif
#endif

       ! Compute the integrals of pfunc*r^2   (not in init_us_1)
       DO nt = 1, ntyp
          IF (tpawp(nt)) THEN
             DO nb = 1, nbeta (nt)
                DO mb = nb, nbeta (nt)
                   aux2(1:msh(nt)) = pfunc_(1:msh(nt), nb, mb, nt) * &
                                     r(1:msh(nt),nt) * r(1:msh(nt),nt)
                   ! add augmentation charge if ps
                   IF (i_what.EQ.2) aux2(1:msh(nt)) = aux2(1:msh(nt)) + &
                                       augfun(1:msh(nt),nb, mb, 0, nt)* &
                                       r(1:msh(nt),nt) * r(1:msh(nt),nt)
                   CALL simpson (msh(nt),aux2,rab(1,nt),int_r2pfunc_(nb,mb,nt))
                END DO
             END DO
          END IF
       END DO
    END DO whattodo
    DEALLOCATE (ylmk0, qtot, besr, aux1, aux)

!    call stop_clock ('init_prad')
    RETURN
  END SUBROUTINE init_prad
!
!----------------------------------------------------------------------------
 SUBROUTINE paw_prod_p
  !----------------------------------------------------------------------------
  !
  ! ... calculates \int \frac{P_{ij}(r)P_{(ij)'}(r')}{|r-r'|}dr dr' in
  ! ... reciprocal space, both for AE and PS P-functions
  !
  USE kinds,         ONLY : DP
  USE gvect,         ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                            ngm, nl, nlm, gg, g, gstart
  USE wvfct,         ONLY : gamma_only
  !  
  USE grid_paw_variables, ONLY : prad, ptrad, tpawp, okpaw, prodp, prodpt,&
                                 prod0p, prod0pt
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,         ONLY : lmaxq, nh, nhm 
  USE wvfct,              ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! local variables
  !
  INTEGER :: gi, ig, na, nt, ih, jh, ijh, ijh2
  ! counters
  !
  REAL(DP), ALLOCATABLE :: qmod (:), ylmk0 (:,:)
  ! the modulus of G
  ! the spherical harmonics
  !
  COMPLEX(DP), ALLOCATABLE ::  qgm(:)
  COMPLEX(DP), ALLOCATABLE ::  pft(:,:,:)
  !
  REAL(DP),    POINTER :: prad_(:,:,:,:)
  COMPLEX(DP), POINTER :: prodp_(:,:,:), prod0p_(:,:,:)
  INTEGER :: i_what
  !
  IF ( .NOT. okpaw ) RETURN
  !
  gi = gstart
  !
  ALLOCATE ( qmod(ngm), qgm(ngm), pft(nhm*(nhm+1)/2,ngm,ntyp), &
             ylmk0(ngm,lmaxq*lmaxq) ) 
  !
  CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  qmod(:) = SQRT(gg(:))
  !
  whattodo: DO i_what=1, 2
     !
     NULLIFY(prad_,prodp_,prod0p_)
     IF (i_what==1) THEN
        prad_ => prad
        prodp_ => prodp
        prod0p_ => prod0p
     ELSE IF (i_what==2) THEN
        prad_ => ptrad
        prodp_ => prodpt
        prod0p_ => prod0pt
     END IF
     !
     prodp_ (:,:,:)  = (0.d0, 0.d0)
     prod0p_ (:,:,:) = (0.d0, 0.d0)
     pft (:,:,:) = (0.d0, 0.d0)
     !
     DO nt = 1, ntyp
        IF (tpawp(nt)) THEN
           ijh = 0
           DO ih = 1, nh (nt)
              DO jh = ih, nh (nt)
                 ijh = ijh + 1
                 CALL pvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0, prad_, &
                      SIZE(prad_,1),SIZE(prad_,2),SIZE(prad_,3),SIZE(prad_,4))
                 pft(ijh,:,nt) =  qgm(:)
              ENDDO ! jh
           ENDDO ! ih
        ENDIF ! tpawp
     ENDDO ! nt
     !
     DO ijh = 1, nhm*(nhm+1)/2
        DO ijh2 = ijh, nhm*(nhm+1)/2
           DO nt = 1, ntyp
              prodp_ (ijh, ijh2, nt) = & 
              SUM( DBLE( CONJG(pft(ijh,gi:,nt))*pft(ijh2,gi:,nt) ) / gg(gi:) )
              !
              prod0p_(ijh, ijh2, nt) = &
              DBLE( CONJG( pft(ijh,1,nt) ) * pft(ijh2,1,nt) )
              !
              prodp_ (ijh2, ijh, nt) = CONJG( prodp_ (ijh, ijh2, nt) )
              prod0p_(ijh2, ijh, nt) = CONJG( prod0p_(ijh, ijh2, nt) )
              !
           ENDDO
        ENDDO
     ENDDO 
     !
  ENDDO whattodo
  !
  DEALLOCATE (qmod, qgm, pft, ylmk0) 
  !
  RETURN
  !
  END SUBROUTINE paw_prod_p

! analogous to compute_onecenter_charges
!
  SUBROUTINE compute_onecenter_charges(becnew, rho1new, rho1tnew)
    !
    USE kinds,                ONLY : DP
    USE ions_base,            ONLY : nat, ntyp => nsp, ityp
    USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                     ngm, nl, nlm, gg, g, eigts1, eigts2, &
                                     eigts3, ig1, ig2, ig3
    USE lsda_mod,             ONLY : nspin
    USE scf,                  ONLY : rho
    USE uspp_param,           ONLY : lmaxq, tvanp, nh, nhm
    USE wvfct,                ONLY : gamma_only
    USE wavefunctions_module, ONLY : psic
    !
    USE cell_base,          ONLY: at
    USE ions_base,          ONLY: tau, atm, ityp
    USE grid_paw_variables, ONLY: prad, ptrad, pp, tpawp, okpaw
    USE us,                 ONLY: qrad
    !
    IMPLICIT NONE
    !
    !first input-output variables
    ! 
    REAL(DP), INTENT(IN) :: becnew (nhm*(nhm+1)/2,nat,nspin)
    REAL(DP), TARGET, INTENT(OUT) :: &
         rho1new(nrxx, nspin, nat), rho1tnew(nrxx,nspin,nat)
    !
    INTEGER :: ig, na, nt, ih, jh, ijh, is ! counters
    !
    REAL(DP), ALLOCATABLE :: qmod (:), & ! the modulus of G
                             ylmk0 (:,:) ! the spherical harmonics
    COMPLEX(DP), ALLOCATABLE :: aux (:,:,:), & ! work space for rho(G,nspin)
                                qgm(:)         ! Fourier transform of q

    REAL(DP), POINTER :: rho1_(:,:,:), prad_(:,:,:,:)
    INTEGER :: i_what

    IF (.NOT.okpaw) RETURN
  
    ALLOCATE (aux(ngm,nspin,nat), qmod(ngm), qgm(ngm), ylmk0(ngm,lmaxq*lmaxq))    
    !  
    CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
    qmod(:) = SQRT(gg(:))

    whattodo: DO i_what=1, 2
       NULLIFY(prad_,rho1_)
       IF (i_what==1) THEN
          prad_ => prad
          rho1_ => rho1new
       ELSE IF (i_what==2) THEN
          prad_ => ptrad
          rho1_ => rho1tnew
       END IF
       aux (:,:,:) = (0.d0, 0.d0)

       DO nt = 1, ntyp
          IF (tpawp (nt) ) THEN
             ijh = 0
             DO ih = 1, nh (nt)
                DO jh = ih, nh (nt)
                   !
                   ijh = ijh + 1
                   CALL pvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0, prad_, &
                        SIZE(prad_,1),SIZE(prad_,2),SIZE(prad_,3),SIZE(prad_,4))
                   DO na = 1, nat
                      !
                      IF (ityp(na).NE.nt) CYCLE
                      DO is = 1, nspin
                         DO ig = 1, ngm
                            aux(ig,is,na) = aux(ig,is,na) +         &
                                            eigts1 (ig1 (ig), na) * &
                                            eigts2 (ig2 (ig), na) * &
                                            eigts3 (ig3 (ig), na) * &
                                            qgm(ig)*becnew(ijh,na,is)
                         ENDDO
                      ENDDO
                      !
                   ENDDO
                   !
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       !
       !     convert aux to real space
       !
       DO na = 1, nat
          IF (tpawp(ityp(na))) THEN
             DO is = 1, nspin
                psic(:) = (0.d0, 0.d0)
                psic( nl(:) ) = aux(:,is,na)
                IF (gamma_only) psic( nlm(:) ) = CONJG(aux(:,is,na))
                CALL cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
                rho1_ (:, is, na) = DBLE (psic (:) )
             ENDDO
          END IF
       END DO
       !
    END DO whattodo
    !
    DEALLOCATE (ylmk0, qgm, qmod, aux)

  END SUBROUTINE compute_onecenter_charges


  ! Analogous to PW/v_of_rho.f90
  ! + evaluation of the spheropole: Int dr r^2 rho(r)
!#define __DEBUG_ONECENTER_POTENTIALS
  SUBROUTINE compute_onecenter_potentials (inp_rho1, inp_rho1t)
    USE kinds,            ONLY : DP
    USE cell_base,        ONLY : at, alat, omega
    USE ions_base,        ONLY: tau, atm, ityp, ntyp=>nsp
    USE ions_base,        ONLY : nat
    USE lsda_mod,         ONLY : nspin
    USE gvect,            ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3,&
                                 nrxx, nl, g, gg
    !
    USE grid_paw_variables, ONLY: vr1, vr1t, radial_distance, & 
         int_r2pfunc, int_r2ptfunc, tpawp, okpaw, ehart1, etxc1, vtxc1, &
         ehart1t, etxc1t, vtxc1t, aerho_core, psrho_core, psvloc_r, aevloc_r
    USE uspp, ONLY: indv, becsum, nhtolm, lpl, ap
    USE uspp_param, ONLY: nh
    USE constants, ONLY: PI
    IMPLICIT NONE
    !
    REAL(DP), TARGET, INTENT(OUT) :: &
         inp_rho1(nrxx, nspin, nat), inp_rho1t(nrxx,nspin,nat)
    !
    REAL(DP), POINTER :: etxc1_(:), vtxc1_(:), ehart1_(:)
    REAL(DP) :: charge, alpha, spheropole
    INTEGER :: na, ih, jh, ijh, nb, mb, is, nt, i
    !
    REAL(DP), POINTER :: rho1_(:,:,:), rho_core_(:,:), vr1_(:,:,:), int_r2pfunc_(:,:,:)
    INTEGER :: i_what
    REAL(DP) :: deltarho, deltav, average_rho, average_pot, rms_rho, rms_pot
    INTEGER :: ir, nir
    !
    IF (.NOT.okpaw) RETURN
    !
    CALL infomsg ('compute_onecenter_potentials','alpha set manually',-1)
    alpha = 0.1_DP
    !
    whattodo: DO i_what=1, 2
       NULLIFY(rho1_,rho_core_,vr1_,int_r2pfunc_,etxc1_,vtxc1_,ehart1_)
       IF (i_what==1) THEN
          rho1_ => inp_rho1 !rho1
          rho_core_ => aerho_core
          vr1_  => vr1
          int_r2pfunc_ => int_r2pfunc
          etxc1_ => etxc1
          vtxc1_ => vtxc1
          ehart1_ => ehart1
       ELSE IF (i_what==2) THEN
          rho1_ => inp_rho1t !rho1t
          rho_core_ => psrho_core
          vr1_  => vr1t
          int_r2pfunc_ => int_r2ptfunc
          etxc1_ => etxc1t
          vtxc1_ => vtxc1t
          ehart1_ => ehart1t
       END IF
       DO na = 1, nat
          IF (tpawp(ityp(na))) THEN
             !
             vr1_(:,:,na)=0._DP
             !
             CALL v_xc( rho1_(:,:,na), rho_core_(:,na), &
                        nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                        nl, ngm, g, nspin, alat, omega, &
                        etxc1_(na), vtxc1_(na), vr1_(:,:,na) )
             !
             spheropole=0.d0
             DO nt = 1, ntyp
                ijh = 0
                DO ih = 1, nh (nt)
                   nb = indv(ih,nt)
                   DO jh = ih, nh (nt)
                      mb = indv(jh,nt)
!                      if (lpl(nhtolm(ih,nt),nhtolm(jh,nt),1).eq.1) &
!                         write (*,*) nhtolm(ih,nt),nhtolm(jh,nt), &
!                                     ap(1,nhtolm(ih,nt),nhtolm(jh,nt))
                      ijh = ijh + 1
                      IF (nhtolm(ih,nt)==nhtolm(jh,nt)) THEN
                         IF (ityp(na)==nt) THEN
                            DO is = 1, nspin
                               spheropole = spheropole + &
                                    int_r2pfunc_(nb,mb,nt) * becsum(ijh,na,is)
                            END DO
                         END IF
                      END IF
                   END DO
                END DO
             END DO
             !
             PRINT *, 'SPHEROPOLE:',na, spheropole
             CALL v_h_grid( rho1_(1,1,na), nr1,nr2,nr3, nrx1,nrx2,nrx3, nrxx, &
                            nl, ngm, gg, gstart, nspin, alat, omega, &
                            ehart1_(na), charge, vr1_(:,:,na), alpha, &
                            spheropole, na)
          END IF
       END DO
    END DO whattodo

  do na=1,nat
     do is=1, nspin
        nir = 0
        average_rho = 0.d0
        rms_rho     = 0.d0
        average_pot = 0.d0
        rms_pot     = 0.d0
        do ir=1,nrxx
           if (radial_distance(ir,na) > 2.0) then
              nir = nir + 1
              !
              deltav = (vr1(ir,is,na)  + aevloc_r(ir,na)) - &
                       (vr1t(ir,is,na) + psvloc_r(ir,na)) 
              average_pot = average_pot + deltav
              rms_pot = rms_pot + deltav*deltav
              !
              deltarho = inp_rho1(ir,is,na) - inp_rho1t(ir,is,na)
              average_rho = average_rho + deltarho
              rms_rho = rms_rho + deltarho*deltarho
           end if
        end do
        if (nir > 0) then
           average_pot = average_pot/nir
           rms_pot = sqrt(rms_pot/nir - average_pot*average_pot)
           average_rho = average_rho/nir
           rms_rho = sqrt(rms_rho/nir - average_rho*average_rho)
           write (*,'(a,2i3,i8,4f18.12)') &
                    ' rho/pot shift for atom ', na, is, nir, &
                      average_rho, rms_rho, average_pot, rms_pot
        else
           write (*,*) 'rho/pot shift for atom failed because nir is ZERO'
        end if
     end do
  end do

  END SUBROUTINE compute_onecenter_potentials
  
  ! Analogous to PW/qvan2.f90
  SUBROUTINE pvan2 (ngy, ih, jh, np, qmod, qg, ylmk0, prad_, s1, s2, s3, s4)
    !
    !#include "f_defs.h"
    USE kinds, ONLY: DP
    USE us, ONLY: dq!, qrad
    USE uspp_param, ONLY: lmaxq, nbrx
    USE uspp, ONLY: nlx, lpl, lpx, ap, indv, nhtolm
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: s1,s2,s3,s4
    REAL(DP), INTENT(IN) :: prad_(s1,s2,s3,s4)

    INTEGER :: ngy, & ! input: the number of G vectors to compute
               ih,  & ! input: the first index of Q
               jh,  & ! input: the second index of Q
               np     ! input: the number of the pseudopotential

    REAL(DP) :: ylmk0 (ngy, lmaxq * lmaxq), & ! the spherical harmonics
                qmod (ngy)         ! input:  moduli of the q+g vectors
    COMPLEX(DP) :: qg (ngy)        ! output: the fourier transform of interest
    !
    !     here the local variables
    !
    COMPLEX(DP) :: sig ! (-i)^L

    INTEGER :: nb,          & ! the atomic index corresponding to ih
               mb,          & ! the atomic index corresponding to jh
               nmb,         & ! combined index (nb,mb)
               ivl,         & ! the lm corresponding to ih
               jvl,         & ! the lm corresponding to jh
               ig,          & ! counter on g vectors
               lp,          & ! the actual LM
               l,           & ! the angular momentum L
               lm,          & ! the possible LM's compatible with ih,j
               i0, i1, i2, i3 ! counters for interpolation table

    REAL(DP) :: sixth,                & ! 1 divided by six
                dqi,                  & ! 1 divided dq
                qm,                   & ! qmod/dq
                px,                   & ! measures for interpolation table
                ux, vx, wx, uvx, pwx, & ! auxiliary variables for intepolation
                work                    ! auxiliary variable
    !
    LOGICAL :: new_qmod
    !
    ! compute the indices which correspond to ih,jh
    !
    sixth = 1.d0 / 6.d0
    dqi = 1 / dq
    nb = indv (ih, np)
    mb = indv (jh, np)
    IF (nb.GE.mb) THEN
       nmb = nb * (nb - 1) / 2 + mb
    ELSE
       nmb = mb * (mb - 1) / 2 + nb
    ENDIF
    ivl = nhtolm(ih, np)
    jvl = nhtolm(jh, np)
    IF (nb.GT.nbrx) CALL errore (' pvan2 ', ' nb.gt.nbrx ', nb)
    IF (mb.GT.nbrx) CALL errore (' pvan2 ', ' mb.gt.nbrx ', mb)
    IF (ivl.GT.nlx) CALL errore (' pvan2 ', ' ivl.gt.nlx  ', ivl)
    IF (jvl.GT.nlx) CALL errore (' pvan2 ', ' jvl.gt.nlx  ', jvl)
    qg(:) = (0.d0, 0.d0)
    !
    !    and make the sum over the non zero LM
    !
    DO lm = 1, lpx (ivl, jvl)
       lp = lpl (ivl, jvl, lm)
       !
       ! extraction of angular momentum l from lp:
       !
       if (lp<1)  CALL errore (' qvan ', ' lp < 1 ', 1)
       l = sqrt(DBLE(lp-1)) + 1
       if (lp>49) CALL errore (' qvan ', ' lp > 49 ', lp)
       !
       sig = (0.d0, -1.d0) ** (l - 1)
       sig = sig * ap (lp, ivl, jvl)
       !
       new_qmod = .true.
       DO ig = 1, ngy
          !
          ! calculate quantites depending on the module of G only when needed
          !
          IF ( ig > 1 ) new_qmod = ABS( qmod(ig) - qmod(ig-1) ) > 1.0D-6
          IF ( new_qmod ) THEN
             qm = qmod (ig) * dqi
             px = qm - INT (qm)
             ux = 1.d0 - px
             vx = 2.d0 - px
             wx = 3.d0 - px
             i0 = INT( qm ) + 1
             i1 = i0 + 1
             i2 = i0 + 2
             i3 = i0 + 3
             uvx = ux * vx * sixth
             pwx = px * wx * 0.5d0
             work = prad_ (i0, nmb, l, np) * uvx * wx + &
                    prad_ (i1, nmb, l, np) * pwx * vx - &
                    prad_ (i2, nmb, l, np) * pwx * ux + &
                    prad_ (i3, nmb, l, np) * px * uvx
          ENDIF
          qg (ig) = qg (ig) + sig * ylmk0 (ig, lp) * work
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE pvan2

  ! From PW/set_rhoc.f90
  SUBROUTINE set_paw_rhoc
    !-----------------------------------------------------------------------
    !
    !    This routine computes the core charge on the real space 3D mesh
    !
    USE io_global, ONLY : stdout
    USE kinds,     ONLY : DP
    USE atom,      ONLY : rho_atc, numeric, msh, r, rab, nlcc
    USE ions_base, ONLY : nat, ityp, ntyp => nsp
    USE cell_base, ONLY : omega, tpiba2, alat
    USE gvect,     ONLY : ngm, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                          nrxx, nl, nlm, ngl, gl, igtongl,      &
                          eigts1, eigts2, eigts3, ig1, ig2, ig3
    USE pseud,     ONLY : a_nlcc, b_nlcc, alpha_nlcc
    USE vlocal,    ONLY : strf
    USE wvfct,     ONLY : gamma_only
    !
    USE grid_paw_variables, ONLY : aerho_atc, psrho_atc, aerho_core, psrho_core
    !
    IMPLICIT NONE
    !
    REAL(DP), POINTER :: rho_atc_(:,:)
    REAL(DP), POINTER :: rho_core_(:,:)
    !
    INTEGER :: i_what
    !
    ! here the local variables
    !
    REAL(DP), PARAMETER :: eps = 1.d-10

    COMPLEX(DP) , ALLOCATABLE :: aux (:) ! used for the fft of the core charge
    REAL(DP) , ALLOCATABLE ::  rhocg(:)  ! the radial fourier trasform
    REAL(DP) ::  rhoima, rhoneg, rhorea  ! used to check the core charge
    REAL(DP) ::  vtxcc                   ! dummy xc energy term
    REAL(DP) , ALLOCATABLE ::  dum(:,:)  ! dummy array containing rho=0
  
    INTEGER :: ir, & ! counter on mesh points
               nt, & ! counter on atomic types
               na, & ! counter on atoms
               ng    ! counter on g vectors
    !
    ALLOCATE (aux(nrxx), rhocg(ngl))
    !    
    whattodo: DO i_what=1, 2
       NULLIFY(rho_atc_,rho_core_)
       IF (i_what==1) THEN
          rho_atc_ => aerho_atc
          rho_core_ => aerho_core
       ELSE IF (i_what==2) THEN
          rho_atc_ => psrho_atc
          rho_core_ => psrho_core
       END IF    
       !
       ! the sum is on atom types
       !
       typ_loop: DO nt = 1, ntyp
          !
          IF ((i_what==2).AND.(.NOT.nlcc(nt))) THEN
             DO na = 1, nat
                IF (ityp (na).EQ.nt) rho_core_(:,na)=0.d0
             END DO
             CYCLE typ_loop
          END IF
          !
          ! drhoc compute the radial fourier transform for each shell of g vec
          !
          CALL drhoc (ngl, gl, omega, tpiba2, numeric (nt), a_nlcc (nt), &
               b_nlcc (nt), alpha_nlcc (nt), msh (nt), r (1, nt), rab (1, nt), &
               rho_atc_ (1, nt), rhocg)
          !
          at_loop: DO na = 1, nat
             at_if: IF (ityp (na) .EQ.nt) THEN
                aux (:) = 0.d0
                aux(nl(:)) = rhocg(igtongl(:)) * &
                             eigts1(ig1(:),na) * &
                             eigts2(ig2(:),na) * &
                             eigts3(ig3(:),na)
                if (gamma_only)  aux(nlm(1:ngm)) = CONJG(aux(nl(1:ngm)))
                !
                ! the core charge in real space
                !
                CALL cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
                !
                ! test on the charge and computation of the core energy
                !
                rhoneg = 0.d0
                rhoima = 0.d0
                DO ir = 1, nrxx
                   rhoneg = rhoneg + min (0.d0,  DBLE (aux (ir) ) )
                   rhoima = rhoima + abs (AIMAG (aux (ir) ) )
                   rho_core_(ir,na) = DBLE (aux(ir))
                ENDDO
                rhoneg = rhoneg / (nr1 * nr2 * nr3)
                rhoima = rhoima / (nr1 * nr2 * nr3)
#ifdef __PARA
                CALL reduce (1, rhoneg)
                CALL reduce (1, rhoima)
#endif
                IF (rhoneg < -1.0d-6 .OR. rhoima > 1.0d-6) then
                   WRITE( stdout, '(/5x,"Check: atom number ", i12)') na
                   WRITE( stdout, '(/12x,"negative/imaginary core charge ", &
                          2f12.6)') rhoneg, rhoima
                ENDIF
                !
             END IF at_if
          END DO at_loop
       END DO typ_loop
    END DO whattodo
    !
    DEALLOCATE (rhocg)
    DEALLOCATE (aux)
    !
    RETURN

  END SUBROUTINE set_paw_rhoc
  
  ! From PW/init_paw_1.f90
  SUBROUTINE step_f(f2,f,r,nrs,nrc,pow,mesh)
    USE kinds , ONLY : dp
    !
    ! This routine apply a function which goes smoothly to zero from rs to rc
    ! 
    IMPLICIT NONE
    INTEGER :: mesh
    REAL(DP), INTENT(out):: f2(mesh)
    REAL(DP), INTENT(in) :: f(mesh), r(mesh)
    REAL(DP), INTENT(in) :: pow
    INTEGER :: nrs, nrc 

    INTEGER :: n,i
    REAL(DP) :: rcp, rsp

    rcp = r(nrc)
    rsp = r(nrs)

    DO i=1,mesh
       IF(r(i).LE.rsp) THEN
          f2(i) = f(i)
       ELSE
          IF(r(i).LE.rcp) THEN
             f2(i)=f(i)* (1.d0-3.d0*((r(i)-rsp)/(rcp-rsp))**2+ &
                  2.d0*((r(i)-rsp)/(rcp-rsp))**3)**pow
          ELSE
             f2(i)=0.d0
          ENDIF
       ENDIF

    END DO

  END SUBROUTINE step_f

! taken from PW/v_of_rho.f90 and adapted to add and substract gaussian charges,
! this gives the hartree potential of isolated atom.
!----------------------------------------------------------------------------
!#define __DEBUG_V_H_GRID
SUBROUTINE v_h_grid( rho, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nl, &
                ngm, gg, gstart, nspin, alat, omega, ehart, charge, v, &
                alpha, spheropole, na )
  !----------------------------------------------------------------------------
  !
  ! ... Hartree potential VH(r) from n(r)
  !
  USE constants, ONLY : fpi, e2
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE wvfct,     ONLY : gamma_only
  USE cell_base, ONLY : tpiba2
  !
  USE constants, ONLY : PI, EPS8
  USE cell_base, ONLY : at
  USE ions_base, ONLY : tau, atm, ityp, ntyp=>nsp
  USE ions_base, ONLY : nat
  USE gvect,     ONLY : eigts1, eigts2, eigts3, ig1, ig2, ig3
  !
  USE grid_paw_variables, ONLY : radial_distance
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nspin, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                         nrxx, ngm, gstart, nl(ngm)
  !
  REAL (DP), INTENT(IN) :: rho(nrxx,nspin), gg(ngm), alat, omega
  !
  REAL(DP), INTENT(IN) :: alpha, spheropole
  INTEGER, INTENT(IN) :: na
  !
  REAL (DP), INTENT(OUT) :: v(nrxx,nspin), ehart, charge
  !
  ! ... local variables
  !
  REAL (DP)              :: fac
  REAL (DP), ALLOCATABLE :: aux(:,:), aux1(:,:)
  INTEGER                     :: ir, is, ig
  !
  REAL(DP) :: dummyx, c(3), r2, gaussrho
  INTEGER :: ir1,ir2,ir3, i,j,k
  COMPLEX(DP) :: skk
  !
  !
  ALLOCATE( aux( 2, nrxx ), aux1( 2, ngm ) )
  !
  ! ... copy total rho in aux
  !
  aux(2,:) = 0.D0
  aux(1,:) = rho(:,1)
  !
  IF ( nspin == 2 ) aux(1,:) = aux(1,:) + rho(:,2)
  !
  ! ... bring rho (aux) to G space
  !
  CALL cft3( aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
  !
  charge = 0.D0
  !
  IF ( gstart == 2 ) charge = omega * aux(1,nl(1))
  !
  CALL reduce( 1, charge )
  !
  PRINT *, 'charge=', charge
  !
  ! ... calculate hartree potential in G-space (NB: only G/=0 )
  !
  ehart     = 0.D0
  aux1(:,:) = 0.D0
  !
  DO ig = gstart, ngm
     !
     ! Define a gaussian distribution of charge, to be substracted
     skk = eigts1(ig1(ig),na) * eigts2(ig2(ig),na) * eigts3(ig3(ig),na)
     gaussrho = charge * EXP(-alpha*gg(ig)*tpiba2) / omega
     !
     fac = 1.D0 / gg(ig)
     !
     ehart = ehart + ( aux(1,nl(ig))**2 + aux(2,nl(ig))**2 - gaussrho**2 ) * fac
     !
     aux1(1,ig) = ( aux(1,nl(ig)) - gaussrho * REAL(skk,DP) ) * fac
     aux1(2,ig) = ( aux(2,nl(ig)) - gaussrho * AIMAG(skk)   ) * fac
     !
  ENDDO
  !
  fac = e2 * fpi / tpiba2
  !
  ehart = ehart * fac
  !
  aux1 = aux1 * fac
  !
  IF (gamma_only) THEN
     !
     ehart = ehart * omega
     !
  ELSE
     !
     ehart = ehart * 0.5D0 * omega
     !
  END IF
  !
#if defined __DEBUG_NEWD_PAW_GRID
!                 WRITE (789,'(8f15.7)') (vr1_(i,1,1), i=1,nrxx)
!                 PRINT '(A,2f20.10)', 'HARTREE1', ehart
#endif
  CALL reduce( 1, ehart )
  !
  ! ... Set G=0 terms
  !
  aux1(1,1) = e2 * (FPI*alpha*charge - 2*PI/3*spheropole) / omega
  ! Add G=0 contribution to Hartree energy
  ehart = ehart + charge * aux1(1,1)
  ! Add analytic contribution from gaussian charges to Hartree energy
  ehart = ehart + 0.5_DP * e2 * charge**2 / SQRT(2*PI*alpha)
  !
#if defined __DEBUG_NEWD_PAW_GRID
!  PRINT '(A,3f20.10)', 'SPHEROPOLE,CHARGE:              ', charge, &
!       spheropole
!  PRINT '(A,3f20.10)', 'V_G=0:              ', aux1(1,1), &
!       e2*(FPI*alpha*charge)/omega, e2*(2*PI/3*spheropole)/omega
!  PRINT '(A,3f20.10)', 'Hartree self-energy:', ehart, &
!     charge * aux1(1,1),    0.5_DP * e2 * charge**2 / SQRT(2*PI*alpha)
#endif
  ! 
  aux(:,:) = 0.D0
  !
  DO ig = 1, ngm
     !
     aux(1,nl(ig)) = aux1(1,ig)
     aux(2,nl(ig)) = aux1(2,ig)
     !
  END DO
  !
  IF ( gamma_only ) THEN
     !
     DO ig = 1, ngm
        !
        aux(1,nlm(ig)) =   aux1(1,ig)
        aux(2,nlm(ig)) = - aux1(2,ig)
        !
     END DO
     !
  END IF
  !
  ! ... transform hartree potential to real space
  !
  CALL cft3( aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
  !
  ! ... Add potential corresponding to gaussian charge (ie, erf(r))
  ! ... in the most inefficient way.
  !
  DO ir=1,nrxx
     r2=radial_distance(ir,na)
     !
     IF (r2 < EPS8) THEN
        aux(1,ir)=aux(1,ir)+e2*charge/SQRT(PI*alpha)
     ELSE
        aux(1,ir)=aux(1,ir)+e2*charge/r2*erf(r2/2/SQRT(alpha))
     END IF
  END DO
  !
  ! ... add hartree potential to the xc potential
  !
  IF ( nspin == 4 ) THEN
     !
     v(:,1) = v(:,1) + aux(1,:)
     !
  ELSE
     !
     DO is = 1, nspin
        !
        v(:,is) = v(:,is) + aux(1,:)
        !
     END DO
     !
  END IF
  !
  DEALLOCATE( aux, aux1 )
  !
  RETURN
  !
END SUBROUTINE v_h_grid

SUBROUTINE set_radial_distance( radial_distance, nr1,nr2,nr3, nrx1,nrx2,nrx3, &
               nrxx, alat)
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  !
  USE constants, ONLY : PI, EPS8
  USE cell_base, ONLY : at, bg
  USE ions_base, ONLY : tau, nat
  !
  IMPLICIT NONE
  !
  INTEGER,   INTENT(IN)  :: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
  REAL (DP), INTENT(IN)  :: alat
  REAL (DP), INTENT(OUT) :: radial_distance(nrxx,nat)
  !
  ! ... local variables
  !
  INTEGER :: i, j, k, ir, idx0, idx, na
  REAL(DP) :: pos_ir(3), pos(3), dist2, dr(3,3)
  idx0 = 0
#if defined (__PARA)
  DO i = 1, me_pool
     idx0 = idx0 + nrx1*nrx2*npp(i)
  END DO
  !
#endif
  dr(:,1) = at(:,1) / DBLE( nr1 )
  dr(:,2) = at(:,2) / DBLE( nr2 )
  dr(:,3) = at(:,3) / DBLE( nr3 )
  do na=1,nat
     do ir=1, nrxx
        ! ... three dimensional indexes
        idx   = idx0 + ir - 1
        k     = idx / (nrx1*nrx2)
        idx   = idx - (nrx1*nrx2)*k
        j     = idx / nrx1
        idx   = idx - nrx1*j
        i     = idx
        !
        pos_ir(:) = i*dr(:,1) + j*dr(:,2) + k*dr(:,3) - tau(:,na)
        ! ... minimum image convenction
        CALL cryst_to_cart( 1, pos_ir, bg, -1 )
        pos_ir(:) = pos_ir(:) - ANINT( pos_ir(:) )
        CALL cryst_to_cart( 1, pos_ir, at, 1 )
        dist2 = pos_ir(1)**2 + pos_ir(2)**2 + pos_ir(3)**2
        do i=-1,1
           do j=-1,1
              do k=-1,1
                 pos(:) = pos_ir(:) + i*at(:,1) + j*at(:,2) + k*at(:,3) 
                 dist2 = min (dist2,pos(1)**2 + pos(2)**2 + pos(3)**2)
              end do
           end do
        end do
        radial_distance(ir,na) = sqrt(dist2) * alat
     end do
  end do
  !
  RETURN
  !
END SUBROUTINE set_radial_distance

! Analogous to PW/newd.f90
!#define __DEBUG_NEWD_PAW_GRID
SUBROUTINE newd_paw_grid
  !
  USE kinds,     ONLY : DP
  !
  USE gvect,              ONLY : nrxx
  USE lsda_mod,           ONLY : nspin
  USE ions_base,          ONLY : nat
  USE grid_paw_variables,   ONLY : okpaw, prad, ptrad, vr1, vr1t, dpaw_ae, dpaw_ps, aevloc_r, psvloc_r, radial_distance
  !
  IMPLICIT NONE
  !
  INTEGER :: i, ir, is, na, nir
  REAL(DP) :: deltav, average, rms
  !
  IF ( .NOT. okpaw ) RETURN
  !
  !PRINT '(A)', 'WARNING newd_paw_grid contains only H+xc potentials'
  !
  CALL integrate_potential_x_charge (prad,  vr1,  aevloc_r, dpaw_ae,  &
       SIZE(prad,1),SIZE(prad,2),SIZE(prad,3),SIZE(prad,4), &
       SIZE(vr1,1),SIZE(vr1,2),SIZE(vr1,3),                 &
       SIZE(aevloc_r,1),SIZE(aevloc_r,2),                   &
       SIZE(dpaw_ae,1),SIZE(dpaw_ae,2),SIZE(dpaw_ae,3),SIZE(dpaw_ae,4))
  CALL integrate_potential_x_charge (ptrad, vr1t, psvloc_r, dpaw_ps,  &
       SIZE(prad,1),SIZE(prad,2),SIZE(prad,3),SIZE(prad,4), &
       SIZE(vr1,1),SIZE(vr1,2),SIZE(vr1,3),                 &
       SIZE(aevloc_r,1),SIZE(aevloc_r,2),                   &
       SIZE(dpaw_ae,1),SIZE(dpaw_ae,2),SIZE(dpaw_ae,3),SIZE(dpaw_ae,4))
  !
CONTAINS
  !
  SUBROUTINE integrate_potential_x_charge (prad_, vr_, vl_, dpaw_, &
       sp1,sp2,sp3,sp4, sv1,sv2,sv3, sl1,sl2, sd1,sd2,sd3,sd4)
    USE kinds,                ONLY : DP
    USE ions_base,            ONLY : nat, ntyp => nsp, ityp
    USE cell_base,            ONLY : omega
    USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                                     g, gg, ngm, gstart, ig1, ig2, ig3, &
                                     eigts1, eigts2, eigts3, nl
    USE lsda_mod,             ONLY : nspin
    USE scf,                  ONLY : vr, vltot
    USE uspp,                 ONLY : deeq, dvan, deeq_nc, dvan_so, okvan
    USE uspp_param,           ONLY : lmaxq, nh, nhm, tvanp
    USE wvfct,                ONLY : gamma_only
    USE wavefunctions_module, ONLY : psic
    USE spin_orb,             ONLY : lspinorb
    USE noncollin_module,     ONLY : noncolin
    !
    USE grid_paw_variables,   ONLY : tpawp
    !
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: sp1,sp2,sp3,sp4, sv1,sv2,sv3, sl1,sl2, sd1,sd2,sd3,sd4
    REAL(DP), INTENT(IN) :: prad_(sp1,sp2,sp3,sp4)
    REAL(DP), INTENT(IN) :: vr_(sv1,sv2,sv3)
    REAL(DP), INTENT(IN) :: vl_(sl1,sl2)
    REAL(DP), INTENT(OUT):: dpaw_(sd1,sd2,sd3,sd4)
    !
    INTEGER :: ig, nt, ih, jh, na, is, nht
    ! counters on g vectors, atom type, beta functions x 2, atoms, spin
    COMPLEX(DP), ALLOCATABLE :: aux(:,:), qgm(:), qgm_na(:)
    ! work space
    REAL(DP), ALLOCATABLE :: ylmk0(:,:), qmod(:)
    ! spherical harmonics, modulus of G
    REAL(DP) :: fact, DDOT
    !
    IF ( gamma_only ) THEN
       fact = 2.D0
    ELSE
       fact = 1.D0
    END IF
    !
    ALLOCATE( aux( ngm, nspin ), qgm_na( ngm ), &
              qgm( ngm ), qmod( ngm ), ylmk0( ngm, lmaxq*lmaxq ) )
    !
    dpaw_(:,:,:,:) = 0._DP
    !
    CALL ylmr2( lmaxq * lmaxq, ngm, g, gg, ylmk0 )
    !
    qmod(1:ngm) = SQRT( gg(1:ngm) )
    !
    DO na=1, nat
       !
       ! ... fourier transform of the total effective potential
       DO is = 1, nspin
          psic(:) = vl_(:,na) + vr_(:,is,na)
          CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1 )
          aux(1:ngm,is) = psic( nl(1:ngm) )
       END DO
       !
       ! ... here we compute the integral Q*V for atom "na",
       ! ...       I = sum_G exp(-iR.G) Q_nm v^*
       !
       ! The type is fixed
       nt = ityp(na)
       !
       IF ( tpawp(nt) ) THEN
          !
          DO ih = 1, nh(nt)
             !
             DO jh = ih, nh(nt)
                !
                ! ... The Q(r) for this atomic species without structure factor
                !
                !CALLqvan2( ngm, ih, jh, nt, qmod, qgm, ylmk0 )
                CALL pvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0, prad_, &
                     sp1,sp2,sp3,sp4)
                !
                ! ... The Q(r) for this specific atom
                !
                qgm_na(1:ngm) = qgm(1:ngm) * eigts1(ig1(1:ngm),na) &
                                           * eigts2(ig2(1:ngm),na) &
                                           * eigts3(ig3(1:ngm),na)
                !
                ! ... and the product with the Q functions
                !
                DO is = 1, nspin
                   !
                   dpaw_(ih,jh,na,is) = fact * omega * &
                        DDOT( 2 * ngm, aux(1,is), 1, qgm_na, 1 )
                   !
                   IF ( gamma_only .AND. gstart == 2 ) &
                        dpaw_(ih,jh,na,is) = dpaw_(ih,jh,na,is) - &
                        omega * DBLE( aux(1,is) * qgm_na(1) )
                   !
                   dpaw_(jh,ih,na,is) = dpaw_(ih,jh,na,is)
                   !
                END DO
                !
             END DO ! jh
             !
          END DO ! ih
          !
       END IF ! tpawp
       !
    END DO ! na
    !
    CALL reduce( nhm * nhm * nat * nspin, dpaw_ )
    !
#if defined __DEBUG_NEWD_PAW_GRID
!    PRINT *, 'D - D1 or D1~'
!    PRINT '(8f15.7)', ((dpaw_(jh,ih,1,1),jh=1,nh(nt)),ih=1,nh(nt))
#endif
    !
    DEALLOCATE( aux, qgm_na, qgm, qmod, ylmk0 )
    !
  END SUBROUTINE integrate_potential_x_charge
    
END SUBROUTINE newd_paw_grid

! Initialize becsum with atomic occupations (for PAW atoms only)
! Notice: requires exact correspondence chi <--> beta in the atom,
! that is that all wavefunctions considered for PAW generation are
! counted in chi (otherwise the array "oc" does not correspond to beta)
!#define __DEBUG_ATOMIC_BECSUM
SUBROUTINE atomic_becsum()
  USE kinds,              ONLY : DP
  USE uspp,               ONLY : becsum, nhtol, indv
  USE uspp_param,         ONLY : nh
  USE ions_base,          ONLY : nat, ityp
  USE lsda_mod,           ONLY : nspin, starting_magnetization
  USE atom,               ONLY : oc
  USE grid_paw_variables, ONLY : tpawp, okpaw
  IMPLICIT NONE
  INTEGER :: ispin, na, nt, ijh, ih, jh, nb, mb
  !
  IF (.NOT. okpaw) RETURN
  !
  if (nspin.GT.2) STOP 'atomic_becsum not implemented'
  !
  na_loop: DO na = 1, nat
     nt = ityp(na)
     is_paw: IF (tpawp(nt)) THEN
        !
        ijh = 1
        ih_loop: DO ih = 1, nh(nt)
           nb = indv(ih,nt)
           !
#if defined __DEBUG_ATOMIC_BECSUM
           PRINT *, ijh,ih,nb,oc(nb,nt),nhtol(ih,nt)
#endif
           IF (nspin==1) THEN
              !
              becsum(ijh,na,1) = oc(nb,nt) / REAL(2*nhtol(ih,nt)+1,DP)
              !
           ELSE IF (nspin==2) THEN
              !
              becsum(ijh,na,1) = 0.5d0 * (1.d0 + starting_magnetization(nt))* &
                                 oc(nb,nt) / REAL(2*nhtol(ih,nt)+1,DP)
              becsum(ijh,na,2) = 0.5d0 * (1.d0 - starting_magnetization(nt))* &
                                 oc(nb,nt) / REAL(2*nhtol(ih,nt)+1,DP)
              !
           END IF
           ijh = ijh + 1
           !
           jh_loop: DO jh = ( ih + 1 ), nh(nt)
              !mb = indv(jh,nt)
              DO ispin = 1, nspin
                 becsum(ijh,na,ispin) = 0._DP
              END DO
              ijh = ijh + 1
              !
           END DO jh_loop
        END DO ih_loop
     END IF is_paw
  END DO na_loop

#if defined __DEBUG_ATOMIC_BECSUM
  PRINT '(1f20.10)', becsum(:,1,1)
#endif

END SUBROUTINE atomic_becsum


  ! Analogous to PW/init_vloc.f90
  SUBROUTINE init_paw_vloc
    !
    USE kinds, ONLY: DP
    USE atom,       ONLY : numeric, msh, mesh, r, rab
    USE ions_base,  ONLY : ntyp => nsp
    USE cell_base,  ONLY : omega, tpiba2
    USE gvect,      ONLY : ngl, gl
    USE pseud,      ONLY : lloc, lmax, cc, nlc, nnl, alpc, alps, aps, zp
    USE grid_paw_variables, ONLY: aevloc_at, psvloc_at, aevloc, psvloc
    !
    USE parameters, ONLY : ndmx
    !
    IMPLICIT NONE
    !
    INTEGER :: nt
    ! counter on atomic types
    !
    REAL(DP), POINTER :: vloc_at_(:,:)
    REAL(DP), POINTER :: vloc_(:,:)
    INTEGER :: i_what
    !  
    whattodo: DO i_what=1, 2
    ! associate a pointer to the AE or PS part
       NULLIFY(vloc_at_,vloc_)
       IF (i_what==1) THEN
          vloc_at_ => aevloc_at
          vloc_    => aevloc
       ELSE IF (i_what==2) THEN
          vloc_at_ => psvloc_at
          vloc_    => psvloc
       END IF
       vloc_(:,:) = 0.d0
       DO nt = 1, ntyp
       !
       ! compute V_loc(G) for a given type of atom
       !
       CALL vloc_of_g_noerf (lloc (nt), lmax (nt), numeric (nt), mesh (nt), &
            msh (nt), rab (1, nt), r (1, nt), vloc_at_ (1, nt), cc (1, &
            nt), alpc (1, nt), nlc (nt), nnl (nt), zp (nt), aps (1, 0, nt), &
            alps (1, 0, nt), tpiba2, ngl, gl, omega, vloc_ (1, nt) )
       END DO
    END DO whattodo 
    !
    RETURN
    !
  END SUBROUTINE init_paw_vloc

! Adapted from vloc_of_g.f90, for bounded potentials: don't add erf(r)/r
! because we don't assume anymore that v(r)\approx 2*e^2*zp/r at large r
! but that it is zero beyond some r_c (see cutoffing in upf_to_internal.f90)
!----------------------------------------------------------------------
subroutine vloc_of_g_noerf (lloc, lmax, numeric, mesh, msh, rab, r, vloc_at, &
     cc, alpc, nlc, nnl, zp, aps, alps, tpiba2, ngl, gl, omega, vloc)
  !----------------------------------------------------------------------
  !
  !    This routine computes the Fourier transform of the local
  !    part of the pseudopotential. Two types of local potentials
  !    are allowed:
  !
  !    a) The pseudopotential is in analytic form and its fourier
  !       transform is computed analytically
  !    b) The pseudopotential is in numeric form and its fourier
  !       transform is computed numerically
  !
  !    The local pseudopotential of the US case is always in
  !    numerical form, expressed in Ry units.
  !
!#include "f_defs.h"
  USE kinds
  implicit none
  !
  !    first the dummy variables
  !
  integer :: nlc, nnl, ngl, lloc, lmax, mesh, msh
  ! input: analytic, number of erf functions
  ! input: analytic, number of gaussian functions
  ! input: the number of shell of G vectors
  ! input: the l taken as local part
  ! input: the maximum non local angular momentum
  ! input: numeric, the dimensions of the mesh
  ! input: numeric, number of mesh points for radial integration

  real(DP) :: cc (2), alpc (2), alps (3, 0:3), aps (6, 0:3), &
       zp, rab (mesh), r (mesh), vloc_at (mesh), tpiba2, omega, gl (ngl), &
       vloc (ngl)
  ! input: analytic, c of the erf functions
  ! input: analytic, alpha of the erf
  ! input: analytic, alpha of the gaussians
  ! input: analytic, a and b of the gaussians
  ! input: valence pseudocharge
  ! input: numeric, the derivative of mesh points
  ! input: numeric, the mesh points
  ! input: numeric, the pseudo on the radial mesh
  ! input: 2 pi / alat
  ! input: the volume of the unit cell
  ! input: the moduli of g vectors for each shell
  ! output: the fourier transform of the potential
  logical :: numeric
  ! input: if true the pseudo is numeric
  !
  real(DP), parameter :: pi = 3.14159265358979d0, fpi= 4.d0 * pi, &
                              e2 = 2.d0, eps= 1.d-8
  !    local variables
  !
  real(DP) :: vlcp, fac, den1, den2, g2a, gx
  real(DP), allocatable :: aux (:), aux1 (:) !  auxiliary variables
  real(DP), external :: erf
  integer :: i, igl, igl0, l, ir
  ! counter on erf functions or gaussians
  ! counter on g shells vectors
  ! first shells with g != 0
  ! the angular momentum
  ! counter on mesh points

  if (.not.numeric) then
     STOP 'vloc_of_g_noerf not implemented'
  else
     !
     ! Pseudopotentials in numerical form (Vloc_at) contain the local part)
!!$NO! in order to perform the Fourier transform, a term erf(r)/r is
!!$NO! subtracted in real space and added again in G space
     !
     allocate ( aux(mesh), aux1(mesh) )
     if (gl (1) < eps) then
        !
        ! first the G=0 term
        !
        do ir = 1, msh
!!$NO      aux (ir) = r (ir) * (r (ir) * vloc_at (ir) + zp * e2)
           aux (ir) = r (ir) * (r (ir) * vloc_at (ir))
        enddo
        call simpson (msh, aux, rab, vlcp)
        vloc (1) = vlcp        
        igl0 = 2
     else
        igl0 = 1
     endif
     !
     !   here the G<>0 terms, we first compute the part of the integrand func
     !   indipendent of |G| in real space
     !
     do ir = 1, msh
!!$NO   aux1 (ir) = r (ir) * vloc_at (ir) + zp * e2 * erf (r (ir) )
        aux1 (ir) = r (ir) * vloc_at (ir)
     enddo
     fac = zp * e2 / tpiba2
     !
     !    and here we perform the integral, after multiplying for the |G|
     !    dependent  part
     !
     do igl = igl0, ngl
        gx = sqrt (gl (igl) * tpiba2)
        do ir = 1, msh
           aux (ir) = aux1 (ir) * sin (gx * r (ir) ) / gx
        enddo
        call simpson (msh, aux, rab, vlcp)
        !
!!$NO   !     here we add the analytic fourier transform of the erf function
!!$NO   !
!!$NO   vlcp = vlcp - fac * exp ( - gl (igl) * tpiba2 * 0.25d0) &
!!$NO        / gl (igl)
        vloc (igl) = vlcp
     enddo
     vloc (:) = vloc(:) * fpi / omega
     deallocate (aux, aux1)
  endif
  return
end subroutine vloc_of_g_noerf

  ! Analogous to setlocal.f90
!----------------------------------------------------------------------
!#define __DEBUG_PAW_GRID_SETLOCAL
subroutine paw_grid_setlocal
  !----------------------------------------------------------------------
  !
  !    This routine computes the local potential in real space vltot(ir)
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : ntyp => nsp
  USE cell_base, ONLY : alat
  USE extfield,  ONLY : tefield, dipfield, etotefield
  USE gvect,     ONLY : igtongl
  USE vlocal,    ONLY : strf !!$ , vloc
  USE wvfct,     ONLY : gamma_only
  USE gvect,     ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nl, nlm, ngm
  USE scf,       ONLY : rho
  !
  USE gvect,     ONLY : eigts1, eigts2, eigts3, ig1, ig2, ig3
  USE grid_paw_variables, ONLY : aevloc, psvloc, aevloc_r, psvloc_r, radial_distance
  USE ions_base, ONLY : nat, ityp
  !
  implicit none
  complex(DP), allocatable :: aux (:)
  ! auxiliary variable
  integer :: nt, ng, ir
  ! counter on atom types
  ! counter on g vectors
  ! counter on r vectors
  !
  REAL(DP), POINTER :: vloc_(:,:)
  REAL(DP), POINTER :: vltot_(:,:)
  INTEGER :: na, i_what
  COMPLEX(DP) :: skk
  !
  allocate (aux( nrxx))    
  !
  whattodo: DO i_what=1, 2
     ! associate a pointer to the AE or PS part
     NULLIFY(vloc_,vltot_)
     IF (i_what==1) THEN
        vloc_   => aevloc
        vltot_  => aevloc_r
     ELSE IF (i_what==2) THEN
        vloc_   => psvloc
        vltot_  => psvloc_r
     END IF

     na_loop: do na=1, nat
        aux(:)=(0.d0,0.d0)
        !
        ! the type is fixed to the specific atom
        nt=ityp(na)
        do ng = 1, ngm
           skk = eigts1(ig1(ng),na) * eigts2(ig2(ng),na) * eigts3(ig3(ng),na)
           aux(nl(ng)) = aux(nl(ng)) + vloc_(igtongl(ng),nt) * skk
        enddo
        if (gamma_only) then
           do ng = 1, ngm
              aux (nlm(ng)) = CONJG(aux (nl(ng)))
           enddo
        end if
        !
        ! aux = potential in G-space . FFT to real space
        !
        call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
        !
        do ir = 1, nrxx
           vltot_ (ir,na) =  DBLE (aux (ir) )
        enddo
        !
        if (tefield.and.(.not.dipfield)) then
           STOP 'paw_grid_setlocal: efield not implemented'
        endif
     end do na_loop
  end DO whattodo
  deallocate(aux)
  ! calculate radial distance from atoms for each grid point
  CALL set_radial_distance( radial_distance, nr1,nr2,nr3, nrx1,nrx2,nrx3, &
               nrxx, alat)
  return
end subroutine paw_grid_setlocal

! Analogous to delta_e in PW/electrons.f90
!-----------------------------------------------------------------------
FUNCTION delta_e_1 ( na, i_what )
!-----------------------------------------------------------------------
  !
  ! ... delta_e = - \int rho(r) V_scf(r)
  !
  USE kinds
  USE gvect,              ONLY : nr1, nr2, nr3
  USE cell_base,          ONLY : omega
  USE grid_paw_variables, ONLY: rho1, rho1t, vr1, vr1t, okpaw
  USE lsda_mod,           ONLY : nspin
  !
  IMPLICIT NONE
  !   
  REAL (DP) :: delta_e_1
  INTEGER :: na, i_what
  !
  REAL(DP), POINTER :: rho_(:,:,:), vr_(:,:,:)
  INTEGER :: ipol
  !
  IF (i_what==1) THEN
     rho_ => rho1
     vr_  => vr1
  ELSE IF (i_what==2) THEN
     rho_ => rho1t
     vr_  => vr1t
  ELSE 
     CALL errore('delta_e_1','wrong i_what',1)
  END IF
  !
  delta_e_1 = 0.D0
  DO ipol = 1, nspin
     delta_e_1 = delta_e_1 - SUM( rho_(:,ipol,na) * vr_(:,ipol,na) )
  END DO
  delta_e_1 = omega * delta_e_1 / ( nr1 * nr2 * nr3 )
  CALL reduce( 1, delta_e_1 )
  !
  RETURN
  !
END FUNCTION delta_e_1

! Analogous to delta_escf in PW/electrons.f90
!-----------------------------------------------------------------------
 FUNCTION delta_e_1scf ( na, i_what )
!-----------------------------------------------------------------------
  !
  ! ... delta_e_1scf = - \int \delta rho(r) V_scf(r)
  ! ... this is the correction needed to have variational energy
  !
  USE kinds
  USE gvect,              ONLY : nr1, nr2, nr3
  USE cell_base,          ONLY : omega
  USE grid_paw_variables, ONLY : rho1, rho1t, vr1, vr1t, rho1new, rho1tnew
  USE lsda_mod,           ONLY : nspin
  !
  IMPLICIT NONE
  !    
  REAL (DP) :: delta_e_1scf
  INTEGER :: na, i_what
  !
  REAL(DP), POINTER :: rho_(:,:,:), rhonew_(:,:,:), vr_(:,:,:)
  INTEGER :: ipol
  !
  IF (i_what==1) THEN
     rho_ => rho1
     rhonew_ => rho1new
     vr_  => vr1
  ELSE IF (i_what==2) THEN
     rho_ => rho1t
     rhonew_ => rho1tnew
     vr_  => vr1t
  ELSE 
     CALL errore('delta_e_1scf','wrong i_what',1)
  END IF
  !
  delta_e_1scf = 0.D0
  DO ipol = 1, nspin
     delta_e_1scf = delta_e_1scf - &
                    SUM( (rhonew_(:,ipol,na)-rho_(:,ipol,na)) * vr_(:,ipol,na) )
  END DO
  delta_e_1scf = omega * delta_e_1scf / ( nr1 * nr2 * nr3 )
  CALL reduce( 1, delta_e_1scf )
  !
  RETURN
  !
END FUNCTION delta_e_1scf

END MODULE grid_paw_routines
