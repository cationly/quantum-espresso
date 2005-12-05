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
    USE gvect,     ONLY : nrxx
    USE lsda_mod,  ONLY : nspin
    USE parameters,       ONLY : nbrx
    USE ions_base,        ONLY : nsp, nat
    USE us,               ONLY : nqxq
    USE uspp_param,       ONLY : lmaxq, nhm
    !
    USE grid_paw_variables, ONLY : pp, ppt, prad, ptrad, rho1, rho1t, &
         vr1, vr1t, int_r2pfunc, int_r2ptfunc
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
    ALLOCATE(rho1(nrxx, nspin, nat))
    ALLOCATE(rho1t(nrxx, nspin, nat))
!!! No more needed since ptfunc already contains the augmentation charge qfunc
!!! ALLOCATE(rho1h(nrxx, nspin, nat))
    !
    ALLOCATE(vr1(nrxx, nspin, nat))
    ALLOCATE(vr1t(nrxx, nspin, nat))
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
    USE uspp,       ONLY : nhtol, nhtoj, nhtolm, dvan, qq, indv, ap, aainit, &
         qq_so, dvan_so, okvan
    USE uspp_param, ONLY : lmaxq, dion, betar, qfunc, qfcoef, rinner, nbeta, &
         kkbeta, nqf, nqlc, lll, jjj, lmaxkb, nh, tvanp, nhm
    USE spin_orb,   ONLY : lspinorb, rot_ylm, fcoef
    !
    USE grid_paw_variables, ONLY: tpawp, pfunc, ptfunc, pp, ppt, prad, ptrad, &
         int_r2pfunc, int_r2ptfunc
    !
    IMPLICIT NONE

    ! NEW
    !
    REAL(DP), POINTER :: pfunc_(:,:,:,:), prad_(:,:,:,:), pp_(:,:,:), int_r2pfunc_(:,:,:)
    !
    INTEGER :: i_what
    REAL(DP) :: aux2(ndmx)
    !
    !     here a few local variables
    !

    INTEGER :: nt, ih, jh, nb, mb, nmb, l, m, ir, iq, is, startq, &
         lastq, ilast, ndm
    ! various counters
    REAL(DP), ALLOCATABLE :: aux (:), aux1 (:), besr (:), qtot (:,:,:)
    ! various work space
    REAL(DP) :: prefr, pref, q, qi
    ! the prefactor of the q functions
    ! the prefactor of the beta functions
    ! the modulus of g for each shell
    ! q-point grid for interpolation
    REAL(DP), ALLOCATABLE :: ylmk0 (:)
    ! the spherical harmonics
    REAL(DP) ::  vll (0:lmaxx), vqint, sqrt2, j
    ! the denominator in KB case
    ! interpolated value
    INTEGER :: n1, m0, m1, n, li, mi, vi, vj, ijs, is1, is2, &
         lk, mk, vk, kh, lh, sph_ind
    COMPLEX(DP) :: coeff, qgm(1)
    REAL(DP) :: spinor, ji, jk

!!$call start_clock ('init_us_1')

    !
    !    Initialization of the variables
    !
    !!kk!! ndm = MAXVAL (kkbeta(1:ntyp))
    ndm = MAXVAL (msh(1:ntyp))
    ALLOCATE (aux ( ndm))    
    ALLOCATE (aux1( ndm))    
    ALLOCATE (besr( ndm))    
    ALLOCATE (qtot( ndm , nbrx , nbrx))
    ALLOCATE (ylmk0( lmaxq * lmaxq))    

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
!!$ [...]
       !
       prefr = fpi / omega
       pp_ (:,:,:)   = 0.d0
       !
!!$ [...]
       !
       !   here for the US types we compute the Fourier transform of the
       !   Q functions.
       !
       CALL divide (nqxq, startq, lastq)
       DO nt = 1, ntyp
          IF (tpawp (nt) ) THEN
             DO l = 0, nqlc (nt) - 1
                !
                !     first we build for each nb,mb,l the total Q(|r|) function
                !     note that l is the true (combined) angular momentum
                !     and that the arrays have dimensions 1..l+1
                !
                DO nb = 1, nbeta (nt)
                   DO mb = nb, nbeta (nt)
                      IF ( (l >= ABS (lll (nb, nt) - lll (mb, nt) ) ) .AND. &
                           (l <= lll (nb, nt) + lll (mb, nt) )        .AND. &
                           (MOD (l + lll (nb, nt) + lll (mb, nt), 2) == 0) ) THEN
                         !!kk!!DO ir = 1, kkbeta (nt)
                         DO ir = 1, msh (nt)
                            !!if (r (ir, nt) >= rinner (l + 1, nt) ) then
!!!                            qtot (ir, nb, mb) = pfunc_ (ir, nb, mb, nt)
                            qtot (ir, nb, mb) = pfunc_ (ir, nb, mb, nt)
                            !!else
                            !!   ilast = ir
                            !!endif
                         ENDDO
                         !!if (rinner (l + 1, nt) > 0.d0) &
                         !!     call setqf(qfcoef (1, l+1, nb, mb, nt), &
                         !!     qtot(1,nb,mb), r(1,nt), nqf(nt),l,ilast)
                      ENDIF
                   ENDDO
                ENDDO
                !
                !     here we compute the spherical bessel function for each |g|
                !
                DO iq = startq, lastq
                   q = (iq - 1) * dq * tpiba
                   !!kk!! CALL sph_bes (kkbeta (nt), r (1, nt), q, l, aux)
                   CALL sph_bes (msh (nt), r (1, nt), q, l, aux)
                   !
                   !   and then we integrate with all the Q functions
                   !
                   DO nb = 1, nbeta (nt)
                      !
                      !    the Q are symmetric with respect to indices
                      !
                      DO mb = nb, nbeta (nt)
                         nmb = mb * (mb - 1) / 2 + nb
                         IF ( (l >= ABS (lll (nb, nt) - lll (mb, nt) ) ) .AND. &
                              (l <= lll (nb, nt) + lll (mb, nt) )        .AND. &
                              (MOD (l + lll(nb, nt) + lll(mb, nt), 2) == 0) ) THEN
                            !!kk!! DO ir = 1, kkbeta (nt)
                            DO ir = 1, msh (nt)
                               aux1 (ir) = aux (ir) * qtot (ir, nb, mb)
                            ENDDO
                            !!kk!! CALL simpson (kkbeta(nt), aux1, rab(1, nt), &
                            CALL simpson (msh(nt), aux1, rab(1, nt), &
                                 prad_(iq,nmb,l + 1, nt) )
                         ENDIF
                      ENDDO
                   ENDDO
                   ! igl
                ENDDO
                ! l
             ENDDO
             prad_ (:, :, :, nt) = prad_ (:, :, :, nt)*prefr
!!$#ifdef __PARA
!!$          call reduce (nqxq * nbrx * (nbrx + 1) / 2 * lmaxq, prad_ (1, 1, 1, nt) )
!!$#endif
          ENDIF
          ! ntyp

       ENDDO
       !
       !   and finally we compute the pp_ coefficients by integrating the Q.
       !   q are the g=0 components of Q.
       !
!!$#ifdef __PARA
!!$    if (gg (1) > 1.0d-8) goto 100
!!$#endif
       CALL ylmr2 (lmaxq * lmaxq, 1, g, gg, ylmk0)
       DO nt = 1, ntyp
          IF (tvanp (nt) ) THEN
!!$          if (lspinorb) then
!!$             ...
!!$          else
             DO ih = 1, nh (nt)
                DO jh = ih, nh (nt)
                   !call qvan2 (1, ih, jh, nt, gg, qgm, ylmk0)
                   CALL pvan2 (1, ih, jh, nt, gg, qgm, ylmk0, prad_, &
                        SIZE(prad_,1),SIZE(prad_,2),SIZE(prad_,3),SIZE(prad_,4))
                   pp_ (ih, jh, nt) = omega *  DBLE (qgm (1) )
                   pp_ (jh, ih, nt) = pp_ (ih, jh, nt)
                ENDDO
             ENDDO
!!$          endif
          ENDIF
       ENDDO
!!$#ifdef __PARA
!!$100 continue
!!$    if (lspinorb) then
!!$       call reduce ( nhm * nhm * ntyp * 8, pp__so )
!!$    else
!!$       call reduce ( nhm * nhm * ntyp, pp_ )
!!$    endif
!!$#endif

       ! Compute the integrals of pfunc*r^2   (not in init_us_1)
       DO nt = 1, ntyp
          IF (tpawp (nt) ) THEN
             DO nb = 1, nbeta (nt)
                DO mb = nb, nbeta (nt)
                   aux2(1:msh(nt)) = pfunc_(1:msh(nt), nb, mb, nt) * &
                        r(1:msh(nt),nt) * r(1:msh(nt),nt)
                   CALL simpson (msh(nt), aux2, rab(1,nt), &
                        int_r2pfunc_(nb,mb,nt))
#if defined __DEBUG_INIT_PRAD
                   WRITE (200000+10000*i_what+100*nb+mb,'(3e20.10)') &
                        (aux2(ir),rab(ir,nt),r(ir,nt),ir=1,mesh(nt))
                   WRITE (*,*)200000+10000*i_what+100*nb+mb, &
                        int_r2pfunc_(nb,mb,nt)
#endif
                END DO
             END DO
          END IF
       END DO
    END DO whattodo

    DEALLOCATE (ylmk0)
    DEALLOCATE (qtot)
    DEALLOCATE (besr)
    DEALLOCATE (aux1)
    DEALLOCATE (aux)

#if defined __DEBUG_INIT_PRAD
    ! Check the integrals
    PRINT *, 'INTEGRALS OF Q, P, Ptilde'
    DO nt = 1, ntyp
       DO ih = 1, nh (nt)
          DO jh = ih, nh (nt)
             PRINT '(3i5,4e15.8)', nt,ih,jh, &
                  qq(ih,jh,nt), pp(ih,jh,nt), ppt(ih,jh,nt), &
                  pp(ih,jh,nt)-ppt(ih,jh,nt)
          END DO
       END DO
    END DO
#endif

!!$    ! Look at the radial fourier transforms
!!$    do nt = 1, ntyp
!!$       do ih = 1, nbrx
!!$          write (10000+nt*100+ih,'(e15.8)') qrad(1:nqxq,ih,1,nt)
!!$          write (20000+nt*100+ih,'(e15.8)') prad(1:nqxq,ih,1,nt)
!!$          write (30000+nt*100+ih,'(e15.8)') ptrad(1:nqxq,ih,1,nt)
!!$       end do
!!$    end do

!!$  call stop_clock ('init_us_1')
#if defined __DEBUG_INIT_PRAD
    PRINT '(A)', 'Writing files with some Q, P1, P1t without structure factor'
    CALL plot_augfun(1,1,1,1000,'QVAN')
    CALL plot_augfun(1,1,1,2000,'P1AE')
    CALL plot_augfun(1,1,1,3000,'P1PS')
#endif
    !
#if defined __DEBUG_INIT_PRAD
    STOP 'STOP __DEBUG_INIT_PRAD'
#endif
    RETURN
  END SUBROUTINE init_prad


  ! Analogous to PW/addusdens.f90
#define __DEBUG_COMPUTE_ONECENTER_CHARGES
  SUBROUTINE compute_onecenter_charges
    USE kinds,                ONLY : DP
    USE ions_base,            ONLY : nat, ntyp => nsp, ityp
    USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
         ngm, nl, nlm, gg, g, eigts1, eigts2, &
         eigts3, ig1, ig2, ig3
    USE lsda_mod,             ONLY : nspin
    USE scf,                  ONLY : rho
    USE uspp,                 ONLY : becsum, okvan
    USE uspp_param,           ONLY : lmaxq, tvanp, nh
    USE wvfct,                ONLY : gamma_only
    USE wavefunctions_module, ONLY : psic
    !
    USE cell_base,          ONLY: alat, at
    USE ions_base,          ONLY: tau, atm
    USE grid_paw_variables, ONLY: prad, ptrad, rho1, rho1t, pp
    USE us,                 ONLY: qrad
    !
    IMPLICIT NONE
    !INTEGER, INTENT(in) :: ih_, jh_, na_, unit_
    !CHARACTER(LEN=4), INTENT(in) :: which_one_
    !
    INTEGER :: ig, na, nt, ih, jh, ijh, is
    ! counters

    REAL(DP), ALLOCATABLE :: qmod (:), ylmk0 (:,:)
    ! the modulus of G
    ! the spherical harmonics

    COMPLEX(DP) :: skk
    COMPLEX(DP), ALLOCATABLE ::  aux (:,:,:), qgm(:)
    ! work space for rho(G,nspin)
    ! Fourier transform of q

    REAL(DP), POINTER :: rho1_(:,:,:), prad_(:,:,:,:)
    INTEGER :: i_what
    REAL(DP) :: charge

    !ALLOCATE (aux ( ngm, nspin))    
    ALLOCATE (aux ( ngm, nspin, nat))    
    ALLOCATE (qmod( ngm))    
    ALLOCATE (qgm( ngm))    
    ALLOCATE (ylmk0( ngm, lmaxq * lmaxq))    
    
    !aux (:,:,:) = (0.d0, 0.d0)
    CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
    DO ig = 1, ngm
       qmod (ig) = SQRT (gg (ig) )
    ENDDO

#if defined  __DEBUG_COMPUTE_ONECENTER_CHARGES
    PRINT '(A)', 'WARNING manually fixing occupation'
    !!2s2 2p4 occupation, for O pseudopotential
    becsum=0.d0
    becsum(1,1,1)=2.d0
    becsum(16,1,1)=4.d0/3
    becsum(22,1,1)=4.d0/3
    becsum(27,1,1)=4.d0/3
#endif

    whattodo: DO i_what=1, 2
       NULLIFY(prad_,rho1_)
       IF (i_what==1) THEN
          prad_ => prad
          rho1_ => rho1
       ELSE IF (i_what==2) THEN
          prad_ => ptrad
          rho1_ => rho1t
!!! No more needed since ptfunc already contains the augmentation charge qfunc
!!!    ELSE IF (i_what==3) THEN
!!!       prad_ => qrad
!!!       rho1_ => rho1h
       END IF
       aux (:,:,:) = (0.d0, 0.d0)
       charge = 0.d0

       DO nt = 1, ntyp
          IF (tvanp (nt) ) THEN
             ijh = 0
             DO ih = 1, nh (nt)
                DO jh = ih, nh (nt)
                   ijh = ijh + 1

                   !CALL qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
                   CALL pvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0, prad_, &
                        SIZE(prad_,1),SIZE(prad_,2),SIZE(prad_,3),SIZE(prad_,4))
                   DO na = 1, nat

                      IF (ityp (na) .EQ.nt) THEN
                         DO is = 1, nspin
#if defined  __DEBUG_COMPUTE_ONECENTER_CHARGES
                            WRITE (3,*) ih,jh,ijh,becsum(ijh,na,is)
                            charge=charge+becsum(ijh,na,is)*pp(ih,jh,nt)
#endif
                            DO ig = 1, ngm
                               skk = eigts1 (ig1 (ig), na) * &
                                    eigts2 (ig2 (ig), na) * &
                                    eigts3 (ig3 (ig), na)
                               !aux(ig,is)=aux(ig,is) + qgm(ig)*skk*becsum(ijh,na,is)
                               aux(ig,is,na)=aux(ig,is,na) + qgm(ig)*skk*becsum(ijh,na,is)
                            ENDDO
                         ENDDO
                      ENDIF

                   ENDDO

                ENDDO
             ENDDO
          ENDIF
       ENDDO
#if defined  __DEBUG_COMPUTE_ONECENTER_CHARGES
       PRINT *, '1c_charges, ', i_what, 'charge=', charge
#endif
       !
       !     convert aux to real space
       !
       DO na = 1, nat
          DO is = 1, nspin
             !!rho1_(:,is,na)=0.d0  
             psic(:) = (0.d0, 0.d0)
             psic( nl(:) ) = aux(:,is,na)
             IF (gamma_only) psic( nlm(:) ) = CONJG(aux(:,is,na))
             CALL cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
             ! Notice: in addusdens augmentation charge is ADDED to rho:
             rho1_ (:, is, na) = DBLE (psic (:) )
          ENDDO
       END DO
       !
#if defined  __DEBUG_COMPUTE_ONECENTER_CHARGES
       PRINT '(A)', 'Writing files with charges of selected atoms'
       CALL xsf_struct (alat, at, nat, tau, atm, ityp, 9870+i_what)
       CALL xsf_fast_datagrid_3d &
            (rho1_(:,1,1), nr1, nr2, nr3, nrx1, nrx2, nrx3, at, alat, 9870+i_what)
#endif
       !

    END DO whattodo
    !
    DEALLOCATE (ylmk0)
    DEALLOCATE (qgm)
    DEALLOCATE (qmod)

    DEALLOCATE (aux)

!!! No more needed since ptfunc already contains the augmentation charge qfunc
!!! ! Add compensation charge to rho1t:
!!! rho1t(:,:,:) = rho1t(:,:,:) + rho1h(:,:,:)

  END SUBROUTINE compute_onecenter_charges


  ! Analogous to PW/v_of_rho.f90
  ! + evaluation of the spheropole: Int dr r^2 rho(r)
#define __DEBUG_ONECENTER_POTENTIALS
  SUBROUTINE compute_onecenter_potentials
    USE kinds,            ONLY : DP
    USE cell_base,        ONLY : at, alat, omega
    USE ions_base,          ONLY: tau, atm, ityp, ntyp=>nsp
    USE ions_base,        ONLY : nat
    USE lsda_mod,         ONLY : nspin
    USE gvect,            ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3,&
         nrxx, nl, g, gg
    !
    USE grid_paw_variables, ONLY: rho1, rho1t, vr1, vr1t, &
         int_r2pfunc, int_r2ptfunc
    USE uspp, ONLY: indv, becsum, nhtolm
    USE uspp_param, ONLY: nh
    USE constants, ONLY: PI
    IMPLICIT NONE
    !
    REAL(DP) :: ehart, charge, alpha, spheropole
    INTEGER :: na, ih, jh, ijh, nb, mb, is, nt
    !
    REAL(DP), POINTER :: rho1_(:,:,:), vr1_(:,:,:), int_r2pfunc_(:,:,:)
    INTEGER :: i_what
    !
    CALL infomsg ('compute_onecenter_potentials','alpha set manually',-1)
    alpha = 0.1d0
    !
    whattodo: DO i_what=1, 2
       NULLIFY(rho1_,vr1_)
       IF (i_what==1) THEN
          rho1_ => rho1
          vr1_  => vr1
          int_r2pfunc_ => int_r2pfunc
       ELSE IF (i_what==2) THEN
          rho1_ => rho1t
          vr1_  => vr1t
          int_r2pfunc_ => int_r2ptfunc
       END IF
       DO na = 1, nat
          !
          vr1_(:,:,na)=0._DP
          !
          !call v_xc....
          !
          spheropole=0.d0
          DO nt = 1, ntyp
             ijh = 0
             DO ih = 1, nh (nt)
                nb = indv(ih,nt)
                DO jh = ih, nh (nt)
                   mb = indv(jh,nt)
                   ijh = ijh + 1
                   IF (nhtolm(ih,nt)==nhtolm(jh,nt)) THEN
                      IF (ityp(na)==nt) THEN
                         DO is = 1, nspin
                            spheropole = spheropole + &
                                 int_r2pfunc_(nb,mb,nt) * becsum(ijh,na,is)
#if defined __DEBUG_ONECENTER_POTENTIALS
                            WRITE (4,'(3i5,f20.10,2i5,f20.10)') ih, jh, ijh, &
                                 becsum(ijh,na,is), nb, mb, int_r2pfunc_(nb,mb,nt)
#endif
                         END DO
                      END IF
                   END IF
                END DO
             END DO
          END DO
          !
#if defined __DEBUG_ONECENTER_POTENTIALS
          PRINT *, 'SPHEROPOLE:',spheropole
#endif
          CALL v_h_grid( rho1_(:,:,na), nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
               nl, ngm, gg, gstart, nspin, alat, omega, ehart, charge, vr1_(:,:,na), &
               alpha, spheropole, na)
       END DO

       !CALL xsf_struct (alat, at, nat, tau, atm, ityp, 93000+i_what)
       !CALL xsf_fast_datagrid_3d &
       !     (vr1_(:,1,1), nr1, nr2, nr3, nrx1, nrx2, nrx3, at, alat, 93000+i_what)

    END DO whattodo
#if defined __DEBUG_ONECENTER_POTENTIALS
    WRITE (93000,'(i5,4f20.10)') (na,rho1(na,1,1),vr1(na,1,1),rho1t(na,1,1),vr1t(na,1,1),na=1,nr1)
    WRITE (93000,*)
#endif

#if defined __DEBUG_ONECENTER_POTENTIALS
    STOP 'STOP __DEBUG_ONECENTER_POTENTIALS'
#endif

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

    INTEGER :: ngy, ih, jh, np
    ! input: the number of G vectors to compute
    ! input: the first index of Q
    ! input: the second index of Q
    ! input: the number of the pseudopotential

    REAL(DP) :: ylmk0 (ngy, lmaxq * lmaxq), qmod (ngy)
    ! the spherical harmonics
    ! input: moduli of the q+g vectors
    COMPLEX(DP) :: qg (ngy)
    ! output: the fourier transform of interest
    !
    !     here the local variables
    !

    COMPLEX(DP) :: sig
    ! (-i)^L

    INTEGER :: nb, mb, nmb, ivl, jvl, ig, lp, l, lm, i0, i1, i2, i3
    ! the atomic index corresponding to ih
    ! the atomic index corresponding to jh
    ! combined index (nb,mb)
    ! the lm corresponding to ih
    ! the lm corresponding to jh
    ! counter on g vectors
    ! the actual LM
    ! the angular momentum L
    ! the possible LM's compatible with ih,j
    ! counters for interpolation table

    REAL(DP) :: sixth, dqi, qm, px, ux, vx, wx, uvx, pwx, work
    ! 1 divided by six
    ! 1 divided dq
    ! qmod/dq
    ! measures for interpolation table
    ! auxiliary variables for intepolation
    ! auxiliary variable
    !
    LOGICAL :: ltest
    !
    !     compute the indices which correspond to ih,jh
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
       !     extraction of angular momentum l from lp:
       !
       IF (lp.EQ.1) THEN
          l = 1
       ELSEIF ( (lp.GE.2) .AND. (lp.LE.4) ) THEN
          l = 2
       ELSEIF ( (lp.GE.5) .AND. (lp.LE.9) ) THEN
          l = 3
       ELSEIF ( (lp.GE.10) .AND. (lp.LE.16) ) THEN
          l = 4
       ELSEIF ( (lp.GE.17) .AND. (lp.LE.25) ) THEN
          l = 5
       ELSEIF ( (lp.GE.26) .AND. (lp.LE.36) ) THEN
          l = 6
       ELSEIF ( (lp.GE.37) .AND. (lp.LE.49) ) THEN
          l = 7
       ELSE
          CALL errore (' qvan ', ' lp > 49 ', lp)
       ENDIF
       sig = (0.d0, -1.d0) ** (l - 1)
       sig = sig * ap (lp, ivl, jvl)
       DO ig = 1, ngy
          !
          ! calculate quantites depending on the module of G only when needed
          !
          IF ( ig > 1 ) ltest = ABS( qmod(ig) - qmod(ig-1) ) > 1.0D-6
          !
          IF ( ig == 1 .OR. ltest ) THEN
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


  ! Analogous to PW/addusdens.f90
  SUBROUTINE plot_augfun(ih_,jh_,na_,unit_,which_one_)
    USE kinds,                ONLY : DP
    USE ions_base,            ONLY : nat, ntyp => nsp, ityp
    USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
         ngm, nl, nlm, gg, g, eigts1, eigts2, &
         eigts3, ig1, ig2, ig3
    USE lsda_mod,             ONLY : nspin
    USE scf,                  ONLY : rho
    USE uspp,                 ONLY : becsum, okvan
    USE uspp_param,           ONLY : lmaxq, tvanp, nh
    USE wvfct,                ONLY : gamma_only
    USE wavefunctions_module, ONLY : psic
    !
    USE cell_base,          ONLY: alat, at
    USE ions_base,          ONLY: tau, atm
    USE grid_paw_variables, ONLY: prad, ptrad
    !
    IMPLICIT NONE
    INTEGER, INTENT(in) :: ih_, jh_, na_, unit_
    CHARACTER(LEN=4), INTENT(in) :: which_one_
    !
    INTEGER :: ig, na, nt, ih, jh, ijh, is
    ! counters

    REAL(DP), ALLOCATABLE :: qmod (:), ylmk0 (:,:)
    ! the modulus of G
    ! the spherical harmonics

    COMPLEX(DP) :: skk
    COMPLEX(DP), ALLOCATABLE ::  aux (:,:), qgm(:)
    ! work space for rho(G,nspin)
    ! Fourier transform of q

    REAL(DP) :: rhor(nrxx)

    IF (.NOT.okvan) RETURN

    ALLOCATE (aux ( ngm, nspin))    
    ALLOCATE (qmod( ngm))    
    ALLOCATE (qgm( ngm))    
    ALLOCATE (ylmk0( ngm, lmaxq * lmaxq))    
    
    aux (:,:) = (0.d0, 0.d0)
    CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
    DO ig = 1, ngm
       qmod (ig) = SQRT (gg (ig) )
    ENDDO

    DO nt = 1, ntyp
       IF (tvanp (nt) ) THEN
          ijh = 0
          DO ih = 1, nh (nt)
             DO jh = ih, nh (nt)
                ijh = ijh + 1
                
                IF ((ih==ih_).AND.(jh==jh_)) THEN  !!!!

                   SELECT CASE (which_one_)
                   CASE ('QVAN')
                      CALL qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
                   CASE ('P1AE')
                      CALL pvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0, prad, &
                           SIZE(prad,1),SIZE(prad,2),SIZE(prad,3),SIZE(prad,4))
                   CASE ('P1PS')
                      CALL pvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0, ptrad, &
                           SIZE(ptrad,1),SIZE(ptrad,2),SIZE(ptrad,3),SIZE(ptrad,4))
                   CASE DEFAULT
                      CALL errore ('plot_augfun','what to plot',-1)
                   END SELECT

                   DO na = 1, nat
                      IF (na==na_) THEN  !!!!
                         IF (ityp (na) .EQ.nt) THEN
                            !write (*,*) 'here we are'
                            DO is = 1, nspin
                               DO ig = 1, ngm
!                                  skk = eigts1 (ig1 (ig), na) * &
!                                       eigts2 (ig2 (ig), na) * &
!                                       eigts3 (ig3 (ig), na)
                                  aux(ig,is)=aux(ig,is) + qgm(ig)!*skk*becsum(ijh,na,is)
                               ENDDO
                            ENDDO
                            !write (*,*)'skk',skk
                            !write (1,*) qgm(1:ngm)
                            !write (2,*) (aux(1:ngm,1)+aux(1:ngm,2))
                         ENDIF
                      END IF
                   ENDDO
                END IF
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    !
    DEALLOCATE (ylmk0)
    DEALLOCATE (qgm)
    DEALLOCATE (qmod)
    !
    !     convert aux to real space
    !
    rhor(:)=0.d0
    DO is = 1, nspin
       psic(:) = (0.d0, 0.d0)
       psic( nl(:) ) = aux(:,is)
       IF (gamma_only) psic( nlm(:) ) = CONJG(aux(:,is))
       CALL cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
       !rho (:, is) = rho (:, is) +  DBLE (psic (:) )
       rhor(:)=rhor(:)+  DBLE (psic (:) )
    ENDDO
    !
    CALL xsf_struct (alat, at, nat, tau, atm, ityp, unit_+999)
    CALL xsf_fast_datagrid_3d &
         (rhor, nr1, nr2, nr3, nrx1, nrx2, nrx3, at, alat, unit_+999)
    !
    DEALLOCATE (aux)
  END SUBROUTINE plot_augfun



! taken from PW/v_of_rho.f90 and adapted to add and substract gaussian charges
!----------------------------------------------------------------------------
#define __DEBUG_V_H_GRID
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
#if defined __DEBUG_V_H_GRID
  PRINT *, 'charge=', charge
#endif
  !
  ! ... calculate hartree potential in G-space (NB: only G/=0 )
  !
  ehart     = 0.D0
  aux1(:,:) = 0.D0
  !
  DO ig = gstart, ngm
     !
     ! Define a gaussian distribution of charge, to be substracted
     skk = eigts1 (ig1 (ig), na) * eigts2 (ig2 (ig), na) * eigts3 (ig3 (ig), na)
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
#if defined __DEBUG_V_H_GRID
  PRINT '(A,3f20.10)', 'V_G=0:              ', aux1(1,1), &
       e2*(FPI*alpha*charge)/omega, e2*(2*PI/3*spheropole)/omega
  PRINT '(A,3f20.10)', 'Hartree self-energy:', ehart, &
     charge * aux1(1,1),    0.5_DP * e2 * charge**2 / SQRT(2*PI*alpha)
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
#if defined __DEBUG_V_H_GRID
  PRINT *, at(:,1)
  PRINT *, at(:,2)
  PRINT *, at(:,3)
  PRINT *, tau(1:3,na)
#endif
  ir=0
  zzz: DO ir3=1,nr3
     c(3)=REAL(ir3-1,DP)/nr3 - tau(3,na)
     IF (c(3)>+0.5d0) c(3)=c(3)-1.d0
     IF (c(3)<-0.5d0) c(3)=c(3)+1.d0
     !
     yyy: DO ir2=1,nr2
        c(2)=REAL(ir2-1,DP)/nr2 - tau(2,na)
        IF (c(2)>+0.5d0) c(2)=c(2)-1.d0
        IF (c(2)<-0.5d0) c(2)=c(2)+1.d0
        !
        xxx: DO ir1=1,nr1
           c(1)=REAL(ir1-1,DP)/nr1 - tau(1,na)
           IF (c(1)>+0.5d0) c(1)=c(1)-1.d0
           IF (c(1)<-0.5d0) c(1)=c(1)+1.d0
           !
           ir=ir+1
           !
           ! \mathbf{r} = Sum_i c_i \mathbf{a}_i
           ! r^2 = Sum_{ij} c_i c_j Sum_k (\mathbf{a}_i)_k (\mathbf{a}_j)_k
           r2=0.d0
           DO i=1,3
              DO j=1,3
                 DO k=1,3
                    r2 = r2 + c(i) * c(j) * at(k,i) * at(k,j)
                 END DO
              END DO
           END DO
           r2=SQRT(r2)*alat
           !
           IF (r2 < EPS8) THEN
              aux(1,ir)=aux(1,ir)+e2*charge/SQRT(PI*alpha)
           ELSE
              aux(1,ir)=aux(1,ir)+e2*charge/r2*erf(r2/2/SQRT(alpha))
           END IF
           !
        END DO xxx
     END DO yyy
  END DO zzz
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








END MODULE grid_paw_routines

!!$==========================================================================

!
! Copyright (C) 2003 Tone Kokalj
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file holds XSF (=Xcrysden Structure File) utilities.
! Routines written by Tone Kokalj on Mon Jan 27 18:51:17 CET 2003
!
! -------------------------------------------------------------------
!   this routine writes the crystal structure in XSF format
! -------------------------------------------------------------------
SUBROUTINE xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER          :: nat, ityp (nat), ounit
  CHARACTER(len=3) :: atm(*)
  REAL(DP)    :: alat, tau (3, nat), at (3, 3)
  ! --
  INTEGER          :: i, j, n
  REAL(DP)    :: at1 (3, 3)
  ! convert lattice vectors to ANGSTROM units ...
  DO i=1,3
     DO j=1,3
        at1(j,i) = at(j,i)*alat*0.529177d0
     ENDDO
  ENDDO

  WRITE(ounit,*) 'CRYSTAL'
  WRITE(ounit,*) 'PRIMVEC'
  WRITE(ounit,'(2(3F15.9/),3f15.9)') at1
  WRITE(ounit,*) 'PRIMCOORD'
  WRITE(ounit,*) nat, 1

  DO n=1,nat
     ! positions are in Angstroms
     WRITE(ounit,'(a3,3x,3f15.9)') atm(ityp(n)), &
          tau(1,n)*alat*0.529177d0, &
          tau(2,n)*alat*0.529177d0, &
          tau(3,n)*alat*0.529177d0
  ENDDO
  RETURN
END SUBROUTINE xsf_struct



! -------------------------------------------------------------------
!   this routine writes the 3D scalar field (i.e. uniform mesh of points)
!   in XSF format using the FFT mesh (i.e. fast write)
! -------------------------------------------------------------------
SUBROUTINE xsf_fast_datagrid_3d &
     (rho, nr1, nr2, nr3, nrx1, nrx2, nrx3, at, alat, ounit)
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER       :: nrx1, nrx2, nrx3, nr1, nr2, nr3, ounit
  REAL(DP) :: alat, at (3, 3), rho(nrx1,nrx2,nrx3)
  ! --
  INTEGER       :: i1, i2, i3, ix, iy, iz, count, i, ii, &
       ind_x(10), ind_y(10),ind_z(10)

  ! XSF scalar-field header
  WRITE(ounit,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
  WRITE(ounit,'(a)') '3D_PWSCF'
  WRITE(ounit,'(a)') 'DATAGRID_3D_UNKNOWN'

  ! number of points in each direction
  WRITE(ounit,*) nr1+1, nr2+1, nr3+1
  ! origin
  WRITE(ounit,'(3f10.6)') 0.0, 0.0, 0.0
  ! 1st spanning (=lattice) vector
  WRITE(ounit,'(3f10.6)') (0.529177d0*alat*at(i,1),i=1,3) ! in ANSTROMS
  ! 2nd spanning (=lattice) vector
  WRITE(ounit,'(3f10.6)') (0.529177d0*alat*at(i,2),i=1,3)
  ! 3rd spanning (=lattice) vector
  WRITE(ounit,'(3f10.6)') (0.529177d0*alat*at(i,3),i=1,3)

  count=0
  DO i3=0,nr3
     !iz = mod(i3,nr3)
     iz = MOD(i3,nr3) + 1

     DO i2=0,nr2
        !iy = mod(i2,nr2)
        iy = MOD(i2,nr2) + 1

        DO i1=0,nr1
           !ix = mod(i1,nr1)
           ix = MOD(i1,nr1) + 1

           !ii = (1+ix) + iy*nrx1 + iz*nrx1*nrx2
           IF (count.LT.6) THEN
              count = count + 1
              !ind(count) = ii
           ELSE
              WRITE(ounit,'(6e13.5)') &
                   (rho(ind_x(i),ind_y(i),ind_z(i)),i=1,6)
              count=1
              !ind(count) = ii
           ENDIF
           ind_x(count) = ix
           ind_y(count) = iy
           ind_z(count) = iz
        ENDDO
     ENDDO
  ENDDO
  WRITE(ounit,'(6e13.5:)') (rho(ind_x(i),ind_y(i),ind_z(i)),i=1,count)
  WRITE(ounit,'(a)') 'END_DATAGRID_3D'
  WRITE(ounit,'(a)') 'END_BLOCK_DATAGRID_3D'
  RETURN
END SUBROUTINE xsf_fast_datagrid_3d
