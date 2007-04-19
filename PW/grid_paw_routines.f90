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
    !
    ALLOCATE(prodp(nhm*(nhm+1)/2,nhm*(nhm+1)/2,ntyp))
    ALLOCATE(prodpt(nhm*(nhm+1)/2,nhm*(nhm+1)/2,ntyp))
    ALLOCATE(prod0p(nhm*(nhm+1)/2,nhm*(nhm+1)/2,ntyp))
    ALLOCATE(prod0pt(nhm*(nhm+1)/2,nhm*(nhm+1)/2,ntyp))
    !
    ALLOCATE(rho1(nrxx, nspin, nat))
    ALLOCATE(rho1t(nrxx, nspin, nat))
!!! No more needed since ptfunc already contains the augmentation charge qfunc
!!! ALLOCATE(rho1h(nrxx, nspin, nat))
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
         kkbeta, nqf, nqlc, lll, jjj, lmaxkb, nh, tvanp, nhm, tvanp
    USE spin_orb,   ONLY : lspinorb, rot_ylm, fcoef
    !
    USE grid_paw_variables, ONLY: tpawp, pfunc, ptfunc, pp, ppt, prad, ptrad, &
         int_r2pfunc, int_r2ptfunc, okpaw, which_paw_augfun, gifpaw, &
         augmom
    !
    IMPLICIT NONE
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
         lk, mk, vk, kh, lh, sph_ind, nnbrx, ll
    COMPLEX(DP) :: coeff, qgm(1)
    REAL(DP) :: spinor, ji, jk
    CHARACTER*20 :: fname 

!!$call start_clock ('init_us_1')
    IF (.NOT.okpaw) RETURN
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
                            prad_ (iq, nmb, l + 1, nt) =     &
                            prad_ (iq, nmb, l + 1, nt) * prefr 
!! NEW-AUG !!
                            IF ((i_what.EQ.2) .AND. &
                            (which_paw_augfun/='DEFAULT')) THEN
                               !
                               prad_ (iq, nmb, l + 1, nt) =      &
                               prad_ (iq, nmb, l + 1, nt) +  &
                               qrad(iq, nmb, l + 1, nt) 
                               !
                            ENDIF
!! NEW-AUG !!
                         ENDIF
                      ENDDO
                   ENDDO
                   ! igl
                ENDDO
                ! l
             ENDDO
!             prad_ (:, :, :, nt) = prad_ (:, :, :, nt)*prefr
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
          IF (tpawp (nt) ) THEN
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
!! NEW-AUG !!
                   IF ( (i_what.EQ.2) .AND. &
                   (which_paw_augfun/='DEFAULT') ) THEN
                   aux2(1:msh(nt)) = aux2(1:msh(nt)) + augmom(nb, mb, 0, nt)* &
                      gifpaw (1:msh(nt), 1) * r(1:msh(nt),nt) * r(1:msh(nt),nt)
                   ENDIF
!! NEW-AUG !!
                   CALL simpson (msh(nt), aux2, rab(1,nt), &
                        int_r2pfunc_(nb,mb,nt))
#if defined __DEBUG_NEWD_PAW_GRID
!                   WRITE (200000+10000*i_what+100*nb+mb,'(3e20.10)') &
!                        (aux2(ir),rab(ir,nt),r(ir,nt),ir=1,mesh(nt))
                   WRITE (777,'(i6,e15.7)') 200000+10000*i_what+100*nb+mb, &
                        int_r2pfunc_(nb,mb,nt)
#endif
                END DO
             END DO
          END IF
       END DO
    END DO whattodo
!STOP
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
    DO nt = 1, ntyp
       DO ih = 1, nbrx
          DO ll = 1, lmaxq
          !
!          WRITE (10000+ll*1000+nt*100+ih,'(e15.6)') ptrad(1:nqxq,ih,ll,nt)
!!$          write (20000+nt*100+ih,'(e15.8)') prad(1:nqxq,ih,1,nt)
!!$          write (30000+nt*100+ih,'(e15.8)') ptrad(1:nqxq,ih,1,nt)
          !CLOSE (UNIT=10000+nt*100+ih,STATUS='KEEP')
          !
          END DO
       END DO
    END DO

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

!! NEW-AUG !!

!#define __DEBUG_COMPUTE_AUGFUN
! compute augmentation functions (Bessel or Gauss)
! this routine is called by init_us_1
  SUBROUTINE compute_pawaugfun
    !
    USE kinds,              ONLY : DP
    USE atom,               ONLY: zmesh, mesh, r, rab, dx, msh
    USE pseud,              ONLY: lmax
    USE parameters,         ONLY : lmaxx, nbrx, lqmax, ndmx
    USE constants,          ONLY : fpi
    USE ions_base,          ONLY : ntyp => nsp
    USE cell_base,          ONLY : omega, tpiba
    USE gvect,              ONLY : g, gg
    USE lsda_mod,           ONLY : nspin
    USE us,                 ONLY : nqxq, dq, nqx, tab, qrad
    USE uspp,               ONLY : qq, qq_so
    USE uspp_param,         ONLY : lmaxq, betar, qfunc, qfcoef, rinner, &
                    nbeta, kkbeta, nqf, nqlc, lll, jjj, lmaxkb, nh, tvanp, nhm
    USE spin_orb,   ONLY : lspinorb, rot_ylm, fcoef
    !
    USE grid_paw_variables, ONLY: tpawp, okpaw, nraug, augmom, which_paw_augfun, gifpaw
    !
    IMPLICIT NONE
    ! 
    !     here a few local variables
    !
    INTEGER :: nt, ih, jh, nb, mb, nmb, l, m, ir, iq, is, startq, &
         lastq, ilast, ndm, iok
    ! various counters
    REAL(DP), ALLOCATABLE :: aux (:)
    ! various work space
    REAL(DP) :: prefr, pref, q, qi, raux
    ! the prefactor of the q functions
    ! the prefactor of the beta functions
    ! the modulus of g for each shell
    ! q-point grid for interpolation
    ! work space
    REAL(DP) :: qc(2), xc(2), b1(2), b2(2)
    REAL(DP), ALLOCATABLE :: j1(:,:)
    INTEGER  :: nc          ! index Bessel funct
    ! 
    !
    IF (.NOT.okpaw) RETURN
    !
    !    Initialization of the variables
    !
    ndm = MAXVAL (msh(1:ntyp))
    ALLOCATE (j1 (ndm,2))    
    ALLOCATE (aux ( ndm))      
    !
    j1(:,:) = 0.d0
    gifpaw(:,:) = 0.d0
    !
    !   here we calculate the gi functions
    !
    CALL divide (nqxq, startq, lastq)
    DO nt = 1, ntyp
       !
       IF (tpawp (nt) ) THEN
          !
          SELECT CASE (which_paw_augfun)
          !
          CASE ('BESSEL')
             ! 
             DO l = 0, nqlc (nt) - 1
                !
                CALL find_aug_qi(qc,nraug(nt),l,2,nt,iok)
                IF (iok.ne.0) & 
                   CALL errore('compute_augfun', 'problems with find_aug_qi',1)                
                DO nc = 1, 2
                   !
                   CALL sph_bes(nraug(nt)+5,r(1,nt),qc(nc),l,j1(1,nc))
                   b1(nc) = j1(nraug(nt),nc)
                   aux(1:nraug(nt)) = j1(1:nraug(nt),nc) * &
                                      r(1:nraug(nt),nt)**(l+2)
                   CALL simpson (nraug(nt), aux, rab(1, nt), b2(nc) )
                   !
                ENDDO
                xc(1) = b1(2) / (b1(2) * b2(1) - b1(1) * b2(2))
                xc(2) = -1._dp * b1(1) * xc(1) / b1(2)
                gifpaw(1:nraug(nt),l+1) = ( xc(1) * j1(1:nraug(nt),1) + &
                    xc(2) * j1(1:nraug(nt),2) ) * r(1:nraug(nt),nt) ** 2
                gifpaw((nraug(nt)+1):msh(nt),l+1) = 0._dp 
                ! l
             ENDDO
             !
          CASE ('GAUSS')
             ! 
             DO l = 0, nqlc (nt) - 1
                !
                gifpaw(1:msh(nt),l+1) = EXP(-(r(1:msh(nt),nt)**2) /  & 
                                 (2._dp*0.25_dp**2)) * r(1:msh(nt),nt) **2  
                CALL simpson (nraug(nt), gifpaw(1:nraug(nt),l+1), rab(1, nt), raux )
                IF (ABS(raux) .LT. 1.d-8) THEN
                   CALL errore('ld1_to_paw','norm of augmentation function too small',nb*100+mb)
                END IF
                gifpaw(1:msh(nt),l+1) = gifpaw(1:msh(nt),l+1) / raux 
                !
             ENDDO
             ! 
          END SELECT
          !
       ENDIF
       ! ntyp
    ENDDO
    DEALLOCATE (j1)    
       DEALLOCATE (aux)    
    !
  END SUBROUTINE compute_pawaugfun

! in this subroutine the zeros of the 1st derivative of j_l are tabulated
! (they are NOT explicitly evaluated)
 !--------------------------------------------------------------------------
SUBROUTINE find_aug_qi(qc,ik,lam,ncn,is,iok)
  !--------------------------------------------------------------------------
  !
  !      This routine finds two values of q such that the
  !      functions f_l have a derivative equal to
  !      0 at the point ik
  !  
  USE kinds,      ONLY : DP
  USE atom,       ONLY : r
  IMPLICIT NONE

  INTEGER ::      &
       ik,    & ! input: the point corresponding to rcaug
       lam,   & ! input: the angular momentum
       ncn,   & ! input: the number of qi to compute
       is,    & ! input: the type of the pseudo
       iok      ! output: if 0 the calculation in this routine is ok

  REAL (DP) :: &
       qc(ncn)  ! output: the values of qi


  REAL (DP) ::   &
       zeroderjl (2,7) ! first two zeros of the first derivative of 
                       ! spherical Bessel function j_l for l = 0,...,6

  INTEGER ::    &
       nc       ! counter on the q found

  data zeroderjl / 0.0_dp,                 4.4934094579614_dp, &
                   2.0815759780862_dp,     5.9403699890844_dp, &
                   3.3420936578747_dp,     7.2899322987026_dp, &
                   4.5140996477983_dp,     8.5837549433127_dp, &
                   5.6467036213923_dp,     9.8404460168549_dp, &
                   6.7564563311363_dp,    11.0702068269176_dp, &
                   7.8510776799611_dp,    12.2793339053177_dp  /
  iok=0
  IF (ncn.gt.2) &
       CALL errore('find_aug_qi','ncn is too large',1)

  IF (lam.gt.6) &
       CALL errore('find_aug_qi','l not programmed',1)
  !
  !    fix deltaq and the maximum step number
  !
  DO nc = 1, ncn
     !
     qc(nc) = zeroderjl (nc, lam + 1) / r (ik, is)
     !
  ENDDO
  RETURN
END SUBROUTINE find_aug_qi

!! NEW-AUG !! 

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
  ALLOCATE (qmod( ngm))    
  ALLOCATE (qgm( ngm))     
  ALLOCATE (pft(nhm*(nhm+1)/2,ngm,ntyp))    
  ALLOCATE (ylmk0( ngm, lmaxq * lmaxq)) 
  !
  CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  DO ig = 1, ngm
     qmod (ig) = SQRT (gg (ig) )
  ENDDO
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
        !
        IF (tpawp (nt) ) THEN
           ! 
           ijh = 0
           DO ih = 1, nh (nt)
              !
              DO jh = ih, nh (nt)
                 !
                 ijh = ijh + 1
                 CALL pvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0, prad_, &
                      SIZE(prad_,1),SIZE(prad_,2),SIZE(prad_,3),SIZE(prad_,4))
                 DO ig = 1, ngm
                    !
                    pft(ijh,ig,nt) =  pft( ijh,ig,nt) + qgm(ig)
                    !
                 ENDDO  
                 !
              ENDDO
              !
           ENDDO  
           !
        ENDIF
        !
     ENDDO
     !
     DO ijh = 1, nhm*(nhm+1)/2
        !
        DO ijh2 = ijh, nhm*(nhm+1)/2
           !
           DO nt = 1, ntyp
              !
              prodp_(ijh, ijh2, nt) = prodp_(ijh, ijh2, nt) + & 
              SUM( DBLE( CONJG( pft(ijh,gi:,nt) ) * pft(ijh2,gi:,nt) )/ &
              gg(gi:) )
              prod0p_(ijh, ijh2, nt) = prod0p_(ijh, ijh2, nt) + &
              DBLE( CONJG( pft(ijh,1,nt) ) * pft(ijh2,1,nt) )
              prodp_(ijh2, ijh, nt) = CONJG( prodp_(ijh, ijh2, nt) )
              prod0p_(ijh2, ijh, nt) = CONJG( prod0p_(ijh, ijh2, nt) )
              !
           ENDDO
           !
        ENDDO
        !
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
    USE grid_paw_variables, ONLY: prad, ptrad, pp, tpawp, okpaw, which_paw_augfun
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
    INTEGER :: ig, na, nt, ih, jh, ijh, is
    ! counters
    !
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

    IF (.NOT.okpaw) RETURN
  
    ALLOCATE (aux ( ngm, nspin, nat))    
    ALLOCATE (qmod( ngm))    
    ALLOCATE (qgm( ngm))    
    ALLOCATE (ylmk0( ngm, lmaxq * lmaxq))    
    !  
    CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
    DO ig = 1, ngm
       qmod (ig) = SQRT (gg (ig) )
    ENDDO

    whattodo: DO i_what=1, 2
       NULLIFY(prad_,rho1_)
       IF (i_what==1) THEN
          prad_ => prad
          rho1_ => rho1new
       ELSE IF (i_what==2) THEN
          prad_ => ptrad
          rho1_ => rho1tnew
!!! No more needed since ptfunc already contains the augmentation charge qfunc
!!!    ELSE IF (i_what==3) THEN
!!!       prad_ => qrad
!!!       rho1_ => rho1newh
       END IF
       aux (:,:,:) = (0.d0, 0.d0)
       charge = 0.d0

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
                      IF (ityp (na) .EQ.nt) THEN
                         DO is = 1, nspin
                            DO ig = 1, ngm
                               skk = eigts1 (ig1 (ig), na) * &
                                    eigts2 (ig2 (ig), na) * &
                                    eigts3 (ig3 (ig), na)
                               aux(ig,is,na)=aux(ig,is,na) + qgm(ig)*skk*becnew(ijh,na,is)
                            ENDDO
                         ENDDO
                      ENDIF
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
    DEALLOCATE (ylmk0)
    DEALLOCATE (qgm)
    DEALLOCATE (qmod)
    DEALLOCATE (aux)

!!! No more needed since ptfunc already contains the augmentation charge qfunc
!!! ! Add compensation charge to rho1tnew:
!!! rho1tnew(:,:,:) = rho1tnew(:,:,:) + rho1newh(:,:,:)

  END SUBROUTINE compute_onecenter_charges


  ! Analogous to PW/v_of_rho.f90
  ! + evaluation of the spheropole: Int dr r^2 rho(r)
!#define __DEBUG_ONECENTER_POTENTIALS
  SUBROUTINE compute_onecenter_potentials (inp_rho1, inp_rho1t)
    USE kinds,            ONLY : DP
    USE cell_base,        ONLY : at, alat, omega
    USE ions_base,          ONLY: tau, atm, ityp, ntyp=>nsp
    USE ions_base,        ONLY : nat
    USE lsda_mod,         ONLY : nspin
    USE gvect,            ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3,&
         nrxx, nl, g, gg
    !
    USE grid_paw_variables, ONLY: vr1, vr1t, & !rho1, rho1t, &
         int_r2pfunc, int_r2ptfunc, tpawp, okpaw, ehart1, etxc1, vtxc1, &
         ehart1t, etxc1t, vtxc1t, aerho_core, psrho_core
    USE uspp, ONLY: indv, becsum, nhtolm
    USE uspp_param, ONLY: nh
    USE constants, ONLY: PI
    IMPLICIT NONE
    !
    REAL(DP), TARGET, INTENT(OUT) :: &
         inp_rho1(nrxx, nspin, nat), inp_rho1t(nrxx,nspin,nat)
    !
    REAL(DP), POINTER :: etxc1_(:), vtxc1_(:), ehart1_(:)
    !REAL(DP) :: rho_core(nrxx)
    REAL(DP) :: charge, alpha, spheropole
    INTEGER :: na, ih, jh, ijh, nb, mb, is, nt, i
    !
    REAL(DP), POINTER :: rho1_(:,:,:), rho_core_(:,:), vr1_(:,:,:), int_r2pfunc_(:,:,:)
    INTEGER :: i_what
    !
    IF (.NOT.okpaw) RETURN
    !
    CALL infomsg ('compute_onecenter_potentials','alpha set manually',-1)
    alpha = 0.1_DP
    !
    !CALL infomsg ('compute_onecenter_potentials','rho_core set to zero',-1)
    !rho_core = 0._DP
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
             CALL v_xc( rho1_, rho_core_, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                  nrxx, nl, ngm, g, nspin, alat, omega, etxc1_(na), vtxc1_(na), vr1_(:,:,na) )
#if defined __DEBUG_NEWD_PAW_GRID
!             IF (i_what==2) THEN
!                 WRITE (456,'(8f15.7)') (vr1_(i,1,1), i=1,nrxx)
!                 PRINT '(A,2f20.10)', 'XC', etxc1_(na), vtxc1_(na)
!             ENDIF
#endif
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
#if defined __DEBUG_NEWD_PAW_GRID
!                               WRITE (4,'(3i5,f20.10,2i5,f20.10)') ih, jh, ijh, &
!                                    becsum(ijh,na,is), nb, mb, int_r2pfunc_(nb,mb,nt)
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
                  nl, ngm, gg, gstart, nspin, alat, omega, ehart1_(na), charge, vr1_(:,:,na), &
                  alpha, spheropole, na)
#if defined __DEBUG_NEWD_PAW_GRID
!             IF (i_what==2) THEN
!                 WRITE (789,'(8f15.7)') (vr1_(i,1,1), i=1,nrxx)
!                 PRINT '(A,2f20.10)', 'HARTREE', ehart1_(na)
!             ENDIF
#endif
          END IF
       END DO

       !CALL xsf_struct (alat, at, nat, tau, atm, ityp, 93000+i_what)
       !CALL xsf_fast_datagrid_3d &
       !     (vr1_(:,1,1), nr1, nr2, nr3, nrx1, nrx2, nrx3, at, alat, 93000+i_what)

    END DO whattodo
#if defined __DEBUG_ONECENTER_POTENTIALS
    WRITE (93000,'(i5,4f20.10)') (na,rho1(na,1,1),vr1(na,1,1),rho1t(na,1,1),vr1t(na,1,1),na=1,nr1)
    WRITE (93000,*)
#endif

!!$#if defined __DEBUG_ONECENTER_POTENTIALS
!!$    STOP 'STOP __DEBUG_ONECENTER_POTENTIALS'
!!$#endif

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

  ! From PW/set_rhoc.f90
!#define __DEBUG_SET_PAW_RHOC
  SUBROUTINE set_paw_rhoc
    !-----------------------------------------------------------------------
    !
    !    This routine computes the core charge on the real space 3D mesh
    !
    !
    USE io_global, ONLY : stdout
    USE kinds,     ONLY : DP
    USE atom,      ONLY : rho_atc, numeric, msh, r, rab, nlcc
    USE ions_base, ONLY : nat, ityp, ntyp => nsp
    USE cell_base, ONLY : omega, tpiba2, alat
    !USE ener,      ONLY : etxcc
    USE gvect,     ONLY : ngm, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                          nrxx, nl, nlm, ngl, gl, igtongl,      &
                          eigts1, eigts2, eigts3, ig1, ig2, ig3
    USE pseud,     ONLY : a_nlcc, b_nlcc, alpha_nlcc
    !USE scf,       ONLY : rho_core
    USE vlocal,    ONLY : strf
    USE wvfct,     ONLY : gamma_only
    !
    USE grid_paw_variables, ONLY : aerho_atc, psrho_atc, aerho_core, psrho_core
    !
    IMPLICIT NONE

    !
    ! NEW
    !
    REAL(DP), POINTER :: rho_atc_(:,:)
    REAL(DP), POINTER :: rho_core_(:,:)
    !
    INTEGER :: i_what
    !
    !     here the local variables
    !
    REAL(DP), PARAMETER :: eps = 1.d-10

    COMPLEX(DP) , ALLOCATABLE :: aux (:)
    ! used for the fft of the core charge

    COMPLEX(DP) :: skk
    ! exp(I*G*R_j) for atom j

    REAL(DP) , ALLOCATABLE ::  rhocg(:)
    ! the radial fourier trasform
    REAL(DP) ::  rhoima, rhoneg, rhorea
    ! used to check the core charge
    REAL(DP) ::  vtxcc
    ! dummy xc energy term
    REAL(DP) , ALLOCATABLE ::  dum(:,:)
    ! dummy array containing rho=0
  
    INTEGER :: ir, nt, na, ng
    ! counter on mesh points
    ! counter on atomic types
    ! counter on atoms
    ! counter on g vectors

    !etxcc = 0.d0
!!$    DO nt = 1, ntyp
!!$       IF (nlcc (nt) ) GOTO 10
!!$    ENDDO
!!$    aerho_core(:,:) = 0.d0
!!$    psrho_core(:,:) = 0.d0
!!$    RETURN
!!$
!!$10  CONTINUE
    !
    ALLOCATE (aux( nrxx))    
    ALLOCATE (rhocg( ngl))
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
       aux (:) = 0.d0
       !
       !    the sum is on atom types
       !
       typ_loop: DO nt = 1, ntyp
          !
          !     drhoc compute the radial fourier transform for each shell of g vec
          !
          CALL drhoc (ngl, gl, omega, tpiba2, numeric (nt), a_nlcc (nt), &
               b_nlcc (nt), alpha_nlcc (nt), msh (nt), r (1, nt), rab (1, nt), &
               rho_atc_ (1, nt), rhocg)
          !
          IF ((i_what==2).AND.(.NOT.nlcc(nt))) THEN
             DO na = 1, nat
                IF (ityp (na) .EQ.nt) THEN
                   rho_core_(:,na)=0.d0
                END IF
             END DO
             CYCLE typ_loop
          END IF
          !
          at_loop: DO na = 1, nat
             at_if: IF (ityp (na) .EQ.nt) THEN
                DO ng = 1, ngm
                   skk = eigts1 (ig1 (ng), na) * &
                        eigts2 (ig2 (ng), na) * &
                        eigts3 (ig3 (ng), na)
                   aux(nl(ng)) = aux(nl(ng)) + skk * rhocg(igtongl(ng))
                ENDDO
!!$       IF (gamma_only) THEN
!!$          DO ng = 1, ngm
!!$             aux(nlm(ng)) = CONJG(aux(nl (ng)))
!!$          END DO
!!$       END IF
                !
                !   the core charge in real space
                !
                CALL cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
                !
                !    test on the charge and computation of the core energy
                !
                rhoneg = 0.d0
                rhoima = 0.d0
                DO ir = 1, nrxx
                   rhoneg = rhoneg + min (0.d0,  DBLE (aux (ir) ) )
                   rhoima = rhoima + abs (AIMAG (aux (ir) ) )
                   rho_core_(ir,na) =  DBLE (aux(ir))
                ENDDO
                rhoneg = rhoneg / (nr1 * nr2 * nr3)
                rhoima = rhoima / (nr1 * nr2 * nr3)
!!$#ifdef __PARA
!!$       CALL reduce (1, rhoneg)
!!$       CALL reduce (1, rhoima)
!!$#endif
                IF (rhoneg < -1.0d-6 .OR. rhoima > 1.0d-6) &
                     WRITE( stdout, '(/5x,"Check: atom number ", i12)') na
                WRITE( stdout, '(/5x,"       negative/imaginary core charge ", 2f12.6)')&
                     rhoneg, rhoima
                !
             END IF at_if
          END DO at_loop
       END DO typ_loop
    END DO whattodo
#if defined __DEBUG_SET_PAW_RHOC
    PRINT '(A)', 'Writing file fort.770'
    WRITE (770,'(3f20.10)') ((ir-1)*alat/nr1,aerho_core(ir,1),psrho_core(ir,1),ir=1,nr1)
    STOP 'STOP __DEBUG_SET_PAW_RHOC'
#endif
    !
    DEALLOCATE (rhocg)
    DEALLOCATE (aux)
    !
    RETURN

  9000 FORMAT (5x,'core-only xc energy         = ',f15.8,' ryd')

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
    USE grid_paw_variables, ONLY: prad, ptrad, tpawp, okpaw
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

    IF (.NOT.okpaw) RETURN

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
       IF (tpawp (nt) ) THEN
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


! Analogous to PW/newd.f90
!#define __DEBUG_NEWD_PAW_GRID
SUBROUTINE newd_paw_grid
  !
  !! NEW-AUG
  USE gvect,              ONLY : nrxx
  !! NEW-AUG
  USE grid_paw_variables,   ONLY : okpaw, prad, ptrad, vr1, vr1t, dpaw_ae, dpaw_ps, aevloc_r, psvloc_r
  !
  IMPLICIT NONE
  !
  !! NEW-AUG
  INTEGER :: i
  !! NEW-AUG
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
#if defined __DEBUG_NEWD_PAW_GRID
!  WRITE (123,'(8f15.7)') (vr1t(i,1,1), i=1,nrxx)
  !WRITE (456,'(8f15.7)') (psvloc_r(i,1), i=1,nrxx)
!  STOP 'STOP __DEBUG_NEWD_PAW_GRID'
#endif
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
                !! loop on atoms is now the outer one:
                !! DO na = 1, nat
                !! IF ( ityp(na) == nt ) THEN
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
                !! END IF
                !! END DO
                !
             END DO
             !
          END DO
          !
       END IF
       !
    END DO
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
!!$    USE uspp_param, ONLY : vloc_at
    USE ions_base,  ONLY : ntyp => nsp
    USE cell_base,  ONLY : omega, tpiba2
!!$    USE vlocal,     ONLY : vloc
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
!!$    open (FILE='../../aevloc_at.dat',UNIT=101,STATUS='unknown')
!!$    open (FILE='../../psvloc_at.dat',UNIT=102,STATUS='unknown')
!!$    do nt = 1, ntyp
!!$       write (101,'(e15.8)') aevloc_at(1:ndmx,nt)
!!$       write (102,'(e15.8)') psvloc_at(1:ndmx,nt)
!!$    end do
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
!!$    open (FILE='../../aevloc.dat',UNIT=103,STATUS='unknown')
!!$    open (FILE='../../psvloc.dat',UNIT=104,STATUS='unknown')
!!$    do nt = 1, ntyp
!!$       write (103,'(e15.8)') aevloc(1:ngl,nt)
!!$       write (104,'(e15.8)') psvloc(1:ngl,nt)
!!$    end do
!!$    ! Look at the fourier transforms
!!$    !CALL plot_augfun(2,2,1,1000,'QVAN')
!!$    !CALL plot_augfun(2,2,1,2000,'P1AE')
!!$    !CALL plot_augfun(2,2,1,3000,'P1PS')
!!$    STOP 'ionpot'
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
  real(DP), allocatable :: aux (:), aux1 (:)
  !  auxiliary variables
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
  USE extfield,  ONLY : tefield, dipfield, etotefield
  USE gvect,     ONLY : igtongl
!!$  USE scf,       ONLY : vltot
  USE vlocal,    ONLY : strf !!$ , vloc
  USE wvfct,     ONLY : gamma_only
  USE gvect,     ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nl, nlm, ngm
  USE scf,       ONLY : rho
  !
  USE gvect,     ONLY : eigts1, eigts2, eigts3, ig1, ig2, ig3
  USE grid_paw_variables, ONLY : aevloc, psvloc, aevloc_r, psvloc_r
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
        ! and not ... do nt = 1, ntyp
        do ng = 1, ngm
           !aux (nl(ng))=aux(nl(ng)) + vloc (igtongl (ng), nt) * strf(ng,nt)
           skk = eigts1 (ig1 (ng), na) * &
                eigts2 (ig2 (ng), na) * &
                eigts3 (ig3 (ng), na)
           aux (nl(ng))=aux(nl(ng)) + vloc_ (igtongl (ng), nt) * skk
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
        !  If required add an electric field to the local potential 
        !
        if (tefield.and.(.not.dipfield)) then
           STOP 'paw_grid_setlocal not implemented'
!!$     call add_efield(rho,vltot,etotefield,0)
!!$! NB rho is not actually used by add_efield in this case ...
!!$!    it should be fixed and removed from this routine
        endif
     end do na_loop
#if defined __DEBUG_PAW_GRID_SETLOCAL
     write (95000+i_what,'(f20.10)') (vltot_(ir,1),ir=1,nr1)
#endif
  end DO whattodo
#if defined __DEBUG_PAW_GRID_SETLOCAL
  STOP 'STOP __DEBUG_PAW_GRID_SETLOCAL'
#endif
  deallocate(aux)
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
  !
  delta_e_1 = 0.D0
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
  DO ipol = 1, nspin
     delta_e_1 = delta_e_1 - SUM( rho_(:,ipol,na) * vr_(:,ipol,na) )
  END DO
  !
  delta_e_1 = omega * delta_e_1 / ( nr1 * nr2 * nr3 )
  !
  CALL reduce( 1, delta_e_1 )
  !
  RETURN
  !
END FUNCTION delta_e_1
!


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
  !
  delta_e_1scf = 0.D0
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
  DO ipol = 1, nspin
     delta_e_1scf = delta_e_1scf - &
                  SUM( ( rhonew_(:,ipol,na) - rho_(:,ipol,na) ) * vr_(:,ipol,na) )
  END DO
  !
  delta_e_1scf = omega * delta_e_1scf / ( nr1 * nr2 * nr3 )
  !
  CALL reduce( 1, delta_e_1scf )
  !
  RETURN
  !
END FUNCTION delta_e_1scf


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
