!
! Copyright (C) 2002-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define ZERO ( 0.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE mix_rho( rhocout, rhocin, becout, becin, nsout, nsin, &
                   alphamix, dr2, tr2_min, iter, n_iter, file_extension, conv )
  !----------------------------------------------------------------------------
  !
  ! ... Modified Broyden's method for charge density mixing
  ! ...         D.D. Johnson PRB 38, 12807 (1988)
  !
  ! ... On output: the mixed density is in rhocin, mixed augmentation
  ! ...            channel occ. is in becin
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE ions_base,            ONLY : nat
  USE gvect,                ONLY : ngm, nl, nlm, gstart
  USE ldaU,                 ONLY : lda_plus_u, Hubbard_lmax
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : imix, tr2
  USE wvfct,                ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE parser,               ONLY : find_free_unit
  USE cell_base,            ONLY : omega
  !
  !!PAW]
  USE uspp_param,           ONLY : nhm
  USE grid_paw_variables,   ONLY : okpaw
  !!PAW]
  !
  IMPLICIT NONE
  !
  ! ... First the I/O variable
  !
  CHARACTER(LEN=256) :: &
    file_extension          !  (in) I/O filename extension for mixing history
                            !  if absent everything is kept in memory
  INTEGER :: &
    iter,                  &!  (in)  counter of the number of iterations
    n_iter                  !  (in)  numb. of iterations used in mixing
  COMPLEX(DP) :: &
    rhocin (ngm,nspin), &
    rhocout(ngm,nspin)
  REAL(DP) :: &
    nsout(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat), &!
    nsin(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat),  &!
    alphamix,              &! (in) mixing factor
    dr2                     ! (out) the estimated errr on the energy
  REAL (DP) :: &
    tr2_min       ! estimated error from diagonalization. If the estimated scf 
                  ! error is smaller than this, exit: a more accurate 
                  ! diagonalization is needed
  LOGICAL :: &
    conv                    ! (out) if true the convergence has been reached
  INTEGER, PARAMETER :: &
    maxmix = 25             ! max number of iterations for charge mixing
  !
  !!PAW[ 
  REAL(DP) :: &
    becin (nhm*(nhm+1)/2,nat,nspin), &
    becout(nhm*(nhm+1)/2,nat,nspin)
  !!PAW] 
  !
  ! ... Here the local variables
  !
  INTEGER ::    &
    iunmix,        &! I/O unit number of charge density file
    iunmix2,       &! I/O unit number of ns file
    iunit,         &! counter on I/O unit numbers
    iter_used,     &! actual number of iterations used
    ipos,          &! index of the present iteration
    inext,         &! index of the next iteration
    i, j,          &! counters on number of iterations
    is,            &! counter on spin component
    ig,            &! counter on G-vectors
    iwork(maxmix), &! dummy array used as output by libr. routines
    info,          &! flag saying if the exec. of libr. routines was ok
    ldim            ! 2 * Hubbard_lmax + 1
  COMPLEX(DP), ALLOCATABLE :: &
    rhoinsave(:,:),     &! rhoinsave(ngm,nspin): work space
    rhoutsave(:,:),     &! rhoutsave(ngm,nspin): work space
    nsinsave(:,:,:,:),  &!
    nsoutsave(:,:,:,:)   !
  REAL(DP) :: &
    betamix(maxmix,maxmix), &
    gamma0,                 &
    work(maxmix),           &
    charge
  LOGICAL :: &
    savetofile,   &! save intermediate steps on file "prefix"."file_extension"
    exst           ! if true the file exists
  !!PAW[  And this for PAW  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(DP), ALLOCATABLE :: &
    becinsave(:,:,:),      &
    becoutsave(:,:,:)
  INTEGER ::    &        
    iunmix3        ! I/O unit number of bec file
  !!PAW] 
  !
  ! ... saved variables and arrays
  !
  INTEGER, SAVE :: &
    mixrho_iter = 0    ! history of mixing
  COMPLEX(DP), ALLOCATABLE, SAVE :: &
    df(:,:,:),        &! information from preceding iterations
    dv(:,:,:)          !     "  "       "     "        "  "
  REAL(DP), ALLOCATABLE, SAVE :: &
    df_ns(:,:,:,:,:), &! idem 
    dv_ns(:,:,:,:,:)   ! idem
  !!PAW[  And this for PAW  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(DP), ALLOCATABLE, SAVE :: &
    df_bec(:,:,:,:), &! idem 
    dv_bec(:,:,:,:)   ! idem
  !!PAW] 
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: rho_dot_product, rho1_dot_product, ns_dot_product
  !
  !
  CALL start_clock( 'mix_rho' )
  !
  mixrho_iter = iter
  !
  IF ( n_iter > maxmix ) CALL errore( 'mix_rho', 'n_iter too big', 1 )
  !
  IF ( lda_plus_u ) ldim = 2 * Hubbard_lmax + 1
  !
  savetofile = ( file_extension /= ' ' )
  !
  rhocout(:,:) = rhocout(:,:) - rhocin(:,:)
  !
  !!PAW[
  IF ( okpaw ) becout(:,:,:) = becout(:,:,:) - becin(:,:,:)
  !!PAW]
  !
  IF ( lda_plus_u ) nsout(:,:,:,:) = nsout(:,:,:,:) - nsin(:,:,:,:)
  !
  !!PAW[
  !!PAW: dr2 contains the contribution of 1-center charges too
  dr2 = rho_dot_product( rhocout, rhocout ) + ns_dot_product( nsout, nsout )+ &
        rho1_dot_product ( becout, becout )
  !!PAW]
  !
  conv = ( dr2 < tr2 )
  !
  ! ... if the self-consistency error (dr2) is smaller than the estimated 
  ! ... error due to diagonalization (tr2_min), exit and leave rhocin and 
  ! ... rhocout unchanged
  !
  IF ( conv .OR. dr2 < tr2_min ) THEN
     !
     IF ( ALLOCATED ( df ) ) DEALLOCATE ( df )
     IF ( ALLOCATED ( dv ) ) DEALLOCATE ( dv )
     IF ( lda_plus_u .AND. ALLOCATED ( df_ns ) ) DEALLOCATE ( df_ns )
     IF ( lda_plus_u .AND. ALLOCATED ( dv_ns ) ) DEALLOCATE ( dv_ns )
     !!PAW[
     IF ( okpaw .AND. ALLOCATED ( df_bec ) ) DEALLOCATE ( df_bec )
     IF ( okpaw .AND. ALLOCATED ( dv_bec ) ) DEALLOCATE ( dv_bec )
     !!PAW]  
     !
     rhocout(:,:) = rhocout(:,:) + rhocin(:,:)
     !
     CALL stop_clock( 'mix_rho' )
     !
     RETURN
     !
  END IF
  !
  !
  IF ( savetofile ) THEN
     !
     iunmix = find_free_unit()
     CALL diropn( iunmix, file_extension, ( 2 * ngm * nspin ), exst )
     !
     IF ( lda_plus_u ) then
        iunmix2 = find_free_unit()
        CALL diropn( iunmix2, TRIM( file_extension ) // '.ns', &
                     ( ldim * ldim * nspin * nat ), exst )
     END IF
     !
     !!PAW[
     IF ( okpaw ) then
        iunmix3 = find_free_unit()
        CALL diropn( iunmix3, TRIM( file_extension ) // '.bec', &
                     ( (nhm * (nhm + 1)/2) * nat * nspin ), exst )
     END IF
     !!PAW]
     !
     IF ( mixrho_iter > 1 .AND. .NOT. exst ) THEN
        !
        CALL infomsg( 'mix_rho','file not found, restarting', -1 )
        mixrho_iter = 1
        !
     END IF
     !
  END IF
  !
  IF ( savetofile .OR. mixrho_iter == 1 ) THEN

     IF ( .NOT. ALLOCATED( df ) ) ALLOCATE( df( ngm, nspin, n_iter ) )
     IF ( .NOT. ALLOCATED( dv ) ) ALLOCATE( dv( ngm, nspin, n_iter ) )
     !
     IF ( lda_plus_u ) THEN
        !
        IF ( .NOT. ALLOCATED( df_ns ) ) &
             ALLOCATE( df_ns( ldim, ldim, nspin, nat, n_iter ) )
        IF ( .NOT. ALLOCATED( dv_ns ) ) &
             ALLOCATE( dv_ns( ldim, ldim, nspin, nat, n_iter ) )
        !
     END IF
     !!PAW[
     IF ( okpaw ) THEN
        IF ( .NOT. ALLOCATED( df_bec ) ) &
             ALLOCATE( df_bec( nhm * (nhm + 1)/2, nat, nspin, n_iter ) )
        IF ( .NOT. ALLOCATED( dv_bec ) ) &
             ALLOCATE( dv_bec( nhm * (nhm + 1)/2, nat, nspin, n_iter ) )
        !
     END IF
     !!PAW]
     !
  END IF
  !
  ! ... iter_used = mixrho_iter-1  if  mixrho_iter <= n_iter
  ! ... iter_used = n_iter         if  mixrho_iter >  n_iter
  !
  iter_used = MIN( ( mixrho_iter - 1 ), n_iter )
  !
  ! ... ipos is the position in which results from the present iteration
  ! ... are stored. ipos=mixrho_iter-1 until ipos=n_iter, then back to 1,2,...
  !
  ipos = mixrho_iter - 1 - ( ( mixrho_iter - 2 ) / n_iter ) * n_iter
  !
  IF ( mixrho_iter > 1 ) THEN
     !
     IF ( savetofile ) THEN
        !
        CALL davcio( df(1,1,ipos), 2*ngm*nspin, iunmix, 1, -1 )
        CALL davcio( dv(1,1,ipos), 2*ngm*nspin, iunmix, 2, -1 )
        !
        IF ( lda_plus_u ) THEN
           !
           CALL davcio( df_ns(1,1,1,1,ipos),ldim*ldim*nspin*nat,iunmix2,1,-1 )
           CALL davcio( dv_ns(1,1,1,1,ipos),ldim*ldim*nspin*nat,iunmix2,2,-1 )
           !
        END IF
        !
        !!PAW[
        IF ( okpaw ) THEN
           !
           CALL davcio( df_bec(1,1,1,ipos),(nhm*(nhm+1)/2)*nat*nspin,iunmix3,1,-1 )
           CALL davcio( dv_bec(1,1,1,ipos),(nhm*(nhm+1)/2)*nat*nspin,iunmix3,2,-1 )
           !
        END IF
        !!PAW]
        !
     END IF
     !
     df(:,:,ipos) = df(:,:,ipos) - rhocout(:,:)
     dv(:,:,ipos) = dv(:,:,ipos) - rhocin (:,:)
     !
     IF ( lda_plus_u ) THEN
        !
        df_ns(:,:,:,:,ipos) = df_ns(:,:,:,:,ipos) - nsout
        dv_ns(:,:,:,:,ipos) = dv_ns(:,:,:,:,ipos) - nsin
        !
     END IF
     !
     !!PAW[
     IF ( okpaw ) THEN
           !
           df_bec(:,:,:,ipos) = df_bec(:,:,:,ipos) - becout
           dv_bec(:,:,:,ipos) = dv_bec(:,:,:,ipos) - becin
           !
     END IF
     !!PAW]
     !
  END IF
  !
  IF ( savetofile ) THEN
     !
     DO i = 1, iter_used
        !
        IF ( i /= ipos ) THEN
           !
           CALL davcio( df(1,1,i), 2*ngm*nspin, iunmix, 2*i+1, -1 )
           CALL davcio( dv(1,1,i), 2*ngm*nspin, iunmix, 2*i+2, -1 )
           !
           IF ( lda_plus_u ) THEN
              !
              CALL davcio(df_ns(1,1,1,1,i),ldim*ldim*nspin*nat,iunmix2,2*i+1,-1)
              CALL davcio(dv_ns(1,1,1,1,i),ldim*ldim*nspin*nat,iunmix2,2*i+2,-1)
              !
           END IF
           !
           !!PAW[
           IF ( okpaw ) THEN
              !
              CALL davcio(df_bec(1,1,1,i),(nhm*(nhm+1)/2)*nat*nspin,iunmix3,2*i+1,-1)
              CALL davcio(dv_bec(1,1,1,i),(nhm*(nhm+1)/2)*nat*nspin,iunmix3,2*i+2,-1)
              !
           END IF
           !!PAW]
           !
        END IF
        !
     END DO
     !
     CALL davcio( rhocout, 2*ngm*nspin, iunmix, 1, 1 )
     CALL davcio( rhocin , 2*ngm*nspin, iunmix, 2, 1 )
     !
     IF ( mixrho_iter > 1 ) THEN
        !
        CALL davcio( df(1,1,ipos), 2*ngm*nspin, iunmix, 2*ipos+1, 1 )
        CALL davcio( dv(1,1,ipos), 2*ngm*nspin, iunmix, 2*ipos+2, 1 )
        !
     END IF
     !
     IF ( lda_plus_u ) THEN
        !
        CALL davcio( nsout, ldim*ldim*nspin*nat, iunmix2, 1, 1 )
        CALL davcio( nsin , ldim*ldim*nspin*nat, iunmix2, 2, 1 )
        !
        IF ( mixrho_iter > 1 ) THEN
           !
           CALL davcio( df_ns(1,1,1,1,ipos), ldim*ldim*nspin*nat, &
                        iunmix2, 2*ipos+1, 1 )
           CALL davcio( dv_ns(1,1,1,1,ipos), ldim*ldim*nspin*nat, &
                        iunmix2, 2*ipos+2, 1 )
        END IF
        !
     END IF
     !
     !!PAW[
     IF ( okpaw ) THEN
        !
        CALL davcio( becout, (nhm*(nhm+1)/2)*nat*nspin, iunmix3, 1, 1 )
        CALL davcio( becin , (nhm*(nhm+1)/2)*nat*nspin, iunmix3, 2, 1 )
        !
        IF ( mixrho_iter > 1 ) THEN
           !
           CALL davcio( df_bec(1,1,1,ipos), (nhm*(nhm+1)/2)*nat*nspin, &
                        iunmix3, 2*ipos+1, 1 )
           CALL davcio( dv_bec(1,1,1,ipos), (nhm*(nhm+1)/2)*nat*nspin, &
                        iunmix3, 2*ipos+2, 1 )
        END IF
        !
     END IF
     !!PAW]
     !
  ELSE
     !
     !
     ALLOCATE( rhoinsave(ngm,nspin), rhoutsave(ngm,nspin) )
     rhoinsave = rhocin (:,:)
     rhoutsave = rhocout(:,:)
     !
     IF ( lda_plus_u ) THEN
        !
        ALLOCATE( nsinsave (ldim,ldim,nspin,nat), &
                  nsoutsave(ldim,ldim,nspin,nat) )
        nsinsave  = nsin
        nsoutsave = nsout
        !
     END IF
     !
     !!PAW[
     IF ( okpaw ) THEN
        !
        ALLOCATE( becinsave (nhm*(nhm+1)/2,nat,nspin), &
                  becoutsave(nhm*(nhm+1)/2,nat,nspin) )
        becinsave  = becin
        becoutsave = becout
        !
     END IF
     !!PAW]
     !
  END IF
  !
  DO i = 1, iter_used
     !
     DO j = i, iter_used
        !
        betamix(i,j) = rho_dot_product( df(1,1,j), df(1,1,i) ) 
        !
        IF ( lda_plus_u ) &
           betamix(i,j) = betamix(i,j) + &
                          ns_dot_product( df_ns(1,1,1,1,j), df_ns(1,1,1,1,i) )
        !
        !!PAW[
        !!PAW: add the contribution of 1-center charges
        IF ( okpaw ) &
           betamix(i,j) = betamix(i,j) + &
                          rho1_dot_product( df_bec(1,1,1,j), df_bec(1,1,1,i) )
        !!PAW]
        !
     END DO
     !
  END DO
  !
  CALL DSYTRF( 'U', iter_used, betamix, maxmix, iwork, work, maxmix, info )
  CALL errore( 'broyden', 'factorization', info )
  !
  CALL DSYTRI( 'U', iter_used, betamix, maxmix, iwork, work, info )
  CALL errore( 'broyden', 'DSYTRI', info )
  !
  FORALL( i = 1 : iter_used, &
          j = 1 : iter_used, j > i ) betamix(j,i) = betamix(i,j)
  !
  DO i = 1, iter_used
     !
     work(i) = rho_dot_product( df(1,1,i), rhocout )
     !
     IF ( lda_plus_u ) &
        work(i) = work(i) + ns_dot_product( df_ns(1,1,1,1,i), nsout )
     !
     !!PAW[
     !!PAW: add the contribution of 1-center charges
     IF ( okpaw ) &
        work(i) = work(i) + rho1_dot_product( df_bec(1,1,1,i), becout )
     !!PAW]
     !
  END DO
  !
  DO i = 1, iter_used
     !
     gamma0 = SUM( betamix(1:iter_used,i) * work(1:iter_used) )
     !
     rhocin (:,:) = rhocin (:,:) - gamma0 * dv(:,:,i)
     rhocout(:,:) = rhocout(:,:) - gamma0 * df(:,:,i)
     !
     IF ( lda_plus_u ) THEN
        !
        nsin  = nsin  - gamma0 * dv_ns(:,:,:,:,i)
        nsout = nsout - gamma0 * df_ns(:,:,:,:,i)
        !
     END IF
     !
     !!PAW[
     IF ( okpaw ) THEN
        !
        becin  = becin  - gamma0 * dv_bec(:,:,:,i)
        becout = becout - gamma0 * df_bec(:,:,:,i)
        !
     END IF
     !!PAW]
     !
  END DO
  !
  ! ... auxiliary vectors dv and df not needed anymore
  !
  IF ( savetofile ) THEN
     !
     IF ( lda_plus_u ) THEN
        !
        CLOSE( iunmix2, STATUS = 'KEEP' )
        !
        DEALLOCATE( df_ns, dv_ns )
        !
     END IF
     !
     !!PAW[
     IF ( okpaw ) THEN
        !
        CLOSE( iunmix3, STATUS = 'KEEP' )
        !
        DEALLOCATE( df_bec, dv_bec )
        !
     END IF
     !!PAW]
     !
     CLOSE( iunmix, STATUS = 'KEEP' )
     !
     DEALLOCATE( df, dv )
     !
  ELSE
     !
     inext = mixrho_iter - ( ( mixrho_iter - 1 ) / n_iter ) * n_iter
     !
     IF ( lda_plus_u ) THEN
        !
        df_ns(:,:,:,:,inext) = nsoutsave
        dv_ns(:,:,:,:,inext) = nsinsave
        !
        DEALLOCATE( nsinsave, nsoutsave )
        !
     END IF
     !
     !!PAW[
     IF ( okpaw ) THEN
        !
        df_bec(:,:,:,inext) = becoutsave
        dv_bec(:,:,:,inext) = becinsave
        !
        DEALLOCATE( becinsave, becoutsave )
        !
     END IF
     !!PAW]
     !
     df(:,:,inext) = rhoutsave(:,:)
     dv(:,:,inext) = rhoinsave(:,:)
     !
     DEALLOCATE( rhoinsave, rhoutsave )
     !
  END IF
  !
  ! ... preconditioning the new search direction
  !
  IF ( imix == 1 ) THEN
     !
     CALL approx_screening( rhocout )
     !
  ELSE IF ( imix == 2 ) THEN
     !
     CALL approx_screening2( rhocout, rhocin )
     !
  END IF
  !
  ! ... set new trial density
  !
  rhocin = rhocin + alphamix * rhocout
  !
  IF ( lda_plus_u ) nsin = nsin + alphamix * nsout
  !
  !!PAW[
  IF ( okpaw ) becin = becin + alphamix * becout
  !!PAW]
  !
  CALL stop_clock( 'mix_rho' )
  !
  RETURN
  !
END SUBROUTINE mix_rho
!
!----------------------------------------------------------------------------
FUNCTION rho_dot_product( rho1, rho2 ) RESULT( rho_ddot )
  !----------------------------------------------------------------------------
  !
  ! ... calculates 4pi/G^2*rho1(-G)*rho2(G) = V1_Hartree(-G)*rho2(G)
  ! ... used as an estimate of the self-consistency error on the energy
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, tpi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, gstart, ngm
  USE lsda_mod,      ONLY : nspin
  USE wvfct,         ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  COMPLEX(DP), INTENT(IN) :: rho1(ngm,nspin), rho2(ngm,nspin)
  !
  REAL(DP)                :: rho_ddot
  !
  ! ... and the local variables
  !
  REAL(DP) :: fac   ! a multiplicative factors
  INTEGER       :: ig, gi
  !
  !
  gi = gstart
  !
  fac = e2 * fpi / tpiba2
  !
  IF ( nspin == 1 ) THEN
     !
     rho_ddot = fac * SUM( DBLE( CONJG( rho1(gi:,1) ) * &
                                        rho2(gi:,1) ) / gg(gi:) )
     !
     IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     ! ... first the charge
     !
     rho_ddot = fac * SUM( DBLE( CONJG( rho1(gi:,1)+rho1(gi:,2) ) * &
                                      ( rho2(gi:,1)+rho2(gi:,2) ) ) / gg(gi:) )
     !
     IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
     !
     ! ... then the magnetization
     !
     fac = e2 * fpi / tpi**2  ! lambda = 1 a.u.
     !
     ! ... G=0 term
     !
     IF ( gstart == 2 ) THEN
        !
        rho_ddot = rho_ddot + &
                   fac * DBLE( CONJG( rho1(1,1) - rho1(1,2) ) * &
                                    ( rho2(1,1) - rho2(1,2) ) )
        !
     END IF
     !
     IF ( gamma_only ) fac = 2.D0 * fac
     !
     rho_ddot = rho_ddot + fac * &
                SUM( DBLE( CONJG( rho1(gi:,1) - rho1(gi:,2) ) * &
                                ( rho2(gi:,1) - rho2(gi:,2) ) ) )
     !

  ELSE IF ( nspin == 4 ) THEN

     rho_ddot = fac * SUM( DBLE( CONJG( rho1(gi:,1) ) * &
                                        rho2(gi:,1) ) / gg(gi:) )
     !
     IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
     !
     fac = e2*fpi / (tpi**2)  ! lambda=1 a.u.
     !
     IF ( gstart == 2 ) THEN
        !
        rho_ddot = rho_ddot + &
             fac *  DBLE( CONJG( rho1(1,2))*(rho2(1,2) )) + &
                    DBLE( CONJG( rho1(1,3))*(rho2(1,3) )) + &
                    DBLE( CONJG( rho1(1,4))*(rho2(1,4) ))
        !
     END IF
     !
     IF ( gamma_only ) fac = 2.D0 * fac
     !
     rho_ddot = rho_ddot + &
          fac * SUM ( DBLE( CONJG( rho1(gi:,2))*(rho2(gi:,2) )) + &
                      DBLE( CONJG( rho1(gi:,3))*(rho2(gi:,3) )) + &
                      DBLE( CONJG( rho1(gi:,4))*(rho2(gi:,4) )) )
     !
  END IF
  !
  rho_ddot = rho_ddot * omega * 0.5D0
  !
  CALL reduce( 1, rho_ddot )
  !
  RETURN
  !
END FUNCTION rho_dot_product
!
!
!----------------------------------------------------------------------------
FUNCTION rho1_dot_product( bec1, bec2 ) RESULT( rho1_ddot )
  !----------------------------------------------------------------------------
  !
  ! ... calculates 4pi/G^2*rho1(-G)*rho2(G) = V1_Hartree(-G)*rho2(G),
  ! ... where rho1 and rho2 are 1-center charges (AE and PS)
  ! ... input variables are the augmentation channel occupations
  ! ... used as an estimate of the self-consistency error on the energy
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, tpi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : ngm, nl, nlm, gg, g, gstart
  USE lsda_mod,      ONLY : nspin
  USE wvfct,         ONLY : gamma_only
  !  
  USE grid_paw_variables, ONLY : pp, okpaw, prodp, prodpt, prod0p, prod0pt
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,         ONLY : lmaxq, nh, nhm 
  USE wvfct,              ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  REAL(DP), INTENT(IN) :: &
     bec1(nhm*(nhm+1)/2,nat,nspin), &
     bec2(nhm*(nhm+1)/2,nat,nspin)
  !
  REAL(DP)                :: rho1_ddot
  !
  ! ... and the local variables
  !
  REAL(DP) :: fac   ! a multiplicative factor
  !
  INTEGER :: gi, ig, na, nt, ih, jh, ijh, ijh2, is
  ! counters
  !
  COMPLEX(DP), POINTER :: prodp_(:,:,:), prod0p_(:,:,:)
  INTEGER :: i_what
  !
  rho1_ddot = 0.D0   
  !
  IF ( .NOT. okpaw ) RETURN
  !
  gi = gstart
  !
  fac = e2 * fpi / tpiba2
  !  
  whattodo: DO i_what=1, 2
     !
     NULLIFY(prodp_,prod0p_)
     IF (i_what==1) THEN
        prodp_ => prodp
        prod0p_ => prod0p
     ELSE IF (i_what==2) THEN
        prodp_ => prodpt
        prod0p_ => prod0pt
     END IF
     !
     DO ijh = 1, nhm*(nhm+1)/2
        !
        DO ijh2 = 1, nhm*(nhm+1)/2
           !
           DO na = 1, nat
              !
              nt = ityp (na)
              IF ( nspin == 1 ) THEN
                 !
                 rho1_ddot = rho1_ddot + DBLE((-1)**(i_what-1)) * fac * &
                 bec1(ijh,na,1) * prodp_(ijh, ijh2, nt) * bec2(ijh2,na,1)
                 !
!!$              gamma_only case not yet implemented
!!$              IF ( gamma_only ) rho1_ddot = 2.D0 * rho1_ddot
                 !
              ELSE IF ( nspin == 2 ) THEN
                 !
                 ! ... first the charge
                 !
                 rho1_ddot = rho1_ddot + DBLE((-1)**(i_what-1)) * fac * &
                  (bec1(ijh,na,1)+bec1(ijh,na,2)) * prodp_(ijh,ijh2,nt) * &
                  (bec2(ijh2,na,1)+bec2(ijh2,na,2))
                 ! 
!!$              IF ( gamma_only ) rho1_ddot = 2.D0 * rho1_ddot
                 !
                 ! ... then the magnetization
                 !
                 fac = e2 * fpi / tpi**2  ! lambda = 1 a.u.
                 !
                 ! ... G=0 term
                 !
                 IF ( gstart == 2 ) THEN
                    !
                    rho1_ddot = rho1_ddot + DBLE((-1)**(i_what-1)) * fac * &
                    (bec1(ijh,na,1)-bec1(ijh,na,2)) * prod0p_(ijh,ijh2,nt) * &
                    (bec2(ijh2,na,1)-bec2(ijh2,na,2))
                    !
                 END IF
                 !
!!$              IF ( gamma_only ) fac = 2.D0 * fac
                 !
                 rho1_ddot = rho1_ddot + DBLE((-1)**(i_what-1)) * fac * &
                   (bec1(ijh,na,1)-bec1(ijh,na,2)) * prodp_(ijh,ijh2,nt) * &
                   (bec2(ijh2,na,1)-bec2(ijh2,na,2))
                 !
!!$           non-collinear case not yet implemented
!!$           ELSE IF ( nspin == 4 ) THEN
!!$              !
!!$              rho1_ddot = rho1_ddot + fac * DBLE((-1)**(i_what-1)) * 
!!$              bec1(ijh,na,1) * prodp_(ijh, ijh2, nt) * bec2(ijh2,na,1)
!!$              !
!!$              IF ( gamma_only ) rho1_ddot = 2.D0 * rho1_ddot
!!$              !
!!$              fac = e2*fpi / (tpi**2)  ! lambda=1 a.u.
!!$              !
!!$              IF ( gstart == 2 ) THEN
!!$                 !
!!$                 rho1_ddot = rho1_ddot + DBLE((-1)**(i_what-1)) * fac * &
!!$                 bec1(ijh,na,2) * prod0p_(ijh,ijh2,nt) * bec2(ijh2,na,2)+ &
!!$                 bec1(ijh,na,3) * prod0p_(ijh,ijh2,nt) * bec2(ijh2,na,3)+ &
!!$                 bec1(ijh,na,4) * prod0p_(ijh,ijh2,nt) * bec2(ijh2,na,4)
!!$                 !
!!$              END IF
!!$              !
!!$              IF ( gamma_only ) fac = 2.D0 * fac
!!$              !
!!$              rho1_ddot = rho1_ddot + DBLE((-1)**(i_what-1)) * fac * &
!!$              ( bec1(ijh,na,2) * prodp_(ijh,ijh2,nt) * bec2(ijh2,na,2)+ &
!!$                bec1(ijh,na,3) * prodp_(ijh,ijh2,nt) * bec2(ijh2,na,3)+ &
!!$                bec1(ijh,na,4) * prodp_(ijh,ijh2,nt) * bec2(ijh2,na,4) )
!!$                 !
              END IF
              !
           END DO
           !
        END DO
        !
     END DO
     !
  END DO whattodo
  !
  rho1_ddot = rho1_ddot * omega * 0.5D0
  !
  CALL reduce( 1, rho1_ddot )
  !
  RETURN
  !
END FUNCTION rho1_dot_product
!
!----------------------------------------------------------------------------
FUNCTION ns_dot_product( ns1, ns2 )
  !----------------------------------------------------------------------------
  !
  ! ... calculates U/2 \sum_i ns1(i)*ns2(i)
  ! ... used as an estimate of the self-consistency error on the 
  ! ... LDA+U correction to the energy
  !
  USE kinds,      ONLY : DP
  USE ldaU,       ONLY : lda_plus_u, Hubbard_lmax, Hubbard_l, Hubbard_U, &
                         Hubbard_alpha
  USE ions_base,  ONLY : nat, ityp
  USE lsda_mod,   ONLY : nspin
  !
  IMPLICIT NONE  
  !
  ! ... I/O variables
  !
  REAL(DP) :: &
    ns_dot_product                           ! (out) the function value
  REAL(DP), INTENT(IN) :: &
    ns1(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat), &
    ns2(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)  ! (in) the two ns 
  !
  ! ... and the local variables
  !
  INTEGER :: na, nt, m1, m2
  !
  !
  ns_dot_product = 0.D0
  !
  IF ( .NOT. lda_plus_u ) RETURN
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     IF ( Hubbard_U(nt) /= 0.D0 .OR. Hubbard_alpha(nt) /= 0.D0 ) THEN
        !
        m1 = 2 * Hubbard_l(nt) + 1
        m2 = 2 * Hubbard_l(nt) + 1
        !
        ns_dot_product = ns_dot_product + 0.5D0 * Hubbard_U(nt) * &
                         SUM( ns1(:m1,:m2,:nspin,na) * ns2(:m1,:m2,:nspin,na) )
        !
     END IF
     !
  END DO
  !
  IF ( nspin == 1 ) ns_dot_product = 2.D0 * ns_dot_product
  !
  RETURN
  !
END FUNCTION ns_dot_product
!
!----------------------------------------------------------------------------
SUBROUTINE approx_screening( drho )
  !----------------------------------------------------------------------------
  !
  ! ... apply an average TF preconditioning to drho
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : e2, pi, fpi
  USE cell_base,      ONLY : omega, tpiba2
  USE gvect,          ONLY : gstart, gg, ngm
  USE klist,          ONLY : nelec
  USE lsda_mod,       ONLY : nspin
  !
  IMPLICIT NONE  
  !
  ! ... I/O variables
  !
  COMPLEX(DP) :: &
    drho(ngm,nspin) ! (in/out)
  !
  ! ... and the local variables
  !
  REAL(DP) :: rrho, rmag, rs, agg0
  INTEGER       :: ig
  !
  !
  rs = ( 3.D0 * omega / fpi / nelec )**( 1.D0 / 3.D0 )
  !
  agg0 = ( 12.D0 / pi )**( 2.D0 / 3.D0 ) / tpiba2 / rs
  !
  IF ( nspin == 1 .OR. nspin == 4 ) THEN
     !
     drho(:,1) =  drho(:,1) * gg(:) / ( gg(:) + agg0 )
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     DO ig = 1, ngm
        !
        rrho = ( drho(ig,1) + drho(ig,2) ) * gg(ig) / ( gg(ig) + agg0 )
        rmag = ( drho(ig,1) - drho(ig,2) )
        !
        drho(ig,1) =  0.5D0 * ( rrho + rmag )
        drho(ig,2) =  0.5D0 * ( rrho - rmag )
        !
     END DO
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE approx_screening
!
!----------------------------------------------------------------------------
SUBROUTINE approx_screening2( drho, rhobest )
  !----------------------------------------------------------------------------
  !
  ! ... apply a local-density dependent TF preconditioning to drho
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : e2, pi, tpi, fpi, eps8, eps32
  USE cell_base,            ONLY : omega, tpiba2
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   nl, nlm, ngm, gg
  USE klist,                ONLY : nelec
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  COMPLEX(DP) :: &
    drho(ngm,nspin), rhobest(ngm,nspin)
  !
  ! ... and the local variables
  !
  INTEGER, PARAMETER :: mmx = 12
  !
  INTEGER :: &
    iwork(mmx), i, j, m, info, nspin_save
  REAL(DP) :: &
    rs, min_rs, max_rs, avg_rsm1, target, dr2_best, ccc, cbest, l2smooth
  REAL(DP) :: &
    aa(mmx,mmx), invaa(mmx,mmx), bb(mmx), work(mmx), vec(mmx), agg0
  COMPLEX(DP), ALLOCATABLE :: &
    v(:,:),     &! v(ngm,mmx)
    w(:,:),     &! w(ngm,mmx)
    dv(:),      &! dv(ngm)
    vbest(:),   &! vbest(ngm)
    wbest(:)     ! wbest(ngm)
  REAL(DP), ALLOCATABLE :: &
    alpha(:)     ! alpha(nrxx)
  !
  COMPLEX(DP) :: rrho, rmag
  INTEGER          :: ir, ig
  !
  REAL(DP), PARAMETER :: one_third = 1.D0 / 3.D0
  !
  REAL(DP), EXTERNAL :: rho_dot_product
  !
  !
  IF ( nspin == 4 ) CALL errore ('approx_screening2', &
       'noncolinear not implemented',1)
  !
  IF ( nspin == 2 ) THEN
     !
     DO ig = 1, ngm
        !
        rrho = drho(ig,1) + drho(ig,2)
        rmag = drho(ig,1) - drho(ig,2)
        !        
        drho(ig,1) = rrho
        drho(ig,2) = rmag
        !
     END DO
     !
  END IF
  !
  nspin_save = nspin
  nspin      = 1
  target     = 0.D0
  !
  IF ( gg(1) < eps8 ) drho(1,1) = ZERO
  !
  ALLOCATE( alpha(nrxx), v(ngm,mmx), w(ngm,mmx), &
            dv(ngm), vbest(ngm), wbest(ngm) )
  !
  v(:,:)   = ZERO
  w(:,:)   = ZERO
  dv(:)    = ZERO
  vbest(:) = ZERO
  wbest(:) = ZERO
  !
  ! ... calculate alpha from density smoothed with a lambda=0 a.u.
  !
  l2smooth = 0.D0
  !
  psic(:) = ZERO
  !
  IF ( nspin == 1 ) THEN
     !
     psic(nl(:)) = rhobest(:,1) * EXP( - 0.5D0 * l2smooth * tpiba2 * gg(:) )
     !
  ELSE
     !
     psic(nl(:)) = ( rhobest(:,1) + rhobest(:,2) ) * &
                   EXP( - 0.5D0 * l2smooth * tpiba2 * gg(:) )
     !
  END IF
  !
  IF ( gamma_only ) psic(nlm(:)) = CONJG( psic(nl(:)) )
  !
  CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
  !
  alpha(:) = DBLE( psic(:) )
  !
  min_rs   = ( 3.D0 * omega / fpi / nelec )**one_third
  max_rs   = min_rs
  avg_rsm1 = 0.D0
  !
  DO ir = 1, nrxx
     !
     alpha(ir) = ABS( alpha(ir) )
     !
     IF ( alpha(ir) > eps32 ) THEN
        !
        rs        = ( 3.D0 / fpi / alpha(ir) )**one_third
        min_rs    = MIN( min_rs, rs )
        avg_rsm1  = avg_rsm1 + 1.D0 / rs
        max_rs    = MAX( max_rs, rs )
        alpha(ir) = rs
        !
     END IF   
     !
  END DO
  !
  CALL reduce( 1, avg_rsm1 )
  !
  CALL extreme( min_rs, -1 )
  CALL extreme( max_rs, +1 )
  !
  alpha = 3.D0 * ( tpi / 3.D0 )**( 5.D0 / 3.D0 ) * alpha
  !
  avg_rsm1 = ( nr1 * nr2 * nr3 ) / avg_rsm1
  rs       = ( 3.D0 * omega / fpi / nelec )**one_third
  agg0     = ( 12.D0 / pi )**( 2.D0 / 3.D0 ) / tpiba2 / avg_rsm1
  !
  ! ... calculate deltaV and the first correction vector
  !
  psic(:) = ZERO
  !
  psic(nl(:)) = drho(:,1)
  !
  IF ( gamma_only ) psic(nlm(:)) = CONJG( psic(nl(:)) )
  !
  CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
  !
  psic(:) = psic(:) * alpha(:)
  !
  CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1 )
  !
  dv(:) = psic(nl(:)) * gg(:) * tpiba2
  v(:,1)= psic(nl(:)) * gg(:) / ( gg(:) + agg0 )
  !
  m       = 1
  ccc     = rho_dot_product( dv, dv )
  aa(:,:) = 0.D0
  bb(:)   = 0.D0
  !
  repeat_loop: DO
     !
     ! ... generate the vector w
     !     
     w(:,m) = fpi * e2 * v(:,m)
     !
     psic(:) = ZERO
     !
     psic(nl(:)) = v(:,m)
     !
     IF ( gamma_only ) psic(nlm(:)) = CONJG( psic(nl(:)) )
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
     !
     psic(:) = psic(:) * alpha(:)
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1 )
     !
     w(:,m) = w(:,m) + gg(:) * tpiba2 * psic(nl(:))
     !
     ! ... build the linear system
     !
     DO i = 1, m
        !
        aa(i,m) = rho_dot_product( w(1,i), w(1,m) )
        !
        aa(m,i) = aa(i,m)
        !
     END DO
     !
     bb(m) = rho_dot_product( w(1,m), dv )
     !
     ! ... solve it -> vec
     !
     invaa = aa
     !
     CALL DSYTRF( 'U', m, invaa, mmx, iwork, work, mmx, info )
     CALL errore( 'broyden', 'factorization', info )
     !
     CALL DSYTRI( 'U', m, invaa, mmx, iwork, work, info )
     CALL errore( 'broyden', 'DSYTRI', info )
     !     
     FORALL( i = 1 : m, j = 1 : m, j > i ) invaa(j,i) = invaa(i,j)
     !
     FORALL( i = 1 : m ) vec(i) = SUM( invaa(i,:) * bb(:) )
     !
     vbest(:) = ZERO
     wbest(:) = dv(:)
     !
     DO i = 1, m
        !
        vbest = vbest + vec(i) * v(:,i)
        wbest = wbest - vec(i) * w(:,i)
        !
     END DO
     !
     cbest = ccc - SUM( bb(:) * vec(:) )
     !
     dr2_best = rho_dot_product( wbest, wbest )
     !
     IF ( target == 0.D0 ) target = 1.D-6 * dr2_best
     !
     IF ( dr2_best < target ) THEN
        !
        drho(:,1) = vbest(:)
        !
        nspin = nspin_save
        !
        IF ( nspin == 2 ) THEN
           !
           DO ig = 1, ngm
              !
              rrho = drho(ig,1)
              rmag = drho(ig,2)
              !
              drho(ig,1) = 0.5D0 * ( rrho + rmag )
              drho(ig,2) = 0.5D0 * ( rrho - rmag )
              !
           END DO
           !
        END IF
        !
        DEALLOCATE( alpha, v, w, dv, vbest, wbest )
        !
        EXIT repeat_loop
        !
     ELSE IF ( m >= mmx ) THEN
        !
        m = 1
        !
        v(:,m)  = vbest(:)
        aa(:,:) = 0.D0
        bb(:)   = 0.D0
        !
        CYCLE repeat_loop
        !
     END IF
     !
     m = m + 1
     !
     v(:,m) = wbest(:) / ( gg(:) + agg0 )
     !
  END DO repeat_loop
  !
  RETURN
  !
END SUBROUTINE approx_screening2
