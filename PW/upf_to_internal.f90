!#define __DEBUG_UPF_TO_INTERNAL
!
! Copyright (C) 2004 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module is USEd, for the time being, as an interface
! between the UPF pseudo type and the pseudo variables internal representation

!=----------------------------------------------------------------------------=!
  MODULE upf_to_internal
!=----------------------------------------------------------------------------=!

  IMPLICIT NONE
  SAVE

!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!

!
!---------------------------------------------------------------------
subroutine set_pseudo_upf (is, upf)
  !---------------------------------------------------------------------
  !
  !   set "is"-th pseudopotential using the Unified Pseudopotential Format
  !   dummy argument ( upf ) - convert and copy to internal PWscf variables
  !
  ! PWSCF modules
  !
  USE parameters, ONLY: ndmx
  USE atom,  ONLY: zmesh, mesh, msh, dx, r, rab, &
       chi, oc, nchi, lchi, jchi, rho_at, rho_atc, nlcc
  USE pseud, ONLY: lloc, lmax, zp
  USE uspp_param, ONLY: vloc_at, dion, betar, qqq, qfcoef, qfunc, nqf, nqlc, &
       rinner, nbeta, kkbeta, lll, jjj, psd, tvanp
  USE funct, ONLY: dft, which_dft, ismeta, ishybrid
  !
  USE ions_base, ONLY: zv
  USE spin_orb, ONLY: lspinorb
  USE pseudo_types
  !
  implicit none
  !
  real(DP), parameter :: rcut = 10.d0
  integer :: is, ir
  !
  !     Local variables
  !
  integer :: nb
  TYPE (pseudo_upf) :: upf
  !
  !
  zp(is)  = upf%zp
  psd (is)= upf%psd
  tvanp(is)=upf%tvanp
  nlcc(is) = upf%nlcc
  dft = upf%dft
  call which_dft( upf%dft )
  !
  IF ( ismeta ) &
    CALL errore( 'upf_to_internal', 'META-GGA not implemented in PWscf', 1 )
#if defined (EXX)
#else
  IF ( ishybrid ) &
    CALL errore( 'upf_to_internal', 'HYBRID XC not implemented in PWscf', 1 )
#endif
  !
  mesh(is) = upf%mesh
  IF ( mesh(is) > ndmx ) &
     CALL errore('upf_to_internal', 'too many grid points', 1)
  !
  nchi(is) = upf%nwfc
  lchi(1:upf%nwfc, is) = upf%lchi(1:upf%nwfc)
  oc(1:upf%nwfc, is) = upf%oc(1:upf%nwfc)
  chi(1:upf%mesh, 1:upf%nwfc, is) = upf%chi(1:upf%mesh, 1:upf%nwfc)
  !
  nbeta(is)= upf%nbeta
  kkbeta(is)=0
  do nb=1,upf%nbeta
     kkbeta(is)=max(upf%kkbeta(nb),kkbeta(is))
  end do
  betar(1:upf%mesh, 1:upf%nbeta, is) = upf%beta(1:upf%mesh, 1:upf%nbeta)
  dion(1:upf%nbeta, 1:upf%nbeta, is) = upf%dion(1:upf%nbeta, 1:upf%nbeta)
  !
  lmax(is) = upf%lmax
  nqlc(is) = upf%nqlc
  nqf (is) = upf%nqf
  lll(1:upf%nbeta,is) = upf%lll(1:upf%nbeta)
  rinner(1:upf%nqlc,is) = upf%rinner(1:upf%nqlc)
  qqq(1:upf%nbeta,1:upf%nbeta,is) = upf%qqq(1:upf%nbeta,1:upf%nbeta)
  qfunc (1:upf%mesh, 1:upf%nbeta, 1:upf%nbeta, is) = &
       upf%qfunc(1:upf%mesh,1:upf%nbeta,1:upf%nbeta)
  qfcoef(1:upf%nqf, 1:upf%nqlc, 1:upf%nbeta, 1:upf%nbeta, is ) = &
       upf%qfcoef( 1:upf%nqf, 1:upf%nqlc, 1:upf%nbeta, 1:upf%nbeta )
  !
  r  (1:upf%mesh, is) = upf%r  (1:upf%mesh)
  rab(1:upf%mesh, is) = upf%rab(1:upf%mesh)

  if (lspinorb.and..not.upf%has_so) &
     call infomsg ('upf_to_internal','At least one non s.o. pseudo', -1)
   
  lspinorb=lspinorb.and.upf%has_so
  if (upf%has_so) then
     jchi(1:upf%nwfc, is) = upf%jchi(1:upf%nwfc)
     jjj(1:upf%nbeta, is) = upf%jjj(1:upf%nbeta)
  else
     jchi(1:upf%nwfc, is) = 0.d0
     jjj(1:upf%nbeta, is) = 0.d0
  endif
  !
  if ( upf%nlcc) then
     rho_atc(1:upf%mesh, is) = upf%rho_atc(1:upf%mesh)
  else
     rho_atc(:,is) = 0.d0
  end if
  rho_at (1:upf%mesh, is) = upf%rho_at (1:upf%mesh)
  !!! TEMP
  lloc(is) = 0
  !!!
  vloc_at(1:upf%mesh,is) = upf%vloc(1:upf%mesh)

  do ir = 1, mesh (is)
    if (r (ir, is) .gt.rcut) then
        msh (is) = ir
        goto 5
    endif
  enddo
  msh (is) = mesh (is)
  !
  ! force msh to be odd for simpson integration
  !
5 msh (is) = 2 * ( (msh (is) + 1) / 2) - 1

  zv(is) = zp(is)  !!! maybe not needed: it is done in setup

#if defined __DEBUG_UPF_TO_INTERNAL
  call write_internals (10000,is)
#endif

end subroutine set_pseudo_upf

!---------------------------------------------------------------------
!#define __DO_NOT_CUTOFF_PAW_FUNC
subroutine set_pseudo_paw (is, pawset)
  !---------------------------------------------------------------------
  !
  !   set "is"-th pseudopotential using the Unified Pseudopotential Format
  !   dummy argument ( upf ) - convert and copy to internal PWscf variables
  !
  ! PWSCF modules
  !
  USE parameters, ONLY: ndmx
  USE atom,  ONLY: zmesh, mesh, msh, dx, r, rab, &
       chi, oc, nchi, lchi, jchi, rho_at, rho_atc, nlcc
  USE pseud, ONLY: lloc, lmax, zp
  USE uspp_param, ONLY: vloc_at, dion, betar, qqq, qfcoef, qfunc, nqf, nqlc, &
       rinner, nbeta, kkbeta, lll, jjj, psd, tvanp
  USE funct, ONLY: dft, which_dft, ismeta, ishybrid
  !
  USE ions_base, ONLY: zv
  USE spin_orb, ONLY: lspinorb
  USE pseudo_types
  USE constants, ONLY: FPI
  !
  USE grid_paw_variables, ONLY : tpawp, pfunc, ptfunc, aevloc_at, psvloc_at, &
                                 aerho_atc, psrho_atc, kdiff, &
                                 augmom, nraug, which_paw_augfun, r2 !!NEW-AUG
  USE grid_paw_routines, ONLY : step_f
  !
  implicit none
  !
  real(DP), parameter :: rcut = 10.d0
  integer :: is, ir
  !
  !     Local variables
  !
  integer :: nb
  TYPE (paw_t) :: pawset
  integer :: i,j, l, nrc, nrs
  real (DP) :: aux (ndmx), pow
  !
#if defined __DO_NOT_CUTOFF_PAW_FUNC
  PRINT '(A)', 'WARNING __DO_NOT_CUTOFF_PAW_FUNC'
#endif
  !
  ! Cutoffing: WARNING: arbitrary right now, for grid calculation
  pow = 1.d0
  nrs =  Count(pawset%r(1:pawset%mesh).le. (pawset%r(pawset%irc)*1.2d0))
  nrc =  Count(pawset%r(1:pawset%mesh).le. (pawset%r(pawset%irc)*1.8d0))
  PRINT *, 'PAW CUTOFF parameters'
  PRINT *, pawset%irc, pawset%r(pawset%irc)
  PRINT *, nrs, pawset%r(nrs)
  PRINT *, nrc, pawset%r(nrc)
  !
  zp(is)  = pawset%zval
  psd (is)= pawset%symbol
  tvanp(is)=.true.
  tpawp(is)=.true.
  nlcc(is) = pawset%nlcc
  dft = pawset%dft
  call which_dft( pawset%dft )
  !
  IF ( ismeta ) &
    CALL errore( 'upf_to_internal', 'META-GGA not implemented in PWscf', 1 )
#if defined (EXX)
#else
  IF ( ishybrid ) &
    CALL errore( 'upf_to_internal', 'HYBRID XC not implemented in PWscf', 1 )
#endif
  !
  mesh(is) = pawset%mesh
  IF ( mesh(is) > ndmx ) &
     CALL errore('upf_to_internal', 'too many grid points', 1)
  !
  ! ... Copy wavefunctions used for PAW construction.
  ! ... Copy also the unoccupied ones, e.g.
  ! ... corresponding to second energy for the same channel
  ! ... (necessary to set starting occupations correctly)
  !
  nchi(is)=0
  do i=1, pawset%nwfc
#if defined __DEBUG_UPF_TO_INTERNAL
     ! take only occupied wfcs (to have exactly the same as for US)
     if (pawset%oc(i)>0._dp) then
#endif
        nchi(is)=nchi(is)+1
        lchi(nchi(is),is)=pawset%l(i)
        oc(nchi(is),is)=MAX(pawset%oc(i),0._DP)
        chi(1:pawset%mesh, nchi(is), is) = pawset%pswfc(1:pawset%mesh, i)
#if defined __DEBUG_UPF_TO_INTERNAL
     end if
#endif
  end do
  !
  nbeta(is)= pawset%nwfc
  kkbeta(is)=0
  do nb=1,pawset%nwfc
     kkbeta(is)=max(pawset%ikk(nb),kkbeta(is))
  end do
  betar(1:pawset%mesh, 1:pawset%nwfc, is) = pawset%proj(1:pawset%mesh, 1:pawset%nwfc)
  dion(1:pawset%nwfc, 1:pawset%nwfc, is) = pawset%dion(1:pawset%nwfc, 1:pawset%nwfc)
  kdiff(1:pawset%nwfc, 1:pawset%nwfc, is) = pawset%kdiff(1:pawset%nwfc, 1:pawset%nwfc)

  !
  lmax(is) = pawset%lmax
  nqlc(is) = 2*pawset%lmax+1
  nqf (is) = 0                   !! no rinner, all numeric
  lll(1:pawset%nwfc,is) = pawset%l(1:pawset%nwfc)
  rinner(1:nqlc(is),is) = 0._dp  !! no rinner, all numeric
  !
  ! integral of augmentation charges vanishes for different values of l
  !
  do i = 1, pawset%nwfc
     do j = 1, pawset%nwfc
        if (pawset%l(i)==pawset%l(j)) then
           qqq(i,j,is) = pawset%augmom(i,j,0) !!gf spherical approximation
        else
           qqq(i,j,is) = 0._dp
        end if
     end do
  end do
  !! NEW-AUG !! 
  nraug(is) = pawset%irc
  dx (is) = pawset%dx
  r2 (1:pawset%mesh, is) = pawset%r2  (1:pawset%mesh)
  do l = 0, 2*pawset%lmax
     do i = 1, pawset%nwfc
        do j = 1, pawset%nwfc
           augmom(i,j,l,is) = pawset%augmom(i,j,l)
        end do 
     end do
  end do
  !! NEW-AUG !!
  qfunc (1:pawset%mesh, 1:pawset%nwfc, 1:pawset%nwfc, is) = &
       pawset%augfun(1:pawset%mesh,1:pawset%nwfc,1:pawset%nwfc,0)
  !
  do i=1,pawset%nwfc
     do j=1,pawset%nwfc
#if defined __DO_NOT_CUTOFF_PAW_FUNC
        pfunc (1:pawset%mesh, i, j, is) = &
             pawset%aewfc(1:pawset%mesh, i) * pawset%aewfc(1:pawset%mesh, j)
        ptfunc (1:pawset%mesh, i, j, is) = &
             pawset%pswfc(1:pawset%mesh, i) * pawset%pswfc(1:pawset%mesh, j)
#else
        aux(1:pawset%mesh) = pawset%aewfc(1:pawset%mesh, i) * &
             pawset%aewfc(1:pawset%mesh, j)
        CALL step_f( pfunc(1:pawset%mesh,i,j,is), aux(1:pawset%mesh), &
             pawset%r(1:pawset%mesh), nrs, nrc, pow, pawset%mesh)
        aux(1:pawset%mesh) = pawset%pswfc(1:pawset%mesh, i) * &
             pawset%pswfc(1:pawset%mesh, j)
        CALL step_f( ptfunc(1:pawset%mesh,i,j,is), aux(1:pawset%mesh), &
             pawset%r(1:pawset%mesh), nrs, nrc, pow, pawset%mesh)
#endif
     end do
  end do
  !
  ! ... Add augmentation charge to ptfunc already here.
  ! ... One should not need \tilde{n}^1 alone in any case.
  !
  !! NEW-AUG !! 
  !which_paw_augfun='DEFAULT'
  !which_paw_augfun='BESSEL'
  which_paw_augfun='GAUSS'
  IF (which_paw_augfun=='DEFAULT') THEN
     ptfunc(1:pawset%mesh, 1:pawset%nwfc, 1:pawset%nwfc, is) =         &
        ptfunc (1:pawset%mesh, 1:pawset%nwfc, 1:pawset%nwfc, is) +   &
        qfunc  (1:pawset%mesh, 1:pawset%nwfc, 1:pawset%nwfc, is)
  ENDIF
  !! NEW-AUG !! 
  !
  !
  ! nqf is always 0 for this PAW format
  ! qfcoef(1:pawset%nqf, 1:pawset%nqlc, 1:pawset%nwfc, 1:pawset%nwfc, is ) = 0._dp
  !
  r  (1:pawset%mesh, is) = pawset%r  (1:pawset%mesh)
  rab(1:pawset%mesh, is) = pawset%r  (1:pawset%mesh)*pawset%dx

  ! NO spin orbit PAW implemented right now (oct 2005)
!!$  if (lspinorb.and..not.pawset%has_so) &
!!$     call infomsg ('pawset_to_internal','At least one non s.o. pseudo', -1)
!!$   
!!$  lspinorb=lspinorb.and.pawset%has_so
!!$  if (pawset%has_so) then
!!$     jchi(1:pawset%nwfc, is) = pawset%jchi(1:pawset%nwfc)
!!$     jjj(1:pawset%nbeta, is) = pawset%jjj(1:pawset%nbeta)
!!$  else
  jchi(1:pawset%nwfc, is) = 0._dp
  jjj(1:pawset%nwfc, is) = 0._dp
!!$  endif
  !
  if ( pawset%nlcc) then
     rho_atc(1:pawset%mesh, is) = pawset%psccharge(1:pawset%mesh) &
          &                       / FPI / pawset%r2(1:pawset%mesh)
  else
     rho_atc(:,is) = 0.d0
  end if

  aerho_atc(1:pawset%mesh, is) = pawset%aeccharge(1:pawset%mesh) &
       &                         / FPI / pawset%r2(1:pawset%mesh)
  if ( pawset%nlcc) then
     psrho_atc(1:pawset%mesh, is) = pawset%psccharge(1:pawset%mesh) &
          &                         / FPI / pawset%r2(1:pawset%mesh)
  else
     psrho_atc(:,is) = 0._dp
  end if
  !
  rho_at (1:pawset%mesh, is) = pawset%pscharge(1:pawset%mesh)

  !!! TEMP       !!! this was already present in set_pseudo_upf. what does it mean?
  lloc(is) = 0
  !!!
  vloc_at(1:pawset%mesh,is) = pawset%psloc(1:pawset%mesh)
#if defined __DO_NOT_CUTOFF_PAW_FUNC
  aevloc_at(1:pawset%mesh,is) = pawset%aeloc(1:pawset%mesh)
  psvloc_at(1:pawset%mesh,is) = pawset%psloc(1:pawset%mesh)
#else
  aux(1:pawset%mesh) = pawset%aeloc(1:pawset%mesh)
  CALL step_f( aevloc_at(1:pawset%mesh,is), aux(1:pawset%mesh), &
       pawset%r(1:pawset%mesh), nrs, nrc, pow, pawset%mesh)
  aux(1:pawset%mesh) = pawset%psloc(1:pawset%mesh)
  CALL step_f( psvloc_at(1:pawset%mesh,is), aux(1:pawset%mesh), &
       pawset%r(1:pawset%mesh), nrs, nrc, pow, pawset%mesh)
#endif

  do ir = 1, mesh (is)
    if (r (ir, is) .gt.rcut) then
        msh (is) = ir
        goto 5
    endif
  enddo
  msh (is) = mesh (is)
  !
  ! force msh to be odd for simpson integration
  !
5 msh (is) = 2 * ( (msh (is) + 1) / 2) - 1

  zv(is) = zp(is)  !!! maybe not needed: it is done in setup

#if defined __DEBUG_UPF_TO_INTERNAL
  call write_internals (20000,is)
#endif

end subroutine set_pseudo_paw

#if defined __DEBUG_UPF_TO_INTERNAL
subroutine write_internals(un,is)
  !
  ! check reading of pseudopotential format
  !

  USE parameters, ONLY: ndmx
  USE atom,  ONLY: zmesh, mesh, msh, dx, r, rab, &
       chi, oc, nchi, lchi, jchi, rho_at, rho_atc, nlcc
  USE pseud, ONLY: lloc, lmax, zp
  USE uspp_param, ONLY: vloc_at, dion, betar, qqq, qfcoef, qfunc, nqf, nqlc, &
       rinner, nbeta, kkbeta, lll, jjj, psd, tvanp
  USE funct, ONLY: dft, which_dft, ismeta, ishybrid
  !
  USE ions_base, ONLY: zv
  USE spin_orb, ONLY: lspinorb

  implicit none
  integer :: un, is
  integer :: i,j

  write (un,*) "zp(is)",zp(is)
  write (un,*) "psd(is)",psd(is)
  write (un,*) "tvanp(is)",tvanp(is)
  write (un,*) "nlcc(is)",nlcc(is)
  write (un,*) "dft",dft
  write (un,*) "mesh(is)",mesh(is)

  write (un,*) "nchi(is)",nchi(is)
  write (un,*) "lchi(1:nchi(is), is)"
  write (un,*) lchi(1:nchi(is), is)
  write (un,*) "oc(1:nchi(is), is)"
  write (un,*) oc(1:nchi(is), is)
  do i=1,nchi(is)
     write (un+1000+i,'(f20.10)') chi(1:mesh(is), i, is)
  end do
  !
  write (un,*) "nbeta(is)",nbeta(is)
  write (un,*) "kkbeta(is)",kkbeta(is)
  do i=1,nbeta(is)
     write (un+2000+i,'(f20.10)') betar(1:mesh(is), i, is)
  end do
  write (un+100,'(f20.10)') dion(1:nbeta(is), 1:nbeta(is), is)

  write (un,*) "lmax(is)",lmax(is)
  write (un,*) "nqlc(is)",nqlc(is)
  write (un,*) "nqf",nqf
  write (un,*) "lll(1:nbeta(is),is)"
  write (un,*) lll(1:nbeta(is),is)
  write (un,*) "rinner(1:nqlc(is),is)"
  write (un,*) rinner(1:nqlc(is),is)
  write (un,*) "qqq(1:nbeta(is),1:nbeta(is),is)"
  write (un,*) qqq(1:nbeta(is),1:nbeta(is),is)
  do i=1,nbeta(is)
     do j=1,nbeta(is)
        write (un+3000+i*10+j,'(f20.10)') qfunc (1:mesh(is), i, j, is)
     end do
  end do
  write (un+200,*) qfcoef(1:nqf(is), 1:nqlc(is), 1:nbeta(is), 1:nbeta(is), is )
  !
  write (un+4000+1,'(f20.10)') r  (1:mesh(is), is)
  write (un+4000+2,'(f20.10)') rab(1:mesh(is), is)

  write (un,*) "lspinorb",lspinorb
  write (un+300+1,'(f20.10)') jchi(1:nchi(is), is)
  write (un+300+2,'(f20.10)') jjj(1:nbeta(is), is)

  write (un+5000,'(f20.10)') rho_atc(1:mesh(is), is)
  write (un+6000,'(f20.10)') rho_at (1:mesh(is), is)

  write (un,*) "lloc(is)",lloc(is)

  write (un+7000,'(f20.10)') vloc_at(1:mesh(is),is)

  write (un,*) "msh(is)",msh(is)

  write (un,*) "zv(is)",zv(is)

end subroutine write_internals
#endif

!=----------------------------------------------------------------------------=!
  END MODULE upf_to_internal
!=----------------------------------------------------------------------------=!
