!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
subroutine compute_chi_tm(lam,ik,ikk_in,phi_in,chi_out,xc,e)
  !--------------------------------------------------------------------------
  !
  !     This routine computes the chi functions:
  !          |chi> = (\epsilon -T -V_{loc)) |psi>
  !      
  use kinds, only : DP
  use ld1inc

  implicit none
  integer :: &
       ik,    & ! the point corresponding to rc
       ikk_in,& ! the point after which the chi should be zero
       lam      ! the angular momentum

  real(DP) :: &
       e,     &       ! input: the energy
       xc(8),       & ! input: the parameters of the fit
       phi_in(ndm), & ! input: pseudo wavefunction
       chi_out(ndm)   ! output: the chi function
  !
  real(DP) :: &
       dpoly           

  real(DP), external :: pr, d2pr, dpr

  integer :: &
       n
  !
  !   Troullier-Martins: use the analytic formula
  !
  do n=1,ik
     dpoly = dpr(xc,xc(7),r(n))
     ! dpr =  first derivate of polynomial pr
     ! d2pr= second derivate of polynomial pr
     chi_out(n) = (e + (2*lam+2)/r(n)*dpoly + &
             d2pr(xc,xc(7),r(n)) + dpoly**2 - vpsloc(n))*phi_in(n)
  enddo
  do n = ik+1,mesh
     chi_out(n) = (vpot(n,1) - vpsloc(n))*phi_in(n)
  enddo
  return
end subroutine compute_chi_tm