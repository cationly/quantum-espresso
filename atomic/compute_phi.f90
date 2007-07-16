!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
subroutine compute_phi(lam,ik,iknorm,nwf0,ns,xc,iflag,iok,occ)
  !--------------------------------------------------------------------------
  !
  !     This routine computes the phi functions by pseudizing the
  !     all_electron chi functions. In input it receives, the point
  !     ik where the cut is done, the angular momentum lam,
  !     and the correspondence with the all electron wavefunction
  !
  !
  !      
  use kinds, only : DP
  use ld1inc
  implicit none

  integer ::    &
       lam,  &   ! the angular momentum
       ik,   &   ! the point corresponding to rc
       iknorm, & ! the point at which set value=0.5
       nwf0, &   ! 
       ns,   &   ! the function to pseudize
       iflag,&   ! if 1 print
       iok       ! output: number of nodes of the pseudized wavefunction

  real(DP) :: &
       xc(8), occ
  !
  real(DP), parameter :: pi=3.14159265358979_DP, eps8=1.d-8
  real(DP) :: &
       fae,    & ! the value of the all-electron function
       f1ae,   & ! its first derivative
       f2ae,   & ! the second derivative
       wmax,   & !
       den,    & ! denominator
       faenor    ! the norm of the function

  integer ::    &
       isign,   &! sign of the max of the ae-wfc
       m, n, nst, nnode, nc, nc1, ij, imax, iq, i

  real(DP) :: &
       chir(ndmx,nwfx), chi_dir(ndmx,2), gi(ndmx), j1(ndmx,4), &
       f1aep1, f1aem1, jnor, psnor, fact(4), &
       cm(10), bm(4), ze2, cn(6), c2, &
       delta, a, b, c, deter, gamma, &
       lamda0, lamda3, lamda4, mu0, mu4, s0, s4, t0, t4

  real(DP), external :: deriv_7pts, deriv2_7pts, int_0_inf_dr, pr
  integer ::  &
       nbes   ! number of Bessel functions to be used

  !
  !   decide whether to use 4 Bessel functions or 3 (default)
  !
  if ( (rho0 /= 0.0_dp) .and. (lam == 0) ) then
     nbes = 4
  else
     nbes = 3
  end if
  !
  nst=(lam+1)*2
  !
  !   compute the reference wavefunction
  !
  if (new(ns)) then
     if (rel == 1) then
        call lschps(3,zed,exp(grid%dx),grid%dx,grid%mesh,grid%mesh,grid%mesh,  &
             1,lam,enls(ns),chir(1,ns),grid%r,vpot)
     elseif (rel == 2) then
        call dir_outward(ndmx,grid%mesh,lam,jjs(ns),enls(ns),grid%dx,chi_dir,grid%r,grid%rab,vpot)
        chir(:,ns)=chi_dir(:,1)
     else
        ze2=-zed*2.0_dp
        call intref(lam,enls(ns),grid%mesh,grid,vpot,ze2,chir(1,ns))
     endif
     !
     !    fix arbitrarily the norm at the cut-off radius equal to 0.5
     !
!!!  Previous version: jnor=chir(ik,ns) 
!!!  Now choose a possibly different point at which the value equals 0.5
     jnor=chir(iknorm,ns)
     do n=1,grid%mesh
        chir(n,ns)=chir(n,ns)*0.5_dp/jnor
    enddo
  else 
     do n=1,grid%mesh
        chir(n,ns)=psi(n,1,nwf0)
     enddo
  endif
  !
  !   save the all-electron function for the PAW setup
  !
  if (lpaw.and.(.not.lnc2paw)) psipaw(1:grid%mesh,ns) = chir(1:grid%mesh,ns)
  !
  !   compute the first and second derivative of all-electron function
  !
  fae=chir(ik,ns)
  f1ae=deriv_7pts(chir(1,ns),ik,grid%r(ik),grid%dx)
  f2ae=deriv2_7pts(chir(1,ns),ik,grid%r(ik),grid%dx)
  !
  !   compute the norm of the all-electron function
  !
  do n=1,ik+1
     gi(n)=chir(n,ns)**2  
  enddo
  faenor=int_0_inf_dr(gi,grid,ik,nst)
  !
  !
  if (tm) then
     !
     ! TM: the pseudo-wavefunction is written as polynomial times exponential
     !
     call find_coefficients &
          (ik, chir(1,ns), enls(ns), grid%r, grid%dx, faenor, vpot, cn, c2, lam, grid%mesh)
      do i=1,ik
         phis(i,ns) =sign( grid%r(i)**(lam+1)*exp(pr(cn,c2,grid%r(i))), &
                          chir(ik,ns) )
      end do
      xc(1:6) = cn(:)
      xc(7) = c2
  else
     !   RRKJ: the pseudo-wavefunction is written as an expansion into (3 or 4) 
     !         spherical bessel functions for r < r_c
     !   find q_i with the correct log derivatives
     
     call find_qi(f1ae/fae,xc(nbes+1),ik,lam,nbes,1,iok)
     if (iok.ne.0) return
     !
     !    compute the bessel functions
     !
     do nc=1,nbes
        call sph_besr(ik+5,grid%r,xc(nbes+nc),lam,j1(1,nc))
        jnor=j1(ik,nc)
        fact(nc)=chir(ik,ns)/jnor
        do n=1,ik+5
           j1(n,nc)=j1(n,nc)*chir(ik,ns)/jnor
        enddo
     enddo
  !
  !    compute the bm functions (second derivative of the j1)
  !    and the integrals for cm
  !
     ij=0
     do nc=1, nbes
        bm(nc)=deriv2_7pts(j1(1,nc),ik,grid%r(ik),grid%dx)
        do nc1=1,nc
           ij=ij+1
           do n=1,ik
              gi(n)=j1(n,nc)*j1(n,nc1)
           enddo
           cm(ij)=int_0_inf_dr(gi,grid,ik,nst)
        enddo
     enddo
     !
     !    solve the second order equation
     !
     if (nbes == 4) then
        wmax=0.0_dp
        do n=1,grid%mesh
           if (abs(chir(n,ns)) > wmax .and. grid%r(n) < 4.0_dp) then
              wmax=abs(chir(n,ns))
              if(chir(n,ns).lt.0.0_dp)then
                 isign=-1
              else
                 isign=+1
              endif
           endif
        enddo
        lamda0=(f2ae-bm(1))/(bm(2)-bm(1))
        lamda3=(bm(3)-bm(1))/(bm(2)-bm(1))
        lamda4=(bm(4)-bm(1))/(bm(2)-bm(1))
        if (new(ns)) then
           den=1.0_dp
        else
           den=abs(occ)
        endif
        mu0=(isign*sqrt(rho0*4.0_dp*pi/den) &
             -fact(1)+fact(1)*lamda0-fact(2)*lamda0)  &
             /(fact(1)*lamda3-fact(1)-fact(2)*lamda3+fact(3))
        mu4=(fact(1)*lamda4-fact(1)-fact(2)*lamda4+fact(4)) &
             /(fact(1)*lamda3-fact(1)-fact(2)*lamda3+fact(3))
        s0=lamda0-mu0*lamda3
        s4=lamda4-mu4*lamda3
        t0=1.0_dp-lamda0+mu0*lamda3-mu0
        t4=1.0_dp+mu4*lamda3-mu4-lamda4

        a=cm(1)*t4**2+cm(3)*s4**2+cm(6)*mu4**2+cm(10)  &
             +2.0_dp*cm(4)*t4*mu4+2.0_dp*cm(2)*t4*s4+2.0_dp*cm(5)*s4*mu4 &
             -2.0_dp*cm(7)*t4-2.0_dp*cm(8)*s4-2.0_dp*cm(9)*mu4
        
        b=-2.0_dp*cm(1)*t0*t4-2.0_dp*cm(3)*s0*s4-2.0_dp*cm(6)*mu0*mu4  &
             -2.0_dp*cm(4)*t0*mu4-2.0_dp*cm(4)*t4*mu0-2.0_dp*cm(2)*t0*s4  &
             -2.0_dp*cm(2)*t4*s0-2.0_dp*cm(5)*s0*mu4-2.0_dp*cm(5)*s4*mu0  &
             +2.0_dp*cm(7)*t0+2.0_dp*cm(8)*s0+2.0_dp*cm(9)*mu0
        
        c=cm(1)*t0**2+cm(3)*s0**2+cm(6)*mu0**2+2.0_dp*cm(4)*t0*mu0 &
             +2.0_dp*cm(2)*t0*s0+2.0_dp*cm(5)*s0*mu0-faenor
        
        deter=b**2-4.0_dp*a*c
        if (deter < 0.0_dp) then
           call infomsg ('compute phi','negative determinant', -1) 
           write(6,100) ns,f1ae/fae, f2ae, faenor
100        format(/5x,i5,' ld= ',f10.6,' f2ae',f10.6,' faenor',f10.6)
           iok=1
           return
        endif
        
        xc(4)=(-b+sqrt(deter))/(2.0_dp*a)
        xc(3)=mu0-mu4*xc(4)
        xc(2)=s0-s4*xc(4)
        xc(1)=t0-t4*xc(4)
        !
        
        do n=1,ik
           phis(n,ns)= xc(1)*j1(n,1) + xc(2)*j1(n,2)  &
                     + xc(3)*j1(n,3) + xc(4)*j1(n,4)
        enddo
     else
        gamma=(bm(3)-bm(1))/(bm(2)-bm(1))
        delta=(f2ae-bm(1))/(bm(2)-bm(1))
       
        a=(gamma-1.0_dp)**2*cm(1)+gamma**2*cm(3)  &
             -2.0_dp*gamma*(gamma-1.0_dp)*cm(2) &
             +2.0_dp*(gamma-1.0_dp)*cm(4)-2.0_dp*gamma*cm(5)+cm(6)
        
        b=2.0_DP*(1.0_dp-delta)*(gamma-1.0_dp)*cm(1) &
             -2.0_dp*gamma*delta*cm(3) &
             -2.0_dp*gamma*(1.0_dp-delta)*cm(2) &
             +2.0_dp*delta*(gamma-1.0_dp)*cm(2) & 
             +2.0_dp*(1.0_dp-delta)*cm(4)+2.0_dp*delta*cm(5)
        
        c=-faenor+(1.0_dp-delta)**2*cm(1)+delta**2*cm(3) &
             +2.0_dp*delta*(1.0_dp-delta)*cm(2)
        
        deter=b**2-4.0_dp*a*c
        if (deter < 0.0_dp) then
           call infomsg ('compute phi','negative determinant', -1) 
           write(6,110) ns,f1ae/fae, f2ae, faenor
110        format (/5x,i5,' ld= ',f10.6,' f2ae',f10.6, ' faenor',f10.6)
           iok=1
           return
        endif
        
        xc(3)=(-b+sqrt(deter))/(2.0_dp*a)
        xc(2)=-xc(3)*gamma+delta
        xc(1)=1.0_dp-xc(2)-xc(3)
        
        do n=1,ik
           phis(n,ns)= xc(1)*j1(n,1) + xc(2)*j1(n,2) &
                     + xc(3)*j1(n,3)
        enddo

     endif
     do nc=1,nbes
        xc(nc)=xc(nc)*fact(nc)
     enddo
  end if
  !
  !      for r > r(ik) the pseudo and all-electron psi(r) coincide
  !
  do n=ik+1,grid%mesh
     phis(n,ns)= chir(n,ns)
  enddo
  !
  !      check for the norm of the pseudo wavefunction
  !
  do n=1,ik
     gi(n)=phis(n,ns)**2
  enddo
  psnor=int_0_inf_dr(gi,grid,ik,nst)

  if (iflag == 1) then
     if (tm) then
        write(6,120) els(ns), rcut(ns)
120     format (/ /5x, ' Wfc  ',a3,'  rcut=',f6.3, &
          '  Using Troullier-Martins method ')
     else
        write(6,130) els(ns),rcut(ns),2.0_dp*xc(6)**2 
130     format (/ /5x, ' Wfc  ',a3,'  rcut=',f6.3, &
          '  Estimated cut-off energy= ', f11.2,' Ry')
        if (nbes == 4) write(6,140) rho0
140        format (5x,' Using 4 Bessel functions for this wfc,', &
                ' rho(0) =',f6.3)
     end if
     if (abs (faenor - psnor) .gt. eps8) then 
        write(6,'(5x," AE norm = ",f12.8,"  PS norm = ",f12.8)') faenor, psnor
        call errore("compute_phi"," AE and PS norm are different",1)
     end if
  end if
  !
  !      check for absence of nodes in the pseudo wavefunction
  !
  nnode=0  
  do n=1,ik+1
     if ( phis(n,ns) .ne. sign(phis(n,ns),phis(n+1,ns)) ) then
        write(6,150) lam,ns,grid%r(n)
150     format (5x,'l=',i4,' ns=',i4,' Node at ',f10.8)
        nnode=nnode+1
     endif
  enddo
  iok=nnode
  if (iflag == 1) write(6,160) nnode,grid%r(ik)
160 format (5x,' This function has ',i4,' nodes', ' for 0 < r < ',f8.3)

  return
end subroutine compute_phi
