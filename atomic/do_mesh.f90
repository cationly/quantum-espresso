!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
      subroutine do_mesh(rmax,zmesh,xmin,dx,ibound,grid)
!---------------------------------------------------------------
!
      use kinds, only : DP
      use radial_grids, only : radial_grid_type, ndmx
      implicit none
      type(radial_grid_type),intent(out) :: grid

      integer, intent(in)::  ibound
      real(DP),intent(in) :: rmax, zmesh, dx
      real(DP),intent(inout):: xmin

      real(DP) :: r(ndmx), r2(ndmx), rab(ndmx), sqr(ndmx), xmax, x
      integer :: mesh, i
!      
      xmax=log(rmax*zmesh)
      mesh=(xmax-xmin)/dx+1
!
!      mesh must be odd for simpson integration. 
!
      mesh=2*(mesh/2)+1
      if(mesh+1 > ndmx) call errore('do_mesh','ndmx is too small',1)
      if(ibound == 1) xmin=xmax-dx*(mesh-1)
!
      do i=1,mesh
         x=xmin+DBLE(i-1)*dx
         r(i)=exp(x)/zmesh
         rab(i)=(r(i)+exp(xmin)/zmesh)*dx
         !!! r(i)=exp(xmin)*(exp((i-1)*dx)-1.0_dp)/zmesh
         !!! rab(i)=r(i)*dx
         r2(i)=r(i)*r(i)
         rab(i)=r(i)*dx
         sqr(i)=sqrt(r(i))
      end do
!
      grid%mesh = mesh
      grid%dx   = dx
      grid%xmin = xmin
      grid%rmax = rmax
      grid%zmesh = zmesh
      grid%r(1:mesh)   = r(1:mesh)
      grid%r2(1:mesh)  = r2(1:mesh)
      grid%rab(1:mesh) = rab(1:mesh)
      grid%sqr(1:mesh) = sqr(1:mesh)

      return
      end subroutine do_mesh
