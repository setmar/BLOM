! ------------------------------------------------------------------------------
! Copyright (C) 2006-2018 Mats Bentsen
!
! This file is part of BLOM.
!
! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

      subroutine diffus(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- diffusion of tracers
c --- ------------------------------------------------------------------
c
      use mod_xc
      use diffus_cuda
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      integer i,j,k,l,kn,km
      real q
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
      real sig
      external sig
c
      real dpeps
      parameter (dpeps=1.e-4)
c
      call xctilr(dp  (1-nbdy,1-nbdy,k1n), 1,kk, 2,2, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1,kk, 2,2, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1,kk, 2,2, halo_ps)
#ifdef TRC
      do nt=1,ntr
#  if defined(TKE) && !defined(TKEIDF)
        if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
        call xctilr(trc(1-nbdy,1-nbdy,k1n,nt), 1,kk, 2,2, halo_ps)
      enddo
#endif
      call xctilr(difiso, 1,kk, 2,2, halo_ps)

c      real, device, allocatable, dimension(:,:) :: usflx_dev, utflx_dev, vsflx_dev, vtflx_dev
c      real, device, allocatable, dimension(:,:) :: usflld_dev, utflld_dev, vsflld_dev, vtflld_dev
      real, device, allocatable, dimension(:,:) :: difiso_dev, scuy_dev, scuxi_dev, scvx_dev, scvyi_dev, dp_dev, scp2_dev
      integer, device, allocatable, dimension(:,:) :: iu_dev
      real, device, allocatable, dimension(:,:) :: saln_dev, temp_dev, sigma_dev
      
c
c     allocate and upload to device
c
c      allocate(usflx_dev(x), utflx_dev(x), vsflx_dev(x), vtflx_dev(x))
c      allocate(usflld_dev(x), utflld_dev(x), vsflld_dev(x), vtflld_dev(x))
      allocate(difiso_dev(x), scuy_dev(x), scuxi_dev(x), scvx_dev(x), scvyi_dev(x), dp_dev(x), scp2_dev(x))
      allocate(iu_dev(x))
      allocate(saln_dev(x), temp_dev(x), sigma_dev(x))
      
c      usflx_dev = usflx
c      utflx_dev = utflx
c      vsflx_dev = vsflx
c      vtflx_dev = vtflx
c      usflld_dev = usflld
c      utflld_dev = utflld
c      vsflld_dev = vsflld
c      vtflld_dev = vtflld
      difiso_dev = difiso
      scuy_dev = scuy
      scuxi_dev = scuxi
      scvx_dev = scvx
      scvyi_dev = scvyi
      dp_dev = dp
      scp2_dev = scp2
      iu_dev = iu
      saln_dev = saln
      temp_dev = temp
      sigma_dev = sigma

c
c     call CUDA kernel    
c
      type(dim3) :: gridsize, blocksize

      dimBlock = dim3(16,16,1)
      dimGrid = dim3(ceiling(real(ii+1)/blocksize%x),ceiling(real(jj+1)/blocksize%y),kk)

c      call diffus_kernel<<dimGrid,dimBlock>>(usflx_dev,utflx_dev,vsflx_dev, 
c     .      vtflx_dev,usflld_dev,utflld_dev,vsflld_dev,vtflld_dev, iu_dev,
c     .      saln_dev,temp_dev,sigma_dev,
c     .      mm,nn,delt1,dpeps)

      call diffus_kernel<<dimGrid,dimBlock>>(iu_dev,
     .      saln_dev,temp_dev,sigma_dev,
     .      mm,nn,delt1,dpeps)     

c
c     download from device and deallocate
c
      saln = saln_dev
      temp = temp_dev
      sigma = sigma_dev

      deallocate(usflx_dev, utflx_dev, vsflx_dev, vtflx_dev)
      deallocate(usflld_dev, utflld_dev, vsflld_dev, vtflld_dev)
      deallocate(difiso_dev, scuy_dev, scuxi_dev, scvx_dev, scvyi_dev, dp_dev, scp2_dev)
      deallocate(iu_dev)
      deallocate(saln_dev, temp_dev, sigma_dev)

      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'diffus:'
        endif
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
#ifdef TRC
        do nt=1,ntr
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
        enddo
#endif
      endif
c
      return
      end
