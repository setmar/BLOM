! ------------------------------------------------------------------------------
! Copyright (C) 2006-2018 Mats Bentsen, Martin Lilleeng Sætra
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

      module diffus_cuda
      attributes(global) subroutine diffus_kernel(usflx_dev,utflx_dev,vsflx_dev, 
     .      vtflx_dev,usflld_dev,utflld_dev,vsflld_dev,vtflld_dev, iu_dev,
     .      saln_dev,temp_dev,sigma_dev,
     .      mm,nn,delt1,dpeps)

      real :: usflx_dev(), utflx_dev(), vsflx_dev(), vtflx_dev()
      real :: usflld_dev(), utflld_dev(), vsflld_dev(), vtflld_dev()
      integer :: iu_dev()
      real :: saln_dev(), temp_dev(), sigma_dev()
      integer :: k,j,i,kn,km

c      real, shared :: usflx_s(), utflx_s(), vsflx_s(), vtflx_s()
c      real, shared :: usflld_s(), utflld_s(), vsflld_s(), vtflld_s()
c      real, shared :: saln_s(), temp_s(), sigma_s()
      
c     Global thread ids
      idx = (blockIdx%x-1) * blockDim%x + threadIdx%x
      idy = (blockIdx%y-1) * blockDim%y + threadIdx%y
      idz = (blockIdx%z-1) * blockDim%z + threadIdx%z

c     Local thread ids
      lidx = (blockIdx%x-1) * blockDim%x + threadIdx%x
      lidy = (blockIdx%y-1) * blockDim%y + threadIdx%y
      lidz = (blockIdx%z-1) * blockDim%z + threadIdx%z

      do k=1,kk
        kn=k+nn
        km=k+mm

        do j=0,jj+1
          do i=0,ii+2
            q=delt1*.5*(difiso(i-1,j,k)+difiso(i,j,k))
     .        *scuy(i,j)*scuxi(i,j)
     .        *max(min(dp(i-1,j,kn),dp(i,j,kn)),dpeps)

            usflld(i,j,km)=q*(saln(i-1,j,kn)-saln(i,j,kn))
            utflld(i,j,km)=q*(temp(i-1,j,kn)-temp(i,j,kn))

            usflx(i,j,km)=usflx(i,j,km)+usflld(i,j,km)
            utflx(i,j,km)=utflx(i,j,km)+utflld(i,j,km)
          enddo
        enddo
 
        do j=0,jj+2
          do i=0,ii+1
            q=delt1*.5*(difiso(i,j-1,k)+difiso(i,j,k))
     .        *scvx(i,j)*scvyi(i,j)
     .        *max(min(dp(i,j-1,kn),dp(i,j,kn)),dpeps)

            vsflld(i,j,km)=q*(saln(i,j-1,kn)-saln(i,j,kn))
            vtflld(i,j,km)=q*(temp(i,j-1,kn)-temp(i,j,kn))

            vsflx(i,j,km)=vsflx(i,j,km)+vsflld(i,j,km)
            vtflx(i,j,km)=vtflx(i,j,km)+vtflld(i,j,km)
          enddo
        enddo
        
        do j=0,jj+1
          do i=0,ii+1
            q=1./(scp2(i,j)*max(dp(i,j,kn),dpeps))

            saln(i,j,kn)=(saln(i,j,kn)
     .                  -q*(usflld(i+1,j,km)-usflld(i,j,km)
     .                     +vsflld(i,j+1,km)-vsflld(i,j,km)))
     .                   *iu(i,j)

            temp(i,j,kn)=(temp(i,j,kn)
     .                  -q*(utflld(i+1,j,km)-utflld(i,j,km)
     .                     +vtflld(i,j+1,km)-vtflld(i,j,km)))
     .                   *iu(i,j)

            sigma(i,j,kn)=sig(temp(i,j,kn),saln(i,j,kn)) *iu(i,j)
          enddo
        enddo
      enddo

      end module
