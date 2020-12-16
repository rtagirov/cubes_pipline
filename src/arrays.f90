Module arrays

 implicit none

 integer numt, numpres

 real(kind=8) tenlog
 parameter ( tenlog = 2.30258509299405d0)


 real(kind=4), allocatable :: buffarr(:), bufffits(:,:,:)
 real(kind=4), allocatable::  temparr(:,:,:), taur(:,:,:),  tau(:, :,:), T(:, :,: ), P(:,:,:), rho(:,:,:)

!--- 
 real(kind=4), allocatable:: newT(:,:,:), newP(:,:,:),newrho(:,:,:)
 

 real(kind=8), allocatable :: tempr(:), tempt(:), tempp(:), kappa(:), taut(:)
 real(kind=8), allocatable:: tabp(:), tabt(:)
 integer(kind=4) ,allocatable :: kappatab(:,:)

! second set to allocate
 real(kind = 8), allocatable:: zgrid(:)
 
! for mu = ... rotated calculations
 real(kind=8), allocatable:: zf(:), xf(:)
 integer, allocatable:: zl(:), xl(:)
 real(kind=8), allocatable:: btemp(:)

! for taugrid mapping 
 real(kind=4), allocatable :: outT(:,:,:), outP(:,:,:), outrho(:,:,:), outz(:,:,:)
! real(kind=8), allocatable ::  ttaugrid(:),taugrid(:), tempa(:) 
 real(kind=8), allocatable ::  ttaugrid(:, :, :),taugrid(:), tempa(:) 
 real(kind=4), allocatable ::  ttaugrid_sp(:, :, :)

! --- these are coefficients that are calculated in paracoe, but used in integ
 real(kind=8), allocatable :: a(:,:), b(:,:), c(:,:)
 
 CONTAINS

 subroutine set_arrays(nx, ny, nz, nycut, nt, np) 

 implicit none
 integer nx, ny, nz
 integer nycut
 integer nt, np
 integer ntot
 ntot = nx*ny*nz 
 allocate(zgrid(nycut))

 allocate(buffarr(ntot))
 allocate(bufffits(nx, nz, ny))
 allocate(T(nx, ny, nz))
 allocate(tau(nx, ny, nycut))
 allocate(taur(nx, ny, nycut))


 allocate(temparr(nx, ny, nz))
 
 allocate(P(nx, ny, nz))
 allocate(rho(nx, ny, nz))

 allocate(tabt(nt))
 allocate(tabp(np))
 allocate(kappatab(nt, np)) 
 T = 0.0d0
 tau = 0.0d0
 taur = 0.0d0 

 buffarr = 0.0
 bufffits = 0.0

 tabt = 0.0d0
 tabp = 0.0d0
 rho = 0.0
 kappatab = 0.0d0
 

 allocate(tempt(nycut))
 allocate(tempp(nycut))
 allocate(tempr(nycut))
 allocate(kappa(nycut))
 allocate(taut(nycut))
 


 tempt = 0.0
 tempp = 0.0
 tempr = 0.0
 taut = 0.0
 kappa = 0.0

 allocate(newT(nx,ny,nycut))
 allocate(newP(nx,ny,nycut))
 allocate(newrho(nx, ny, nycut))

 allocate(a(nycut, 2))
 allocate(b(nycut, 2))
 allocate(c(nycut, 2))


 end subroutine set_arrays 

 subroutine set_muarray(nycut)
  implicit none

  integer nycut

  allocate(zf(nycut))
  allocate(xf(nycut))
  allocate(zl(nycut))
  allocate(xl(nycut))
 
 end subroutine set_muarray

 subroutine set_tarrays(nx, ny, nz)
   implicit none 
   integer, intent(in):: nx, ny, nz

   allocate(taugrid(nz))
!   allocate(ttaugrid(nz))
   allocate(ttaugrid(nx, ny, nz))
   allocate(ttaugrid_sp(nx, ny, nz))

   allocate(tempa(nz))
   allocate(outz(nx, ny, nz))
   allocate(outT(nx, ny, nz))
   allocate(outP(nx, ny, nz))
   allocate(outrho(nx, ny, nz))
  

 end subroutine set_tarrays


 subroutine close_arrays 
  implicit none

  deallocate(zgrid)
  deallocate(buffarr)
  deallocate(bufffits)
  deallocate(T)
  deallocate(tau)
  deallocate(taur)


  deallocate(temparr)
 
  deallocate(rho)
  deallocate(P)
  deallocate(tabp)
  deallocate(tabt)
  deallocate(kappatab)

  deallocate(tempp)
  deallocate(tempr)
  deallocate(tempt)
  deallocate(taut)
  deallocate(kappa)


  deallocate(newT)
  deallocate(newP)
  deallocate(newrho) 

  deallocate(a)
  deallocate(b)
  deallocate(c) 


 end subroutine close_arrays 

 subroutine close_muarrays
   implicit none

  deallocate(xf)
  deallocate(zf)
  deallocate(xl)
  deallocate(zl)


 end subroutine close_muarrays

 subroutine close_tarrays
  implicit none

  deallocate(outT)
  deallocate(outP) 
  deallocate(outrho)
  deallocate(outz)
  deallocate(taugrid)
  deallocate(ttaugrid)
  deallocate(ttaugrid_sp)

  deallocate(tempa)
  


 end subroutine close_tarrays

 end module arrays

