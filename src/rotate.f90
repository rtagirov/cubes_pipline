  module rotate

  use arrays


  CONTAINS

  subroutine rotate_cube(mu, pivot, nx, ny, nz, dx, dy, dz)

   implicit none 
   integer, intent(in) :: nx, ny, nz

   integer  Nzcut 
   integer  i,j,jj,k, xn, xnp 
   real(kind=8), intent(in) :: mu, pivot, dx, dy, dz 
   real(kind=8)  theta, pivotdx,  newdx
    



   theta = acos(mu)
   pivotdx = tan(theta) * pivot
   newdx = dz * sin(theta)

   Nzcut = Nz * int((1.0d0/mu))

!   ---  allocate rotation arrays

   call set_muarray(Nzcut)
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!
       zgrid(1)=1.0d0

        do k=2,Nzcut
         zgrid(k)=zgrid(k-1)+DZ*mu

        end do

! if mu not equal 1 we need new z indices! 

        do k =1,Nzcut
          zf(k) = zgrid(k)/dz
          zl(k) = int(zf(k))
          zf(k) = zf(k) - dble(zl(k))
!  this is for new x gird for Nx = 1, for any other Nx one has to add nx to x_low 

          xf(k) = ((k-1)*newdx - pivotdx)/dx
          xl(k) = int(xf(k))
          xf(k) = abs( xf(k) - dble(xl(k)) )
        end do


!---- start to do the whole cube 
!--- over the x-y plane
      do i = 1, Nx
        do k = 1, Ny
           do j = 1, Nzcut
             jj = (Nzcut -j +1)

             xn = xl(j)+(k-1)

             if (xn .lt. 1) then
              xn = xn + Ny
             else if (xn .gt. Ny) then
              xn = mod(xn, Ny) 
             end if

             xnp = xn+1
             if(xnp .gt. Ny) xnp = mod(xnp, Ny)
             
          newT(i,k,jj) = T(i, xn, (Nz-zl(j))) * (1.0d0-xf(j)) + xf(j)*T(i, xnp,  (Nz-zl(j)))
          newT(i,k,jj) = newT(i,k,jj)*(1.0d0 -xf(j)) +   xf(j)*(T(i, xn, (Nz-zl(j)-1))*(1.0d0-xf(j))+ xf(j)*T(i, xnp,(Nz-zl(j)-1)))
!

          newP(i,k,jj) = P(i,xn,(Nz-zl(j))) * (1.0d0-xf(j)) +   xf(j)*P(i, xnp, (Nz-zl(j)))
          newP(i,k,jj) = newP(i,k,jj)*(1.0d0 -xf(j)) + xf(j)*(P(i, xn, (Nz-zl(j)-1))*(1.0d0-xf(j))+   xf(j)*P(i, xnp, (Nz-zl(j)-1)))
!

        newrho(i,k,jj) =rho(i, xn, (Nz-zl(j)))*(1.0d0-xf(j)) +   xf(j)*rho(i, xnp,  (Nz-zl(j)))
        newrho(i,k,jj) = newrho(i,k,jj)*(1.0d0 -xf(j)) +  xf(j)*(rho(i, xn, (Nz-zl(j)-1))*(1.0d0-xf(j))+  xf(j)*rho(i, xnp, (Nz-zl(j)-1))  )
!

           end do
         end do 
       end do 










  end subroutine 





  end module 
