  module rotate

  use arrays


  CONTAINS

  subroutine rotate_cube(mu, pivot, nx, ny, nz, dx, dy, dz)

   implicit none 
   integer, intent(in) :: nx, ny, nz

   integer  Nzcut 
   integer  i,j,jj,k, xn, xnp, x_l, z_l 
   real(kind=8), intent(in) :: mu, pivot, dx, dy, dz 
   real(kind=8)  theta, pivotdx,  newdx, newx, newz 
   real(kind=8)   z_f,  x_f 
   real(kind=8)   zold, znew, Tnewd 

   zold = 0.0d0
   znew = 0.0d0 
  
   theta = acos(mu)
   pivotdx = tan(theta) * pivot
   newdx = dz * sin(theta)
!   Nzcut = min(int(Nz * (0.9d0/mu)), int(3.5*Nz))
   Nzcut = int(Nz * (0.9d0 / mu))

!   ---  allocate rotation arrays
!   call set_muarray(Nzcut)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!
   zgrid(1)=0.0d0
   do k=2,Nzcut
      zgrid(k)=zgrid(k-1)+dz
   end do

! if mu not equal 1 we need new z indices! 

 

!---- start to do the whole cube 
!--- over the x-y plan


      do i = 1, Nx
        do k = 1, Ny
           do j = 1, Nzcut 
             newx = (k-1)*dx - pivotdx +(j-1)*newdx 
             newz = (j-1)*(dz)*mu;
             if (newx .lt. 0.0) then 
                 newx = newx+(Ny*dx)
             elseif (newx .ge. (Ny*dx)) then 
                 newx = newx - (floor(newx/(dble(Ny)*dx))*(Ny*dx))   
             endif 

             x_f = newx/dx
             x_l = int(x_f)
             x_f = x_f-dble(x_l)
             z_f = newz/dz
             z_l = int(z_f)
             z_f = z_f-dble(z_l)  
             z_l = z_l +1
             xn = x_l + 1
             xnp = xn+1 
             if (xnp .gt. Ny) xnp = mod(xnp, Ny)  

             newT(i,k,j)  = ((T(i, xn, z_l)) * (1.0d0-x_f) + x_f*(T(i, xnp,  z_l))) *(1.0d0- z_f)
             newT(i,k,j)  = newT(i,k,j) +  z_f*((T(i, xn, z_l+1))*(1.0d0-x_f)+ x_f*(T(i, xnp,z_l+1)))

             newP(i,k,j)  = ((P(i, xn, z_l)) * (1.0d0-x_f) + x_f*(P(i, xnp,  z_l))) *(1.0d0- z_f)
             newP(i,k,j)  = newP(i,k,j) +  z_f*((P(i, xn, z_l+1))*(1.0d0-x_f)+ x_f*(P(i, xnp,z_l+1)))


             newrho(i,k,j)  = ((rho(i, xn, z_l)) * (1.0d0-x_f) + x_f*(rho(i, xnp,  z_l))) *(1.0d0- z_f)
             newrho(i,k,j)  = newrho(i,k,j) +  z_f*((rho(i, xn, z_l+1))*(1.0d0-x_f)+ x_f*(rho(i, xnp,z_l+1)))

            znew = zold + dz*(j-1)*1.0d-5 

           end do
         end do 
       end do 










  end subroutine 





  end module 
