 program cube_taucalc 

  use initcalc  
  use arrays
  use io_cubes
  use rotate 

  implicit none
!--------------------------------
! nc routines can be potentially used in parallel:
  integer sizee, myrank
!--- cube dimensions 
  integer  Nzz, Ni, Nt, Nx, Ny, Nz, Nzcut 

! --- for writting nc files
  integer ier, ncid, fnx, fny, fnz 
! --- loop integers
  integer i,j,k,m

 
!--- reading cube:
  character(len=80)  filename1, filename2, filename3, filename 
  character(len=80)  folder
  character(len=6)   filenumber
  character(len=50)  numberx

! --- for zgrid 
  real(kind=8)  dx, dz, dy

! --- for rotation 
  integer indum
  real(kind=8) meanzt, onepoint 
  real(kind=8) mu, theta, pivot, pivotdx,  newdx
  real(kind=8) summean 

! --- for setting up tau grid on which to  map :
  integer Ngrid 
  integer Ngridmax 
  parameter( Ngridmax = 82)
   
  real(kind=8) maxz, tau1lg, step, tau2lg
 
  
!--- functions -----------------------------------------------!  
  real(kind=8) introssk
  integer      map1



!       filename='result_0.120000.fits'
!       folder='/scratch/yeo/SATIRE3D/D000/snapshots/'
    read(*,*) folder
    read(*,*) filenumber
    sizee = 1
    myrank = 0  

!----------------------------------------------------------------
!      initialize
! ---- read control file and set all logicals ------------------
 
    call init_calc(mu, tau1lg, step, tau2lg)
    if (gettaug) then 
     Ngrid  = int((tau2lg - tau1lg)/step) +1 
     if (Ngrid .gt. Ngridmax) then 
        print*,' Tau grid configuration results in too many points; Ngrid =', Ngrid 
        print*, ' Application will be aborted ' 
        stop
     endif
    endif 

! get cube dimensions :

  if (mpiread .or. binread) then 

    filename=trim(folder)//'Header.'//trim(filenumber)
    open (unit = 1, file= filename, form='formatted',status ='old')
    read(1,*) Nz, Nx, Ny, dz, dx, dy
    print*, ' Nx, Ny, Nz', Nx, Ny, Nz
    close(unit=1)

  else

    print*, ' file name for dimensions nx, ny, nz and dz, dx, dy'
    read(*,*) filename
    open(unit = 2, file=filename, form='formatted', status='old')
    read(2,*) Nz, Nx, Ny, dz, dx, dy
    close(unit=2) 

  endif 



!---- for tau - integration we do not need the whole depth of the cube
!--- set depth of temporary arrays
   if (ifmu ) then  
     Nzcut = Nz * int((1.0d0/mu))
   else
     mu = 1.0
     Nzcut = int(Nz/2)
   endif 

!----------------------------------------------------------------------------
! ---           read kappa table -------------------------------------------! 

! kappa-table: includes kappa
  if(ifmu .or. gettaug) then 
    open(unit = 2, file='ross_table.dat', form='formatted', status='old')  
  else 
    open(unit =2, file='kappa_table.dat', form='formatted', status='old')
  endif 

    read(2,*) numt, numpres

!----- allocate  arrays

    call set_arrays(Nx, Ny, Nz, Nzcut, numt, numpres)
    if (gettaug) call set_tarrays(Nx, Ny, Ngrid)

!---- read in the kappa table 1 :

    read(2,*) (tabt(i), i = 1 , numt)
    read(2,*) (tabp(i), i = 1, numpres)
    do i = 1, numt
       read(2,*) (kappatab(i,j), j = 1, numpres )
    end do

    close(unit=2)

!--------------------------------------------------------------------------
!------     READ in the CUBE ----------------------------------------------
!--------------------------------------------------------------------------

   filename1=trim(folder)//'eosT.'//trim(filenumber)
   filename2=trim(folder)//'eosP.'//trim(filenumber)
   filename3=trim(folder)//'result_prim_0.'//trim(filenumber)

   if (mpiread) then 
    !    Temperature  
    Nt = nx*ny*nz

    call init_comm

    call read_cube_mpi(filename1, Nt, nx, ny, nz, T)
    call read_cube_mpi(filename2, Nt, nx, ny, nz, P)
    call read_cube_mpi(filename3, Nt, nx, ny, nz, rho)

    call fin_comm
!   --- if binary read call binary_read
   else if (binread) then 

    call read_cube_bin(filename1, nx, ny, nz, T)
    call read_cube_bin(filename2,  nx, ny, nz, P)
    call read_cube_bin(filename3,  nx, ny, nz, rho)

   else if (fitsread) then 
! ---- if fits call fits read 
     call read_cube_fits(filename3, nx, ny, nz, rho)
     call read_cube_fits(filename1, nx, ny, nz, T)
     call read_cube_fits(filename2, nx, ny, nz, P)
   else 

    print*,' ERROR, it is not specified in which format to read the cubes' 
    stop

   endif 

! z-grid needs to be set for the cube as it is read in:

   zgrid(1)=1.0d0

   do k=2,Nzcut
     zgrid(k)=zgrid(k-1)+dz
   end do


!----- if the cube needs to be rotated, we need the pivot point
    
   if (ifmu) then 
! this is because we need to calculate the tau rosseland on the actual cube, and Nz < Nzcut 
     Nzcut = Nz 
   endif 

!----  ------------------------------------------------------------! 
!---  GET rays from top inwards and calculate TAU -----------------! 

! ----- NOTE:       in this place the tau is calculated, now if we rotate, then it is the tau Rosseland
!-------            if there is no rotation, then it is the higher tau opacity in the UV.
!-------            if the first is true, and the second needs to be done, the another kappa-table need to be read later and applied. 

    summean = 0.0d0 
    onepoint = 1.0d0 

!--- for debugging: 
    Nx = 2

      do i = 1, Nx
        do k = 1, Ny
            do j = 1, Nzcut

             tempt(j) = T(i,k, Nz-j+1)
             tempp(j) = P(i,k, Nz-j+1)
             tempr(j) = rho(i,k, Nz-j+1)
!    get kappa* rho
             kappa(j) = introssk(tempt(j), tempp(j))
             kappa(j) = kappa(j)* tempr(j)

            end do
! inegrate to get tau
         call integ(zgrid,kappa,taut,Nzcut,(kappa(1)*zgrid(1)))
!--- note that the cube has its top at Nz, and we integrated from top inwards, 
! --- but for convenience we store the tau array having its top at 1
         tau(i, k, 1:Nzcut) = taut(1:Nzcut)
         indum = map1(taut, zgrid, Nzcut, onepoint, meanzt, 1 ) 
         summean = summean + meanzt
        end do 
      end do 

     summean = summean/(Nx*Ny)
     print*, ' Summean = ', summean 
  
!-------------------------------------------------------------------------!
!    --------   DO the ROTATION ------------------------------------------!

   if (ifmu) then 

     Nzcut = Nz * int((1.0d0/mu)) 
     pivot = summean 
     call rotate_cube(mu, pivot, nx, ny, nz, dx, dy, dz)  
!---- after rotation was performed, the arrays are stored in newT, newP, newrho!
! okey so now we actually have Nzcut points in the vertical direction
       Nzcut = Nz * int((1.0d0/mu))
!
     print*,  ' Finished rotation' 

     if(tau200) then 
       print*, ' Read kappa table for tau - 200 calculation after rotation '
!---- get the other kappa table: 
       open(unit =2, file='kappa_table.dat', form='formatted', status='old')
       read(2,*) numt, numpres
! ---- make sure that the numt, numpress in this table are <= than in the tau table because 
!      of allocation, IN CASE THIS IS NOT TURE: 
!      - deallocate kappatab, tabt, tabp, 
!      - allocate with new dimensions, and read in !
 
       read(2,*) (tabt(i), i = 1 , numt)
       read(2,*) (tabp(i), i = 1, numpres)
       do i = 1, numt
          read(2,*) (kappatab(i,j), j = 1, numpres )
       end do

       close(unit=2)


!------------ again calc tau, this time on the rotated cube:

      do i = 1, Nx
        do k = 1, Ny
            do j = 1, Nzcut

             tempt(j) = newT(i,k, Nzcut-j+1)
             tempp(j) = newP(i,k, Nzcut-j+1)
             tempr(j) = newrho(i,k, Nzcut-j+1)
!    get kappa* rho
             kappa(j) = introssk(tempt(j), tempp(j))
             kappa(j) = kappa(j)* tempr(j)

            end do
! inegrate to get tau
         call integ(zgrid,kappa,taut,Nzcut,(kappa(1)*zgrid(1)))
!--- note that the cube has its top at Nz, and we integrated from top inwards, 
! --- but for convenience we store the tau array having its top at 1
         tau(i, k, 1:Nzcut) = taut(1:Nzcut)
        end do
      end do

         print*, ' Finished calculating tau-200 after rotation' 
 
     endif 


   endif 
!-------------------------------------------------------------------------!
!   ----- Do we need to map onto a different tau grid?  ------------------!

!--- Note that gettau and tau200 should not be used together!!!! 

   if (gettaug)  then 

   print*, ' Start to set up tau-grid onto which to interpolate' 
     
!   --- get up taugrid 
    do i = 1, ngrid
     taugrid(i) = tau1lg+(i-1)*step
     taugrid(i) = 10**(taugrid(i)) 
    end do 

!--- since after rotation our arrays are called differently, and I did not come up 
!  --- with an easier solution there are two possibilities: 

!   1) ----:

     if (ifmu) then 
      print*, ' get temporary arrays, get tau Rosseland and interpolate'  
      do i = 1, Nx
        do k = 1, Ny
            do j = 1, Nzcut

             tempt(j) = newT(i,k, Nzcut-j+1)
             tempp(j) = newP(i,k, Nzcut-j+1)
             tempr(j) = newrho(i,k, Nzcut-j+1)
!    get kappa* rho
             kappa(j) = introssk(tempt(j), tempp(j))
             kappa(j) = kappa(j)* tempr(j)

            end do
! debug: 

! inegrate to get tau
         call integ(zgrid,kappa,taut,Nzcut,(kappa(1)*zgrid(1)))

         indum = map1(taut, tempt, Nzcut, taugrid, tempa, Ngrid)
         outT(i,k,1:Ngrid) = tempa(1:Ngrid)

         indum = map1(taut, tempp, Nzcut, taugrid, tempa, Ngrid)
         outP(i,k,1:Ngrid) = tempa(1:Ngrid)

         indum = map1(taut, tempr, Nzcut, taugrid, tempa, Ngrid)
         outrho(i,k,1:Ngrid) = tempa(1:Ngrid)

         indum = map1(taut, zgrid, Nzcut, taugrid, tempa, Ngrid)
! note z starts from top pointing inwards, Rinat needs z pointinh outwards! 
         maxz = maxval(tempa)
         do  j = 1, Ngrid  
           outz(i,k,j) = abs(tempa(j) - maxz)
         end do  
          
         

        end do
      end do

!   2) -----: note if there was no rotation then, tau holds already the rosseland tau, so it does not need to be re-calculated. 
     else
       do i = 1, Nx
          do k = 1, Ny
            do j = 1, Nzcut

             tempt(j) = T(i,k, Nzcut-j+1)
             tempp(j) = P(i,k, Nzcut-j+1)
             tempr(j) = rho(i,k, Nzcut-j+1)
             taut(j) =  tau(i,k, j)
            end do 
            
            indum = map1(taut, tempt, Nzcut, taugrid, tempa, Ngrid)
            outT(i,k,1:Ngrid) = tempa(1:Ngrid)

            indum = map1(taut, tempp, Nzcut, taugrid, tempa, Ngrid)
            outP(i,k,1:Ngrid) = tempa(1:Ngrid)

            indum = map1(taut, tempr, Nzcut, taugrid, tempa, Ngrid)
            outrho(i,k,1:Ngrid) = tempa(1:Ngrid)

            indum = map1(taut, zgrid, Nzcut, taugrid, tempa, Ngrid)
! note z starts from top pointing inwards, Rinat needs z pointinh outwards! 
            maxz = maxval(tempa)

            do  j = 1, Nzcut
              outz(i,k,j) = abs(tempa(j) - maxz)
            end do


          end do  
       end do 
     endif 




   endif 

!-------------------------------------------------------------------------!
! ---             Finally , PRINT OUT the results ------------------------!
!---- the if structure is as follows, since if tau200 is set, then gettaug can not be set. But rotation and tau200 can be set simulaneously. For gettau the end array will always be in the nx, ny, ngrid with and without rotation


    print*, ' Finished all claculations let us write stuff out!' 
    
    if (gettaug) then 

      Nzz = Ngrid
      !---- need to get mu as a number for the file name:
       call str(int(mu*10), numberx)
! since the cube is rotated, the dz is changed, write out a new Headerfile:
       filename='Header_mu_'//trim(numberx)//'.'//trim(filenumber)
       open (unit = 1, file= filename, form='formatted',status ='new')
       write (1,*) ' tau-grid points, start lgtau, step, finish lgtau,  Nx, Ny, dx,  dy' 
       write(1,*) Ngrid, tau1lg, step , tau2lg, Nx, Ny, dx, dy 
       close(unit=1)


!      Temperature 
       filename='T_onTau.'//trim(filenumber)//'.nc'
        call  create_netcdf(ncid, filename, 'T',  nx, ny, nzz, ier)
        call write_netcdf(ncid, myrank, sizee, 'T', outT, nx, nx, ny, nzz, comm, ier)
        call close_netcdf(ncid, ier)

!      Presssure  
       filename='P_onTau.'//trim(filenumber)//'.nc'
        call  create_netcdf(ncid, filename, 'P',  nx, ny, nzz, ier)
        call write_netcdf(ncid, myrank, sizee, 'P', outP, nx, nx, ny, nzz, comm, ier)
        call close_netcdf(ncid, ier)
!      density 
       filename='rho_onTau.'//trim(filenumber)//'.nc'
        call  create_netcdf(ncid, filename, 'R',  nx, ny, nzz, ier)
        call write_netcdf(ncid, myrank, sizee, 'R', outrho, nx, nx, ny, nzz, comm, ier)
        call close_netcdf(ncid, ier)
!---- zgrid 

       filename='Z_onTau.'//trim(filenumber)//'.nc'
        call  create_netcdf(ncid, filename, 'Z',  nx, ny, nzz, ier)
        call write_netcdf(ncid, myrank, sizee, 'Z', outz, nx, nx, ny, nzz, comm, ier)
        call close_netcdf(ncid, ier)




    else if (tau200 .or. ifmu) then 
! write new NetCDF file for bv and tau_ross
        Nzz = Nz
        if(ifmu) Nzz = Nzcut

        if (tau200) then 
          filename='result_ttau200.'//trim(filenumber)//'.nc'
          call  create_netcdf(ncid, filename, 'tau',  nx, ny, nzz, ier)
          call write_netcdf(ncid, myrank, sizee, 'tau', tau, nx, nx, ny, nzz, comm, ier)
          call close_netcdf(ncid, ier)
        endif 



!---- need to get mu as a number for the file name:
       call str(int(mu*10), numberx)
! since the cube is rotated, the dz is changed, write out a new Headerfile:
       filename='Header_mu_'//trim(numberx)//'.'//trim(filenumber)
       open (unit = 1, file= filename, form='formatted',status ='new')
       dz = dz*mu
       write (1,*) Nzz, Nx, Ny, dz, dx, dy
       close(unit=1)


!      Temperature 
       filename='T_mu_'//trim(numberx)//'.'//trim(filenumber)//'.nc'
        call  create_netcdf(ncid, filename, 'T',  nx, ny, nzz, ier)
        call write_netcdf(ncid, myrank, sizee, 'T', newT, nx, nx, ny, nzz, comm, ier)
        call close_netcdf(ncid, ier)

!      Presssure  
       filename='P_mu_'//trim(numberx)//'.'//trim(filenumber)//'.nc'
        call  create_netcdf(ncid, filename, 'P',  nx, ny, nzz, ier)
        call write_netcdf(ncid, myrank, sizee, 'P', newP, nx, nx, ny, nzz, comm, ier)
        call close_netcdf(ncid, ier)
!      density 
       filename='rho_mu_'//trim(numberx)//'.'//trim(filenumber)//'.nc'
        call  create_netcdf(ncid, filename, 'R',  nx, ny, nzz, ier)
        call write_netcdf(ncid, myrank, sizee, 'R', newrho, nx, nx, ny, nzz, comm, ier)
        call close_netcdf(ncid, ier)


     end if 



!------ tidy up ! ---------- never leave a mess ;) ---------------------!
 
       call close_arrays 
       if (ifmu) call close_muarrays
       if (gettaug) call close_tarrays


   end
!-----        END of main routine -----------------------------------------

!--------------------------------------------------------------------------

!----------------------------------------------------------------------------
! writing numbers to strings!
   subroutine str(i, s)
      implicit none
      character*50, intent(out) :: s
      integer, intent(in) :: i

      write (s, *) i
      s = adjustl(s)
      s = trim(s)
   end subroutine




      subroutine integ (x, f, fint, n, start)
 
      use arrays
      implicit none
!
!--------------------------- DUMMY ARGUMENTS ---------------------------
!
      double precision  f(*), fint(*), start, x(*)
      integer  n
!
!--------------------------- LOCAL VARIABLES --------------------------
!
      double precision  db, df, wt
      integer           i, nn
      double precision a1(n), b1(n), c1(n)

!
!-------------------------------  EXTERNALS ----------------------------
!
      external parcoe 
!
!------------------------------- EXECUTION -----------------------------
!

      nn = n
      call parcoe (f, x, (nn))
      fint(1) = start
      df      = x(2) - x(1)
      fint(2) = fint(1) + (a(1,1) + (b(1,1) * 0.5d0 + c(1,1) / & 
     &          3.0d0 * df) * df) * df
!
      do i = 2, nn - 2
         wt        = 1.0d0
         if (c(i,1) .ne. 0.0d0 .or.  & 
     &       c(i,2) .ne. 0.0d0 ) wt = abs(c(i,1)) / (abs(c(i,1)) + & 
     &      abs(c(i,2)))
         df        = (x(i) - x(i+1))
         db        = -df
         fint(i+1) = fint(i) + wt * (a(i,2) + (b(i,2) * 0.5d0 + c(i,2)/  & 
     &      3.0d0 * db) * db) * db - (1.0d0 - wt) * (a(i,1) + (b(i,1) *  & 
     &      0.5d0 + c(i,1) / 3.0d0 * df) * df) * df
      end do
!
      df       = x(nn) - x(nn - 1)
      fint(nn) = fint(nn - 1) + (a(nn,1) + (b(nn,1) * 0.5d0 + c(nn,1) /  & 
     &         3.0d0 * df) * df) * df


      end
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   I N T E G !!!!!!!!!!!!!!
      subroutine parcoe (f, x, n)

      use arrays
      implicit none
!
!.... NOW ALWAYS IN DOUBLE PRECISION
!.... THIS HAS BEEN REWRITTEN IN THE SPIRIT OF THE U LONDON SUGGESTION
!.....  VAR(M,1 OR 2)  1=FORWARD INTERPOLATION, 2=BACKWARD INTERPOLATION
!
!--------------------------- DUMMY ARGUMENTS ---------------------------
!
      double precision  f(*), x(*)
      integer           n
!--------------------------- LOCAL VARIABLES --------------------------
!
      double precision  d, dx, wt
      integer           j, n1, jp2, jp1, jm1
!
!------------------------------- EXECUTION -----------------------------
!

      n1      = n - 1
      c(1, 1) = 0.0d0
      c(1, 2) = c(1, 1)

      b(1, 1) = (f(2) - f(1)) / (x(2) - x(1))
      b(1, 2) = b(1, 1)

      a(1, 1) = f(1)
      a(1, 2) = a(1, 1)

      c(n, 1) = 0.0d0
      c(n, 2) = c(n, 1)

      b(n, 1) = (f(n) - f(n1)) / (x(n) - x(n1))
      b(n, 2) = b(n, 1)

      a(n, 1) = f(n1)
      a(n, 2) = a(n, 1)
!
!--------------------
      if (n .gt. 2) then
!
         do j = 2, n1 - 1
            jp2     = j + 2
            jp1     = j + 1
            jm1     = j - 1
            d       = (f(jp1) - f(j)) / (x(jp1) - x(j))
            dx      = x(jp2) - x(j)
            c(j, 1) = (f(jp2) - f(jp1)) / (x(jp2) - x(jp1)) / dx - d / dx
            b(j, 1) = d + c(j,1) * (x(jp1) - x(j))
            a(j, 1) = f(jp1)
!
            if (j .eq. 2) then
!.... MUST COMPUTE THE BACKWARD COEFFICIENT WHEN J = 2
               d       = (f(2) - f(1)) / (x(2) - x(1))
               dx      = x(3) - x(1)
               c(j, 2) = (f(3) - f(2)) / (x(3) - x(2)) / dx - d / dx
               b(j, 2) = d + c(j,2) * (x(2) - x(1))
               a(j, 2) = f(2)
!
            else
               a(j, 2) = a(jm1, 1)
               b(j, 2) = b(jm1, 1)
               c(j, 2) = c(jm1, 1)
            end if
!
         end do
!
         a(n1, 1) = f(n)
         b(n1, 1) = -d
         c(n1, 1) = 0.0d0
         a(n1, 2) = f(n1)
         b(n1, 2) = b(n1-1, 1)
         c(n1, 2) = c(n1-1, 1)
      end if
!
      end




! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 double precision function introssk( tv, pv)


!--- common arrays; tabrosskap,  tabt , tabp 
! ---- common constants tenlog 
      use arrays 
      implicit none

!.... ASSUMES THAT VTURB IS CONSTANT AND THAT THE Rosstab has been calculated
!.... ONLY FOR THAT VTURB
!
!-------------------------------- COMMONS -----------------------------
!
!
!--------------------------------- CONSTANTS ---------------------------
!
      double precision  ttenlg
      parameter ( ttenlg = 0.001d0 * tenlog )
!
!------------------------------- LOCAL VARIABLES -----------------------
!
      double precision tv, pv

      double precision  resulta, co1, co2, co3, co4
      double precision  plog, tl,  x, y

      integer    j, ip, it
!
!---------------------------------- EXECUTION --------------------------
!
     
            tl = min(max(log10(tv) , tabt(1)), tabt(numt) )
            it = 2
!
            do while (tabt(it) .le. tl .and. it .lt. numt)
               it = it + 1
            end do
!
            plog = min(max(log10(pv), tabp(1)), tabp(numpres))
            ip   = 2
!
            do while (tabp(ip) .le. plog .and. ip .lt. numpres)
               ip = ip + 1
            end do

!
            x      = (tl   - tabt(it-1)) / (tabt(it) - tabt(it-1))
            y      = (plog - tabp(ip-1)) / (tabp(ip) - tabp(ip-1))
            co1 = (1.0d0 - x) * (1.0d0 - y) * ttenlg
            co2 = (1.0d0 - x) * y * ttenlg
            co3 = x * (1.0d0 - y) * ttenlg
            co4 = x * y * ttenlg
!.... THE STEPS HAVE BEEN SCALED BY 1000
!
         resulta  = exp( co1 * dble(kappatab( it-1, ip-1)) + &
     &             co2 * dble(kappatab( it-1  , ip)) + &
     &             co3 * dble(kappatab(it, ip-1)) + &
     &             co4 * dble(kappatab( it  , ip)) )
         introssk= resulta
 end function 


!--------------------------- Map function -----------------------------------

 integer function map1 (xold, fold, nold, xnew, fnew, nnew)
      implicit none
!
!--------------------------- DUMMY ARGUMENTS ---------------------------
!
      double precision  fnew(*), fold(*), xnew(*), xold(*)
      integer  nnew, nold
!
!--------------------------- LOCAL VARIABLES --------------------------
!
      double precision  ab, af, bb, bf, cb, cf, d, db, df, wt
      integer  l, ll, k, lm1, lm2, lp1
!
!------------------------------- EXECUTION -----------------------------
!
      l = 2
      ll = 0

      do k = 1, nnew

         do while (l .le. nold .and. xnew(k) .ge. xold(l))
            l = l + 1
         end do

         if (l .gt. nold) l = nold

         if (l .gt. 2 .and. l .lt. nold) then
!
!.... PARABOLIC CASE
!
            if (l .ne. ll) then

               if (l .gt. 3 .and. l .eq. ll+1) then
                  ab = af
                  bb = bf
                  cb = cf

               else
!
!.... MUST COMPUTE THE BACKWARD COEFFICIENTS
!
                  lm1 = l - 1
                  lm2 = l - 2
                  d = (fold(lm1) - fold(lm2)) / (xold(lm1) - xold(lm2))

                  cb = ((fold(l) - fold(lm1)) / (xold(l) - xold(lm1)) - d) / (xold(l) - xold(lm2))

                  bb = d + cb * (xold(lm1) - xold(lm2))
                  ab = fold(lm1)
               end if
!
               lp1 = l + 1
               lm1 = l - 1

               d = (fold(l) - fold(lm1)) / (xold(l) - xold(lm1))
               cf = ((fold(lp1) - fold(l)) / (xold(lp1) - xold(l)) - d) / (xold(lp1) - xold(lm1))
               bf = d + cf * (xold(l) - xold(lm1))
               af = fold(l)
               wt = 0.0d0
               if (cf .ne. 0.0d0) wt = abs(cf) / (abs(cf) + abs(cb))
               ll = l
            end if
!
            df = xnew(k) - xold(l)
            db = xnew(k) - xold(lm1)
            fnew(k) = (1.0d0 - wt) * (af + (bf + cf * df) * df) + wt *  (ab + (bb + cb * db) * db)
!
         else
!
            if (l .ne. ll) then
               ll = l
               lm1 = l - 1
               af = fold(lm1)
               bf = (fold(l) - fold(lm1)) / (xold(l) - xold(lm1))
            end if
!
            fnew(k) = af + bf * (xnew(k) - xold(lm1))
         end if
!
      end do
!
      map1 = ll - 1
   end
!
!*********** E N D   O F   F U N C T I O N   M A P 1 ********************


