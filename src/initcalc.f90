module initcalc

 implicit none 

 logical ifmu, gettaug 
 logical  tau200, mpiread, fitsread, binread, sliceread 
 

contains

 subroutine init_calc(mu, tau1lg, step, tau2lg, pivot_in ) 
 
   implicit none 
   real(kind=8), intent(out) :: mu, step, tau1lg, tau2lg, pivot_in 

   integer iread,  rot, interpol, tau
  
   mpiread = .false.
   fitsread = .false.
   binread = .false. 
   sliceread = .false. 
   ifmu = .false. 
   tau200 = .false. 
   gettaug = .false. 

   mu = 1.0d0 

   open(unit=2, file='INPUT', form='formatted', status='old' )
    read(2,*) 
    read(2,*) iread
    read(2,*) rot
    read(2,*) mu
    read(2,*) tau  
    read(2,*) interpol
    read(2,*) pivot_in
    read(2,*) tau1lg, step, tau2lg 
   close(unit=2) 
!------------------------------
  if (iread .eq. 1) binread =.true.
  if (iread .eq. 2) mpiread = .true. 
  if (iread .eq. 3) fitsread = .true.
  if (iread .eq. 4) sliceread = .true. 
 
  if (rot .eq. 1) ifmu = .true.
  if (tau .eq. 1) tau200 = .true. 
  if (interpol .eq. 1)  gettaug = .true. 
  if (gettaug) tau200 = .true. 
! ----------------------------------------

    

!--- write out what is going to be calculated: ------!
   print*,'---------------------------------------------------'
   print*,'Routine to convert MURam cubes to Rays             '
   print*,'by V. Witkze                                       '
   print*,'                                                   '
   print*,'---------------------------------------------------'
   print*,' Summary of the calculation:                       '
   print*,'---------------------------------------------------' 
   print*,'                                                   '    
   if (binread) print*, ' - Cubes are read from binary files             ' 
   if (mpiread) print*, '  - Cubes are read from MURam mpiwrite files     ' 
   if (fitsread) print*, ' - Cubes are read form Fits files              ' 
   if (sliceread) print*, ' --- Slice of the cube is beeing read ---      '
    print*,'                                                   '
   if (ifmu) print*, ' the cubes will be rotated to the view angle = ', mu
   if (tau200) print*, ' the highest line opacity in the UV is calculated   '
   if (interpol .eq. 1) print*, ' the obtained rays are interpolated on a coarse tau grid ' 
 




 end subroutine init_calc





end module initcalc
