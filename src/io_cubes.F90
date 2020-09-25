  module io_cubes 

   use arrays

   implicit none
#ifdef MPIF
   include 'mpif.h'
#endif
   integer comm, ierror

  CONTAINS 
  
  subroutine read_cube_mpi(filename, Nt, nx, ny, nz, outarr)

   implicit none
   integer, intent(in):: Nt, nx, ny, nz
   character(*), intent(in):: filename

   real(kind=4), intent(out) :: outarr(nx, ny, nz)

   integer mfh, count, datatype, status, amode, mpi_info
   integer i,j,k, Ni

!---------- New CUbes are written with MPI_write

#ifdef MPIF
      call  MPI_File_open(comm,trim(filename),MPI_MODE_RDONLY,MPI_INFO_NULL,mfh, ierror)

      call MPI_File_read_all(mfh, buffarr,  Nt,  MPI_FLOAT, status, ierror)

      call  MPI_File_close(mfh, ierror)

#else
     print*, ' your are trying to read mpi-read without mpi, we will abort'
     stop

#endif
       Ni = 0
       do i = 1, Nx
         do j = 1, Ny
           do k = 1, Nz
             Ni = Ni +1
             outarr(i,j,k) = (buffarr(Ni))
           end do
         end do
       end do

      return 


  end subroutine read_cube_mpi

  subroutine fin_comm
   implicit none

#ifdef MPIF  
   call MPI_Finalize(ierror) 
#endif 

  end subroutine fin_comm



  subroutine init_comm
   implicit none

#ifdef MPIF
    call MPI_Init ( ierror )
    comm = MPI_COMM_WORLD   
#else
    print*, ' You are trying to initiat MPI, without MPIF flag, abort!' 
    stop
#endif 

  end subroutine init_comm
!---------------------------------------------------------------------------------
!   Read binary 



  subroutine read_cube_bin(filename, nx, ny, nz, outarr)
 
   implicit none 
   
   integer, intent(in) :: nx, ny, nz
   real(kind=4), intent(out):: outarr(nx, ny, nz)
   character(*) :: filename
   integer i, j, k

 
   open (unit = 1, file= filename, form='unformatted',status ='old')
       do i = 1, Nx
         do j = 1, Ny
           do k = 1, Nz
             read(1) outarr(i,j,k)

           end do
         end do
       end do

   close(unit=1)


   return 


 end subroutine read_cube_bin





! --------------- READ FITS files + fits routines ---------------------------------
#ifdef FITS

  subroutine read_cube_fits(filename, nx, ny, nz, outarr)
   implicit none 
   
   integer, intent(in) :: nx, ny, nz
   real(kind=4), intent(out):: outarr(nx, ny, nz)
   character(*) :: filename
   integer i, j, k
   


   call readheader(trim(filename))
   call readimage(trim(filename), nx, ny, nz)
       do i = 1, Nx
         do j = 1, Ny
           do k = 1, Nz
           outarr(i,j,k) = bufffits(i,k,j)
           end do
         end do
       end do
   return 

  end subroutine read_cube_fits
 

  subroutine readheader(filename)
!  Print out all the header keywords in all extensions of a FITS file
    integer status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i,j
    character filename*80,record*80
!  The STATUS parameter must always be initialized.
    status=0
!  Get an unused Logical Unit Number to use to open the FITS file.
    call ftgiou(unit,status)
!     name of FITS file 
!      filename='ATESTFILEZ.FITS'
!     open the FITS file, with read-only access.  The returned BLOCKSIZE
!     parameter is obsolete and should be ignored. 
    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
      j = 0
100   continue
      j = j + 1
      print *,'Header listing for HDU', j
!  The FTGHSP subroutine returns the number of existing keywords in the
!  current header data unit (CHDU), not counting the required END keyword,
      call ftghsp(unit,nkeys,nspace,status)
!  Read each 80-character keyword record, and print it out.
      do i = 1, nkeys
          call ftgrec(unit,i,record,status)
          print *,record
      end do
!  Print out an END record, and a blank line to mark the end of the header.
      if (status .eq. 0)then
          print *,'END'
          print *,' '
      end if
!  Try moving to the next extension in the FITS file, if it exists.
!  The FTMRHD subroutine attempts to move to the next HDU, as specified by
!  the second parameter.   This subroutine moves by a relative number of
!  HDUs from the current HDU.  The related FTMAHD routine may be used to
!  move to an absolute HDU number in the FITS file.  If the end-of-file is
!  encountered when trying to move to the specified extension, then a
!  status = 107 is returned.
      call ftmrhd(unit,1,hdutype,status)
      if (status .eq. 0)then
!         success, so jump back and print out keywords in this extension
          go to 100
      else if (status .eq. 107)then
!         hit end of file, so quit
          status=0
      end if
!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)
!  Check for any error, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      if (status .gt. 0)call printerror(status)
  end

!--------------------------------------------------------------------------

!C *************************************************************************
      subroutine printerror(status)
!  This subroutine prints out the descriptive text corresponding to the
!  error status value and prints out the contents of the internal
!  error message stack generated by FITSIO whenever an error occurs.
      integer status
      character errtext*30,errmessage*80
!  Check if status is OK (no error); if so, simply return
      if (status .le. 0)return
!  The FTGERR subroutine returns a descriptive 30-character text string that
!  corresponds to the integer error status number.  A complete list of all
!  the error numbers can be found in the back of the FITSIO User's Guide.
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

!  FITSIO usually generates an internal stack of error messages whenever
!  an error occurs.  These messages provide much more information on the
!  cause of the problem than can be provided by the single integer error
!  status value.  The FTGMSG subroutine retrieves the oldest message from
!  the stack and shifts any remaining messages on the stack down one
!  position.  FTGMSG is called repeatedly until a blank message is
!  returned, which indicates that the stack is empty.  Each error message
!  may be up to 80 characters in length.  Another subroutine, called
!  FTCMSG, is available to simply clear the whole error message stack in
!  cases where one is not interested in the contents.
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do
      end
      subroutine readimage(filename, nx, nz, ny)
!  Read a FITS image and determine the minimum and maximum pixel value.
!  Rather than reading the entire image in
!  at once (which could require a very large array), the image is read
!  in pieces, 100 pixels at a time.  
      use arrays
      implicit none
      integer status,unit,readwrite,blocksize,naxes(3),nfound
      integer dim1, dim2, dim3
      integer, intent(in):: Nx, Ny, Nz
!      parameter( Nx = 512,  Ny = 324, Nz = 512)


      integer count2, count3
      integer group,firstpix,nbuffer,npixels,i
      real datamin,datamax,nullval,buffer(nx)
      logical anynull
      character filename*80

!  The STATUS parameter must always be initialized.
      status=0
!  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)
!  Open the FITS file previously created by WRITEIMAGE
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)
!  Determine the size of the image.
      call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
!  Check that it found both NAXIS1 and NAXIS2 keywords.
      if (nfound .ne. 3)then
          print *,'READIMAGE failed to read the NAXISn keywords.'
          return
       end if
!  Initialize variables
      npixels=naxes(1)*naxes(2)*naxes(3)
      group=1

     firstpix=1
      nullval=-999
      datamin=1.0E30
      datamax=-1.0E30
      count2 = 0
      count3 = 1
      do while (npixels .gt. 0)
          count2 = count2 + 1
          if ((count2) .gt. 324 ) then
              count3 = count3 +1
              count2 = 1
          endif

!         read up to 100 pixels at a time 
          nbuffer=min(512,npixels)

          call ftgpve(unit,group,firstpix,nbuffer,nullval, buffer,anynull,status)
!         find the min and max values
          do i=1,nbuffer
              datamin=min(datamin,buffer(i))
              datamax=max(datamax,buffer(i))
              bufffits(i, count2, count3) = buffer(i)
          end do
!         increment pointers and loop back to read the next group of pixels
          npixels=npixels-nbuffer
          firstpix=firstpix+nbuffer
      end do
      print *
      print *,'Min and max image pixels = ',datamin,datamax
!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)
!  Check for any error, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      if (status .gt. 0)call printerror(status)
      end



! *************************************************************************
!C *************************************************************************
      subroutine deletefile(filename,status)

!  A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*(*) filename

!  Simply return if status is greater than zero
      if (status .gt. 0)return

!  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)

!  Try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
!         file was opened;  so now delete it 
          call ftdelt(unit,status)
      else if (status .eq. 103)then
!         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
      else
!         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if
!  Free the unit number for later reuse
      call ftfiou(unit, status)
      end

! *************************************************************************

#endif






  end module 
