 subroutine create_netcdf(ncid, name, varname,  n1, n2, n3, ier)
 implicit none
 include 'netcdf.inc'

 character(80) name
 character(*) varname
 integer :: ier, ncid
 integer :: n1, n2, n3

 integer :: status

 integer :: dimx, dimy, dimz
 integer :: dims(3), varid 
 integer :: numberdims

 
 ier = 0 
 status = nf_create( name, NF_CLOBBER, ncid )


  
 status = nf_def_dim( ncid, 'nx', n1, dimx )
 if ( status .ne. NF_NOERR ) call  handle_error(status) 

 status = nf_def_dim( ncid, 'ny', n2, dimy )
 if ( status .ne. NF_NOERR ) call  handle_error(status)     

 if (n3 .eq. 1 ) then 
   numberdims = 2
 else 
   status = nf_def_dim( ncid, 'nz', n3, dimz )
   if ( status .ne. NF_NOERR ) call  handle_error(status)     
   numberdims = 3
 endif 

 dims(1) = dimx
 dims(2) = dimy
 if (numberdims .eq. 3)  dims(3) = dimz

 ! create variables

 status = nf_def_var( ncid, trim(varname), NF_DOUBLE, numberdims ,dims, varid )
 if ( status .ne. NF_NOERR ) call  handle_error(status)
 
! nf_enddef puts named file out of "define" mode
      status = nf_enddef(ncid)
      if ( status .ne. NF_NOERR ) call  handle_error(status)

! nf_sync synchronises disk writes with memory buffers (immediately 
! available after writing)
      status = nf_sync(ncid)
      if ( status .ne. NF_NOERR ) call  handle_error(status)


 end subroutine 



 subroutine close_netcdf(ncid, ier)
 implicit none

 include 'netcdf.inc'
 integer ier, ncid, status

 ier = 0

      status = nf_sync(ncid)
      if ( status .ne. NF_NOERR ) call  handle_error(status)

! Closes off named NetCDF file
      status = nf_close(ncid)
      if ( status .ne. NF_NOERR ) call handle_error(status)


      return
 

 end subroutine


 subroutine write_netcdf(ncid, myrank, nproc, varname,  var , nlayer,  n1, n2, n3, comm,  ier)

  implicit none

  include 'netcdf.inc'
#ifdef MPI
  include 'mpif.h'
  real(kind=8) dwk(nlayer,n2, n3)
  integer ncount, vartag

#endif
   integer myrank, comm, nproc
   integer nlayer
   integer ier, ncid, n1, n2, n3
   integer status
   integer varid
   integer icount3(3), istart3(3)

   integer i, j
   real(kind = 8) var(nlayer,n2,  n3)
   character(*) varname

   ier = 0

#ifdef MPI
   vartag = 11
   call MPI_Barrier(comm, ier) 

   if (myrank .eq. 0) then

     call gatherwrite(ncid, comm, nproc, var, dwk, trim(varname), nlayer,  n1, n2, n3, vartag, ier )
     print*, 'Finish gather write NetCDF'

   else
     ncount = n2*nlayer*n3

     call MPI_Send(var, ncount, MPI_double_precision, 0, vartag, comm, ier)

   endif
#else
   status = nf_inq_varid( ncid,  trim(varname) , varid )
   if ( status .ne. NF_NOERR) call handle_error(status)

   do i = 1, n3
      do j = 1, n2

         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = n1
         icount3(2) = 1
         icount3(3) = 1

          status = nf_put_vara_double(ncid,varid,istart3, icount3, var(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

         end do
      end do

   return
#endif


 end subroutine


 subroutine open_netcdf(ncid, name, n1, n2, n3, ier ) 

!
      implicit none
      include 'netcdf.inc'
      integer          ier
      character(len=*)   name
      integer::        n1, n2, n3 
      integer::        n1id, n2id, n3id
      integer::        status, varid, ncid
      integer::        numberdims

      ier = -1
      ncid = -1

! NF_OPEN reads named file, NF_WRITE implies file is writeable
! ncid is returned NetCDF file ID.

      status =  NF_OPEN( trim(name), NF_WRITE, ncid )
      if ( status .ne. NF_NOERR ) then
         print*,'ERROR: Unable to open file ', name
         call  handle_error(status)
         return
      end if



!.......Read the dimension id's
! nf_inq_dimid returns dimension ID (integer) of "string" in ncid

      status = nf_inq_dimid( ncid, 'nx', n1id )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimid( ncid, 'ny', n2id )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      if (n3 .eq. 1 ) then 
        numberdims = 2
      else 
        numberdims = 3
        status = nf_inq_dimid( ncid, 'nz', n3id )
        if ( status .ne. NF_NOERR ) call  handle_error(status)
      endif 

!.......Read the actual dimensions
! nf_inq_dimlen returns length of dimension
!
      status = nf_inq_dimlen( ncid, n1id, n1 )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimlen( ncid, n2id, n2 )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      if (numberdims .eq. 3) then 
      status = nf_inq_dimlen( ncid, n3id, n3 )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      end if 

      ier = 0

      return


 end subroutine 

 subroutine read_netcdf(ncid, myrank, nproc, varname,  var , nlayer,  n1, n2, n3, comm,  ier)
 
  implicit none

  include 'netcdf.inc'
#ifdef MPI
  include 'mpif.h'
  real(kind=8) dwk(nlayer,n2, n3)
  integer ncount, vartag
  integer mstat(MPI_STATUS_SIZE,1)
  integer ireq

#endif
   integer myrank, comm, nproc
   integer nlayer
   integer ier, ncid, n1, n2, n3
   integer n1id, n2id, n3id
   integer numberdims
   integer status
   integer varid
   integer icount3(3), istart3(3)

   integer i, j
   real(kind = 8) var(nlayer,n2,  n3)
   character(*) varname
 
   numberdims = 2
   ier = 0




   if (myrank .eq. 0) then 
!.......Read the dimension id's
! nf_inq_dimid returns dimension ID (integer) of "string" in ncid

      status = nf_inq_dimid( ncid, 'nx', n1id )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimid( ncid, 'ny', n2id )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      if (n3 .gt. 1) then 
        numberdims = 3 
        status = nf_inq_dimid( ncid, 'nz', n3id )
        if ( status .ne. NF_NOERR ) call  handle_error(status)
      endif 

      status = nf_inq_dimlen( ncid, n1id, n1 )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimlen( ncid, n2id, n2 )

      if (numberdims .eq. 3 )  then 
        if ( status .ne. NF_NOERR ) call  handle_error(status)
        status = nf_inq_dimlen( ncid, n3id, n3 )
      endif 


   end if 



#ifdef MPI 

   vartag = 11

   call MPI_BARRIER(comm, ier )

   if ( myrank .eq. 0 ) then

    call ReadBcast (ncid, comm, nproc, var, dwk, trim(varname), nlayer,  n1, n2, n3, vartag, ier)
    print*, 'Finished readin NetCDF in parallel'

   else 
    ncount = nlayer*n2*n3
  
    call  MPI_Irecv(var(1,1,1), ncount, MPI_double_precision, 0, vartag, comm, ireq, ier)
    call MPI_Waitall(1, ireq, mstat, ier)

   endif 


   return

#else



      status = nf_inq_varid( ncid, trim(varname), varid )
      if ( status .ne. NF_NOERR) call handle_error(status)




   do i = 1, n3
      do j = 1, n2

         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = n1
         icount3(2) = 1
         icount3(3) = 1


          status = nf_get_vara_double(ncid,varid,istart3, icount3, var(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)


         end do
      end do

    return

#endif 



 end subroutine 


 
 subroutine handle_error(status)
      implicit none

      integer  status
      include 'netcdf.inc'

      if ( status .ne. NF_NOERR ) then
        print*, NF_STRERROR(STATUS)
      end if
end  subroutine handle_error

#ifdef MPI
 subroutine gatherwrite(ncid, comm, nproc, var, dwk, varname, nlayer,  n1, n2, n3, vartag, ier)
  implicit none
  include 'mpif.h'
  include 'netcdf.inc'

  integer ncid, comm, nrpoc, n1, n2, nlayer, n3, ier
  integer vartag, varid, nproc
  integer status(MPI_STATUS_SIZE)

  real(kind=8) var( nlayer,n2, n3)
  real(kind=8) dwk( nlayer,n2, n3)

  character(*) varname

  integer i,j,k,l
  integer istart3(3), icount3(3), iproc, ncount
  integer istt

  ier = nf_inq_varid(ncid, trim(varname), varid)
  if ( ier .ne. NF_NOERR ) call  handle_error(ier)

! loop over all procs

  do iproc = 0, nproc-1

    istt = (iproc * n1)/nproc
    ncount = n2*nlayer*n3

   if(iproc .eq. 0) then
     do i = 1, nlayer
      do j = 1, n2
       do k = 1, n3
        dwk(i,j,k) = var(i,j,k)
       end do
      end do
     end do
   else
     call MPI_Recv( dwk, ncount, MPI_DOUBLE_PRECISION, iproc, vartag, comm, status, ier )
   end if

    do i = 1, n3
      do j = 1, n2

         istart3(1) = istt+1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = nlayer
         icount3(2) = 1
         icount3(3) = 1

          ier  = nf_put_vara_double(ncid,varid,istart3, icount3, dwk(1,j,i) )
          if ( ier .ne. NF_NOERR) call handle_error(ier)


         end do
      end do
! close loop over all procs
  end do

 end subroutine gatherwrite

 subroutine ReadBcast (ncid, comm, nproc, var, dwk, varname, nlayer,  n1, n2, n3, vartag, ier)

  implicit none
  include 'mpif.h'
  include 'netcdf.inc'

  integer ncid, comm, nrpoc, n1, n2, nlayer, n3, ier
  integer vartag, varid, nproc
  integer status(MPI_STATUS_SIZE)

  real(kind=8) var( nlayer,n2, n3)
  real(kind=8) dwk( nlayer,n2, n3)

  character(*) varname

  integer i,j,k,l
  integer istart3(3), icount3(3), iproc, ncount
  integer istt



     ier = nf_inq_varid( ncid, trim(varname), varid )
     if ( ier .ne. NF_NOERR ) call  handle_error(ier )

     do iproc=0,nproc-1
       istt = ((iproc)*n1)/nproc +1

       ncount = nlayer*n2*n3

       do k= 1, n3 
         do j=1,n2
           istart3(1) = istt
           istart3(2) = j
           istart3(3) = k
           icount3(1) = nlayer
           icount3(2) = 1
           icount3(3) = 1

             ier  = nf_get_vara_double(ncid, varid, istart3, icount3, dwk(1,j,k))
         end do
       end do



       if ( iproc .eq. 0 ) then
         do i= 1, nlayer
           do j=1,n2
             do k=1,n3
               var(i,j,k) = dwk(i,j,k)
             end do
           end do
         end do
       else
         call MPI_Send( dwk, ncount, MPI_DOUBLE_PRECISION, iproc, vartag, comm, ier )
       end if
     end do
 
   ier = 0 
   return



 end subroutine ReadBcast



#endif

 
