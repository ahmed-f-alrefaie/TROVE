module MPI_module
#if ((usempi .and. .not.(intelmpi)) .or. (usempi .and. intelmpi .and. .not.(mpi64bitinteger)))
   !The MPI library is not Intel MPI or it is Intel MPI but using 32bit integers.
   use mpi
#endif
   use accuracy, only: out, redirect_stdout_to_file
   implicit none
#if ((usempi .and. (intelmpi) .and. mpi64bitinteger))
!We're using Intel MPI using 64bit integers so we have to take MPI constants from the mpif.h file rather than from the mpi module (see https://software.intel.com/en-us/node/528844).
#include <mpif.h>
#endif   


   private

   !> This can be set only once at the beginning of the program run by a call to MPI_MOD_START. The value of this variable is accessed by MPI_MOD_CHECK which aborts the program if MPI is not running.
   logical, public, protected :: mpi_running = .false.

   !> Total number of processes. Set by mpi_mod_start.
   integer(kind=mpiint), public, protected :: nprocs = -1

   !> The local process number (rank). Set by mpi_mod_start.
   integer(kind=mpiint), public, protected :: myrank = -1

   !> Total number of processes on the node this task is bound to. Set by mpi_mod_start.
   integer(kind=mpiint), public, protected :: local_nprocs = 1

   !> The local process number (rank) on the node this task is bound to. Set by mpi_mod_start.
   integer(kind=mpiint), public, protected :: local_rank = 0

   !> The routine mpi_start creates the shared_communicator which groups all tasks sitting on the same node. If we're using MPI 3.0 standard the shared_communicator enables creation of
   !> shared memory areas. In this case shared_enabled is set to .true. If we're not using MPI 3.0 then shared_enabled is set to .false.
   logical, public, protected :: shared_enabled = .false.

   !> ID of the master process.
   !> \warning This must be set to 0 so don't change it.
   integer(kind=mpiint), parameter, public :: master = 0

   !> This variable is used only in max_data_count below but it MUST be 32 bit integer.
   integer(kind=ik), parameter :: dummy_32bit_integer = 1
#if usempi
   !> Type value corresponding to the default integer type. Set by mpi_mod_start.
   integer(kind=mpiint) :: mpi_mod_int = -1

   !> Type value corresponding to the MPI integer type. Set by mpi_mod_start.
   integer(kind=mpiint) :: mpi_mod_mpiint = -1

   !> MPI type value corresponding to the default (cfp) float. Set by mpi_mod_start.
   integer(kind=mpiint) :: mpi_mod_cfp = -1

   !> MPI type value corresponding to the double precision float. Set by mpi_mod_start.
   integer(kind=mpiint) :: mpi_mod_wp = -1

   !> MPI type value corresponding to the quad precision float. Set by mpi_mod_start.
   integer(kind=mpiint) :: mpi_mod_ep = -1

   !> Name of the processor on which the current process is running.
   character(len=MPI_MAX_PROCESSOR_NAME), public, protected :: procname

   !> Intra-node communicator created by mpi_mod_start.
   integer(kind=mpiint), public, protected :: shared_communicator

   !> Integer identifying the node the MPI task belongs to.
   integer(kind=mpiint) :: node_colour
  
   !> The largest allowed data type count (i.e. the number of array elements) that can be passed to the MPI routines. This limitation exists in Intel MPI but perhaps also in other MPI libraries?
   integer, parameter :: max_data_count = huge(dummy_32bit_integer)-2
#else
   character(len=*), parameter :: procname = "N/A"
#endif

contains



   !> Interface to the routine MPI_BARRIER.
   subroutine mpi_mod_barrier(error)
      implicit none
      integer :: error
#if usempi
      integer(kind=mpiint) :: ierr

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         error = ierr
#endif

   end subroutine mpi_mod_barrier

   !> This is a lightweight routine which aborts the program if MPI is not found running.
   subroutine check_mpi_running
      implicit none

#if usempi
         if (.not.(mpi_running)) stop "MPI not running, Abort"
#endif

   end subroutine check_mpi_running


   !> Initializes MPI, assigns the rank for each process and finds out how many processors we're using. It also maps the current floating point precision (kind=cfp) to the corresponding MPI numeric type.
   !> This routine also sets the unit for standard output that each rank will use (the value of stdout). In case of serial run the value stdout is not set here and is left to the default value input_unit
   !> as specified in the module const.
   !> \warning This must the first statement in every level3 program. 
   subroutine mpi_mod_start(do_stdout)
      implicit none
      logical,optional,intent(in)	::	do_stdout
#if usempi
      integer(kind=mpiint) :: ierr=0, error=0, isize=0, zero=0
      integer :: myint
      character(len=MPI_MAX_ERROR_STRING) :: estring
      integer(kind=mpiint) :: error_class,error2,error_length
      logical :: do_stdout_
        do_stdout_ = .false.
	if(present(do_stdout)) do_stdout_ = do_stdout


         if (mpi_running) then
   
            stop "Attempt to start MPI while it has been initialized before. The program will abort."
            call MPI_ABORT(MPI_COMM_WORLD,ierr,error)
   
         else

            !Initialize MPI without threading
            call MPI_INIT(error)

            if (error .ne. MPI_SUCCESS) then
               stop "MPI initialization has failed. The program will abort."
               call MPI_Error_class(error, error_class,error2);
               call MPI_Error_string(error, estring,error_length,error2);
               write(*,*) estring(1:error_length)
               write(*,*) 'Error class',error_class
               call MPI_ABORT(MPI_COMM_WORLD,ierr,error)
            else

               !Determine basic properties of the parallelism: number of processes and rank of my MPI process.
               call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,error)
               if (error .ne. MPI_SUCCESS) then
                  stop "MPI_COMM_SIZE has failed. The program will abort."
                  call MPI_ABORT(MPI_COMM_WORLD,ierr,error)
               endif

               call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,error)
               if (error .ne. MPI_SUCCESS) then
                  stop "MPI_COMM_RANK has failed. The program will abort."
                  call MPI_ABORT(MPI_COMM_WORLD,ierr,error)
               endif

               !This is where we redirect all standard output to a file. By standard output we mean all output using the write statement where the unit number used is 'stdout'.
               !Following the call to redirect_stdout_to_file the value of the variable const/stdout will be associated with a file intended for
               !standard output for the process with rank=myrank. Finally, remember that redirect_stdout_to_file can be called only once!
               call redirect_stdout_to_file(myrank,do_stdout_)

               call MPI_GET_PROCESSOR_NAME(procname,isize,error)
               if (error .ne. MPI_SUCCESS) then
                  stop "MPI_GET_PROCESSOR_NAME has failed. The program will abort."
                  call MPI_ABORT(MPI_COMM_WORLD,ierr,error)
               endif

               mpi_mod_wp = MPI_DOUBLE_PRECISION
   
               !Determine the MPI type corresponding to the default real floating point numbers.

               mpi_mod_cfp = MPI_DOUBLE_PRECISION

   
               !Type of the default integer.
               isize = bit_size(myint)/8
               call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER,isize,mpi_mod_int,ierr)

               !Type of the MPI integer.
               mpi_mod_mpiint = MPI_INTEGER
   
               mpi_running = .true.

               write(out,'(10X,/,"MPI running with ",i," processes.")') nprocs
               write(out,'(10X,"My rank is ",i)') myrank
               write(out,'(10X,"I am running on the processor with name: ",a,/)') procname
#if (mpithree)
               shared_enabled = .true.

               CALL MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,   MPI_INFO_NULL, shared_communicator,error)
               call MPI_COMM_SIZE(shared_communicator,local_nprocs,error)
               call MPI_COMM_RANK(shared_communicator,local_rank,error)
               write(out,'(10X,"MPI-3.0 standard used: Using shared memory")')
               write(out,'(10X,"Number of MPI tasks on this node ",i)') local_nprocs
               write(out,'(10X,"My rank on this node is ",i)') local_rank
#else
               shared_enabled = .false.
               shared_communicator = MPI_COMM_WORLD
               local_nprocs = nprocs
               local_rank = myrank
#endif

            endif
   
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#else

         nprocs = 1
         myrank = 0

         !This is where we redirect all standard output to a file. Remember that redirect_stdout_to_file can be called only once!
         call redirect_stdout_to_file(myrank)
#endif
         write(stdout,'("--------->","done:mpi_mod:mpi_mod_start")')

   end subroutine mpi_mod_start



end module
