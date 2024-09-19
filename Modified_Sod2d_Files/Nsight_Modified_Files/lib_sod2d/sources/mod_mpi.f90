module mod_mpi
   use mpi
   use mod_constants
   use mod_nvtx
#ifndef NOACC
   use openacc
#endif

   implicit none

   integer :: mpi_rank, mpi_size, mpi_err
   integer :: mpi_world_rank, mpi_world_size
   integer :: mpi_integer_size,mpi_int4_size,mpi_int8_size
   integer :: mpi_real_size,mpi_real4_size,mpi_real8_size
   integer :: mpi_datatype_int,mpi_datatype_int4,mpi_datatype_int8
   integer :: mpi_datatype_real,mpi_datatype_real4,mpi_datatype_real8
   integer :: mpi_datatype_real_int

   integer :: smNode_comm,smNode_rank,smNode_size
   integer :: num_devices, id_device
   integer :: app_comm, app_color, world_comm

contains

   subroutine init_mpi()
      implicit none

      integer(mpi_address_kind) :: color_ptr
      logical :: mpi_app_num_flag
      integer(4) :: sizeCustomType
      integer,dimension(:),allocatable :: blockLengthCustom, dataTypesCustom
      integer(mpi_address_kind),dimension(:),allocatable :: displacementsCustom

      ! Init MPI_COMM_WORLD communicator (shared accros all applications launched with mpirun MPMD)
      call nvtxStartRange("init_mpi: mpi_init")
      call mpi_init(mpi_err)
      call nvtxEndRange()
      
      world_comm = MPI_COMM_WORLD

      call nvtxStartRange("init_mpi: mpi_comm_rank")
      call mpi_comm_rank(world_comm, mpi_world_rank, mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("init_mpi: mpi_comm_size")
      call mpi_comm_size(world_comm, mpi_world_size, mpi_err)
      call nvtxEndRange()

      ! Get the app number (color)
      call nvtxStartRange("init_mpi: mpi_comm_get_attr")
      call MPI_Comm_get_attr(world_comm, MPI_APPNUM, color_ptr, mpi_app_num_flag, mpi_err)
      call nvtxEndRange()
      app_color = color_ptr ! necessary to get integer from pointer

      ! Split world_comm and create a communicator for this app only (color must be unique for each application)
      if (mpi_app_num_flag) then
         call nvtxStartRange("init_mpi: mpi_comm_split")
         call MPI_Comm_split(world_comm, app_color, mpi_world_rank, app_comm, mpi_err)
         call nvtxEndRange()

         call nvtxStartRange("init_mpi: mpi_comm_rank")
         call MPI_Comm_rank(app_comm, mpi_rank, mpi_err)
         call nvtxEndRange()

         call nvtxStartRange("init_mpi: mpi_comm_size")
         call MPI_Comm_size(app_comm, mpi_size, mpi_err)
         call nvtxEndRange()
      else
         write(*,*) 'Fatal error in init_mpi()! Cannot find MPI_APPNUM.'
         call nvtxStartRange("init_mpi: mpi_abort")
         call MPI_Abort(world_comm,-1,mpi_err)
         call nvtxEndRange()
      end if
      !write(*,982), app_color, mpi_rank, mpi_world_rank, mpi_app_num_flag
      !982 format ('app_color = ',I2,'; mpi_world_rank = ',I4,'; mpi_rank = ',I4,'; mpi_app_num_flag = ',L)

      mpi_datatype_int = MPI_INTEGER
      mpi_datatype_int4 = MPI_INTEGER4
      mpi_datatype_int8 = MPI_INTEGER8
      mpi_datatype_real4 = MPI_REAL4
      mpi_datatype_real8 = MPI_REAL8

      if(rp.eq.4) then
         mpi_datatype_real = MPI_REAL4
      else if(rp.eq.8) then
         mpi_datatype_real = MPI_REAL8
      else
         write(*,*) 'Fatal error in init_mpi()! rp is not 4 or 8 >> CRASH!'
         call nvtxStartRange("init_mpi: mpi_abort")
         call MPI_Abort(app_comm,-1,mpi_err)
         call nvtxEndRange()
      end if

      call nvtxStartRange("init_mpi: mpi_type_size")
      call MPI_Type_size(mpi_datatype_int,mpi_integer_size, mpi_err)
      call MPI_Type_size(mpi_datatype_int4,mpi_int4_size, mpi_err)
      call MPI_Type_size(mpi_datatype_int8,mpi_int8_size, mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("init_mpi: mpi_type_size")
      call MPI_Type_size(mpi_datatype_real,mpi_real_size,mpi_err)
      call MPI_Type_size(mpi_datatype_real4,mpi_real4_size,mpi_err)
      call MPI_Type_size(mpi_datatype_real8,mpi_real8_size,mpi_err)
      call nvtxEndRange()

      !Create MPI structure with [real(rp), integer]
      call nvtxStartRange("init_mpi: CPU-DATA")
      sizeCustomType = 2
      allocate(blockLengthCustom(sizeCustomType),dataTypesCustom(sizeCustomType),displacementsCustom(sizeCustomType))

      blockLengthCustom(:)  = 1
      displacementsCustom(1) = 0
      displacementsCustom(2) = mpi_real_size
      dataTypesCustom(1) = mpi_datatype_real
      dataTypesCustom(2) = mpi_datatype_int
      call nvtxEndRange()

      call nvtxStartRange("init_mpi: MPI_Type_create_struct")
      call MPI_Type_create_struct(sizeCustomType,blockLengthCustom,displacementsCustom,dataTypesCustom,mpi_datatype_real_int,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("init_mpi: MPI_Type_commit")
      call MPI_Type_commit(mpi_datatype_real_int,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("init_mpi: CPU-DATA")
      deallocate(blockLengthCustom,dataTypesCustom,displacementsCustom)
      call nvtxEndRange()

      call nvtxStartRange("init_mpi: init_sharedMemoryNode_comm")
      call init_sharedMemoryNode_comm()
      call nvtxEndRange()

   end subroutine init_mpi

   subroutine init_sharedMemoryNode_comm()
      implicit none
#ifndef NOACC
      call MPI_Comm_split_type(world_comm,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,smNode_comm,mpi_err)
      call mpi_comm_rank(smNode_comm, smNode_rank, mpi_err)
      call mpi_comm_size(smNode_comm, smNode_size, mpi_err)

      call nvtxStartRange("init_sharedMemoryNode_comm: CPU-DATA")
      num_devices = acc_get_num_devices(acc_device_nvidia)
      if(num_devices.ge.1) then

         id_device   = mod(smNode_rank,num_devices)
         call acc_set_device_num(id_device, acc_device_nvidia)
         write(*,736) mpi_world_rank,mpi_world_size,mpi_rank,mpi_size,smNode_rank,smNode_size,num_devices,id_device
736      format ('r(w) = ',I4,'; s(w) = ',I4,'; r(app) = ',I4,'; s(app) = ',I4,'; r(sm) = ',I4,'; s(sm) = ',I4,'; num_devices = ',I4,'; id_device = ' I4)

      else
         id_device = 0
         write(*,*) 'NO GPU FOUND IN THIS NODE!'
      end if
#else
      num_devices = 0
      id_device = 0
#endif   
      call nvtxEndRange()
   end subroutine init_sharedMemoryNode_comm

   subroutine end_mpi()
      call nvtxStartRange("end_mpi: mpi_finalize")
      call mpi_finalize(mpi_err)
      call nvtxEndRange()
   end subroutine end_mpi

end module mod_mpi
