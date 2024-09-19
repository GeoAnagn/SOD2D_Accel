module mod_comms
   use mod_mpi_mesh
   use mod_nvtx
#ifdef NCCL_COMMS
   use nccl
#endif

!-- Select type of communication for mpi_atomic_updates
#define _SENDRCV_ 0
#define _ISENDIRCV_ 1
#define _PUTFENCE_ 0

   implicit none

   !---- for the comms---
   integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
   integer(KIND=MPI_ADDRESS_KIND) :: memPos_t
   integer :: worldGroup,commGroup

   integer(4),dimension(:),allocatable :: aux_intField_s,  aux_intField_r
   real(rp),dimension(:),allocatable   :: aux_realField_s, aux_realField_r

   integer :: window_id_int,window_id_real
   integer :: window_id_sm
   integer :: beginFence=0,endFence=0
   integer :: startAssert=0,postAssert=0

   logical :: isInt,isReal
   logical :: isLockBarrier,isPSCWBarrier

#ifdef NCCL_COMMS
   integer                        :: cuda_stat
   type(ncclResult)               :: nccl_stat
   type(ncclUniqueId)             :: nccl_uid
   type(ncclComm)                 :: nccl_comm
   integer(kind=cuda_stream_kind) :: nccl_stream
#endif

contains

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
   subroutine init_comms(useInt,useReal)
      implicit none
      logical, intent(in) :: useInt,useReal
      logical :: useFenceFlags,useAssertNoCheckFlags,useLockBarrier

#if _SENDRCV_
      if(mpi_rank.eq.0) write(111,*) "--| Comm. scheme: Send-Recv"
#endif
#if _ISENDIRCV_
      if(mpi_rank.eq.0) write(111,*) "--| Comm. scheme: iSend-iRecv"
#endif
#if _PUTFENCE_
      if(mpi_rank.eq.0) write(111,*) "--| Comm. scheme: Put(Fence)"
#endif
#ifdef NCCL_COMMS
      if(mpi_rank.eq.0) write(111,*) "--| Comm. scheme: NCCL"
#endif

      isInt=.false.
      isReal=.false.

      if(useInt) then
         isInt = .true.

         call nvtxStartRange("init_comms: CPU-DATA")
         allocate(aux_intField_s(numNodesToComm))
         allocate(aux_intField_r(numNodesToComm))
         call nvtxEndRange()

         call nvtxStartRange("init_comms: GPU-DATA")
         !$acc enter data create(aux_intField_s(:))
         !$acc enter data create(aux_intField_r(:))
         call nvtxEndRange()

         call nvtxStartRange("init_comms: init_window_intField")
         call init_window_intField()
         call nvtxEndRange()
      end if

      if(useReal) then
         isReal = .true.

         call nvtxStartRange("init_comms: CPU-DATA")
         allocate(aux_realField_s(numNodesToComm))
         allocate(aux_realField_r(numNodesToComm))
         call nvtxEndRange()

         call nvtxStartRange("init_comms: GPU-DATA")
         !$acc enter data create(aux_realField_s(:))
         !$acc enter data create(aux_realField_r(:))
         call nvtxEndRange()

         call nvtxStartRange("init_comms: init_window_realField")
         call init_window_realField()
         call nvtxEndRange()
      end if

      call nvtxStartRange("init_comms: MPI_Comm_group")
      call MPI_Comm_group(app_comm,worldGroup,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("init_comms: MPI_Group_incl")
      call MPI_Group_incl(worldGroup,numRanksWithComms,ranksToComm,commGroup,mpi_err);
      call nvtxEndRange()

      useFenceFlags=.false. !by default
      useAssertNoCheckFlags=.true. !by default
      isPSCWBarrier=.true.!si faig molts loops amb aquesta opció a false la comm queda bloquejada
      useLockBarrier=.true.!.false. with false it fails!

      call nvtxStartRange("init_comms: setFenceFlags")
      call setFenceFlags(useFenceFlags)
      call nvtxEndRange()

      call nvtxStartRange("init_comms: setPSCWAssertNoCheckFlags")
      call setPSCWAssertNoCheckFlags(useAssertNoCheckFlags)
      call nvtxEndRange()

      call nvtxStartRange("init_comms: setLockBarrier")
      call setLockBarrier(useLockBarrier)
      call nvtxEndRange()

#ifdef NCCL_COMMS
      if (mpi_rank == 0) then
         call nvtxStartRange("init_comms: ncclGetUniqueId")
         nccl_stat = ncclGetUniqueId(nccl_uid)
         call nvtxEndRange()
      end if
      call nvtxStartRange("init_comms: MPI_Bcast")
      call MPI_Bcast(nccl_uid, int( sizeof(ncclUniqueId), kind = 4 ), MPI_BYTE, 0, app_comm, mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("init_comms: ncclCommInitRank")
      nccl_stat = ncclCommInitRank(nccl_comm, mpi_size, nccl_uid, mpi_rank);
      call nvtxEndRange()

      call nvtxStartRange("init_comms: cudaStreamCreate")
      cuda_stat = cudaStreamCreate(nccl_stream);
      call nvtxEndRange()
#endif

   end subroutine init_comms

   subroutine end_comms()
      implicit none

      if(isInt) then
         call nvtxStartRange("end_comms: GPU-DATA")
         !$acc exit data delete(aux_intField_s(:))
         !$acc exit data delete(aux_intField_r(:))
         call nvtxEndRange()

         call nvtxStartRange("end_comms: CPU-DATA")
         deallocate(aux_intField_s)
         deallocate(aux_intField_r)
         call nvtxEndRange()

         call nvtxStartRange("end_comms: close_window_intField")
         call close_window_intField()
         call nvtxEndRange()
      end if

      if(isreal) then
         call nvtxStartRange("end_comms: GPU-DATA")
         !$acc exit data delete(aux_realField_s(:))
         !$acc exit data delete(aux_realField_r(:))
         call nvtxEndRange()

         call nvtxStartRange("end_comms: CPU-DATA")
         deallocate(aux_realField_s)
         deallocate(aux_realField_r)
         call nvtxEndRange()

         call nvtxStartRange("end_comms: close_window_realField")
         call close_window_realField()
         call nvtxEndRange()
      end if

#ifdef NCCL_COMMS
      call nvtxStartRange("end_comms: ncclCommDestroy")
      nccl_stat = ncclCommDestroy(nccl_comm)
      call nvtxEndRange()

      call nvtxStartRange("end_comms: cudaStreamDestroy")
      cuda_stat = cudaStreamDestroy(nccl_stream)
      call nvtxEndRange()
#endif
   end subroutine end_comms

   subroutine setFenceFlags(useFenceFlags)
      implicit none
      logical,intent(in) :: useFenceFlags

      call nvtxStartRange("setFenceFlags: CPU-DATA")
      if(useFenceFlags) then
         beginFence = MPI_MODE_NOPRECEDE
         endFence   = IOR(MPI_MODE_NOSTORE,IOR(MPI_MODE_NOPUT,MPI_MODE_NOSUCCEED))
         !endFence   = IOR(MPI_MODE_NOSTORE,MPI_MODE_NOPUT)
      else
         beginFence = 0
         endFence   = 0
      end if
      call nvtxEndRange()
   end subroutine

   subroutine setPSCWAssertNoCheckFlags(useAssertNoCheckFlags)
      implicit none
      logical,intent(in) :: useAssertNoCheckFlags

      call nvtxStartRange("setPSCWAssertNoCheckFlags: CPU-DATA")
      if(useAssertNoCheckFlags) then
         postAssert  = MPI_MODE_NOCHECK
         startAssert = MPI_MODE_NOCHECK
      else
         postAssert  = 0
         startAssert = 0
      endif
      call nvtxEndRange()
   end subroutine

   subroutine setLockBarrier(useLockBarrier)
      implicit none
      logical,intent(in) :: useLockBarrier
      
      call nvtxStartRange("setLockBarrier: CPU-DATA")
      isLockBarrier = useLockBarrier
      call nvtxEndRange()
   end subroutine
!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
   subroutine init_window_intField()
      implicit none

      call nvtxStartRange("init_window_intField: CPU-COMPUTE")
      window_buffer_size = mpi_integer_size*numNodesToComm
      call nvtxEndRange()

      call nvtxStartRange("init_window_intField: MPI_Win_create")
      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      call MPI_Win_create(aux_intField_r,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id_int,mpi_err)
      !$acc end host_data
      call nvtxEndRange()
   end subroutine init_window_intField

   subroutine close_window_intField()
      implicit none
      call nvtxStartRange("close_window_intField: MPI_Win_free")
      call MPI_Win_free(window_id_int,mpi_err)
      call nvtxEndRange()
   end subroutine close_window_intField
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
   subroutine init_window_realField()
      implicit none

      call nvtxStartRange("init_window_realField: CPU-COMPUTE")
      window_buffer_size = mpi_real_size*numNodesToComm
      call nvtxEndRange()

      call nvtxStartRange("init_window_realField: MPI_Win_create")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      call MPI_Win_create(aux_realField_r,window_buffer_size,mpi_real_size,MPI_INFO_NULL,app_comm,window_id_real,mpi_err)
      !$acc end host_data
      call nvtxEndRange()
   end subroutine init_window_realField

   subroutine close_window_realField()
      implicit none

      call nvtxStartRange("close_window_realField: MPI_Win_free")
      call MPI_Win_free(window_id_real,mpi_err)
      call nvtxEndRange()
   end subroutine close_window_realField

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
!    SUBROUTINES COPY/FROM SEND/RCV BUFFERS
!-----------------------------------------------------------------------------------------------------------------------
   subroutine fill_sendBuffer_int(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,iNodeL

      call nvtxStartRange("fill_sendBuffer_int: GPU-COMPUTE")
      !$acc parallel loop
      do i=1,numNodesToComm
         iNodeL = nodesToComm(i)
         aux_intField_s(i) = intField(iNodeL)
      end do
      !$acc end parallel loop
      call nvtxEndRange()

      call nvtxStartRange("fill_sendBuffer_int: GPU-DATA")
      !$acc kernels
      aux_intField_r(:)=0
      !$acc end kernels
      call nvtxEndRange()
   end subroutine fill_sendBuffer_int
!-------------------------------------------------------------------------
   subroutine fill_sendBuffer_real(realField)
      implicit none
      real(rp),intent(inout) :: realField(:)
      integer(4) :: i,iNodeL

      call nvtxStartRange("fill_sendBuffer_real: GPU-COMPUTE")
      !$acc parallel loop
      do i=1,numNodesToComm
         iNodeL = nodesToComm(i)
         aux_realField_s(i) = realField(iNodeL)
      end do
      !$acc end parallel loop
      call nvtxEndRange()

      call nvtxStartRange("fill_sendBuffer_real: GPU-DATA")
      !$acc kernels
      aux_realField_r(:)=0.
      !$acc end kernels
      call nvtxEndRange()
   end subroutine fill_sendBuffer_real
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
   subroutine fill_sendBuffer_get_int(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,iNodeL

      call nvtxStartRange("fill_sendBuffer_get_int: GPU-COMPUTE")
      !$acc parallel loop
      do i=1,numNodesToComm
         iNodeL = nodesToComm(i)
         aux_intField_r(i) = intField(iNodeL)
      end do
      !$acc end parallel loop
      call nvtxEndRange()

      call nvtxStartRange("fill_sendBuffer_get_int: GPU-DATA")
      !$acc kernels
      aux_intField_s(:)=0
      !$acc end kernels
      call nvtxEndRange()
   end subroutine fill_sendBuffer_get_int
!-------------------------------------------------------------------------
   subroutine fill_sendBuffer_get_real(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,iNodeL

      call nvtxStartRange("fill_sendBuffer_get_real: GPU-COMPUTE")
      !$acc parallel loop
      do i=1,numNodesToComm
         iNodeL = nodesToComm(i)
         aux_realField_r(i) = realField(iNodeL)
      end do
      !$acc end parallel loop
      call nvtxEndRange()

      call nvtxStartRange("fill_sendBuffer_get_real: GPU-DATA")
      !$acc kernels
      aux_realField_s(:)=0.
      !$acc end kernels
      call nvtxEndRange()
   end subroutine fill_sendBuffer_get_real
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
   subroutine copy_from_rcvBuffer_int(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,iNodeL

      call nvtxStartRange("copy_from_rcvBuffer_int: GPU-COMPUTE")
      !$acc parallel loop
      do i=1,numNodesToComm
         iNodeL = nodesToComm(i)
         !$acc atomic update
         intField(iNodeL) = intField(iNodeL) + aux_intField_r(i)
         !$acc end atomic
      end do
      !$acc end parallel loop
      call nvtxEndRange()
   end subroutine copy_from_rcvBuffer_int
!-------------------------------------------------------------------------
   subroutine copy_from_rcvBuffer_real(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,iNodeL

      call nvtxStartRange("copy_from_rcvBuffer_real: GPU-COMPUTE")
      !$acc parallel loop
      do i=1,numNodesToComm
         iNodeL = nodesToComm(i)
         !$acc atomic update
         realField(iNodeL) = realField(iNodeL) + aux_realField_r(i)
         !$acc end atomic
      end do
      !$acc end parallel loop
      call nvtxEndRange()
   end subroutine copy_from_rcvBuffer_real
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
   subroutine copy_from_min_rcvBuffer_real(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,iNodeL

      call nvtxStartRange("copy_from_min_rcvBuffer_real: GPU-COMPUTE")
      !$acc parallel loop
      do i=1,numNodesToComm
         iNodeL = nodesToComm(i)
         !$acc atomic update
         realField(iNodeL) = min(realField(iNodeL) , aux_realField_r(i))
         !$acc end atomic
      end do
      !$acc end parallel loop
      call nvtxEndRange()
   end subroutine copy_from_min_rcvBuffer_real
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
   subroutine copy_from_rcvBuffer_get_int(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,iNodeL

      call nvtxStartRange("copy_from_rcvBuffer_get_int: GPU-COMPUTE")
      !$acc parallel loop
      do i=1,numNodesToComm
         iNodeL = nodesToComm(i)
         !$acc atomic update
         intField(iNodeL) = intField(iNodeL) + aux_intField_s(i)
         !$acc end atomic
      end do
      !$acc end parallel loop
      call nvtxEndRange()
   end subroutine copy_from_rcvBuffer_get_int
!-------------------------------------------------------------------------
   subroutine copy_from_rcvBuffer_get_real(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,iNodeL

      call nvtxStartRange("copy_from_rcvBuffer_get_real: GPU-COMPUTE")
      !$acc parallel loop
      do i=1,numNodesToComm
         iNodeL = nodesToComm(i)
         !$acc atomic update
         realField(iNodeL) = realField(iNodeL) + aux_realField_s(i)
         !$acc end atomic
      end do
      !$acc end parallel loop
      call nvtxEndRange()
   end subroutine copy_from_rcvBuffer_get_real
!-----------------------------------------------------------------------------------------------------------------------
   subroutine copy_from_conditional_ave_rcvBuffer_real(cond,realField)
      implicit none
      real(rp),intent(in) :: cond
      real(rp),intent(inout) :: realField(:)
      integer(4) :: i,iNodeL

      call nvtxStartRange("copy_from_conditional_ave_rcvBuffer_real: GPU-COMPUTE")
      !$acc parallel loop
      do i=1,numNodesToComm
         iNodeL = nodesToComm(i)
         if(abs(aux_realField_r(i)).gt. cond) then
            if(abs(realField(iNodeL)).gt. cond) then
               !$acc atomic update
               realField(iNodeL) = realField(iNodeL)+aux_realField_r(i)
               !$acc end atomic
            else
               !$acc atomic write
               realField(iNodeL) = aux_realField_r(i)
               !$acc end atomic
            end if
         end if
      end do
      !$acc end parallel loop
      call nvtxEndRange()
   end subroutine copy_from_conditional_ave_rcvBuffer_real
!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------

   subroutine mpi_halo_atomic_update_int(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
#if _SENDRCV_
      call nvtxStartRange("mpi_halo_atomic_update_int: mpi_halo_atomic_update_int_sendRcv")
      call mpi_halo_atomic_update_int_sendRcv(intField)
      call nvtxEndRange()
#endif
#if _ISENDIRCV_
      call nvtxStartRange("mpi_halo_atomic_update_int: mpi_halo_atomic_update_int_iSendiRcv")
      call mpi_halo_atomic_update_int_iSendiRcv(intField)
      call nvtxEndRange()
#endif
#if _PUTFENCE_
      call nvtxStartRange("mpi_halo_atomic_update_int: mpi_halo_atomic_update_int_put_fence")
      call mpi_halo_atomic_update_int_put_fence(intField)
      call nvtxEndRange()
#endif
   end subroutine mpi_halo_atomic_update_int

   subroutine mpi_halo_atomic_update_real(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)

#if _SENDRCV_
      call nvtxStartRange("mpi_halo_atomic_update_real: mpi_halo_atomic_update_real_sendRcv")
      call mpi_halo_atomic_update_real_sendRcv(realField)
      call nvtxEndRange()
#endif
#if _ISENDIRCV_
      call nvtxStartRange("mpi_halo_atomic_update_real: mpi_halo_atomic_update_real_iSendiRcv")
      call mpi_halo_atomic_update_real_iSendiRcv(realField)
      call nvtxEndRange()
#endif
#if _PUTFENCE_
      call nvtxStartRange("mpi_halo_atomic_update_real: mpi_halo_atomic_update_real_put_fence")
      call mpi_halo_atomic_update_real_put_fence(realField)
      call nvtxEndRange()
#endif

   end subroutine mpi_halo_atomic_update_real

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
!------------- SEND/RECV -------------------------------------------
   ! INTEGER ---------------------------------------------------
   subroutine mpi_halo_atomic_update_int_sendRcv(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,ngbRank,tagComm
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_int_sendRcv: fill_sendBuffer_int")
      call fill_sendBuffer_int(intField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_sendRcv: MPI_Sendrecv")
      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         tagComm  = 0
         memPos_l = commsMemPosInLoc(i)
         memSize  = commsMemSize(i)

         call MPI_Sendrecv(aux_intField_s(mempos_l), memSize, mpi_datatype_int, ngbRank, tagComm, &
            aux_intField_r(mempos_l), memSize, mpi_datatype_int, ngbRank, tagComm, &
            app_comm, MPI_STATUS_IGNORE, mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_sendRcv: copy_from_rcvBuffer_int")
      call copy_from_rcvBuffer_int(intField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_int_sendRcv
   ! REAL ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_sendRcv(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,ngbRank,tagComm
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_real_sendRcv: fill_sendBuffer_real")
      call fill_sendBuffer_real(realField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_sendRcv: MPI_Sendrecv")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         tagComm  = 0
         memPos_l = commsMemPosInLoc(i)
         memSize  = commsMemSize(i)

         call MPI_Sendrecv(aux_realField_s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
            aux_realField_r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
            app_comm, MPI_STATUS_IGNORE, mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_sendRcv: copy_from_rcvBuffer_real")
      call copy_from_rcvBuffer_real(realField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_real_sendRcv

   ! REAL ---------------------------------------------------
   subroutine mpi_halo_atomic_min_update_real_sendRcv(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,ngbRank,tagComm
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_min_update_real_sendRcv: fill_sendBuffer_real")
      call fill_sendBuffer_real(realField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_min_update_real_sendRcv: MPI_Sendrecv")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         tagComm  = 0
         memPos_l = commsMemPosInLoc(i)
         memSize  = commsMemSize(i)

         call MPI_Sendrecv(aux_realField_s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
            aux_realField_r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
            app_comm, MPI_STATUS_IGNORE, mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_min_update_real_sendRcv: copy_from_min_rcvBuffer_real")
      call copy_from_min_rcvBuffer_real(realField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_min_update_real_sendRcv


!------------- ISEND/IRECV -------------------------------------------
   !INTEGER ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_iSendiRcv(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,ireq,ngbRank,tagComm
      integer(4) :: memPos_l,memSize
      integer(4) :: requests(2*numRanksWithComms)

      call nvtxStartRange("mpi_halo_atomic_update_int_iSendiRcv: fill_sendBuffer_int")
      call fill_sendBuffer_int(intField)
      call nvtxEndRange()

      ireq=0
      call nvtxStartRange("mpi_halo_atomic_update_int_iSendiRcv: MPI_Irecv_MPI_ISend")
      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         tagComm  = 0
         memPos_l = commsMemPosInLoc(i)
         memSize  = commsMemSize(i)

         ireq = ireq+1
         call MPI_Irecv(aux_intField_r(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
         ireq = ireq+1
         call MPI_ISend(aux_intField_s(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_iSendiRcv: MPI_Waitall")
      call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_iSendiRcv: copy_from_rcvBuffer_int")
      call copy_from_rcvBuffer_int(intField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_int_iSendiRcv
   !REAL ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_iSendiRcv(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,ireq,ngbRank,tagComm
      integer(4) :: memPos_l,memSize
      integer(4) :: requests(2*numRanksWithComms)

      call nvtxStartRange("mpi_halo_atomic_update_real_iSendiRcv: fill_sendBuffer_real")
      call fill_sendBuffer_real(realField)
      call nvtxEndRange()

#if NCCL_COMMS
      call nvtxStartRange("mpi_halo_atomic_update_real_iSendiRcv: ncclGroupStart")
      nccl_stat = ncclGroupStart()
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_iSendiRcv: ncclRecv_ncclSend")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = commsMemPosInLoc(i)
         memSize  = commsMemSize(i)

         nccl_stat = ncclRecv(aux_realField_r(mempos_l), memSize, ncclFloat, ngbRank, nccl_comm, nccl_stream)
         nccl_stat = ncclSend(aux_realField_s(mempos_l), memSize, ncclFloat, ngbRank, nccl_comm, nccl_stream)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_iSendiRcv: ncclGroupEnd")
      nccl_stat = ncclGroupEnd()
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_iSendiRcv: cudaStreamSynchronize")
      cuda_stat = cudaStreamSynchronize(nccl_stream)
      call nvtxEndRange()
#else
      ireq=0
      call nvtxStartRange("mpi_halo_atomic_update_real_iSendiRcv: MPI_Irecv_MPI_ISend")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         tagComm  = 0
         memPos_l = commsMemPosInLoc(i)
         memSize  = commsMemSize(i)

         ireq = ireq+1
         call MPI_Irecv(aux_realField_r(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
         ireq = ireq+1
         call MPI_ISend(aux_realField_s(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_iSendiRcv: MPI_Waitall")
      call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)
      call nvtxEndRange()
#endif

      call nvtxStartRange("mpi_halo_atomic_update_real_iSendiRcv: copy_from_rcvBuffer_real")
      call copy_from_rcvBuffer_real(realField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_real_iSendiRcv

!------------- conditional average ISEND/IRECV -------------------------------------------
   ! REAL ---------------------------------------------------
   subroutine mpi_halo_conditional_ave_update_real_iSendiRcv(cond,realField)
      implicit none
      real(rp),intent(in) :: cond
      real(rp),intent(inout) :: realField(:)
      integer(4) :: i,ireq,ngbRank,tagComm
      integer(4) :: memPos_l,memSize
      integer(4) :: requests(2*numRanksWithComms)

      call nvtxStartRange("mpi_halo_conditional_ave_update_real_iSendiRcv: fill_sendBuffer_real")
      call fill_sendBuffer_real(realField)
      call nvtxEndRange()

      ireq=0
      call nvtxStartRange("mpi_halo_conditional_ave_update_real_iSendiRcv: MPI_Irecv_MPI_ISend")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         tagComm  = 0
         memPos_l = commsMemPosInLoc(i)
         memSize  = commsMemSize(i)

         ireq = ireq+1
         call MPI_Irecv(aux_realfield_r(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
         ireq = ireq+1
         call MPI_ISend(aux_realfield_s(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_conditional_ave_update_real_iSendiRcv: MPI_Waitall")
      call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_conditional_ave_update_real_iSendiRcv: copy_from_conditional_ave_rcvBuffer_real")
      call copy_from_conditional_ave_rcvBuffer_real(cond,realField)
      call nvtxEndRange()
   end subroutine mpi_halo_conditional_ave_update_real_iSendiRcv
!------------- PUT FENCE -------------------------------------------
   !INT
   subroutine mpi_halo_atomic_update_int_put_fence(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,ngbRank
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_int_put_fence: fill_sendBuffer_int")
      call fill_sendBuffer_int(intField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_fence: MPI_Win_fence")
      call MPI_Win_fence(beginFence,window_id_int,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_fence: MPI_Put")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = commsMemPosInLoc(i)
         memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
         memSize  = commsMemSize(i)

         call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_fence: MPI_Win_fence")
      call MPI_Win_fence(endFence,window_id_int,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_fence: copy_from_rcvBuffer_int")
      call copy_from_rcvBuffer_int(intField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_int_put_fence
   !REAL
   subroutine mpi_halo_atomic_update_real_put_fence(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,ngbRank
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_real_put_fence: fill_sendBuffer_real")
      call fill_sendBuffer_real(realField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_fence: MPI_Win_fence")
      call MPI_Win_fence(beginFence,window_id_real,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_fence: MPI_Put")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = commsMemPosInLoc(i)
         memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
         memSize  = commsMemSize(i)

         call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_fence: MPI_Win_fence")
      call MPI_Win_fence(endFence,window_id_real,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_fence: copy_from_rcvBuffer_real")
      call copy_from_rcvBuffer_real(realField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_real_put_fence
!------------- PUT PSCW -------------------------------------------
   !INTEGER--------------------------------------------------
   subroutine mpi_halo_atomic_update_int_put_pscw(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,ngbRank
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_int_put_pscw: fill_sendBuffer_int")
      call fill_sendBuffer_int(intField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_pscw: MPI_Barrier")
      if(isPSCWBarrier) call MPI_Barrier(app_comm,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_pscw: MPI_Win_post")
      call MPI_Win_post(commGroup,postAssert,window_id_int,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_pscw: MPI_Win_start")
      call MPI_Win_start(commGroup,startAssert,window_id_int,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_pscw: MPI_Put")
      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = commsMemPosInLoc(i)
         memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
         memSize  = commsMemSize(i)

         call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_pscw: MPI_Win_complete")
      call MPI_Win_complete(window_id_int,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_pscw: MPI_Win_wait")
      call MPI_Win_wait(window_id_int,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_pscw: copy_from_rcvBuffer_int")
      call copy_from_rcvBuffer_int(intField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_int_put_pscw

   !REAL-------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_put_pscw(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,ngbRank
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_real_put_pscw: fill_sendBuffer_real")
      call fill_sendBuffer_real(realField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_pscw: MPI_Barrier")
      if(isPSCWBarrier) call MPI_Barrier(app_comm,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_pscw: MPI_Win_post")
      call MPI_Win_post(commGroup,postAssert,window_id_real,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_pscw: MPI_Win_start")
      call MPI_Win_start(commGroup,startAssert,window_id_real,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_pscw: MPI_Put")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = commsMemPosInLoc(i)
         memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
         memSize  = commsMemSize(i)

         call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_pscw: MPI_Win_complete")
      call MPI_Win_complete(window_id_real,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_pscw: MPI_Win_wait")
      call MPI_Win_wait(window_id_real,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_pscw: copy_from_rcvBuffer_real")
      call copy_from_rcvBuffer_real(realField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_real_put_pscw

!------------- PUT LOCK -------------------------------------------
   !INTEGER-------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_put_lock(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,ngbRank
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_int_put_lock: fill_sendBuffer_int")
      call fill_sendBuffer_int(intField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_lock: MPI_Win_lock_put_unlock")
      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = commsMemPosInLoc(i)
         memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
         memSize  = commsMemSize(i)

         call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_int,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
         call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
         call MPI_Win_unlock(ngbRank,window_id_int,mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_lock: MPI_Barrier")
      if(isLockBarrier) call MPI_Barrier(app_comm,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_put_lock: copy_from_rcvBuffer_int")
      call copy_from_rcvBuffer_int(intField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_int_put_lock

   !REAL-------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_put_lock(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,ngbRank
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_real_put_lock: fill_sendBuffer_real")
      call fill_sendBuffer_real(realField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_lock: MPI_Win_lock_put_unlock")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = commsMemPosInLoc(i)
         memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
         memSize  = commsMemSize(i)

         call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_real,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
         call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
         call MPI_Win_unlock(ngbRank,window_id_real,mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_lock: MPI_Barrier")
      if(isLockBarrier) call MPI_Barrier(app_comm,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_put_lock: copy_from_rcvBuffer_real")
      call copy_from_rcvBuffer_real(realField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_real_put_lock

!------------- GET FENCE -------------------------------------------
   !INTEGER-------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_get_fence(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,ngbRank
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_int_get_fence: fill_sendBuffer_get_int")
      call fill_sendBuffer_get_int(intField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_fence: MPI_Win_fence")
      call MPI_Win_fence(beginFence,window_id_int,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_fence: MPI_Get")
      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = commsMemPosInLoc(i)
         memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
         memSize  = commsMemSize(i)

         call MPI_Get(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_fence: MPI_Win_fence")
      call MPI_Win_fence(endFence,window_id_int,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_fence: copy_from_rcvBuffer_get_int")
      call copy_from_rcvBuffer_get_int(intField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_int_get_fence

   !REAL-------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_get_fence(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,ngbRank
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_real_get_fence: fill_sendBuffer_get_real")
      call fill_sendBuffer_get_real(realField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_fence: MPI_Win_fence")
      call MPI_Win_fence(beginFence,window_id_real,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_fence: MPI_Get")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = commsMemPosInLoc(i)
         memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
         memSize  = commsMemSize(i)

         call MPI_Get(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_fence: MPI_Win_fence")
      call MPI_Win_fence(endFence,window_id_real,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_fence: copy_from_rcvBuffer_get_real")
      call copy_from_rcvBuffer_get_real(realField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_real_get_fence

!------------- GET PSCW -------------------------------------------
   !INTEGER-------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_get_pscw(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,ngbRank
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_int_get_pscw: fill_sendBuffer_get_int")
      call fill_sendBuffer_get_int(intField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_pscw: MPI_Barrier")
      if(isPSCWBarrier) call MPI_Barrier(app_comm,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_pscw: MPI_Win_post")
      call MPI_Win_post(commGroup,postAssert,window_id_int,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_pscw: MPI_Win_start")
      call MPI_Win_start(commGroup,startAssert,window_id_int,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_pscw: MPI_Get")
      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank   = ranksToComm(i)
         memPos_l  = commsMemPosInLoc(i)
         memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
         memSize   = commsMemSize(i)

         call MPI_Get(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_pscw: MPI_Win_complete")
      call MPI_Win_complete(window_id_int,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_pscw: MPI_Win_wait")
      call MPI_Win_wait(window_id_int,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_pscw: copy_from_rcvBuffer_get_int")
      call copy_from_rcvBuffer_get_int(intField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_int_get_pscw
   !REAL-------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_get_pscw(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,ngbRank
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_real_get_pscw: fill_sendBuffer_get_real")
      call fill_sendBuffer_get_real(realField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_pscw: MPI_Barrier")
      if(isPSCWBarrier) call MPI_Barrier(app_comm,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_pscw: MPI_Win_post")
      call MPI_Win_post(commGroup,postAssert,window_id_real,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_pscw: MPI_Win_start")
      call MPI_Win_start(commGroup,startAssert,window_id_real,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_pscw: MPI_Get")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank   = ranksToComm(i)
         memPos_l  = commsMemPosInLoc(i)
         memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
         memSize   = commsMemSize(i)

         call MPI_Get(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_pscw: MPI_Win_complete")
      call MPI_Win_complete(window_id_real,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_pscw: MPI_Win_wait")
      call MPI_Win_wait(window_id_real,mpi_err);
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_pscw: copy_from_rcvBuffer_get_real")
      call copy_from_rcvBuffer_get_real(realField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_real_get_pscw

!------------- GET LOCK -------------------------------------------
   !INTEGER-------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_get_lock(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,ngbRank
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_int_get_lock: fill_sendBuffer_get_int")
      call fill_sendBuffer_get_int(intField)
      call nvtxEndRange()
      
      call nvtxStartRange("mpi_halo_atomic_update_int_get_lock: MPI_Win_lock_get_unlock")
      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = commsMemPosInLoc(i)
         memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
         memSize  = commsMemSize(i)

         call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_int,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
         call MPI_Get(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
         call MPI_Win_unlock(ngbRank,window_id_int,mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_lock: MPI_Barrier")
      if(isLockBarrier) call MPI_Barrier(app_comm,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_get_lock: copy_from_rcvBuffer_get_int")
      call copy_from_rcvBuffer_get_int(intField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_int_get_lock
   !REAL-------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_get_lock(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,ngbRank
      integer(4) :: memPos_l,memSize

      call nvtxStartRange("mpi_halo_atomic_update_real_get_lock: fill_sendBuffer_get_real")
      call fill_sendBuffer_get_real(realField)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_lock: MPI_Win_lock_get_unlock")
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = commsMemPosInLoc(i)
         memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
         memSize  = commsMemSize(i)

         call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_real,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
         call MPI_Get(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
         call MPI_Win_unlock(ngbRank,window_id_real,mpi_err)
      end do
      !$acc end host_data
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_lock: MPI_Barrier")
      if(isLockBarrier) call MPI_Barrier(app_comm,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_get_lock: copy_from_rcvBuffer_get_real")
      call copy_from_rcvBuffer_get_real(realField)
      call nvtxEndRange()
   end subroutine mpi_halo_atomic_update_real_get_lock

!------------- ONLY BUFFERS -------------------------------------------
!for testing and devel stuff
   !INTEGER ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_onlybuffers(intfield)
      implicit none
      integer(4), intent(inout) :: intfield(:)
      integer(4) :: i,ngbrank,tagcomm
      integer(4) :: mempos_l,memsize

      call nvtxStartRange("mpi_halo_atomic_update_int_onlybuffers: fill_sendbuffer_int")
      call fill_sendbuffer_int(intfield)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_int_onlybuffers: copy_from_rcvbuffer_int")
      call copy_from_rcvbuffer_int(intfield)
      call nvtxEndRange()

   end subroutine mpi_halo_atomic_update_int_onlybuffers
   !REAL ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_onlybuffers(realfield)
      implicit none
      real(rp), intent(inout) :: realfield(:)
      integer(4) :: i,ngbrank,tagcomm
      integer(4) :: mempos_l,memsize

      call nvtxStartRange("mpi_halo_atomic_update_real_onlybuffers: fill_sendbuffer_real")
      call fill_sendbuffer_real(realfield)
      call nvtxEndRange()

      call nvtxStartRange("mpi_halo_atomic_update_real_onlybuffers: copy_from_rcvbuffer_real")
      call copy_from_rcvbuffer_real(realfield)
      call nvtxEndRange()

   end subroutine mpi_halo_atomic_update_real_onlybuffers

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------

end module mod_comms
