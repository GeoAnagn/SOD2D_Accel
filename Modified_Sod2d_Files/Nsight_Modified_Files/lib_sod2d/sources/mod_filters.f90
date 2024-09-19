module mod_filters
   use mod_constants
   use mod_nvtx

   implicit none

   ! integer ---------------------------------------------------
   integer(4), allocatable :: convertIJK(:)

   ! real ------------------------------------------------------
   real(rp), allocatable :: al_weights(:),am_weights(:),an_weights(:)

contains

   subroutine init_filters()
      implicit none
      integer(4) :: ii

      call nvtxStartRange('init_filters: CPU-DATA')
      allocate(al_weights(-1:1))
      call nvtxEndRange()

      call nvtxStartRange('init_filters: GPU-DATA')
      !$acc enter data create(al_weights(:))
      al_weights(-1) = 1.0_rp/4.0_rp
      al_weights(0)  = 2.0_rp/4.0_rp
      al_weights(1)  = 1.0_rp/4.0_rp
      !$acc update device(al_weights(:))
      call nvtxEndRange()

      call nvtxStartRange('init_filters: CPU-DATA')
      allocate(am_weights(-1:1))
      call nvtxEndRange()

      call nvtxStartRange('init_filters: GPU-DATA')
      !$acc enter data create(am_weights(:))
      am_weights(-1) = 1.0_rp/4.0_rp
      am_weights(0)  = 2.0_rp/4.0_rp
      am_weights(1)  = 1.0_rp/4.0_rp
      !$acc update device(am_weights(:))
      call nvtxEndRange()

      call nvtxStartRange('init_filters: CPU-DATA')
      allocate(an_weights(-1:1))
      call nvtxEndRange()

      call nvtxStartRange('init_filters: GPU-DATA')
      !$acc enter data create(an_weights(:))
      an_weights(-1) = 1.0_rp/4.0_rp
      an_weights(0)  = 2.0_rp/4.0_rp
      an_weights(1)  = 1.0_rp/4.0_rp
      !$acc update device(an_weights(:))
      call nvtxEndRange()

      call nvtxStartRange('init_filters: CPU-DATA')
      allocate(convertIJK(0:porder+2))
      call nvtxEndRange()

      call nvtxStartRange('init_filters: GPU-DATA')
      !$acc enter data create(convertIJK(:))
      do ii=3,porder+1
         convertIJK(ii-1) = ii
      end do
      convertIJK(0) = 3
      convertIJK(1) = 1
      convertIJK(porder+1) = 2
      convertIJK(porder+2) = porder
      !$acc update device(convertIJK(:))
      call nvtxEndRange()
   end subroutine init_filters

   subroutine deallocate_filters()

      call nvtxStartRange('deallocate_filters: GPU-DATA')
      !$acc exit data delete(al_weights(:))
      call nvtxEndRange()

      call nvtxStartRange('deallocate_filters: CPU-DATA')
      deallocate(al_weights)
      call nvtxEndRange()

      call nvtxStartRange('deallocate_filters: GPU-DATA')
      !$acc exit data delete(am_weights(:))
      call nvtxEndRange()

      call nvtxStartRange('deallocate_filters: CPU-DATA')
      deallocate(am_weights)
      call nvtxEndRange()

      call nvtxStartRange('deallocate_filters: GPU-DATA')
      !$acc exit data delete(an_weights(:))
      call nvtxEndRange()

      call nvtxStartRange('deallocate_filters: CPU-DATA')
      deallocate(an_weights)
      call nvtxEndRange()

      call nvtxStartRange('deallocate_filters: GPU-DATA')
      !$acc exit data delete(convertIJK(:))    
      call nvtxEndRange()

      call nvtxStartRange('deallocate_filters: CPU-DATA')
      deallocate(convertIJK)
      call nvtxEndRange()

   end subroutine deallocate_filters
end module mod_filters