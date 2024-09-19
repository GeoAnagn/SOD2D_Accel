module elem_qua

   use mod_numerical_params
   use mod_maths
   use mod_nvtx

   implicit none
   integer(4), allocatable :: quad_order_edges(:,:)

contains

   subroutine init_quad_info()
      implicit none

      call nvtxStartRange('init_quad_info: CPU-DATA')
      allocate(quad_order_edges(4,2))
      call nvtxEndRange()

      call nvtxStartRange('init_quad_info: CPU-COMPUTE')
      quad_order_edges = transpose(reshape([1,2,2,3,3,4,4,1],(/2,4/)))
      call nvtxEndRange()

      call nvtxStartRange('init_quad_info: GPU-DATA')
      !$acc enter data copyin(quad_order_edges)
      call nvtxEndRange()
   end subroutine init_quad_info

   subroutine quad_highorder(mporder,mnpbou,xi,eta,atoIJ,N,dN) ! QUA16 element
      implicit none
      integer(4),intent(in) :: mporder,mnpbou
      real(rp),intent(in)   :: xi,eta
      integer(4),intent(in) :: atoIJ(mnpbou)
      real(rp),intent(out)  :: N(mnpbou), dN(2,mnpbou)
      real(rp)              :: xi_grid(mporder+1)

      call nvtxStartRange('quad_highorder: getGaussLobattoLegendre_roots')
      call getGaussLobattoLegendre_roots(mporder,xi_grid)
      call nvtxEndRange()

      call nvtxStartRange('quad_highorder: DoubleTensorProduct')
      call DoubleTensorProduct(mporder,mnpbou,xi_grid,xi,eta,atoIJ,N,dN)
      call nvtxEndRange()
   end subroutine quad_highorder

end module
