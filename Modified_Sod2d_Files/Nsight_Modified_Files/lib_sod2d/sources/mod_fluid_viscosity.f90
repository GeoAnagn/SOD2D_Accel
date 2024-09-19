module mod_fluid_viscosity

   use mod_numerical_params
   use mod_nvtx
   contains

   subroutine constant_viscosity(npoin,val,mu_fluid)

      implicit none

      integer(4), intent(in)    :: npoin
      real(rp)  , intent(in)    :: val
      real(rp)  , intent(inout) :: mu_fluid(npoin)

      call nvtxStartRange('constant_viscosity: GPU-DATA')
      !$acc kernels
      mu_fluid(:) = val
      !$acc end kernels
      call nvtxEndRange()
   end subroutine constant_viscosity

   subroutine sutherland_viscosity(npoin,Tem,mu_factor,mu_fluid)

      implicit none
      
      integer(4), intent(in)    :: npoin
      real(rp),   intent(in)    :: Tem(npoin)
      real(rp),   intent(in)    :: mu_factor(npoin)
      real(rp),   intent(inout) :: mu_fluid(npoin)
      integer(4)                :: iNodeL
      real(rp)                  :: C1, S

      call nvtxStartRange('sutherland_viscosity: CPU-DATA')
      C1 = 0.000001458_rp
      S = 110.4_rp
      call nvtxEndRange()

      call nvtxStartRange('sutherland_viscosity: GPU-COMPUTE')
      !$acc parallel loop
      do iNodeL = 1,npoin
         mu_fluid(iNodeL) = mu_factor(iNodeL)*C1*(Tem(iNodeL)**1.5_rp)/(Tem(iNodeL)+S)
      end do
      !$acc end parallel loop
      call nvtxEndRange()
   end subroutine sutherland_viscosity

end module mod_fluid_viscosity
