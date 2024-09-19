
module time_integ

   use mod_nvtx
   use elem_convec
   use elem_diffu
   use elem_source
   use mod_solver
   use mod_entropy_viscosity
   use mod_numerical_params
   use mod_fluid_viscosity
   use mod_sgs_viscosity
   use mod_sgs_ilsa_viscosity
   use mod_bc_routines
   use mod_wall_model

   implicit none

   real(rp), allocatable, dimension(:)     :: Rmass,Rener,Reta
   real(rp), allocatable, dimension(:,:)   :: Rmom
   real(rp), allocatable, dimension(:)   :: aux_rho, aux_pr, aux_E, aux_Tem, aux_e_int,aux_eta,aux_h
   real(rp), allocatable, dimension(:,:) :: aux_u,aux_q
   real(rp), allocatable, dimension(:)   :: Rmass_sum,Rener_sum,Reta_sum,Rdiff_mass,Rdiff_ener
   real(rp), allocatable, dimension(:,:) :: Rmom_sum,Rdiff_mom,f_eta

   real(rp), allocatable, dimension(:)   :: a_i, b_i, c_i,b_i2
   real(rp), allocatable, dimension(:,:) :: a_ij

   contains

   subroutine init_rk4_solver(npoin)
      implicit none
      integer(4),intent(in) :: npoin
      integer(4) :: numSteps

      call nvtxStartRange('init_rk4_solver: CPU-DATA')
      allocate(Rmass(npoin),Rener(npoin),Reta(npoin),Rmom(npoin,ndime))
      call nvtxEndRange()

      call nvtxStartRange('init_rk4_solver: GPU-DATA')
      !$acc enter data create(Rmass(:))
      !$acc enter data create(Rener(:))
      !$acc enter data create(Reta(:))
      !$acc enter data create(Rmom(:,:))
      call nvtxEndRange()

      call nvtxStartRange('init_rk4_solver: CPU-DATA')
      allocate(aux_rho(npoin),aux_pr(npoin),aux_E(npoin),aux_Tem(npoin),aux_e_int(npoin),aux_eta(npoin),aux_h(npoin))
      call nvtxEndRange()

      call nvtxStartRange('init_rk4_solver: GPU-DATA')
      !$acc enter data create(aux_rho(:))
      !$acc enter data create(aux_pr(:))
      !$acc enter data create(aux_E(:))
      !$acc enter data create(aux_Tem(:))
      !$acc enter data create(aux_e_int(:))
      !$acc enter data create(aux_eta(:))
      !$acc enter data create(aux_h(:))
      call nvtxEndRange()

      call nvtxStartRange('init_rk4_solver: CPU-DATA')
      allocate(aux_u(npoin,ndime),aux_q(npoin,ndime))
      call nvtxEndRange()

      call nvtxStartRange('init_rk4_solver: GPU-DATA')
      !$acc enter data create(aux_u(:,:))
      !$acc enter data create(aux_q(:,:))
      call nvtxEndRange()

      call nvtxStartRange('init_rk4_solver: CPU-DATA')
      allocate(Rmass_sum(npoin),Rener_sum(npoin),Reta_sum(npoin),Rdiff_mass(npoin),Rdiff_ener(npoin))
      call nvtxEndRange()

      call nvtxStartRange('init_rk4_solver: GPU-DATA')
      !$acc enter data create(Rmass_sum(:))
      !$acc enter data create(Rener_sum(:))
      !$acc enter data create(Reta_sum(:))
      !$acc enter data create(Rdiff_mass(:))
      !$acc enter data create(Rdiff_ener(:))
      call nvtxEndRange()

      call nvtxStartRange('init_rk4_solver: CPU-DATA')
      allocate(Rmom_sum(npoin,ndime),Rdiff_mom(npoin,ndime),f_eta(npoin,ndime))
      call nvtxEndRange()

      call nvtxStartRange('init_rk4_solver: GPU-DATA')
      !$acc enter data create(Rmom_sum(:,:))
      !$acc enter data create(Rdiff_mom(:,:))
      !$acc enter data create(f_eta(:,:))
      call nvtxEndRange()

      call nvtxStartRange('init_rk4_solver: CPU-DATA')
      allocate(a_i(4),b_i(4),c_i(4))
      call nvtxEndRange()

      call nvtxStartRange('init_rk4_solver: GPU-DATA')
      !$acc enter data create(a_i(:))
      !$acc enter data create(b_i(:))
      !$acc enter data create(c_i(:))
      call nvtxEndRange()

      if (flag_rk_order == 1) then
         a_i = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
         c_i = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
         b_i = [1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
      else if (flag_rk_order == 2) then
         a_i = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
         c_i = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
         b_i = [0.5_rp, 0.5_rp, 0.0_rp, 0.0_rp]
      else if (flag_rk_order == 3) then
         write(1,*) "--| NOT CODED FOR RK3 YET!"
         stop 1
      else if (flag_rk_order == 4) then
         a_i = [0.0_rp, 0.5_rp, 0.5_rp, 1.0_rp]
         c_i = [0.0_rp, 0.5_rp, 0.5_rp, 1.0_rp]
         b_i = [1.0_rp/6.0_rp, 1.0_rp/3.0_rp, 1.0_rp/3.0_rp, 1.0_rp/6.0_rp]
      else
         write(1,*) "--| NOT CODED FOR RK > 4 YET!"
         stop 1
      end if

      call nvtxStartRange('init_rk4_solver: GPU-DATA')
      !$acc update device(a_i(:))
      !$acc update device(b_i(:))
      !$acc update device(c_i(:))
      call nvtxEndRange()

   end subroutine init_rk4_solver

   subroutine end_rk4_solver()
      implicit none

      call nvtxStartRange('end_rk4_solver: GPU-DATA')
      !$acc exit data delete(Rmass(:))
      !$acc exit data delete(Rener(:))
      !$acc exit data delete(Reta(:))
      !$acc exit data delete(Rmom(:,:))
      call nvtxEndRange()

      call nvtxStartRange('end_rk4_solver: CPU-DATA')
      deallocate(Rmass,Rener,Reta,Rmom)
      call nvtxEndRange()

      call nvtxStartRange('end_rk4_solver: GPU-DATA')
      !$acc exit data delete(aux_rho(:))
      !$acc exit data delete(aux_pr(:))
      !$acc exit data delete(aux_E(:))
      !$acc exit data delete(aux_Tem(:))
      !$acc exit data delete(aux_e_int(:))
      !$acc exit data delete(aux_eta(:))
      !$acc exit data delete(aux_h(:))
      call nvtxEndRange()

      call nvtxStartRange('end_rk4_solver: CPU-DATA')
      deallocate(aux_rho,aux_pr,aux_E,aux_Tem,aux_e_int,aux_eta,aux_h)
      call nvtxEndRange()

      call nvtxStartRange('end_rk4_solver: GPU-DATA')
      !$acc exit data delete(aux_u(:,:))
      !$acc exit data delete(aux_q(:,:))
      call nvtxEndRange()

      call nvtxStartRange('end_rk4_solver: CPU-DATA')
      deallocate(aux_u,aux_q)
      call nvtxEndRange()

      call nvtxStartRange('end_rk4_solver: GPU-DATA')
      !$acc exit data delete(Rmass_sum(:))
      !$acc exit data delete(Rener_sum(:))
      !$acc exit data delete(Reta_sum(:))
      !$acc exit data delete(Rdiff_mass(:))
      !$acc exit data delete(Rdiff_ener(:))
      call nvtxEndRange()

      call nvtxStartRange('end_rk4_solver: CPU-DATA')
      deallocate(Rmass_sum,Rener_sum,Reta_sum,Rdiff_mass,Rdiff_ener)
      call nvtxEndRange()

      call nvtxStartRange('end_rk4_solver: GPU-DATA')
      !$acc exit data delete(Rmom_sum(:,:))
      !$acc exit data delete(Rdiff_mom(:,:))
      !$acc exit data delete(f_eta(:,:))
      call nvtxEndRange()

      call nvtxStartRange('end_rk4_solver: CPU-DATA')
      deallocate(Rmom_sum,Rdiff_mom,f_eta)
      call nvtxEndRange()

      call nvtxStartRange('end_rk4_solver: GPU-DATA')
      !$acc exit data delete(a_i(:))
      !$acc exit data delete(b_i(:))
      !$acc exit data delete(c_i(:))
      call nvtxEndRange()

      call nvtxStartRange('end_rk4_solver: CPU-DATA')
      deallocate(a_i,b_i,c_i)
      call nvtxEndRange()
   end subroutine end_rk4_solver

         subroutine rk_4_main(noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor,mue_l, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term,walave_u,zo)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn_nodes(npoin)
            integer(4),           intent(in)    :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode),gmshAtoJ(nnode),gmshAtoK(nnode)
            integer(4),           intent(in)    :: ppow
            real(rp),             intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),             intent(in)    :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime)
            real(rp),             intent(in)    :: gpvol(1,ngaus,nelem)
            real(rp),             intent(in)    :: dt, helem(nelem)
            real(rp),             intent(in)    :: helem_l(nelem,nnode)
            real(rp),             intent(in)    :: Ml(npoin)
            real(rp),             intent(in)    :: mu_factor(npoin)
            real(rp),             intent(in)    :: Rgas, gamma_gas, Cp, Prt
            real(rp),             intent(inout) :: rho(npoin,3)
            real(rp),             intent(inout) :: u(npoin,ndime,3)
            real(rp),             intent(inout) :: q(npoin,ndime,3)
            real(rp),             intent(inout) :: pr(npoin,3)
            real(rp),             intent(inout) :: E(npoin,3)
            real(rp),             intent(inout) :: Tem(npoin,2)
            real(rp),             intent(inout) :: e_int(npoin,2)
            real(rp),             intent(inout) :: eta(npoin,3)
            real(rp),             intent(inout) :: mu_fluid(npoin)
            real(rp),             intent(inout) :: csound(npoin)
            real(rp),             intent(inout) :: machno(npoin)
            real(rp),             intent(inout) :: mu_e(nelem,ngaus)
            real(rp),             intent(inout) :: mu_sgs(nelem,ngaus)
            real(rp),             intent(inout) :: kres(npoin)
            real(rp),             intent(inout) :: etot(npoin)
            real(rp),             intent(inout) :: au(npoin,ndime)
            real(rp),             intent(inout) :: ax1(npoin)
            real(rp),             intent(inout) :: ax2(npoin)
            real(rp),             intent(inout) :: ax3(npoin)
            real(rp),             intent(inout) :: mue_l(nelem,nnode)
            real(rp),             intent(in)    :: coord(npoin,ndime)
            real(rp),             intent(in)  ::  wgp(ngaus)
            integer(4),            intent(in)    :: numBoundsWM
            integer(4), optional, intent(in)    :: ndof, nbnodes, ldof(*), lbnodes(*)
            integer(4), optional, intent(in)    :: bound(nboun,npbou), bou_codes(nboun), bou_codes_nodes(npoin)
            integer(4), optional, intent(in)    :: listBoundsWM(*)
            real(rp), optional, intent(in)      :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
            real(rp), optional,   intent(in)    :: u_buffer(npoin,ndime)
            real(rp), optional,   intent(inout) :: tauw(npoin,ndime)
            real(rp), optional, intent(in)      :: source_term(npoin,ndime)
            real(rp), optional, intent(in)      :: walave_u(npoin,ndime)
            real(rp), optional, intent(in)      :: zo(npoin)
            integer(4)                          :: pos
            integer(4)                          :: istep, ipoin, idime,icode
            real(rp),    dimension(npoin)       :: Rrho
            real(rp)                            :: umag


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! New version of RK4 using loops                 !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            pos = 2 ! Set correction as default value

            !
            ! Initialize variables to zero
            !
            call nvtxStartRange('rk_4_main: GPU-DATA')
            !$acc kernels
            aux_rho(1:npoin) = 0.0_rp
            aux_u(1:npoin,1:ndime) = 0.0_rp
            aux_q(1:npoin,1:ndime) = 0.0_rp
            aux_pr(1:npoin) = 0.0_rp
            aux_E(1:npoin) = 0.0_rp
            aux_Tem(1:npoin) = 0.0_rp
            aux_e_int(1:npoin) = 0.0_rp
            aux_eta(1:npoin) = 0.0_rp
            Rdiff_mass(1:npoin) = 0.0_rp
            Rdiff_mom(1:npoin,1:ndime) = 0.0_rp
            Rdiff_ener(1:npoin) = 0.0_rp
            Rmass(1:npoin) = 0.0_rp
            Rmom(1:npoin,1:ndime) = 0.0_rp
            Rener(1:npoin) = 0.0_rp
            Reta(1:npoin) = 0.0_rp
            Rmass_sum(1:npoin) = 0.0_rp
            Rener_sum(1:npoin) = 0.0_rp
            Reta_sum(1:npoin) = 0.0_rp
            Rmom_sum(1:npoin,1:ndime) = 0.0_rp
            !$acc end kernels
            call nvtxEndRange()

            !
            ! Loop over all RK steps
            !

            do istep = 1,flag_rk_order
               !
               ! Compute variable at substep (y_i = y_n+dt*A_ij*R_j)
               !
               call nvtxStartRange('rk_4_main: GPU-COMPUTE')
               !$acc parallel loop
               do ipoin = 1,npoin
                  aux_rho(ipoin) = rho(ipoin,pos) - dt*a_i(istep)*Rmass(ipoin)
                  aux_E(ipoin)   = E(ipoin,pos)   - dt*a_i(istep)*Rener(ipoin)
                  !$acc loop seq
                  do idime = 1,ndime
                     aux_q(ipoin,idime) = q(ipoin,idime,pos) - dt*a_i(istep)*Rmom(ipoin,idime)
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange()

               call nvtxStartRange('rk_4_main: updateBuffer')
               if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,aux_rho,aux_q,aux_E,u_buffer)
               call nvtxEndRange()

               !
               ! Apply bcs after update
               !
               call nvtxStartRange('rk_4_main: temporary_bc_routine_dirichlet_prim')
               if (noBoundaries .eqv. .false.) then
                  call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,aux_rho(:),aux_q(:,:),aux_u(:,:),aux_pr(:),aux_E(:),u_buffer)  
               end if
               call nvtxEndRange()

               !
               ! Update velocity and equations of state
               !
               call nvtxStartRange('rk_4_main: GPU-COMPUTE')
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  umag = 0.0_rp
                  !$acc loop seq
                  do idime = 1,ndime
                     aux_u(lpoin_w(ipoin),idime) = aux_q(lpoin_w(ipoin),idime)/aux_rho(lpoin_w(ipoin))
                     umag = umag + (aux_u(lpoin_w(ipoin),idime)*aux_u(lpoin_w(ipoin),idime))
                  end do
                  aux_e_int(lpoin_w(ipoin)) = (aux_E(lpoin_w(ipoin))/aux_rho(lpoin_w(ipoin)))-0.5_rp*umag
                  aux_pr(lpoin_w(ipoin)) = aux_rho(lpoin_w(ipoin))*(gamma_gas-1.0_rp)*aux_e_int(lpoin_w(ipoin))
                  aux_h(lpoin_w(ipoin)) = (gamma_gas/(gamma_gas-1.0_rp))*aux_pr(lpoin_w(ipoin))/aux_rho(lpoin_w(ipoin))
                  aux_Tem(lpoin_w(ipoin)) = aux_pr(lpoin_w(ipoin))/(aux_rho(lpoin_w(ipoin))*Rgas)
                  aux_eta(lpoin_w(ipoin)) = (aux_rho(lpoin_w(ipoin))/(gamma_gas-1.0_rp))* &
                  log(aux_pr(lpoin_w(ipoin))/(aux_rho(lpoin_w(ipoin))**gamma_gas))
                  !$acc loop seq
                  do idime = 1,ndime
                     f_eta(lpoin_w(ipoin),idime) = aux_u(lpoin_w(ipoin),idime)*aux_eta(lpoin_w(ipoin))
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange()

               call nvtxStartRange('rk_4_main: generic_scalar_convec_ijk')
               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                  gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,aux_eta(:),aux_u(:,:),Reta(:))
               call nvtxEndRange()

               call nvtxStartRange('rk_4_main: mpi_halo_atomic_update_real')
               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(Reta(:))
               end if
               call nvtxEndRange()

               call nvtxStartRange('rk_4_main: bc_fix_dirichlet_residual_entropy')
               if (noBoundaries .eqv. .false.) then
                  call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,Reta)
               end if
               call nvtxEndRange()

               call nvtxStartRange('rk_4_main: lumped_solver_scal')
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta)
               call nvtxEndRange()

               !
               ! Compute viscosities and diffusion
               !
               !
               ! Update viscosity if Sutherland's law is active
               !
               call nvtxStartRange('rk_4_main: sutherland_viscosity')
               if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
                  call sutherland_viscosity(npoin,aux_Tem,mu_factor,mu_fluid)
               end if
               call nvtxEndRange()
               ! 
               ! Compute diffusion terms with values at current substep
               !
               call nvtxStartRange('rk_4_main: full_diffusion_ijk')
               call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,aux_rho,aux_rho,aux_u,aux_Tem,mu_fluid,mu_e,mu_sgs,Ml,Rdiff_mass,Rdiff_mom,Rdiff_ener)
               call nvtxEndRange()
               !
               ! Call source term if applicable
               !
               call nvtxStartRange('rk_4_main: mom_source_const_vect')
               if(present(source_term)) then
                  call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_u,source_term,Rdiff_mom)
               end if
               call nvtxEndRange()
               !
               ! Evaluate wall models

               if((isWallModelOn) .and. (numBoundsWM .ne. 0)) then
                  if(flag_walave) then
                     if(flag_type_wmles == wmles_type_reichardt) then
                        call nvtxStartRange('rk_4_main: evalWallModelReichardt')
                        call evalWallModelReichardt(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,aux_rho(:),walave_u(:,:),tauw,Rdiff_mom)
                        call nvtxEndRange()
                     else if (flag_type_wmles == wmles_type_abl) then
                        call nvtxStartRange('rk_4_main: evalWallModelABL')
                        call evalWallModelABL(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                                             bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                                             aux_rho(:),walave_u(:,:),zo,tauw,Rdiff_mom)
                        call nvtxEndRange()
                     end if   
                  else
                     if(flag_type_wmles == wmles_type_reichardt) then
                        call nvtxStartRange('rk_4_main: evalWallModelReichardt')
                        call evalWallModelReichardt(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,aux_rho(:),aux_u(:,:),tauw,Rdiff_mom)
                        call nvtxEndRange()
                     else
                        write(1,*) "--| Only Reichardt wall model can work without time filtering!"
                        stop 1
                     end if
                  end if
               end if
               !
               !
               ! Compute convective terms
               !
               if(flag_total_enthalpy .eqv. .true.) then
                  call nvtxStartRange('rk_4_main: full_convec_ijk_H')
                  call full_convec_ijk_H(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,aux_u,aux_q,aux_rho,aux_pr,aux_E,Rmass(:),Rmom(:,:),Rener(:))
                  call nvtxEndRange()
               else 
                  call nvtxStartRange('rk_4_main: full_convec_ijk')
                  call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,aux_u,aux_q,aux_rho,aux_pr,aux_h,Rmass(:),Rmom(:,:),Rener(:))
                  call nvtxEndRange()
               end if

               call nvtxStartRange('rk_4_main: GPU-DATA')
               !$acc kernels
               Rmass(:) = Rmass(:) + Rdiff_mass(:)
               Rener(:) = Rener(:) + Rdiff_ener(:)
               Rmom(:,:) = Rmom(:,:) + Rdiff_mom(:,:)
               !$acc end kernels
               call nvtxEndRange()

               !TESTING NEW LOCATION FOR MPICOMMS
               call nvtxStartRange('rk_4_main: mpi_halo_atomic_update_real')
               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(Rmass(:))
                  call mpi_halo_atomic_update_real(Rener(:))
                  do idime = 1,ndime
                     call mpi_halo_atomic_update_real(Rmom(:,idime))
                  end do
               end if
               call nvtxEndRange()

               !
               ! Call lumped mass matrix solver
               !
               call nvtxStartRange('rk_4_main: lumped_solver_scal')
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass(:))
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener(:))
               call nvtxEndRange()

               call nvtxStartRange('rk_4_main: lumped_solver_vect')
               call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,Rmom(:,:))
               call nvtxEndRange()
               !
               ! Accumulate the residuals
               !
               call nvtxStartRange('rk_4_main: GPU-COMPUTE')
               !$acc parallel loop
               do ipoin = 1,npoin
                  Rmass_sum(ipoin) = Rmass_sum(ipoin) + b_i(istep)*Rmass(ipoin)
                  Rener_sum(ipoin) = Rener_sum(ipoin) + b_i(istep)*Rener(ipoin)
                  Reta_sum(ipoin)  = Reta_sum(ipoin)  + b_i(istep)*Reta(ipoin)
                  !$acc loop seq
                  do idime = 1,ndime
                     Rmom_sum(ipoin,idime) = Rmom_sum(ipoin,idime) + b_i(istep)*Rmom(ipoin,idime)
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange()
            end do
            !
            ! RK update to variables
            !
            call nvtxStartRange('rk_4_main: GPU-COMPUTE')
            !$acc parallel loop
            do ipoin = 1,npoin
               rho(ipoin,pos) = rho(ipoin,pos)-dt*Rmass_sum(ipoin)
               E(ipoin,pos) = E(ipoin,pos)-dt*Rener_sum(ipoin)
               !$acc loop seq
               do idime = 1,ndime
                  q(ipoin,idime,pos) = q(ipoin,idime,pos)-dt*Rmom_sum(ipoin,idime)
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange()

            call nvtxStartRange('rk_4_main: updateBuffer')
            if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,rho(:,pos),q(:,:,pos),E(:,pos),u_buffer)
            call nvtxEndRange()
            !
            ! Apply bcs after update
            !
            call nvtxStartRange('rk_4_main: temporary_bc_routine_dirichlet_prim')
            if (noBoundaries .eqv. .false.) then
               call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,rho(:,pos),q(:,:,pos),u(:,:,pos),pr(:,pos),E(:,pos),u_buffer)
            end if
            call nvtxEndRange()
            !
            ! Update velocity and equations of state
            !

            call nvtxStartRange('rk_4_main: GPU-COMPUTE')
            !$acc parallel loop
            do ipoin = 1,npoin_w
               umag = 0.0_rp
               !$acc loop seq
               do idime = 1,ndime
                  u(lpoin_w(ipoin),idime,pos) = q(lpoin_w(ipoin),idime,pos)/rho(lpoin_w(ipoin),pos)
                  umag = umag + u(lpoin_w(ipoin),idime,pos)**2
               end do
               e_int(lpoin_w(ipoin),pos) = (E(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos))- &
                  0.5_rp*umag
               pr(lpoin_w(ipoin),pos) = rho(lpoin_w(ipoin),pos)*(gamma_gas-1.0_rp)*e_int(lpoin_w(ipoin),pos)
               csound(lpoin_w(ipoin)) = sqrt(gamma_gas*pr(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos))
               umag = sqrt(umag)
               machno(lpoin_w(ipoin)) = umag/csound(lpoin_w(ipoin))
               Tem(lpoin_w(ipoin),pos) = pr(lpoin_w(ipoin),pos)/(rho(lpoin_w(ipoin),pos)*Rgas)
               eta(lpoin_w(ipoin),pos) = (rho(lpoin_w(ipoin),pos)/(gamma_gas-1.0_rp))* &
                  log(pr(lpoin_w(ipoin),pos)/(rho(lpoin_w(ipoin),pos)**gamma_gas))
            end do
            !$acc end parallel loop

            !$acc parallel loop
            do ipoin = 1,npoin_w
               Reta(lpoin_w(ipoin)) =  -Reta_sum(lpoin_w(ipoin))!-(eta(lpoin_w(ipoin),2)-eta(lpoin_w(ipoin),1))/dt
            end do
            !$acc end parallel loop
            call nvtxEndRange()

            call nvtxStartRange('rk_4_main: bc_fix_dirichlet_residual_entropy')
            if (noBoundaries .eqv. .false.) then
               call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,Reta)
            end if
            call nvtxEndRange()
            !
            ! If using Sutherland viscosity model:
            !
            call nvtxStartRange('rk_4_main: sutherland_viscosity')
            if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
               call sutherland_viscosity(npoin,Tem(:,pos),mu_factor,mu_fluid)
            end if
            call nvtxEndRange()
            !
            ! Compute entropy viscosity
            !
            call nvtxStartRange('rk_4_main: smart_visc_spectral')
            call smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,Reta,Ngp,coord,dNgp,gpvol,wgp, &
               gamma_gas,rho(:,pos),u(:,:,pos),csound,Tem(:,pos),eta(:,pos),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
            call nvtxEndRange()
            !
            ! Compute subgrid viscosity if active
            !
            if(flag_les == 1) then
               if(flag_les_ilsa == 1) then
                  call nvtxStartRange('rk_4_main: sgs_ilsa_visc')
                  call sgs_ilsa_visc(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dt,rho(:,pos),u(:,:,pos),mu_sgs,mu_fluid,mu_e,kres,etot,au,ax1,ax2,ax3,mue_l) 
                  call nvtxEndRange()
               else
                  call nvtxStartRange('rk_4_main: sgs_visc')
                  call sgs_visc(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,pos),u(:,:,pos),Ml,mu_sgs,mue_l)
                  call nvtxEndRange()
               end if
            end if

         end subroutine rk_4_main

         subroutine updateBuffer(npoin,npoin_w,coord,lpoin_w,rho,q,E,u_buffer)
            integer(4),           intent(in)    :: npoin
            integer(4),           intent(in)    :: npoin_w
            real(rp),             intent(in)    :: coord(npoin,ndime)
            integer(4),           intent(in)    :: lpoin_w(npoin_w)
            real(rp),             intent(inout) :: rho(npoin)
            real(rp),             intent(inout) :: E(npoin)
            real(rp),             intent(inout) :: q(npoin,ndime)
            real(rp),             intent(in) :: u_buffer(npoin,ndime)
            integer(4) :: ipoin
            real(rp)   :: xs,xb,xi,c1,c2, E_inf

            call nvtxStartRange('updateBuffer: CPU-COMPUTE')
            c1 = 0.01_rp
            c2 = 10.0_rp
            E_inf = (nscbc_rho_inf*0.5_rp*nscbc_u_inf**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0_rp))
            call nvtxEndRange()
            
            call nvtxStartRange('updateBuffer: GPU-COMPUTE')
            !$acc parallel loop
            do ipoin = 1,npoin_w
               xi = 1.0_rp
               !east
               if(flag_buffer_on_east .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),1)
                  if(xs>flag_buffer_e_min) then
                     xb = (xs-flag_buffer_e_min)/flag_buffer_e_size
                     xi = min((1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))),xi)
                  end if
               end if
               !west
               if(flag_buffer_on_west .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),1)
                  if(xs<flag_buffer_w_min) then
                     xb = (flag_buffer_w_min-xs)/flag_buffer_w_size
                     xi = min((1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))),xi)
                  end if
               end if
               !north
               if(flag_buffer_on_north .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),2)
                  if(xs>flag_buffer_n_min) then
                     xb = (xs-flag_buffer_n_min)/flag_buffer_n_size
                     xi = min((1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))),xi)
                  end if
               end if
               !south
               if(flag_buffer_on_south .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),2)
                  if(xs<flag_buffer_s_min) then
                     xb = (flag_buffer_s_min-xs)/flag_buffer_s_size
                     xi = min((1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))),xi)
                  end if
               end if
               !north
               if(flag_buffer_on_top .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),3)
                  if(xs>flag_buffer_t_min) then
                     xb = (xs-flag_buffer_t_min)/flag_buffer_t_size
                     xi = min((1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))),xi)
                  end if
               end if
               !bottom
               if(flag_buffer_on_bottom .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),3)
                  if(xs<flag_buffer_b_min) then
                     xb = (flag_buffer_b_min-xs)/flag_buffer_b_size
                     xi = min((1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))),xi)
                  end if
               end if

               q(lpoin_w(ipoin),1) = u_buffer(lpoin_w(ipoin),1)*rho(lpoin_w(ipoin)) + xi*(q(lpoin_w(ipoin),1)-u_buffer(lpoin_w(ipoin),1)*rho(lpoin_w(ipoin)))
               q(lpoin_w(ipoin),2) = u_buffer(lpoin_w(ipoin),2)*rho(lpoin_w(ipoin)) + xi*(q(lpoin_w(ipoin),2)-u_buffer(lpoin_w(ipoin),2)*rho(lpoin_w(ipoin)))
               q(lpoin_w(ipoin),3) = u_buffer(lpoin_w(ipoin),3)*rho(lpoin_w(ipoin)) + xi*(q(lpoin_w(ipoin),3)-u_buffer(lpoin_w(ipoin),3)*rho(lpoin_w(ipoin)))

               E(lpoin_w(ipoin)) = E_inf + xi*(E(lpoin_w(ipoin))-E_inf)
               rho(lpoin_w(ipoin)) = nscbc_rho_inf + xi*(rho(lpoin_w(ipoin))-nscbc_rho_inf)
            end do
            !$acc end parallel loop
            call nvtxEndRange()
         end subroutine updateBuffer

      end module time_integ
