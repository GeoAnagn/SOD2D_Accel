module elem_diffu

   use mod_nvtx
   use mod_constants
   use mod_veclen
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use mod_comms

   ! TODO: Create unit tests for all subroutines

contains

   subroutine mass_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,rho,mu_e,Rmass)

      ! TODO: Add stab. viscosity

      implicit none

      integer(4), intent(in)  :: nelem, npoin
      integer(4), intent(in)  :: connec(nelem,nnode)
      real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
      real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
      real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
      real(rp),    intent(in)  :: rho(npoin), mu_e(nelem,ngaus)
      real(rp),    intent(out) :: Rmass(npoin)
      integer(4)              :: ielem, igaus, inode, idime, jdime
      real(rp)                 :: Re(nnode), nu_e
      real(rp)                 :: tmp1, gpcar(ndime,nnode), tmp2

      call nvtxStartRange("Mass diffusion")
      !$acc kernels
      Rmass(:) = 0.0_rp
      !$acc end kernels
      !$acc parallel loop gang  private(gpcar,Re) vector_length(vecLength)
      do ielem = 1,nelem
         !$acc loop vector
         do inode = 1,nnode
            Re(inode) = 0.0_rp
         end do
         !$acc loop seq
         do igaus = 1,ngaus
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop vector
               do inode = 1,nnode
                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
               end do
            end do
         !$acc loop seq
         do idime = 1,ndime
            tmp1 = 0.0_rp
            tmp2 = 0.0_rp
            !$acc loop vector reduction(+:tmp1,tmp2)
            do inode = 1,nnode
               tmp1 = tmp1+gpcar(idime,inode)*rho(connec(ielem,inode))
               tmp2 = tmp2+Ngp(igaus,inode)*rho(connec(ielem,inode)) 
            end do
            nu_e = c_rho*mu_e(ielem,igaus)/tmp2
            !$acc loop vector
            do inode = 1,nnode
               Re(inode) = Re(inode)+nu_e*gpvol(1,igaus,ielem)* &
                           gpcar(idime,inode)*tmp1
            end do
         end do
         end do
         !$acc loop vector
         do inode = 1,nnode
            !$acc atomic update
            Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode))+Re(inode)
            !$acc end atomic
         end do
      end do
      !$acc end parallel loop
      call nvtxEndRange

   end subroutine mass_diffusion

   subroutine mom_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u,mu_fluid,mu_e,mu_sgs,Rmom)

      ! TODO: Add. stab. viscosity

      implicit none

      integer(4), intent(in)  :: nelem, npoin
      integer(4), intent(in)  :: connec(nelem,nnode)
      real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
      real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
      real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
      real(rp),    intent(in)  :: u(npoin,ndime),mu_e(nelem,ngaus),mu_sgs(nelem,ngaus)
      real(rp),    intent(in)  :: mu_fluid(npoin)
      real(rp),    intent(out) :: Rmom(npoin,ndime)
      integer(4)              :: ielem, igaus, idime, jdime, inode, jnode
      real(rp)                 :: Re(nnode,ndime), twoThirds, gpcar(ndime,nnode), tau(ndime,ndime)
      real(rp)                 :: divU, gradU(ndime,ndime), mu_fgp, aux

      twoThirds = 2.0_rp/3.0_rp
      call nvtxStartRange("Momentum diffusion")
      !$acc kernels
      Rmom(:,:) = 0.0_rp
      !$acc end kernels
      !$acc parallel loop gang  private(Re,gpcar,tau,gradU) vector_length(vecLength)
      do ielem = 1,nelem
         !$acc loop vector collapse(2)
         do inode = 1,nnode
            do idime = 1,ndime
               Re(inode,idime) = 0.0_rp
            end do
         end do
         !
         ! Gradient structure:
         !
         !         | u1,1 u1,2 u1,3 |
         ! u_i,j = | u2,1 u2,2 u2,3 |
         !         | u3,1 u3,2 u3,3 |
         !
         !$acc loop seq
         do igaus = 1,ngaus
            !
            ! Compute gpcar
            !
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop vector
               do inode = 1,nnode
                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
               end do
            end do
            !
            ! Compute combined viscosity at Gauss point
            !
            mu_fgp = 0.0_rp
            !$acc loop vector reduction(+:mu_fgp)
            do inode = 1,nnode
               mu_fgp = mu_fgp+(Ngp(igaus,inode)*mu_fluid(connec(ielem,inode)))
            end do
            mu_fgp = mu_fgp+mu_e(ielem,igaus)+mu_sgs(ielem,igaus)
            !
            ! Compute grad(u)
            !
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  aux = 0.0_rp
                  !$acc loop vector reduction(+:aux)
                  do inode = 1,nnode
                     aux = aux+gpcar(jdime,inode)*u(connec(ielem,inode),idime)
                  end do
                  gradU(idime,jdime) = aux
               end do
            end do
            !
            ! Compute div(u)
            !
            divU = 0.0_rp
            ! divU = gradU(1,1)+gradU(2,2)+gradU(3,3)
            !$acc loop seq
            do idime = 1,ndime
               divU = divU+gradU(idime,idime)
            end do

            ! TODO: Remove this
            !divU = 0.0_rp
            !!$acc loop vector collapse(2) reduction(+:divU)
            !do idime = 1,ndime
            !   do inode = 1,nnode
            !      divU = divU+gpcar(idime,inode)*u(connec(ielem,inode),idime)
            !   end do
            !end do

            !
            ! Finish computing tau_ij = u_i,j + u_j,i - (2/3)*div(u)*d_ij
            !
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  tau(idime,jdime) = gradU(idime,jdime)+gradU(jdime,idime)
               end do
            end do
            !$acc loop seq
            do idime = 1,ndime
               tau(idime,idime) = tau(idime,idime)-twoThirds*divU
            end do
            !
            ! Compute div(tau) at the Gauss point
            !
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  !$acc loop vector
                  do inode = 1,nnode
                     Re(inode,idime) = Re(inode,idime)+gpvol(1,igaus,ielem)* &
                        gpcar(jdime,inode)*mu_fgp*tau(idime,jdime)
                  end do
               end do
            end do
         end do
         !
         ! Assembly
         !
         !$acc loop vector collapse(2)
         do inode = 1,nnode
            do idime = 1,ndime
               !$acc atomic update
               Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+1.0_rp*Re(inode,idime)
               !$acc end atomic
            end do
         end do
      end do
      !$acc end parallel loop
      call nvtxEndRange

   end subroutine mom_diffusion

   subroutine ener_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u,Tem,mu_fluid,mu_e,mu_sgs,Rener)

      implicit none

      integer(4), intent(in)  :: nelem, npoin
      integer(4), intent(in)  :: connec(nelem,nnode)
      real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
      real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
      real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
      real(rp),    intent(in)  :: u(npoin,ndime), Tem(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus)
      real(rp),    intent(in)  :: mu_fluid(npoin)
      real(rp),    intent(out) :: Rener(npoin)
      integer(4)              :: ielem, igaus, inode, idime, jdime
      real(rp)                 :: Re(nnode), kappa_e, mu_fgp, aux, divU, tauU(ndime), twoThirds
      real(rp)                 :: gpcar(ndime,nnode), gradU(ndime,ndime), gradT(ndime)

      call nvtxStartRange("Energy diffusion")
      twoThirds = 2.0_rp/3.0_rp
      !$acc kernels
      Rener(:) = 0.0_rp
      !$acc end kernels
      !$acc parallel loop gang  private(Re,gpcar,gradU,gradT,tauU)  vector_length(vecLength)
      do ielem = 1,nelem
         !$acc loop vector
         do inode = 1,nnode
            Re(inode) = 0.0_rp
         end do
         !$acc loop seq
         do igaus = 1,ngaus
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop vector
               do inode = 1,nnode
                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
               end do
            end do
            !
            ! Compute viscosity and conductivity at Gauss point
            !
            mu_fgp = 0.0_rp
            !$acc loop vector reduction(+:mu_fgp)
            do inode = 1,nnode
               mu_fgp = mu_fgp+(Ngp(igaus,inode)*mu_fluid(connec(ielem,inode)))
            end do
            kappa_e =mu_fgp*1004.0_rp/0.71_rp+c_ener*mu_e(ielem,igaus)/0.4_rp + mu_sgs(ielem,igaus)/0.9_rp
            mu_fgp = mu_fgp+mu_e(ielem,igaus)
            !
            ! Compute grad(U)
            !
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  aux = 0.0_rp
                  !$acc loop vector reduction(+:aux)
                  do inode = 1,nnode
                     aux = aux+gpcar(jdime,inode)*u(connec(ielem,inode),idime)
                  end do
                  gradU(idime,jdime) = aux
               end do
            end do
            !
            ! Compute div(u)
            divU = 0.0_rp
            !$acc loop seq
            do idime = 1,ndime
               divU = divU+gradU(idime,idime)
            end do

            ! TODO: Remove this lines
            !divU = 0.0_rp
            !!$acc loop vector collapse(2) reduction(+:divU)
            !do idime = 1,ndime
            !   do inode = 1,nnode
            !      divU = divU+gpcar(idime,inode)*u(connec(ielem,inode),idime)
            !   end do
            !end do

            !
            ! Compute tau*u
            !
            !$acc loop seq
            do idime = 1,ndime
               tauU(idime) = 0.0_rp
               !$acc loop seq
               do jdime = 1,ndime
                  aux = 0.0_rp
                  !$acc loop vector reduction(+:aux)
                  do inode = 1,nnode
                     aux = aux+Ngp(igaus,inode)*u(connec(ielem,inode),jdime)
                  end do
                  tauU(idime) = tauU(idime) + &
                     gradU(idime,jdime)*aux + gradU(jdime,idime)*aux
               end do
               aux = 0.0_rp
               !$acc loop vector reduction(+:aux)
               do inode = 1,nnode
                  aux = aux+Ngp(igaus,inode)*u(connec(ielem,inode),idime)
               end do
               tauU(idime) = tauU(idime)-twoThirds*divU*aux
            end do
            !
            ! Compute grad(T)
            !
            !$acc loop seq
            do idime = 1,ndime
               aux = 0.0_rp
               !$acc loop vector reduction(+:aux)
               do inode = 1,nnode
                  aux = aux+gpcar(idime,inode)*Tem(connec(ielem,inode))
               end do
               gradT(idime) = aux
            end do
            !
            ! Gaussian integ
            !
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop vector
               do inode = 1,nnode
                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*gpcar(idime,inode)* &
                     (mu_fgp*tauU(idime)+kappa_e*gradT(idime))
               end do
            end do
         end do
         !$acc loop vector
         do inode = 1,nnode
            !$acc atomic update
            Rener(connec(ielem,inode)) = Rener(connec(ielem,inode))+1.0_rp*Re(inode)
            !$acc end atomic
         end do
      end do
      !$acc end parallel loop
      call nvtxEndRange

   end subroutine ener_diffusion

   subroutine full_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,Cp,Pr,rho,u,Tem,mu_fluid,mu_e,mu_sgs,Ml,Rmass,Rmom,Rener)
      implicit none

      integer(4), intent(in)  :: nelem, npoin
      integer(4), intent(in)  :: connec(nelem,nnode)
      real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
      real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
      real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
      real(rp),    intent(in)  :: Cp,Pr,rho(npoin),u(npoin,ndime), Tem(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus),Ml(npoin)
      real(rp),    intent(in)  :: mu_fluid(npoin)
      real(rp),    intent(out) :: Rmass(npoin)
      real(rp),    intent(out) :: Rmom(npoin,ndime)
      real(rp),    intent(out) :: Rener(npoin)
      integer(4)              :: ielem, igaus, inode, idime, jdime
      real(rp)                 :: Re_mass(nnode),Re_mom(nnode,ndime),Re_ener(nnode)
      real(rp)                 :: kappa_e, mu_fgp, mu_egp,aux, divU, tauU(ndime), twoThirds,nu_e,tau(ndime,ndime)
      real(rp)                 :: gpcar(ndime,nnode), gradU(ndime,ndime), gradT(ndime),tmp1,ugp(ndime),vol,arho
      real(rp)                 :: ul(nnode,ndime), rhol(nnode),Teml(nnode),mufluidl(nnode), aux2

      call nvtxStartRange("Full diffusion")
      twoThirds = 2.0_rp/3.0_rp
      !$acc kernels
      Rmass(:) = 0.0_rp
      Rmom(:,:) = 0.0_rp
      Rener(:) = 0.0_rp
      !$acc end kernels

      !$acc parallel loop gang  private(ul,Teml,rhol,mufluidl,Re_mass,Re_mom,Re_ener) !!vector_length(vecLength)
      do ielem = 1,nelem
         !$acc loop vector
         do inode = 1,nnode
            Re_mass(inode) = 0.0_rp
            Re_ener(inode) = 0.0_rp
            rhol(inode) = rho(connec(ielem,inode))
            Teml(inode) = Tem(connec(ielem,inode))
            mufluidl(inode) = mu_fluid(connec(ielem,inode))
         end do
         !$acc loop vector collapse(2)
         do inode = 1,nnode
            do idime = 1,ndime
               Re_mom(inode,idime) = 0.0_rp
               ul(inode,idime) = u(connec(ielem,inode),idime)
            end do
         end do
         !$acc loop worker private(gpcar,tau,gradU,gradT,tauU,ugp)
         do igaus = 1,ngaus
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop vector
               do inode = 1,nnode
                  aux =  dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                  gpcar(idime,inode) = aux
               end do
               ugp(idime) = ul(igaus,idime)
            end do
            nu_e = c_rho*mu_e(ielem,igaus)/rhol(igaus)
            mu_fgp = mufluidl(igaus)+rhol(igaus)*mu_sgs(ielem,igaus)
            mu_egp = mu_e(ielem,igaus)
            kappa_e =mufluidl(igaus)*Cp/Pr+c_ener*mu_e(ielem,igaus)/0.4_rp + rhol(igaus)*mu_sgs(ielem,igaus)/0.9_rp
            !
            ! Compute grad(u)
            !
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  aux = 0.0_rp
                  !$acc loop vector reduction(+:aux)
                  do inode = 1,nnode
                     aux = aux+gpcar(jdime,inode)*ul(inode,idime)
                  end do
                  gradU(idime,jdime) = aux
               end do
            end do
            !
            ! Compute div(u)
            !
            divU = gradU(1,1)+gradU(2,2)+gradU(3,3)
            !
            ! Finish computing tau_ij = u_i,j + u_j,i - (2/3)*div(u)*d_ij
            ! Compute tau*u
            !
            !$acc loop seq
            do idime = 1,ndime
               tauU(idime) = 0.0_rp
               !$acc loop seq
               do jdime = 1,ndime
                  tauU(idime) = tauU(idime) + &
                     (mu_fgp+mu_egp)*(gradU(idime,jdime)+ gradU(jdime,idime))*ugp(jdime)
                  tau(idime,jdime) = (mu_fgp+mu_egp)*(gradU(idime,jdime)+gradU(jdime,idime))
               end do
               tauU(idime) = tauU(idime)-(mu_fgp)*twoThirds*divU*ugp(idime)
               tau(idime,idime) = tau(idime,idime)-(mu_fgp)*twoThirds*divU
            end do

            ! Dif rho
            ! Dif T
            ! Compute div(tau) at the Gauss point
            !$acc loop seq
            do idime = 1,ndime
               tmp1 = 0.0_rp
               aux = 0.0_rp
               !$acc loop vector reduction(+:tmp1,aux)
               do inode = 1,nnode
                  tmp1 = tmp1+gpcar(idime,inode)*rhol(inode)
                  aux = aux+gpcar(idime,inode)*Teml(inode)
               end do
               !$acc loop vector
               do inode = 1,nnode
                  !$acc atomic update
                  Re_mass(inode) = Re_mass(inode)+nu_e*gpvol(1,igaus,ielem)* &
                     gpcar(idime,inode)*tmp1
                  !$acc end atomic
                  !$acc atomic update
                  Re_ener(inode) = Re_ener(inode)+gpvol(1,igaus,ielem)*gpcar(idime,inode)* &
                     (tauU(idime)+kappa_e*aux) 
                  !$acc end atomic
               end do
               !$acc loop vector collapse(2)
               do jdime = 1,ndime
                  do inode = 1,nnode
                     !$acc atomic update
                     Re_mom(inode,idime) = Re_mom(inode,idime)+gpvol(1,igaus,ielem)* &
                        gpcar(jdime,inode)*tau(idime,jdime)
                     !$acc end atomic
                  end do
               end do
            end do
         end do
         !$acc loop vector
         do inode = 1,nnode
            !$acc atomic update
            Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode))+Re_mass(inode)
            !$acc end atomic
            !$acc atomic update
            Rener(connec(ielem,inode)) = Rener(connec(ielem,inode))+Re_ener(inode)
            !$acc end atomic
         end do
         !$acc loop vector collapse(2)
         do inode = 1,nnode
            do idime = 1,ndime
               !$acc atomic update
               Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re_mom(inode,idime)
               !$acc end atomic
            end do
         end do
      end do
      !$acc end parallel loop
      call nvtxEndRange

   end subroutine full_diffusion
   
   subroutine StripSpaces_2(string)
      character(len=*) :: string
      integer :: stringLen
      integer :: last, actual

      stringLen = len (string)
      last = 1
      actual = 1

      do while (actual < stringLen)
         if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
         else
            last = last + 1
            if (actual < last) &
               actual = last
         endif
      end do

   end subroutine

   subroutine full_diffusion_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Pr,rho,u,Tem,mu_fluid,mu_e,mu_sgs,Ml,Rmass,Rmom,Rener)
      implicit none

      integer(4), intent(in)  :: nelem, npoin
      integer(4), intent(in)  :: connec(nelem,nnode)
      real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
      real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
      real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
      integer(4), intent(in)  :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),    intent(in)  :: Cp,Pr,rho(npoin),u(npoin,ndime), Tem(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus),Ml(npoin)
      real(rp),    intent(in)  :: mu_fluid(npoin)
      real(rp),    intent(out) :: Rmass(npoin)
      real(rp),    intent(out) :: Rmom(npoin,ndime)
      real(rp),    intent(out) :: Rener(npoin)
      integer(4)              :: ielem, igaus, inode, idime, jdime, isoI, isoJ, isoK,kdime,ii
      real(rp)                 :: kappa_e, mu_fgp, mu_egp,divU, tauU(ndime), twoThirds,nu_e,tau(ndime,ndime)
      real(rp)                 :: gradU(ndime,ndime), gradT(ndime),tmp1,vol,arho
      real(rp)                 :: gradIsoRho(ndime),gradIsoT(ndime),gradIsoU(ndime,ndime)
      real(rp)                 :: gradRho(ndime),divDm(ndime),divDr,divDe
      real(rp)                 :: ul(nnode,ndime), rhol(nnode),Teml(nnode),mufluidl(nnode)
      real(rp)                 :: tauXl(nnode,ndime), tauYl(nnode,ndime), tauZl(nnode,ndime)
      real(rp)                 :: gradTl(nnode,ndime),gradRhol(nnode,ndime),tauUl(nnode,ndime)

      !---------------------------------------------------------------------------------------------------------------------------------------------------------

      character(100) :: file_num_char
      character(200) :: path
      character(300) :: filepath
      character(300) :: mkdir_command

      integer :: file_num
      integer :: info

      ! Create appropriate folder for call data.
      call get_environment_variable("file_num", file_num_char)

      ! Create appropriate folder for call data.
      call get_environment_variable("file_num", file_num_char)
      path = '/home/apps/Examples/Clean_code/archive/data' // '/data_' // file_num_char
      call StripSpaces_2(path)
      mkdir_command = "mkdir " // path
      call system(mkdir_command)

      path = '/home/apps/Examples/Clean_code/archive/data' // '/data_' // file_num_char // '/elem_diffu'
      call StripSpaces_2(path)
      mkdir_command = "mkdir " // path
      call system(mkdir_command)

      ! Write nelem variable to nelem.bin file.
      filepath = path // '/nelem.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) nelem ! check
      close(1)
      
      ! Write npoin variable to npoin.bin file.
      filepath = path // '/npoin.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) npoin ! check
      close(1)

      ! Write nnode variable to nnode.bin file.
      filepath = path // '/nnode.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) nnode ! check
      close(1)
      
      ! Write ngaus variable to ngaus.bin file.
      filepath = path // '/ngaus.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) ngaus ! check
      close(1)

      ! Write ndime variable to ndime.bin file.
      filepath = path // '/ndime.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) ndime ! check
      close(1)

      ! Write porder variable to porder.bin file.
      filepath = path // '/porder.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) porder ! check
      close(1)
      
      ! Write connec variable to connec.bin file.
      filepath = path // '/connec.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) connec
      close(1)

      ! Write Ngp variable to Ngp.bin file.
      filepath = path // '/Ngp.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) Ngp! check
      close(1)

      ! Write dNgp variable to dNgp.bin file.
      filepath = path // '/dNgp.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) dNgp ! check
      close(1)

      ! Write He variable to He.bin file.
      filepath = path // '/He.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) He ! check
      close(1)

      ! Write xgp variable to xgp.bin file.
      filepath = path // '/xgp.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) xgp ! check
      close(1)

      ! Write dlxigp_ip variable to dlxigp_ip.bin file.
      filepath = path // '/dlxigp_ip.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) dlxigp_ip ! check
      close(1)

      ! Write gpvol variable to gpvol.bin file.
      filepath = path // '/gpvol.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) gpvol ! check
      close(1)

      ! Write atoIJK variable to atoIJK.bin file.
      filepath = path // '/atoIJK.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) atoIJK ! check
      close(1)

      ! Write invAtoIJK variable to invAtoIJK.bin file.
      filepath = path // '/invAtoIJK.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) invAtoIJK! check
      close(1)

      ! Write gmshAtoI variable to gmshAtoI.bin file.
      filepath = path // '/gmshAtoI.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) gmshAtoI ! check
      close(1)

      ! Write gmshAtoJ variable to gmshAtoJ.bin file.
      filepath = path // '/gmshAtoJ.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) gmshAtoJ ! check
      close(1)

      ! Write gmshAtoK variable to gmshAtoK.bin file.
      filepath = path // '/gmshAtoK.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) gmshAtoK ! check
      close(1)

      ! Write Cp variable to Cp.bin file.
      filepath = path // '/Cp.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) Cp ! check
      close(1)

      ! Write Pr variable to Pr.bin file.
      filepath = path // '/Pr.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) Pr ! check
      close(1)

      ! Write rho variable to rho.bin file.
      filepath = path // '/rho.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) rho ! check
      close(1)

      ! Write u variable to u.bin file.
      filepath = path // '/u.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) u ! check
      close(1)

      ! Write Tem variable to Tem.bin file.
      filepath = path // '/Tem.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) Tem ! check
      close(1)

      ! Write mu_e variable to mu_e.bin file.
      filepath = path // '/mu_e.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) mu_e ! check
      close(1)

      ! Write mu_sgs variable to mu_sgs.bin file.
      filepath = path // '/mu_sgs.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) mu_sgs ! check
      close(1)

      ! Write Ml variable to Ml.bin file.
      filepath = path // '/Ml.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) Ml ! check
      close(1)

      ! Write mu_fluid variable to mu_fluid.bin file.
      filepath = path // '/mu_fluid.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) mu_fluid ! check
      close(1)

      !---------------------------------------------------------------------------------------------------------------------------------------------------------

      call nvtxStartRange("Full diffusion")
      twoThirds = 2.0_rp/3.0_rp
      !$acc kernels
      Rmass(:) = 0.0_rp
      Rmom(:,:) = 0.0_rp
      Rener(:) = 0.0_rp
      !$acc end kernels

      !$acc parallel loop gang private(ul,Teml,rhol,mufluidl,gradRhol,gradTl,tauUl,tauXl,tauYl,tauZl)
      do ielem = 1,nelem
         !$acc loop vector
         do inode = 1,nnode
            rhol(inode) = rho(connec(ielem,inode))
            Teml(inode) = Tem(connec(ielem,inode))
            mufluidl(inode) = mu_fluid(connec(ielem,inode))
         end do
         !$acc loop vector collapse(2)
         do inode = 1,nnode
            do idime = 1,ndime
               ul(inode,idime) = u(connec(ielem,inode),idime)
            end do
         end do
         tauXl(:,:) = 0.0_rp
         tauYl(:,:) = 0.0_rp
         tauZl(:,:) = 0.0_rp
         gradTl(:,:) = 0.0_rp
         gradRhol(:,:) = 0.0_rp
         tauUl(:,:) = 0.0_rp

         !$acc loop vector private(tau,gradU,gradT,tauU,gradIsoRho,gradIsoT,gradIsoU,gradRho,divU)
         do igaus = 1,ngaus
            nu_e = c_rho*mu_e(ielem,igaus)/rhol(igaus)
            mu_fgp = mufluidl(igaus)+rhol(igaus)*mu_sgs(ielem,igaus)
            mu_egp = mu_e(ielem,igaus)
            kappa_e =mufluidl(igaus)*Cp/Pr+c_ener*mu_e(ielem,igaus)/0.4_rp + rhol(igaus)*mu_sgs(ielem,igaus)/0.9_rp

            isoI = gmshAtoI(igaus) 
            isoJ = gmshAtoJ(igaus) 
            isoK = gmshAtoK(igaus) 

            gradIsoRho(:) = 0.0_rp
            gradIsoT(:) = 0.0_rp
            gradIsoU(:,:) = 0.0_rp
            !$acc loop seq
            do ii=1,porder+1
               gradIsoRho(1) = gradIsoRho(1) + dlxigp_ip(igaus,1,ii)*rhol(invAtoIJK(ii,isoJ,isoK))
               gradIsoRho(2) = gradIsoRho(2) + dlxigp_ip(igaus,2,ii)*rhol(invAtoIJK(isoI,ii,isoK))
               gradIsoRho(3) = gradIsoRho(3) + dlxigp_ip(igaus,3,ii)*rhol(invAtoIJK(isoI,isoJ,ii))

               gradIsoT(1) = gradIsoT(1) + dlxigp_ip(igaus,1,ii)*Teml(invAtoIJK(ii,isoJ,isoK))
               gradIsoT(2) = gradIsoT(2) + dlxigp_ip(igaus,2,ii)*Teml(invAtoIJK(isoI,ii,isoK))
               gradIsoT(3) = gradIsoT(3) + dlxigp_ip(igaus,3,ii)*Teml(invAtoIJK(isoI,isoJ,ii))

               !$acc loop seq
               do idime=1,ndime
                  gradIsoU(idime,1) = gradIsoU(idime,1) + dlxigp_ip(igaus,1,ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                  gradIsoU(idime,2) = gradIsoU(idime,2) + dlxigp_ip(igaus,2,ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                  gradIsoU(idime,3) = gradIsoU(idime,3) + dlxigp_ip(igaus,3,ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)
               end do
            end do

            gradRho(:) = 0.0_rp
            gradT(:) = 0.0_rp
            gradU(:,:) = 0.0_rp
            !$acc loop seq
            do idime=1, ndime
               !$acc loop seq
               do jdime=1, ndime
                  gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                  gradT(idime)   = gradT(idime)   + He(idime,jdime,igaus,ielem) * gradIsoT(jdime)
                  !$acc loop seq
                  do kdime=1,ndime
                     gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                  end do
               end do
            end do

            divU = gradU(1,1)+gradU(2,2)+gradU(3,3)

            tauU(:) = 0.0_rp
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  tauU(idime) = tauU(idime) + &
                     (mu_fgp+mu_egp)*(gradU(idime,jdime)+ gradU(jdime,idime))*ul(igaus,jdime)
                  tau(idime,jdime) = (mu_fgp+mu_egp)*(gradU(idime,jdime)+gradU(jdime,idime))
               end do
               tauU(idime) = tauU(idime)-(mu_fgp)*twoThirds*divU*ul(igaus,idime)
               tau(idime,idime) = tau(idime,idime)-(mu_fgp)*twoThirds*divU
            end do

            !$acc loop seq
            do idime = 1,ndime
               tauXl(igaus,idime) =  tau(1,idime)
               tauYl(igaus,idime) =  tau(2,idime)
               tauZl(igaus,idime) =  tau(3,idime)
               gradTl(igaus,idime) =  gradT(idime)
               gradRhol(igaus,idime) =  gradRho(idime)
               tauUl(igaus,idime) =  tauU(idime)
            end do
         end do

         !$acc loop vector private(divDm,divDr,divDe) 
         do igaus = 1,ngaus
            nu_e = c_rho*mu_e(ielem,igaus)/rhol(igaus)
            mu_fgp = mufluidl(igaus)+rhol(igaus)*mu_sgs(ielem,igaus)
            mu_egp = mu_e(ielem,igaus)
            kappa_e =mufluidl(igaus)*Cp/Pr+c_ener*mu_e(ielem,igaus)*Cp/Pr + Cp*rhol(igaus)*mu_sgs(ielem,igaus)/0.9_rp

            isoI = gmshAtoI(igaus) 
            isoJ = gmshAtoJ(igaus) 
            isoK = gmshAtoK(igaus) 

            divDe = 0.0_rp
            divDr = 0.0_rp
            divDm(:) = 0.0_rp
            
            !$acc loop seq
            do ii=1,porder+1
               !$acc loop seq
               do idime=1,ndime
                  divDr = divDr + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*gradRhol(invAtoIJK(ii,isoJ,isoK),idime)
                  divDr = divDr + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*gradRhol(invAtoIJK(isoI,ii,isoK),idime)
                  divDr = divDr + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*gradRhol(invAtoIJK(isoI,isoJ,ii),idime)

                  divDe = divDe + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*(tauUl(invAtoIJK(ii,isoJ,isoK),idime)+kappa_e*gradTl(invAtoIJK(ii,isoJ,isoK),idime))
                  divDe = divDe + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*(tauUl(invAtoIJK(isoI,ii,isoK),idime)+kappa_e*gradTl(invAtoIJK(isoI,ii,isoK),idime))
                  divDe = divDe + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*(tauUl(invAtoIJK(isoI,isoJ,ii),idime)+kappa_e*gradTl(invAtoIJK(isoI,isoJ,ii),idime))
                  
                  divDm(1) = divDm(1) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauXl(invAtoIJK(ii,isoJ,isoK),idime)
                  divDm(1) = divDm(1) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauXl(invAtoIJK(isoI,ii,isoK),idime)
                  divDm(1) = divDm(1) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauXl(invAtoIJK(isoI,isoJ,ii),idime)

                  divDm(2) = divDm(2) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauYl(invAtoIJK(ii,isoJ,isoK),idime)
                  divDm(2) = divDm(2) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauYl(invAtoIJK(isoI,ii,isoK),idime)
                  divDm(2) = divDm(2) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauYl(invAtoIJK(isoI,isoJ,ii),idime)

                  divDm(3) = divDm(3) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauZl(invAtoIJK(ii,isoJ,isoK),idime)
                  divDm(3) = divDm(3) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauZl(invAtoIJK(isoI,ii,isoK),idime)
                  divDm(3) = divDm(3) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauZl(invAtoIJK(isoI,isoJ,ii),idime)
               end do
            end do

            !$acc atomic update
            Rmass(connec(ielem,igaus)) = Rmass(connec(ielem,igaus))+nu_e*divDr
            !$acc end atomic
            !$acc atomic update
            Rener(connec(ielem,igaus)) = Rener(connec(ielem,igaus))+divDe
            !$acc end atomic
            do idime = 1,ndime
               !$acc atomic update
               Rmom(connec(ielem,igaus),idime) = Rmom(connec(ielem,igaus),idime)+divDm(idime)
               !$acc end atomic
            end do
         end do
      end do
      !$acc end parallel loop
      call nvtxEndRange
      
      !---------------------------------------------------------------------------------------------------------------------------------------------------------

      ! Write Rmass variable to Rmass.bin file.
      filepath = path // '/Rmass.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) Rmass ! check
      close(1)

      ! Write Rmom variable to Rmom.bin file.
      filepath = path // '/Rmom.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) Rmom ! check
      close(1)

      ! Write Rener variable to Rener.bin file.
      filepath = path // '/Rener.bin'
      call StripSpaces_2(filepath)
      open(1, file = filepath, form="unformatted")
      write(1) Rener ! check
      close(1)
      
      !---------------------------------------------------------------------------------------------------------------------------------------------------------

   end subroutine full_diffusion_ijk

end module elem_diffu
