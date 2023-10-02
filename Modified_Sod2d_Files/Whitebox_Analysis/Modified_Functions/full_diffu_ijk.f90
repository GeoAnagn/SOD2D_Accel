! DONT FORGET TO INPUT APPROPRIATE FILEPATH NAME IN <<path>> variable at line 28
! This is the static version for signle bin file use.

program elem_diffu_solo
    use mod_nvtx
    use mod_constants

    implicit none

    ! elem_convec variables
    integer(4) :: nelem, npoin
    
    ! Read variables
    character(200) :: path
    character(300) :: filepath

    ! Create appropriate folder for call data.
    call get_environment_variable("folder_path", path)
    call StripSpaces(path)

    ! Write nelem variable to nelem.bin file.
    filepath = path // '/nelem.bin'
    call StripSpaces(filepath)
    open(1, file = filepath, form="unformatted", action='read')
    read(1) nelem ! check
    close(1)
    
    ! Write npoin variable to npoin.bin file.
    filepath = path // '/npoin.bin'
    call StripSpaces(filepath)
    open(1, file = filepath, form="unformatted", action='read')
    read(1) npoin ! check
    close(1)
    
    call read_arrays(nelem, npoin, path)

    call exit(0)

contains

    subroutine StripSpaces(string)
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

    subroutine read_arrays(nelem, npoin, path)
        ! Function Variables
        integer(4), intent(in) :: nelem, npoin
        integer(4) :: connec(nelem,nnode)
        real(rp) :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
        real(rp) :: He(ndime,ndime,ngaus,nelem), xgp(ngaus,ndime), dlxigp_ip(ngaus,ndime,porder+1)
        real(rp) :: gpvol(1,ngaus,nelem)
        integer(4) :: atoIJK(nnode), invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
        real(rp)  :: Cp, Pr, rho(npoin), u(npoin,ndime), Tem(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus), Ml(npoin)
        real(rp)  :: mu_fluid(npoin)
        real(rp) :: Rmass(npoin)
        real(rp) :: Rmom(npoin,ndime)
        real(rp) :: Rener(npoin)

        ! Folder and File path variables.
        character(200), intent(in) :: path
        character(300) :: filepath  
    
        ! Write connec variable to connec.bin file.
        filepath = path // '/connec.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) connec
        close(1)

        ! Write Ngp variable to Ngp.bin file.
        filepath = path // '/Ngp.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) Ngp! check
        close(1)

        ! Write dNgp variable to dNgp.bin file.
        filepath = path // '/dNgp.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) dNgp ! check
        close(1)

        ! Write He variable to He.bin file.
        filepath = path // '/He.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) He ! check
        close(1)

        ! Write xgp variable to xgp.bin file.
        filepath = path // '/xgp.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) xgp ! check
        close(1)

        ! Write dlxigp_ip variable to dlxigp_ip.bin file.
        filepath = path // '/dlxigp_ip.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) dlxigp_ip ! check
        close(1)

        ! Write gpvol variable to gpvol.bin file.
        filepath = path // '/gpvol.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) gpvol ! check
        close(1)

        ! Write atoIJK variable to atoIJK.bin file.
        filepath = path // '/atoIJK.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) atoIJK ! check
        close(1)

        ! Write invAtoIJK variable to invAtoIJK.bin file.
        filepath = path // '/invAtoIJK.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) invAtoIJK! check
        close(1)

        ! Write gmshAtoI variable to gmshAtoI.bin file.
        filepath = path // '/gmshAtoI.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) gmshAtoI ! check
        close(1)

        ! Write gmshAtoJ variable to gmshAtoJ.bin file.
        filepath = path // '/gmshAtoJ.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) gmshAtoJ ! check
        close(1)

        ! Write gmshAtoK variable to gmshAtoK.bin file.
        filepath = path // '/gmshAtoK.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) gmshAtoK ! check
        close(1)

        ! Write Cp variable to Cp.bin file.
        filepath = path // '/Cp.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) Cp ! check
        close(1)

        ! Write Pr variable to Pr.bin file.
        filepath = path // '/Pr.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) Pr ! check
        close(1)

        ! Write rho variable to rho.bin file.
        filepath = path // '/rho.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) rho ! check
        close(1)

        ! Write u variable to u.bin file.
        filepath = path // '/u.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) u ! check
        close(1)

        ! Write Tem variable to Tem.bin file.
        filepath = path // '/Tem.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) Tem ! check
        close(1)

        ! Write mu_e variable to mu_e.bin file.
        filepath = path // '/mu_e.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) mu_e ! check
        close(1)

        ! Write mu_sgs variable to mu_sgs.bin file.
        filepath = path // '/mu_sgs.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) mu_sgs ! check
        close(1)

        ! Write Ml variable to Ml.bin file.
        filepath = path // '/Ml.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) Ml ! check
        close(1)

        ! Write mu_fluid variable to mu_fluid.bin file.
        filepath = path // '/mu_fluid.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) mu_fluid ! check
        close(1)

        call full_diffusion_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Pr,rho,u,Tem,mu_fluid,mu_e,mu_sgs,Ml,Rmass,Rmom,Rener)

        open(1, file = "./Rmass.bin", form="unformatted")
        write(1) Rmass
        close(1)

        open(1, file = "./Rmom.bin", form="unformatted")
        write(1) Rmom
        close(1)

        open(1, file = "./Rener.bin", form="unformatted")
        write(1) Rener
        close(1)

    end subroutine read_arrays

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

        !Declare Environment Variables for OpenTuner
        CHARACTER(100) :: gangchar
        CHARACTER(100) :: workerchar
        CHARACTER(100) :: vectorchar
        integer :: gang_num
        integer :: worker_num
        integer :: vector_num
        
        !Read Environment variables for OpenTuner
        CALL get_environment_variable("gang_num_fdijk", gangchar)
        CALL get_environment_variable("worker_num_fdijk", workerchar)
        CALL get_environment_variable("vector_num_fdijk", vectorchar)
        READ(gangchar,*)gang_num
        READ(workerchar,*)worker_num
        READ(vectorchar,*)vector_num

        call nvtxStartRange("Full diffusion")
        twoThirds = 2.0_rp/3.0_rp
        !$acc kernels
        Rmass(:) = 0.0_rp
        Rmom(:,:) = 0.0_rp
        Rener(:) = 0.0_rp
        !$acc end kernels

        !$acc parallel loop gang private(ul,Teml,rhol,mufluidl,gradRhol,gradTl,tauUl,tauXl,tauYl,tauZl) num_gangs(gang_num) vector_length(vector_num)
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

    end subroutine full_diffusion_ijk

end program elem_diffu_solo