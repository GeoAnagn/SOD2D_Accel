! DONT FORGET TO INPUT APPROPRIATE FILEPATH NAME IN <<path>> variable at line 28
! This is the static version for signle bin file use.

!/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/latest/ompi/bin/mpif90  -I/home/apps/libraries/hdf5/1.12.0/gnu/include -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel -gpu=cuda11.7,managed,lineinfo -cuda -acc -fast -c /home/apps/sod2d/src/lib_sod2d/sources/mod_constants.f90 -o ./mod_constants.f90.o
!/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/latest/ompi/bin/mpif90  -I/home/apps/libraries/hdf5/1.12.0/gnu/include -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel -gpu=cuda11.7,managed,lineinfo -cuda -acc -fast -c /home/apps/sod2d/src/lib_sod2d/sources/mod_nvtx.f90 -o ./mod_nvtx.f90.o
!/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/latest/ompi/bin/mpif90  -I/home/apps/libraries/hdf5/1.12.0/gnu/include -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel -gpu=cuda11.7,managed,lineinfo -cuda -acc -fast -c /home/apps/sod2d/src/lib_sod2d/sources/elem_convec.f90 -o CMakeFiles/sod2d.dir/__/lib_sod2d/sources/elem_convec.f90.o

program elem_convec_solo
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
        real(rp) :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
        real(rp) :: gpvol(1,ngaus,nelem)
        integer(4) :: atoIJK(nnode), invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
        real(rp) :: q(npoin,ndime), u(npoin,ndime), rho(npoin),pr(npoin), E(npoin)
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

        ! Write q variable to q.bin file.
        filepath = path // '/q.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) q ! check
        close(1)

        ! Write u variable to u.bin file.
        filepath = path // '/u.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) u ! check
        close(1)

        ! Write rho variable to rho.bin file.
        filepath = path // '/rho.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) rho ! check
        close(1)
        
        ! Write pr variable to pr.bin file.
        filepath = path // '/pr.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) pr ! check
        close(1)

        ! Write E variable to E.bin file.
        filepath = path // '/E.bin'
        call StripSpaces(filepath)
        open(1, file = filepath, form="unformatted", action='read', status='old')
        read(1) E ! check
        close(1)

        call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,pr,E,Rmass,Rmom,Rener)

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

    subroutine full_convec_ijk(nelem, npoin, connec, Ngp, dNgp, He, gpvol, dlxigp_ip, xgp, atoIJK, invAtoIJK, gmshAtoI, gmshAtoJ, gmshAtoK, u, q, rho, pr, E, Rmass, Rmom, Rener)

        implicit none

        integer(4), intent(in)  :: nelem, npoin
        integer(4), intent(in)  :: connec(nelem,nnode)
        real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
        real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
        real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
        integer(4), intent(in)  :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
        real(rp),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), rho(npoin),pr(npoin), E(npoin)
        real(rp),    intent(out) :: Rmass(npoin)
        real(rp),    intent(out) :: Rmom(npoin,ndime)
        real(rp),    intent(out) :: Rener(npoin)
        integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,kdime,ii
        real(rp)                 :: Re_mom(nnode,ndime)
        real(rp)                 :: Re_mass(nnode), Re_ener(nnode)
        real(rp)                 :: gradIsoRho(ndime),gradIsoP(ndime), gradIsoE(ndime),gradIsoU(ndime,ndime), gradIsoF(ndime,ndime,ndime), gradIsoQ(ndime,ndime), gradIsoFe(ndime,ndime)
        real(rp)                 :: gradRho(ndime),gradP(ndime),gradE(ndime),gradU(ndime,ndime),divF(ndime),divU,divFe,divQ
        real(rp)                 :: ul(nnode,ndime), ql(nnode,ndime), rhol(nnode), prl(nnode),El(nnode),fel(nnode,ndime),fl(nnode,ndime,ndime)
        real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

        !Declare Environment Variables for OpenTuner
        CHARACTER(100) :: gangchar
        CHARACTER(100) :: workerchar
        CHARACTER(100) :: vectorchar
        integer :: gang_num
        integer :: worker_num
        integer :: vector_num
        
        !Read Environment variables for OpenTuner
        CALL get_environment_variable("gang_num_fcijk", gangchar)
        CALL get_environment_variable("worker_num_fcijk", workerchar)
        CALL get_environment_variable("vector_num_fcijk", vectorchar)
        READ(gangchar,*)gang_num
        READ(workerchar,*)worker_num
        READ(vectorchar,*)vector_num

        call nvtxStartRange("Full convection")
        !$acc kernels
        Rmom(:,:) = 0.0_rp
        Rmass(:) = 0.0_rp
        Rener(:) = 0.0_rp
        !$acc end kernels

        !$acc parallel loop gang private(Re_ener,Re_mass,Re_mom,ul,ql,rhol,prl,El,fl,fel) num_gangs(gang_num) vector_length(vector_num)
        do ielem = 1,nelem
            !$acc loop vector collapse(2)
            do idime = 1,ndime
                do inode = 1,nnode
                    ul(inode,idime) = u(connec(ielem,inode),idime)
                    ql(inode,idime) = q(connec(ielem,inode),idime)
                    fel(inode,idime) = (E(connec(ielem,inode))+pr(connec(ielem,inode)))*u(connec(ielem,inode),idime)
                end do
            end do
            !$acc loop vector collapse(3)
            do idime = 1,ndime
                do jdime = 1,ndime
                    do inode = 1,nnode
                        fl(inode,idime,jdime)  = q(connec(ielem,inode),idime)*u(connec(ielem,inode),jdime)
                    end do
                end do
            end do
            !$acc loop vector
            do inode = 1,nnode
                rhol(inode) = rho(connec(ielem,inode))
                El(inode) = E(connec(ielem,inode))
                prl(inode) = pr(connec(ielem,inode))
            end do
            !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip, gradIsoRho,gradIsoP, gradIsoE,gradIsoU, gradIsoF, gradIsoQ, gradIsoFe,gradRho,gradP,gradE,gradU,divF,divU,divQ,divFe)
            do igaus = 1,ngaus
                !$acc loop seq
                do ii=1,porder+1
                    dlxi_ip(ii) = dlxigp_ip(igaus,1,ii)
                    dleta_ip(ii) = dlxigp_ip(igaus,2,ii)
                    dlzeta_ip(ii) = dlxigp_ip(igaus,3,ii)
                end do
                isoI = gmshAtoI(igaus)
                isoJ = gmshAtoJ(igaus)
                isoK = gmshAtoK(igaus)

                gradIsoRho(:) = 0.0_rp
                gradIsoP(:) = 0.0_rp
                gradIsoE(:) = 0.0_rp
                gradIsoU(:,:) = 0.0_rp
                gradIsoF(:,:,:) = 0._rp
                gradIsoQ(:,:) = 0._rp
                gradIsoFe(:,:) = 0._rp
                !$acc loop seq
                do ii=1,porder+1
                    gradIsoRho(1) = gradIsoRho(1) + dlxi_ip(ii)*rhol(invAtoIJK(ii,isoJ,isoK))
                    gradIsoRho(2) = gradIsoRho(2) + dleta_ip(ii)*rhol(invAtoIJK(isoI,ii,isoK))
                    gradIsoRho(3) = gradIsoRho(3) + dlzeta_ip(ii)*rhol(invAtoIJK(isoI,isoJ,ii))

                    gradIsoP(1) = gradIsoP(1) + dlxi_ip(ii)*prl(invAtoIJK(ii,isoJ,isoK))
                    gradIsoP(2) = gradIsoP(2) + dleta_ip(ii)*prl(invAtoIJK(isoI,ii,isoK))
                    gradIsoP(3) = gradIsoP(3) + dlzeta_ip(ii)*prl(invAtoIJK(isoI,isoJ,ii))

                    gradIsoE(1) = gradIsoE(1) + dlxi_ip(ii)*(prl(invAtoIJK(ii,isoJ,isoK))  + El(invAtoIJK(ii,isoJ,isoK)))
                    gradIsoE(2) = gradIsoE(2) + dleta_ip(ii)*(prl(invAtoIJK(isoI,ii,isoK)) + El(invAtoIJK(isoI,ii,isoK)))
                    gradIsoE(3) = gradIsoE(3) + dlzeta_ip(ii)*(prl(invAtoIJK(isoI,isoJ,ii))+ El(invAtoIJK(isoI,isoJ,ii)))

                    !$acc loop seq
                    do idime=1,ndime
                        gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoQ(idime,1) = gradIsoQ(idime,1) + dlxi_ip(ii)*ql(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoQ(idime,2) = gradIsoQ(idime,2) + dleta_ip(ii)*ql(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoQ(idime,3) = gradIsoQ(idime,3) + dlzeta_ip(ii)*ql(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoFe(idime,1) = gradIsoFe(idime,1) + dlxi_ip(ii)*fel(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFe(idime,2) = gradIsoFe(idime,2) + dleta_ip(ii)*fel(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFe(idime,3) = gradIsoFe(idime,3) + dlzeta_ip(ii)*fel(invAtoIJK(isoI,isoJ,ii),idime)

                        !$acc loop seq
                        do jdime=1, ndime
                            gradIsoF(idime,jdime,1) = gradIsoF(idime,jdime,1) + dlxi_ip(ii)*fl(invAtoIJK(ii,isoJ,isoK),idime,jdime)
                            gradIsoF(idime,jdime,2) = gradIsoF(idime,jdime,2) + dleta_ip(ii)*fl(invAtoIJK(isoI,ii,isoK),idime,jdime)
                            gradIsoF(idime,jdime,3) = gradIsoF(idime,jdime,3) + dlzeta_ip(ii)*fl(invAtoIJK(isoI,isoJ,ii),idime,jdime)
                        end do
                    end do
                end do

                gradRho(:) = 0.0_rp
                gradP(:) = 0.0_rp
                gradE(:) = 0.0_rp
                gradU(:,:) = 0.0_rp
                divF(:) = 0.0_rp
                divQ = 0.0_rp
                divFe = 0.0_rp
                !$acc loop seq
                do idime=1, ndime
                    !$acc loop seq
                    do jdime=1, ndime
                        gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                        gradP(idime)   = gradP(idime) + He(idime,jdime,igaus,ielem) * gradIsoP(jdime)
                        gradE(idime)   = gradE(idime) + He(idime,jdime,igaus,ielem) * gradIsoE(jdime)
                        divQ = divQ + He(idime,jdime,igaus,ielem) * gradIsoQ(idime,jdime)
                        divFe = divFe + He(idime,jdime,igaus,ielem) * gradIsoFe(idime,jdime)
                        !$acc loop seq
                        do kdime=1,ndime
                            gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                            divF(idime) = divF(idime) + He(jdime,kdime,igaus,ielem)*gradIsoF(idime,jdime,kdime)
                        end do
                    end do
                end do
            
                divU  = gradU(1,1)  + gradU(2,2)  + gradU(3,3)
                Re_mass(igaus) = 0.5_rp*(divQ+rhol(igaus)*divU)
                Re_ener(igaus) = 0.5_rp*(divFe+(El(igaus)+prl(igaus))*divU)
                !$acc loop seq
                do idime=1, ndime
                    Re_mom(igaus,idime) = 0.5_rp*(divF(idime)+ql(igaus,idime)*divU)+ gradP(idime)
                    Re_mass(igaus) = Re_mass(igaus) + 0.5_rp*(ul(igaus,idime)*gradRho(idime))
                    Re_ener(igaus) = Re_ener(igaus) + 0.5_rp*(ul(igaus,idime)*gradE(idime))
                    !$acc loop seq
                    do jdime=1, ndime
                        Re_mom(igaus,idime) = Re_mom(igaus,idime) + 0.5_rp*(ul(igaus,idime)*ul(igaus,jdime)*gradRho(jdime) &
                        + ql(igaus,jdime)*gradU(idime,jdime))
                    end do
                    Re_mom(igaus,idime) = gpvol(1,igaus,ielem)*Re_mom(igaus,idime)
                end do
                Re_mass(igaus) = gpvol(1,igaus,ielem)*Re_mass(igaus)
                Re_ener(igaus) = gpvol(1,igaus,ielem)*Re_ener(igaus)
            end do

            !
            ! Final assembly
            !
            !$acc loop vector collapse(2)
            do idime = 1,ndime
                do inode = 1,nnode
                    !$acc atomic update
                    Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re_mom(inode,idime)
                    !$acc end atomic
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
        end do
        !$acc end parallel loop
        
        call nvtxEndRange

    end subroutine full_convec_ijk

end program elem_convec_solo