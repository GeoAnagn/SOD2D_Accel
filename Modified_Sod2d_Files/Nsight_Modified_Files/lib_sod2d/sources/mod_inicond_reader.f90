module mod_inicond_reader

   use mod_numerical_params
   use mod_utils
   use mod_maths
   use mod_mpi
   use mod_comms
   use mod_nvtx
   contains

   subroutine order_matrix_globalIdSrl(numNodesInRank,globalIdsArray,matGidSrlOrdered)
      implicit none
      integer(4), intent(in)  :: numNodesInRank
      integer(8), intent(in)  :: globalIdsArray(numNodesInRank)
      integer(8), intent(out) :: matGidSrlOrdered(numNodesInRank,2)
      integer(4) :: iNodeL

      call nvtxStartRange("order_matrix_globalIdSrl: CPU-DATA")
      do iNodeL=1,numNodesInRank
         matGidSrlOrdered(iNodeL,1) = iNodeL
         matGidSrlOrdered(iNodeL,2) = globalIdsArray(iNodeL)
      end do
      call nvtxEndRange()

      call nvtxStartRange("order_matrix_globalIdSrl: quicksort_matrix_int8")
      call quicksort_matrix_int8(matGidSrlOrdered,2)
      call nvtxEndRange()

   end subroutine order_matrix_globalIdSrl

   subroutine avg_randomField_in_sharedNodes_Par(numNodesInRank,realField)
      implicit none
      integer, intent(in) :: numNodesInRank
      real(rp), intent(inout) :: realField(numNodesInRank)
      integer :: numRanksNodeCnt(numNodesInRank)
      integer :: i,iNodeL

      call nvtxStartRange("avg_randomField_in_sharedNodes_Par: CPU-DATA")
      numRanksNodeCnt(:)=1
      do i= 1,numNodesToComm
         iNodeL = nodesToComm(i)
         numRanksNodeCnt(iNodeL) = numRanksNodeCnt(iNodeL) + 1
      end do 
      call nvtxEndRange()

      call nvtxStartRange("avg_randomField_in_sharedNodes_Par: mpi_halo_atomic_update_real")
      call mpi_halo_atomic_update_real(realField)
      call nvtxEndRange()

      call nvtxStartRange("avg_randomField_in_sharedNodes_Par: CPU-COMPUTE")
      do iNodeL = 1,numNodesInRank
         realField(iNodeL) = realField(iNodeL) / real(numRanksNodeCnt(iNodeL),rp)
      end do
      call nvtxEndRange()
   end subroutine avg_randomField_in_sharedNodes_Par

   subroutine read_densi_from_file_Par(numElemsRank,numNodesRank,totalNumNodesInSerial,file_path,rho,connecPar,Ngp_l,matGidSrlOrdered)
      implicit none
      character(500), intent(in) :: file_path
      integer(4), intent(in)     :: numElemsRank,numNodesRank
      integer(8), intent(in)     :: totalNumNodesInSerial
      real(rp), intent(out)      :: rho(numNodesRank)
      integer(4), intent(in)     :: connecPar(numElemsRank,nnode)
      real(rp), intent(in)       :: Ngp_l(ngaus,nnode)
      integer(8), intent(in)     :: matGidSrlOrdered(numNodesRank,2)
      character(500)             :: file_type, file_name
      integer(4)                 :: iNodeL,iElem,igp,idime,auxCnt,readInd
      integer(8)                 :: iLine,iNodeGSrl
      real(8)                    :: readValue
      real(rp)                   :: aux_2(numNodesRank)

      call nvtxStartRange("read_densi_from_file_Par: CPU-DATA")
      write(file_type,*) ".alya"
      write(file_name,*) "DENSI"
      if(mpi_rank.eq.0) write(*,*) "--| Reading file DENSI.alya in parallel..."

      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
      
      auxCnt = 1
      do iLine = 1,totalNumNodesInSerial
         read(99,*) readInd, readValue
         if(iLine.eq.matGidSrlOrdered(auxCnt,2)) then
            iNodeL = matGidSrlOrdered(auxCnt,1)
            auxCnt=auxCnt+1
            if(rp.eq.4) then
               rho(iNodeL) = real(readValue,rp)
            else
               rho(iNodeL) = readValue
            end if
         end if         
      end do

      close(99)
      call nvtxEndRange()

      if(mpi_rank.eq.0) write(*,*) "--| Interpolating density from file coords to new mesh coords..."
      aux_2(:) = rho(:)
      do iElem = 1,numElemsRank
         do igp = (2**ndime)+1,nnode
            call nvtxStartRange("read_densi_from_file_Par: var_interpolate")
            call var_interpolate(nnode,aux_2(connecPar(iElem,:)),Ngp_l(igp,:),rho(connecPar(iElem,igp)))
            call nvtxEndRange()
         end do
      end do

      call nvtxStartRange("read_densi_from_file_Par: GPU-DATA")
      !$acc update device(rho(:))
      call nvtxEndRange()
   end subroutine read_densi_from_file_Par

   subroutine read_veloc_from_file_Par(numElemsRank,numNodesRank,totalNumNodesInSerial,file_path,u,connecPar,Ngp_l,matGidSrlOrdered)
      implicit none
      character(500), intent(in) :: file_path
      integer(4), intent(in)     :: numElemsRank,numNodesRank
      integer(8), intent(in)     :: totalNumNodesInSerial
      real(rp), intent(out)      :: u(numNodesRank,ndime)
      integer(4), intent(in)     :: connecPar(numElemsRank,nnode)
      real(rp), intent(in)       :: Ngp_l(ngaus,nnode)
      integer(8), intent(in)     :: matGidSrlOrdered(numNodesRank,2)
      character(500)             :: file_type, file_name
      integer(4)                 :: iNodeL,iElem,igp,idime,auxCnt,readInd
      integer(8)                 :: iLine,iNodeGSrl
      real(8)                    :: readValue_ux,readValue_uy,readValue_uz
      real(rp)                   :: aux_1(numNodesRank,ndime)

      call nvtxStartRange("read_veloc_from_file_Par: CPU-DATA")
      write(file_type,*) ".alya"
      write(file_name,*) "VELOC"
      if(mpi_rank.eq.0) write(*,*) "--| Reading file VELOC.alya in parallel..."

      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
      
      auxCnt = 1
      serialLoop: do iLine = 1,totalNumNodesInSerial
         read(99,*) readInd, readValue_ux, readValue_uy, readValue_uz
         if(iLine.eq.matGidSrlOrdered(auxCnt,2)) then
            iNodeL = matGidSrlOrdered(auxCnt,1)
            auxCnt=auxCnt+1
            if(rp.eq.4) then
               u(iNodeL,1)=real(readValue_ux,rp)
               u(iNodeL,2)=real(readValue_uy,rp)
               u(iNodeL,3)=real(readValue_uz,rp)
            else
               u(iNodeL,1)=readValue_ux
               u(iNodeL,2)=readValue_uy
               u(iNodeL,3)=readValue_uz
            end if
         end if         
         if(auxCnt.gt.numNodesRank) then
            exit serialLoop
         end if
      end do serialLoop
      
      close(99)
      call nvtxEndRange()

      if(mpi_rank.eq.0) write(*,*) "--| Interpolating velocity from file coords to new mesh coords..."
      aux_1(:,:) = u(:,:)
      do iElem = 1,numElemsRank
         do igp = (2**ndime)+1,nnode
            do idime = 1,ndime
               call nvtxStartRange("read_veloc_from_file_Par: var_interpolate")
               call var_interpolate(nnode,aux_1(connecPar(iElem,:),idime),Ngp_l(igp,:),u(connecPar(iElem,igp),idime))
               call nvtxEndRange()
            end do
         end do
      end do

      !$acc update device(u(:,:))

   end subroutine read_veloc_from_file_Par

   subroutine read_press_from_file_Par(numElemsRank,numNodesRank,totalNumNodesInSerial,file_path,pr,connecPar,Ngp_l,matGidSrlOrdered)
      implicit none
      character(500), intent(in) :: file_path
      integer(4), intent(in)     :: numElemsRank,numNodesRank
      integer(8), intent(in)     :: totalNumNodesInSerial
      real(rp), intent(out)      :: pr(numNodesRank)
      integer(4), intent(in)     :: connecPar(numElemsRank,nnode)
      real(rp), intent(in)       :: Ngp_l(ngaus,nnode)
      integer(8), intent(in)     :: matGidSrlOrdered(numNodesRank,2)
      character(500)             :: file_type, file_name
      integer(4)                 :: iNodeL,iElem,igp,idime,auxCnt,readInd
      integer(8)                 :: iLine,iNodeGSrl
      real(8)                    :: readValue
      real(rp)                   :: aux_2(numNodesRank)

      call nvtxStartRange("read_press_from_file_Par: CPU-DATA")
      write(file_type,*) ".alya"
      write(file_name,*) "PRESS"
      if(mpi_rank.eq.0) write(*,*) "--| Reading file PRESS.alya in parallel..."

      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
      
      auxCnt = 1
      do iLine = 1,totalNumNodesInSerial
         read(99,*) readInd, readValue
         if(iLine.eq.matGidSrlOrdered(auxCnt,2)) then
            iNodeL = matGidSrlOrdered(auxCnt,1)
            auxCnt=auxCnt+1
            if(rp.eq.4) then
               pr(iNodeL) = real(readValue,rp)
            else
               pr(iNodeL) = readValue
            end if
         end if         
      end do
      
      close(99)
      call nvtxEndRange()


      if(mpi_rank.eq.0) write(*,*) "--| Interpolating pressure from file coords to new mesh coords..."
      aux_2(:) = pr(:)
      do iElem = 1,numElemsRank
         do igp = (2**ndime)+1,nnode
            call nvtxStartRange("read_press_from_file_Par: var_interpolate")
            call var_interpolate(nnode,aux_2(connecPar(iElem,:)),Ngp_l(igp,:),pr(connecPar(iElem,igp)))
            call nvtxEndRange()
         end do
      end do

      call nvtxStartRange("read_press_from_file_Par: GPU-DATA")
      !$acc update device(pr(:))
      call nvtxEndRange()

   end subroutine read_press_from_file_Par

   subroutine read_temper_from_file_Par(numElemsRank,numNodesRank,totalNumNodesInSerial,file_path,temp,connecPar,Ngp_l,matGidSrlOrdered)
      implicit none
      character(500), intent(in) :: file_path
      integer(4), intent(in)     :: numElemsRank,numNodesRank
      integer(8), intent(in)     :: totalNumNodesInSerial
      real(rp), intent(out)      :: temp(numNodesRank)
      integer(4), intent(in)     :: connecPar(numElemsRank,nnode)
      real(rp), intent(in)       :: Ngp_l(ngaus,nnode)
      integer(8), intent(in)     :: matGidSrlOrdered(numNodesRank,2)
      character(500)             :: file_type, file_name
      integer(4)                 :: iNodeL,iElem,igp,idime,auxCnt,readInd
      integer(8)                 :: iLine,iNodeGSrl
      real(8)                    :: readValue
      real(rp)                   :: aux_2(numNodesRank)

      call nvtxStartRange("read_temper_from_file_Par: CPU-DATA")
      write(file_type,*) ".alya"
      write(file_name,*) "TEMPE"
      if(mpi_rank.eq.0) write(*,*) "--| Reading file TEMPE.alya in parallel..."

      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
      
      auxCnt = 1
      do iLine = 1,totalNumNodesInSerial
         read(99,*) readInd, readValue
         if(iLine.eq.matGidSrlOrdered(auxCnt,2)) then
            iNodeL = matGidSrlOrdered(auxCnt,1)
            auxCnt=auxCnt+1
            if(rp.eq.4) then
               temp(iNodeL) = real(readValue,rp)
            else
               temp(iNodeL) = readValue
            end if
         end if         
      end do
      
      close(99)
      call nvtxEndRange()

      if(mpi_rank.eq.0) write(*,*) "--| Interpolating temperature from file coords to new mesh coords..."
      aux_2(:) = temp(:)
      do iElem = 1,numElemsRank
         do igp = (2**ndime)+1,nnode
            call nvtxStartRange("read_temper_from_file_Par: var_interpolate")
            call var_interpolate(nnode,aux_2(connecPar(iElem,:)),Ngp_l(igp,:),temp(connecPar(iElem,igp)))
            call nvtxEndRange()
         end do
      end do

      call nvtxStartRange("read_temper_from_file_Par: GPU-DATA")
      !$acc update device(temp(:))
      call nvtxEndRange()

   end subroutine read_temper_from_file_Par

end module mod_inicond_reader
