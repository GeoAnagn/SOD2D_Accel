module mod_hdf5
   use hdf5
   use mod_constants
   use mod_mpi
   use mod_nvtx
   use mod_mpi_mesh
   use mod_comms
   use mod_comms_boundaries
   use mod_custom_types
   implicit none

   character(256) :: meshFile_h5_name,surface_meshFile_h5_name
   character(256) :: base_resultsFile_h5_name,base_avgResultsFile_h5_name,base_restartFile_h5_name

   integer(hid_t) :: h5_datatype_uint1,h5_datatype_int1,h5_datatype_int4,h5_datatype_int8
   integer(hid_t) :: h5_datatype_real4,h5_datatype_real8

   real(rp_vtk), allocatable :: auxInterpNodeScalarField(:),auxInterpNodeVectorField(:,:)
   real(rp_vtk), allocatable :: auxNodeScalarField_vtk(:),auxNodeVectorField_vtk(:,:)

contains

   subroutine init_hdf5_interface()
      implicit none
      integer :: h5err

      call nvtxStartRange('init_hdf5_interface: CPU-DATA')
      !.init h5 interface
      call h5open_f(h5err)

      h5_datatype_uint1 = H5T_STD_U8LE
      h5_datatype_int1 = H5T_STD_I8LE
      h5_datatype_int4 = H5T_NATIVE_INTEGER
      h5_datatype_int8 = H5T_STD_I64LE

      h5_datatype_real4 = H5T_NATIVE_REAL
      h5_datatype_real8 = H5T_NATIVE_DOUBLE
      call nvtxEndRange()
   end subroutine init_hdf5_interface

   subroutine end_hdf5_interface()
      implicit none
      integer :: h5err

      !close h5 interface
      call nvtxStartRange('end_hdf5_interface: CPU-DATA')
      call h5close_f(h5err)
      call nvtxEndRange()
   end subroutine end_hdf5_interface

   subroutine init_hdf5_auxiliar_saving_arrays()
      implicit none

      call nvtxStartRange('init_hdf5_auxiliar_saving_arrays: CPU-DATA')
      if(mesh_isLoaded .eqv. .false.) then
         write(*,*) 'FATAL ERROR! Mesh not loaded when calling init_hdf5_auxiliar_saving_arrays()! CRASHING!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      allocate(auxInterpNodeScalarField(numNodesRankPar))
      allocate(auxInterpNodeVectorField(numNodesRankPar,ndime))
      call nvtxEndRange()

      call nvtxStartRange('init_hdf5_auxiliar_saving_arrays: GPU-DATA')
      !$acc enter data create(auxInterpNodeScalarField(:))
      !$acc enter data create(auxInterpNodeVectorField(:,:))
      call nvtxEndRange()

      call nvtxStartRange('init_hdf5_auxiliar_saving_arrays: CPU-DATA')
      allocate(auxNodeScalarField_vtk(numNodesRankPar))
      call nvtxEndRange()

      call nvtxStartRange('init_hdf5_auxiliar_saving_arrays: GPU-DATA')
      !$acc enter data create(auxNodeScalarField_vtk(:))
      call nvtxEndRange()

      call nvtxStartRange('init_hdf5_auxiliar_saving_arrays: CPU-DATA')
      allocate(auxNodeVectorField_vtk(numNodesRankPar,ndime))
      call nvtxEndRange()

      call nvtxStartRange('init_hdf5_auxiliar_saving_arrays: GPU-DATA')
      !$acc enter data create(auxNodeVectorField_vtk(:,:))
      call nvtxEndRange()

   end subroutine

   subroutine end_hdf5_auxiliar_saving_arrays()
      implicit none

      call nvtxStartRange('end_hdf5_auxiliar_saving_arrays: GPU-DATA')
      !$acc exit data delete(auxInterpNodeScalarField(:))
      !$acc exit data delete(auxInterpNodeVectorField(:,:))
      call nvtxEndRange()

      call nvtxStartRange('end_hdf5_auxiliar_saving_arrays: CPU-DATA')
      deallocate(auxInterpNodeScalarField)
      deallocate(auxInterpNodeVectorField)
      call nvtxEndRange()

      call nvtxStartRange('end_hdf5_auxiliar_saving_arrays: GPU-DATA')
      !$acc exit data delete(auxNodeScalarField_vtk(:))
      !$acc exit data delete(auxNodeVectorField_vtk(:,:))
      call nvtxEndRange()

      call nvtxStartRange('end_hdf5_auxiliar_saving_arrays: CPU-DATA')
      deallocate(auxNodeScalarField_vtk)
      deallocate(auxNodeVectorField_vtk)
      call nvtxEndRange()

   end subroutine


   subroutine set_hdf5_meshFile_name(file_path,file_name,numRanks)
      implicit none
      character(len=*), intent(in) :: file_path,file_name
      integer,intent(in) :: numRanks
      character(len=12) :: aux_numRanks

      call nvtxStartRange('set_hdf5_meshFile_name: CPU-DATA')
      write(aux_numRanks,'(I0)') numRanks
      meshFile_h5_name = trim(adjustl(file_path))//trim(adjustl(file_name))//'-'//trim(aux_numRanks)//'.hdf'
      call nvtxEndRange()
   end subroutine set_hdf5_meshFile_name

   subroutine set_hdf5_surface_meshFile_name()
      implicit none

      call nvtxStartRange('set_hdf5_surface_meshFile_name: CPU-DATA')
      surface_meshFile_h5_name = 'surface_'//trim(adjustl(meshFile_h5_name))
      call nvtxEndRange()
   end subroutine set_hdf5_surface_meshFile_name

   subroutine set_hdf5_baseResultsFile_name(res_filePath,res_fileName,mesh_fileName,numRanks)
      implicit none
      character(len=*), intent(in) :: res_filePath,res_fileName,mesh_fileName
      integer,intent(in) :: numRanks
      character(len=12) :: aux_numRanks

      call nvtxStartRange('set_hdf5_baseResultsFile_name: CPU-DATA')
      write(aux_numRanks,'(I0)') numRanks
      base_resultsFile_h5_name = trim(adjustl(res_filePath))//trim(adjustl(res_fileName))//'_'&
         //trim(adjustl(mesh_fileName))//'-'//trim(aux_numRanks)//'_'
      base_avgResultsFile_h5_name = trim(adjustl(res_filePath))//trim(adjustl(res_fileName))//'_AVG_'&
         //trim(adjustl(mesh_fileName))//'-'//trim(aux_numRanks)//'_'

      base_restartFile_h5_name = trim(adjustl(res_filePath))//'restart_'&
         //trim(adjustl(mesh_fileName))//'-'//trim(aux_numRanks)//'_'
      call nvtxEndRange()
   end subroutine set_hdf5_baseResultsFile_name

!---------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------

   subroutine create_hdf5_file(full_hdf5FileName,hdf5_file_id)
      implicit none
      character(*), intent(in) :: full_hdf5FileName
      integer(hid_t),intent(inout) :: hdf5_file_id
      integer(hid_t) :: plist_id
      integer :: h5err

      call nvtxStartRange('create_hdf5_file: CPU-DATA')
      !-----------------------------------------------------------------------------------------------
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_hdf5FileName,H5F_ACC_TRUNC_F,hdf5_file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create results file ',trim(adjustl(full_hdf5FileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)
      call nvtxEndRange()
   end subroutine create_hdf5_file

   subroutine open_hdf5_file(full_hdf5FileName,hdf5_file_id)
      implicit none
      character(*), intent(in) :: full_hdf5FileName
      integer(hid_t),intent(inout) :: hdf5_file_id
      integer(hid_t) :: plist_id
      integer :: h5err

      call nvtxStartRange('open_hdf5_file: CPU-DATA')
      !-----------------------------------------------------------------------------------------------
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_hdf5FileName,H5F_ACC_RDWR_F,hdf5_file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load hdf5 file ',trim(adjustl(full_hdf5FileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)
      call nvtxEndRange()
   end subroutine open_hdf5_file

   subroutine close_hdf5_file(hdf5_file_id)
      implicit none
      integer(hid_t),intent(inout) :: hdf5_file_id
      integer :: h5err

      call nvtxStartRange('close_hdf5_file: CPU-DATA')
      !close the file.
      call h5fclose_f(hdf5_file_id,h5err)
      call nvtxEndRange()
   end subroutine close_hdf5_file

!---------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------

   subroutine create_hdf5_groups_datasets_in_meshFile_from_tool(mnnode,mnpbou,hdf5_file_id,isPeriodic,isBoundaries,isLinealOutput,numMshRanks2Part,numElemsGmsh,numNodesParTotal_i8,&
      vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank,&
      vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank)
      implicit none
      integer(4),intent(in) :: mnnode,mnpbou
      integer(hid_t),intent(in) :: hdf5_file_id
      logical,intent(in) :: isPeriodic,isBoundaries,isLinealOutput
      integer(4),intent(in) :: numMshRanks2Part,numElemsGmsh
      integer(8),intent(in) :: numNodesParTotal_i8
      integer(4),intent(in),dimension(0:numMshRanks2Part-1) :: vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank
      integer(4),intent(in),dimension(0:numMshRanks2Part-1) :: vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank
      integer(4),intent(in),dimension(0:numMshRanks2Part-1) :: vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank
      integer(4),intent(in),dimension(0:numMshRanks2Part-1) :: vecNumPerNodesMshRank
      integer(hid_t) :: dset_id,dspace_id,group_id
      integer(hid_t) :: dtype
      integer(hsize_t) :: ds_dims(1),ds_dims2d(2)
      integer(4) :: ds_rank,h5err
      character(256) :: groupname,dsetname

      !--------------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = numNodesParTotal_i8
      dtype = h5_datatype_int8

      call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_group_hdf5')
      groupname = '/globalIds'
      call create_group_hdf5(hdf5_file_id,groupname)
      call nvtxEndRange()

      call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_dataspace_hdf5')
      dsetname = '/globalIds/globalIdSrl'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_dataspace_hdf5')
      dsetname = '/globalIds/globalIdPar'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_dataspace_hdf5')
      dtype = h5_datatype_int4
      ds_dims(1) = numElemsGmsh
      dsetname = '/globalIds/elemGid'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      !--------------------------------------------------------------------------------
      if(.not.(isLinealOutput)) then
         call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: select_dtype_rp')
         call select_dtype_rp(dtype)
         call nvtxEndRange()

         ds_rank = 2
         ds_dims2d(1) = ndime
         ds_dims2d(2) = numNodesParTotal_i8

         call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_group_hdf5')
         groupname = '/Coords'
         call create_group_hdf5(hdf5_file_id,groupname)
         call nvtxEndRange()

         call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_dataspace_hdf5')
         dsetname = '/Coords/Points'
         call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims2d,dtype)
         call nvtxEndRange()
      end if
      !--------------------------------------------------------------------------------

      call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_groups_datasets_connectivity_workingNodes_hdf5')
      call create_groups_datasets_connectivity_workingNodes_hdf5(mnnode,hdf5_file_id,numMshRanks2Part,numElemsGmsh,vecNumWorkingNodes)
      call nvtxEndRange()

      if(isBoundaries) then
         call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_groups_datasets_boundary_data_hdf5')
         call create_groups_datasets_boundary_data_hdf5(mnpbou,hdf5_file_id,numMshRanks2Part,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank)
         call nvtxEndRange()
      end if

      if(numMshRanks2Part.ge.2) then
         call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_groups_datasets_parallel_data_hdf5')
         call create_groups_datasets_parallel_data_hdf5(hdf5_file_id,numMshRanks2Part,vecNumMshRanksWithComms,vecNumNodesToCommMshRank)
         call nvtxEndRange()
         if(isBoundaries) then
            call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_groups_datasets_parallel_data_boundary_hdf5')
            call create_groups_datasets_parallel_data_boundary_hdf5(hdf5_file_id,numMshRanks2Part,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank)
            call nvtxEndRange()
         end if
      end if

      if(isPeriodic) then
         call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_groups_datasets_periodic_data_hdf5')
         call create_groups_datasets_periodic_data_hdf5(hdf5_file_id,numMshRanks2Part,vecNumPerNodesMshRank)
         call nvtxEndRange()
      end if

      !---------------------------------------------------------
      ! order

      dtype = h5_datatype_int4

      call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_group_hdf5')
      groupname = '/order'
      call create_group_hdf5(hdf5_file_id,groupname)
      call nvtxEndRange()

      call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_dataspace_hdf5')
      ds_rank = 1
      ds_dims(1) = 1
      dsetname = '/order/porder'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = mnnode
      dsetname = '/order/a2ijk'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/order/gmsh2ijk'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/order/vtk2ijk'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = mnpbou
      dsetname = '/order/a2ij'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/order/gmsh2ij'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/order/vtk2ij'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      !---------------------------------------------------------
      ! meshOutputInfo

      call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_group_hdf5')
      groupname = '/meshOutputInfo'
      call create_group_hdf5(hdf5_file_id,groupname)
      call nvtxEndRange()

      call nvtxStartRange('create_hdf5_groups_datasets_in_meshFile_from_tool: create_dataspace_hdf5')
      dtype = h5_datatype_uint1
      ds_rank = 1
      ds_dims(1) = 1
      dsetname = '/meshOutputInfo/isLinealOutput'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)

      dtype = h5_datatype_int4
      dsetname = '/meshOutputInfo/nnodeVTK'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/meshOutputInfo/numVTKElemsPerMshElem'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = numMshRanks2Part
      dsetname = '/meshOutputInfo/numElemsVTKMshRank'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/meshOutputInfo/sizeConnecVTKMshRank'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

   end subroutine create_hdf5_groups_datasets_in_meshFile_from_tool

   subroutine create_groups_datasets_connectivity_workingNodes_hdf5(mnnode,file_id,numMshRanks2Part,numElemsGmsh,vecNumWorkingNodes)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: mnnode,numMshRanks2Part,numElemsGmsh
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumWorkingNodes
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims
      integer(hid_t) :: dtype
      integer(4) :: ds_rank,mshRank,accumVal

      call nvtxStartRange('create_groups_datasets_connectivity_workingNodes_hdf5: create_group_hdf5')
      groupname = trim('/Connectivity')
      call create_group_hdf5(file_id,groupname)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_connectivity_workingNodes_hdf5: CPU-COMPUTE')
      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = int(numElemsGmsh,hsize_t)*int(mnnode,hsize_t)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_connectivity_workingNodes_hdf5: create_dataspace_hdf5')
      dsetname = '/Connectivity/connecParOrig'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Connectivity/connecParWork'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = numMshRanks2Part
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_connectivity_workingNodes_hdf5: CPU-COMPUTE')
      ds_dims(1)=0
      do mshRank=0,numMshRanks2Part-1
         ds_dims(1)=ds_dims(1)+int(vecNumWorkingNodes(mshRank),hsize_t)
      end do
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_connectivity_workingNodes_hdf5: create_dataspace_hdf5')
      dsetname = '/Connectivity/workingNodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

   end subroutine create_groups_datasets_connectivity_workingNodes_hdf5

   subroutine create_groups_datasets_parallel_data_hdf5(file_id,numMshRanks2Part,vecNumMshRanksWithComms,vecNumNodesToCommMshRank)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: numMshRanks2Part
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumMshRanksWithComms,vecNumNodesToCommMshRank
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,mshRank,accumVal

      call nvtxStartRange('create_groups_datasets_parallel_data_hdf5: create_group_hdf5')
      groupname = trim('/Parallel_data')
      call create_group_hdf5(file_id,groupname)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_parallel_data_hdf5: create_dataspace_hdf5')
      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = numMshRanks2Part
      dsetname = '/Parallel_data/rankNodeStart'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/rankNodeEnd'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/rankElemStart'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/rankElemEnd'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/numRanksWithComms'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/numNodesToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_parallel_data_hdf5: CPU-COMPUTE')
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumMshRanksWithComms(mshRank)
      end do
      ds_dims(1) = accumVal
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_parallel_data_hdf5: create_dataspace_hdf5')
      dsetname = '/Parallel_data/ranksToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/commsMemPosInLoc'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/commsMemSize'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/commsMemPosInNgb'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_parallel_data_hdf5: CPU-COMPUTE')
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumNodesToCommMshRank(mshRank)
      end do
      ds_dims(1) = accumVal
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_parallel_data_hdf5: create_dataspace_hdf5')
      dsetname = '/Parallel_data/nodesToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

   end subroutine create_groups_datasets_parallel_data_hdf5

   subroutine create_groups_datasets_parallel_data_boundary_hdf5(file_id,numMshRanks2Part,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: numMshRanks2Part
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,mshRank,accumVal

      call nvtxStartRange('create_groups_datasets_parallel_data_boundary_hdf5: create_group_hdf5')
      groupname = trim('/Parallel_data_boundary')
      call create_group_hdf5(file_id,groupname)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_parallel_data_boundary_hdf5: create_dataspace_hdf5')
      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = numMshRanks2Part
      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_parallel_data_boundary_hdf5: CPU-COMPUTE')
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecBndNumMshRanksWithComms(mshRank)
      end do
      ds_dims(1) = accumVal
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_parallel_data_boundary_hdf5: create_dataspace_hdf5')
      dsetname = '/Parallel_data_boundary/ranksToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data_boundary/commsMemSize'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_parallel_data_boundary_hdf5: CPU-COMPUTE')
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecBndNumNodesToCommMshRank(mshRank)
      end do
      ds_dims(1) = accumVal
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_parallel_data_boundary_hdf5: create_dataspace_hdf5')
      dsetname = '/Parallel_data_boundary/nodesToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

   end subroutine create_groups_datasets_parallel_data_boundary_hdf5

   subroutine create_groups_datasets_periodic_data_hdf5(file_id,numMshRanks2Part,vecNumPerNodesMshRank)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: numMshRanks2Part
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumPerNodesMshRank
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,mshRank,accumVal

      call nvtxStartRange('create_groups_datasets_periodic_data_hdf5: create_group_hdf5')
      groupname = trim('/Periodic_data')
      call create_group_hdf5(file_id,groupname)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_periodic_data_hdf5: create_dataspace_hdf5')
      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = numMshRanks2Part
      dsetname = '/Periodic_data/nPerRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_periodic_data_hdf5: CPU-COMPUTE')
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumPerNodesMshRank(mshRank)
      end do
      ds_dims(1) = accumVal
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_periodic_data_hdf5: create_dataspace_hdf5')
      dsetname = '/Periodic_data/masSlaRankPar1'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Periodic_data/masSlaRankPar2'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

   end subroutine create_groups_datasets_periodic_data_hdf5

   subroutine create_groups_datasets_boundary_data_hdf5(mnpbou,file_id,numMshRanks2Part,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: mnpbou,numMshRanks2Part
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,h5err
      integer :: mshRank,accumVal

      call nvtxStartRange('create_groups_datasets_boundary_data_hdf5: create_group_hdf5')
      groupname = trim('/Boundary_data')
      call create_group_hdf5(file_id,groupname)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_boundary_data_hdf5: create_dataspace_hdf5')
      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = 1
      dsetname = '/Boundary_data/numBoundCodes'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = numMshRanks2Part
      dsetname = '/Boundary_data/numBoundsRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Boundary_data/ndofRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Boundary_data/numBoundaryNodesRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_boundary_data_hdf5: CPU-COMPUTE')
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumBoundFacesMshRank(mshRank)
      end do
      ds_dims(1) = accumVal
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_boundary_data_hdf5: create_dataspace_hdf5')
      dsetname = '/Boundary_data/bouCodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = ds_dims(1)*mnpbou
      dsetname = '/Boundary_data/boundPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Boundary_data/boundParOrig'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_boundary_data_hdf5: CPU-COMPUTE')
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumDoFMshRank(mshRank)
      end do
      ds_dims(1) = accumVal
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_boundary_data_hdf5: create_dataspace_hdf5')
      dsetname = '/Boundary_data/ldofPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_boundary_data_hdf5: CPU-COMPUTE')
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumBoundaryNodesMshRank(mshRank)
      end do
      ds_dims(1) = accumVal
      call nvtxEndRange()

      call nvtxStartRange('create_groups_datasets_boundary_data_hdf5: create_dataspace_hdf5')
      dsetname = '/Boundary_data/lbnodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()
      
   end subroutine create_groups_datasets_boundary_data_hdf5

   subroutine write_mshRank_data_in_hdf5_meshFile_from_tool(mporder,mnnode,mnpbou,hdf5_file_id,mshRank,numMshRanks2Part,isPeriodic,isBoundaries,isLinealOutput,numElemsGmsh,numBoundFacesGmsh,&
      numElemsMshRank,mshRankElemStart,mshRankElemEnd,mshRankNodeStart_i8,mshRankNodeEnd_i8,numNodesMshRank,numWorkingNodesMshRank,numBoundFacesMshRank,numBoundaryNodesMshRank,numDoFMshRank,maxBoundCode,&
      numElemsVTKMshRank,sizeConnecVTKMshRank,mnnodeVTK,numVTKElemsPerMshElem,&
      a2ijk,a2ij,gmsh2ijk,gmsh2ij,vtk2ijk,vtk2ij,&
      elemGidMshRank,globalIdSrlMshRank_i8,globalIdParMshRank_i8,connecParOrigMshRank,connecParWorkMshRank,coordParMshRank,workingNodesMshRank,&
      boundaryNodesMshRank,dofNodesMshRank,boundFacesCodesMshRank,boundFacesOrigMshRank,boundFacesMshRank,numPerNodesMshRank,masSlaNodesMshRank,&
      numNodesToCommMshRank,numMshRanksWithComms,nodesToCommMshRank,commsMemPosInLocMshRank,commsMemSizeMshRank,commsMemPosInNgbMshRank,ranksToCommMshRank,&
      bnd_numNodesToCommMshRank,bnd_numMshRanksWithComms,bnd_nodesToCommMshRank,bnd_commsMemPosInLocMshRank,bnd_commsMemSizeMshRank,bnd_commsMemPosInNgbMshRank,bnd_ranksToCommMshRank,&
      vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank)
      implicit none
      integer(hid_t),intent(in) :: hdf5_file_id
      logical,intent(in) :: isPeriodic,isBoundaries,isLinealOutput
      integer(4),intent(in) :: mporder,mnnode,mnpbou,mshRank,numMshRanks2Part,numElemsGmsh,numBoundFacesGmsh
      integer(4),intent(in) :: numElemsMshRank,mshRankElemStart,mshRankElemEnd
      integer(8),intent(in) :: mshRankNodeStart_i8,mshRankNodeEnd_i8
      integer(4),intent(in) :: numNodesMshRank,numWorkingNodesMshRank,numPerNodesMshRank
      integer(4),intent(in) :: numBoundFacesMshRank,numBoundaryNodesMshRank,numDoFMshRank
      integer(4),intent(in) :: maxBoundCode,numElemsVTKMshRank,sizeConnecVTKMshRank,mnnodeVTK,numVTKElemsPerMshElem

      integer(4),intent(in),dimension(mnnode) :: a2ijk,gmsh2ijk,vtk2ijk
      integer(4),intent(in),dimension(mnpbou) :: a2ij,gmsh2ij,vtk2ij
      integer(4),intent(in) :: elemGidMshRank(numElemsMshRank),workingNodesMshRank(numWorkingNodesMshRank)
      integer(8),intent(in) :: globalIdSrlMshRank_i8(numNodesMshRank),globalIdParMshRank_i8(numNodesMshRank)
      integer(4),intent(in) :: connecParOrigMshRank(numElemsMshRank,mnnode),connecParWorkMshRank(numElemsMshRank,mnnode)
      real(rp),intent(in)   :: coordParMshRank(numNodesMshRank,3)

      integer(4),intent(in) :: boundaryNodesMshRank(numBoundaryNodesMshRank),dofNodesMshRank(numDoFMshRank)
      integer(4),intent(in) :: boundFacesCodesMshRank(numBoundFacesMshRank),boundFacesOrigMshRank(numBoundFacesMshRank,mnpbou),boundFacesMshRank(numBoundFacesMshRank,mnpbou)
      integer(4),intent(in) :: numNodesToCommMshRank,numMshRanksWithComms
      integer(4),intent(in) :: masSlaNodesMshRank(numPerNodesMshRank,2)
      integer(4),intent(in),dimension(numMshRanksWithComms) :: nodesToCommMshRank,commsMemPosInLocMshRank,commsMemSizeMshRank,commsMemPosInNgbMshRank,ranksToCommMshRank
      integer(4),intent(in) :: bnd_numNodesToCommMshRank,bnd_numMshRanksWithComms
      integer(4),intent(in),dimension(bnd_numMshRanksWithComms) :: bnd_nodesToCommMshRank,bnd_commsMemPosInLocMshRank,bnd_commsMemSizeMshRank,bnd_commsMemPosInNgbMshRank,bnd_ranksToCommMshRank

      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank
      
      character(128) :: dsetname
      integer(hsize_t) :: ms_dims(1),ms_dims2d(2)
      integer(4) :: i,m,iBound,iElemL
      integer(hssize_t) :: ms_offset(1),ms_offset2d(2)
      integer(1),allocatable :: aux_array_i1(:)
      integer(4),allocatable :: aux_array_i4(:)
      integer(8),allocatable :: aux_array_i8(:)
      
      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int8_hyperslab_parallel')
      ms_dims(1) = int(numNodesMshRank,hsize_t)
      ms_offset(1) = int(mshRankNodeStart_i8,hssize_t)-1
      dsetname = '/globalIds/globalIdSrl'
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,globalIdSrlMshRank_i8)

      dsetname = '/globalIds/globalIdPar'
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,globalIdParMshRank_i8)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      ms_dims(1) = int(numElemsMshRank,hsize_t)
      ms_offset(1) = int(mshRankElemStart,hssize_t)-1
      dsetname = '/globalIds/elemGid'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,elemGidMshRank)
      call nvtxEndRange()

      !---------------------------------------------------------------------------------------------------------------------
      if(.not.(isLinealOutput)) then
         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_2d_tr_real_rp_hyperslab_parallel')
         ms_dims2d(1) = ndime
         ms_dims2d(2) = int(numNodesMshRank,hsize_t)
         ms_offset2d(1) = 0
         ms_offset2d(2) = int(mshRankNodeStart_i8,hssize_t)-1
         dsetname = '/Coords/Points'
         call write_dataspace_2d_tr_real_rp_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims2d,ms_offset2d,coordParMshRank)
         call nvtxEndRange()
      end if
      !---------------------------------------------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------------------------------------------
      if(numMshRanks2Part.ge.2) then
         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
         ms_dims(1) = 1
         ms_offset(1) = int(mshRank,hssize_t)
         allocate(aux_array_i8(1))
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int8_hyperslab_parallel')
         dsetname = '/Parallel_data/rankNodeStart'
         aux_array_i8(1)=mshRankNodeStart_i8
         call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

         dsetname = '/Parallel_data/rankNodeEnd'
         aux_array_i8(1)=mshRankNodeEnd_i8
         call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
         deallocate(aux_array_i8)
         allocate(aux_array_i4(1))
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Parallel_data/rankElemStart'
         aux_array_i4(1)=mshRankElemStart
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)

         dsetname = '/Parallel_data/rankElemEnd'
         aux_array_i4(1)=mshRankElemEnd
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)

         dsetname = '/Parallel_data/numRanksWithComms'
         aux_array_i4(1)=numMshRanksWithComms
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)

         dsetname = '/Parallel_data/numNodesToComm'
         aux_array_i4(1)=numNodesToCommMshRank
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
         deallocate(aux_array_i4)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-COMPUTE')
         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumMshRanksWithComms(i),hssize_t)
         end do
         ms_dims(1)=int(numMshRanksWithComms,hsize_t)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Parallel_data/ranksToComm'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,ranksToCommMshRank)

         dsetname = '/Parallel_data/commsMemPosInLoc'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,commsMemPosInLocMshRank)

         dsetname = '/Parallel_data/commsMemSize'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,commsMemSizeMshRank)

         dsetname = '/Parallel_data/commsMemPosInNgb'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,commsMemPosInNgbMshRank)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-COMPUTE')
         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumNodesToCommMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numNodesToCommMshRank,hsize_t)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Parallel_data/nodesToComm'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,nodesToCommMshRank)
         call nvtxEndRange()
      end if
      !------------------------------------------------------------------------------------------------------------------------

      if(isBoundaries) then

         dsetname = '/Boundary_data/numBoundCodes'
         if(mshRank.eq.0) then
            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
            ms_dims(1) = 1
            ms_offset(1) = int(mshRank,hssize_t)
            allocate(aux_array_i4(1))
            aux_array_i4(1)=maxBoundCode
            call nvtxEndRange()

            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
            call nvtxEndRange()

            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
            deallocate(aux_array_i4)
            call nvtxEndRange()
         else
            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
            ms_dims(1) = 0
            ms_offset(1) = 0
            allocate(aux_array_i4(0))
            call nvtxEndRange()

            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
            call nvtxEndRange()

            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
            deallocate(aux_array_i4)
            call nvtxEndRange()
         end if

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
         ms_dims(1) = 1
         ms_offset(1) = int(mshRank,hssize_t)
         allocate(aux_array_i4(1))
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Boundary_data/numBoundsRankPar'
         aux_array_i4(1)=numBoundFacesMshRank
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)

         dsetname = '/Boundary_data/ndofRankPar'
         aux_array_i4(1)=numDoFMshRank
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)

         dsetname = '/Boundary_data/numBoundaryNodesRankPar'
         aux_array_i4(1)=numBoundaryNodesMshRank
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
         deallocate(aux_array_i4)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-COMPUTE')
         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumDoFMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numDoFMshRank,hsize_t)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Boundary_data/ldofPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,dofNodesMshRank)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-COMPUTE')
         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumBoundaryNodesMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numBoundaryNodesMshRank,hsize_t)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Boundary_data/lbnodesPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,boundaryNodesMshRank)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-COMPUTE')
         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumBoundFacesMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numBoundFacesMshRank,hsize_t)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Boundary_data/bouCodesPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,boundFacesCodesMshRank)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
         allocate(aux_array_i4(numBoundFacesMshRank*mnpbou))

         i=1
         do iBound=1,numBoundFacesMshRank
            do m=1,mnpbou
               aux_array_i4(i)=boundFacesMshRank(iBound,m)
               i=i+1
            end do
         end do
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-COMPUTE')
         ms_dims(1) = int(ms_dims(1),hsize_t)*int(mnpbou,hsize_t)
         ms_offset(1) = ms_offset(1)*int(mnpbou,hssize_t)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Boundary_data/boundPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
         i=1
         do iBound=1,numBoundFacesMshRank
            do m=1,mnpbou
               aux_array_i4(i)=boundFacesOrigMshRank(iBound,m)
               i=i+1
            end do
         end do
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Boundary_data/boundParOrig'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
         deallocate(aux_array_i4)
         call nvtxEndRange()

         if(numMshRanks2Part.ge.2) then
            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
            ms_dims(1) = 1
            ms_offset(1) = int(mshRank,hssize_t)
            allocate(aux_array_i4(1))
            call nvtxEndRange()

            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
            dsetname = '/Parallel_data_boundary/numRanksWithComms'
            aux_array_i4(1)=bnd_numMshRanksWithComms
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)

            dsetname = '/Parallel_data_boundary/numNodesToComm'
            aux_array_i4(1)=bnd_numNodesToCommMshRank
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
            call nvtxEndRange()

            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
            deallocate(aux_array_i4)
            call nvtxEndRange()

            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-COMPUTE')
            ms_offset(1)=0
            do i=0,mshRank-1 !from rank 0 mpi_rank-1
               ms_offset(1)=ms_offset(1)+int(vecBndNumMshRanksWithComms(i),hssize_t)
            end do
            ms_dims(1)=int(bnd_numMshRanksWithComms,hsize_t)
            call nvtxEndRange()

            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
            dsetname = '/Parallel_data_boundary/ranksToComm'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,bnd_ranksToCommMshRank)

            dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,bnd_commsMemPosInLocMshRank)

            dsetname = '/Parallel_data_boundary/commsMemSize'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,bnd_commsMemSizeMshRank)

            dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,bnd_commsMemPosInNgbMshRank)
            call nvtxEndRange()

            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-COMPUTE')
            ms_offset(1)=0
            do i=0,mshRank-1 !from rank 0 mpi_rank-1
               ms_offset(1)=ms_offset(1)+int(vecBndNumNodesToCommMshRank(i),hssize_t)
            end do
            ms_dims(1)=int(bnd_numNodesToCommMshRank,hsize_t)
            call nvtxEndRange()

            call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
            dsetname = '/Parallel_data_boundary/nodesToComm'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,bnd_nodesToCommMshRank)
            call nvtxEndRange()
         end if
      end if

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      ms_dims(1) = 1
      ms_offset(1) = int(mshRank,hssize_t)
      allocate(aux_array_i4(1))
      call nvtxEndRange()
      
      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      aux_array_i4(1)=numWorkingNodesMshRank
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      deallocate(aux_array_i4)
      allocate(aux_array_i4(numElemsMshRank*mnnode)) ! SAVING connecParOrig(:,:)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-COMPUTE')
      ms_dims(1) = int(numElemsMshRank,hsize_t)*int(mnnode,hsize_t)
      ms_offset(1) = int((mshRankElemStart-1),hssize_t)*int(mnnode,hssize_t)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      i=1
      do iElemL=1,numElemsMshRank
         do m=1,mnnode
            aux_array_i4(i)=connecParOrigMshRank(iElemL,m)
            i=i+1
         end do
      end do
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Connectivity/connecParOrig'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      call nvtxEndRange()

      !  SAVING connecParWork(:,:)
      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')      
      i=1
      do iElemL=1,numElemsMshRank
         do m=1,mnnode
            aux_array_i4(i)=connecParWorkMshRank(iElemL,m)
            i=i+1
         end do
      end do
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Connectivity/connecParWork'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      deallocate(aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-COMPUTE')
      ms_offset(1)=0
      do i=0,mshRank-1 !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(vecNumWorkingNodes(i),hssize_t)
      end do
      ms_dims(1)=int(numWorkingNodesMshRank,hsize_t)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Connectivity/workingNodesPar'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,workingNodesMshRank)
      call nvtxEndRange()

      if(isPeriodic) then
         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
         ms_dims(1) = 1
         ms_offset(1) = int(mshRank,hssize_t)
         allocate(aux_array_i4(1))
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Periodic_data/nPerRankPar'
         aux_array_i4(1)= numPerNodesMshRank
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
         deallocate(aux_array_i4)
         allocate(aux_array_i4(numPerNodesMshRank*2))
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-COMPUTE') 
         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumPerNodesMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numPerNodesMshRank,hsize_t)
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Periodic_data/masSlaRankPar1'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,masSlaNodesMshRank(:,1))
         call nvtxEndRange()

         call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Periodic_data/masSlaRankPar2'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,masSlaNodesMshRank(:,2))
         call nvtxEndRange()
      end if

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      dsetname = '/order/porder'
      ms_dims(1) = 0
      ms_offset(1) = 0
      if(mshRank.eq.0) then
         ms_dims(1) = 1
      endif

      allocate(aux_array_i4(ms_dims(1)))

      if(mshRank.eq.0) then
         aux_array_i4(1) = mporder
      end if
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      deallocate(aux_array_i4)
   
      ms_dims(1) = 0
      ms_offset(1) = 0
      if(mshRank.eq.0) then
         ms_dims(1) = mnnode
      endif

      allocate(aux_array_i4(ms_dims(1)))

      dsetname = '/order/a2ijk'
      if(mshRank.eq.0) then
         aux_array_i4(:) = a2ijk(:)
      end if
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      dsetname = '/order/gmsh2ijk'
      if(mshRank.eq.0) then
         aux_array_i4(:) = gmsh2ijk(:)
      end if
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      dsetname = '/order/vtk2ijk'
      if(mshRank.eq.0) then
         aux_array_i4(:) = vtk2ijk(:)
      end if
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      deallocate(aux_array_i4)

      ms_dims(1) = 0
      ms_offset(1) = 0
      if(mshRank.eq.0) then
         ms_dims(1) = mnpbou
      endif

      allocate(aux_array_i4(ms_dims(1)))

      dsetname = '/order/a2ij'
      if(mshRank.eq.0) then
         aux_array_i4(:) = a2ij(:)
      end if
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      dsetname = '/order/gmsh2ij'
      if(mshRank.eq.0) then
         aux_array_i4(:) = gmsh2ij(:)
      end if
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      dsetname = '/order/vtk2ij'
      if(mshRank.eq.0) then
         aux_array_i4(:) = vtk2ij(:)
      end if
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      deallocate(aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      ms_dims(1) = 0
      ms_offset(1) = 0
      if(mshRank.eq.0) then
         ms_dims(1) = 1
      endif

      allocate(aux_array_i1(ms_dims(1)))

      dsetname = '/meshOutputInfo/isLinealOutput'
      if(mshRank.eq.0) then
         aux_array_i1(1) = 0
         if(isLinealOutput) aux_array_i1(1) = 1
      end if
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_uint1_hyperslab_parallel')
      call write_dataspace_1d_uint1_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i1)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      deallocate(aux_array_i1)
      allocate(aux_array_i4(ms_dims(1)))

      dsetname = '/meshOutputInfo/nnodeVTK'
      if(mshRank.eq.0) then
         aux_array_i4(1) = mnnodeVTK
      end if
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      dsetname = '/meshOutputInfo/numVTKElemsPerMshElem'
      if(mshRank.eq.0) then
         aux_array_i4(1) = numVTKElemsPerMshElem
      end if
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)

      ms_dims(1) = 1
      ms_offset(1) = int(mshRank,hssize_t)
      dsetname = '/meshOutputInfo/numElemsVTKMshRank'
      aux_array_i4(1) = numElemsVTKMshRank
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)

      dsetname = '/meshOutputInfo/sizeConnecVTKMshRank'
      aux_array_i4(1) = sizeConnecVTKMshRank
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      deallocate(aux_array_i4)
      call nvtxEndRange()

   end subroutine write_mshRank_data_in_hdf5_meshFile_from_tool

   subroutine dummy_write_mshRank_data_in_hdf5_meshFile_from_tool(hdf5_file_id,numMshRanks2Part,isPeriodic,isBoundaries,isLinealOutput)
      implicit none
      integer(hid_t),intent(in) :: hdf5_file_id
      logical,intent(in) :: isPeriodic,isBoundaries,isLinealOutput
      integer(4),intent(in) :: numMshRanks2Part
      
      character(128) :: dsetname
      integer(hsize_t) :: ms_dims(1),ms_dims2d(2)
      integer(hssize_t) :: ms_offset(1),ms_offset2d(2)
      integer(1),allocatable :: empty_array_i1(:)
      integer(4),allocatable :: empty_array_i4(:)
      integer(8),allocatable :: empty_array_i8(:)
      real(rp),allocatable :: empty_array_rp(:),empty_array2d_rp(:,:)

      call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      allocate(empty_array_i1(0))
      allocate(empty_array_i4(0))
      allocate(empty_array_i8(0))
      allocate(empty_array_rp(0))
      allocate(empty_array2d_rp(0,0))

      ms_dims(1) = 0
      ms_offset(1) = 0
      ms_dims2d(1) = 0
      ms_dims2d(2) = 0
      ms_offset2d(1) = 0
      ms_offset2d(2) = 0
      call nvtxEndRange()

      call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int8_hyperslab_parallel')
      dsetname = '/globalIds/globalIdSrl'
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

      dsetname = '/globalIds/globalIdPar'
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/globalIds/elemGid'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
      call nvtxEndRange()

      if(.not.(isLinealOutput)) then
         call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_2d_tr_real_rp_hyperslab_parallel')
         dsetname = '/Coords/Points'
         call write_dataspace_2d_tr_real_rp_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims2d,ms_offset2d,empty_array2d_rp)
         call nvtxEndRange()
      end if

      if(numMshRanks2Part.ge.2) then
         call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int8_hyperslab_parallel')
         dsetname = '/Parallel_data/rankNodeStart'
         call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

         dsetname = '/Parallel_data/rankNodeEnd'
         call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i8)
         call nvtxEndRange()

         call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Parallel_data/rankElemStart'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/rankElemEnd'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/numRanksWithComms'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/numNodesToComm'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/ranksToComm'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/commsMemPosInLoc'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/commsMemSize'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/commsMemPosInNgb'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/nodesToComm'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
         call nvtxEndRange()
      end if

      if(isBoundaries) then
         call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Boundary_data/numBoundCodes'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/numBoundsRankPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/ndofRankPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/numBoundaryNodesRankPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/ldofPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/lbnodesPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/bouCodesPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/boundPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/boundParOrig'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
         call nvtxEndRange()

         if(numMshRanks2Part.ge.2) then
            call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
            dsetname = '/Parallel_data_boundary/numRanksWithComms'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/numNodesToComm'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/ranksToComm'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/commsMemSize'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/nodesToComm'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
            call nvtxEndRange()
         end if
      end if

      call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
      
      dsetname = '/Connectivity/connecParOrig'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
      
      dsetname = '/Connectivity/connecParWork'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
      
      dsetname = '/Connectivity/workingNodesPar'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
      call nvtxEndRange()

      if(isPeriodic) then
         call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Periodic_data/nPerRankPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Periodic_data/masSlaRankPar1'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Periodic_data/masSlaRankPar2'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
         call nvtxEndRange()
      end if

      call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/order/porder'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/order/a2ijk'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/order/a2ij'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/order/gmsh2ijk'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/order/gmsh2ij'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/order/vtk2ijk'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/order/vtk2ij'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_uint1_hyperslab_parallel')
      dsetname = '/meshOutputInfo/isLinealOutput'
      call write_dataspace_1d_uint1_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i1)
      call nvtxEndRange()

      call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: write_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/meshOutputInfo/nnodeVTK'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/meshOutputInfo/numVTKElemsPerMshElem'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/meshOutputInfo/numElemsVTKMshRank'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/meshOutputInfo/sizeConnecVTKMshRank'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
      call nvtxEndRange()

      call nvtxStartRange('dummy_write_mshRank_data_in_hdf5_meshFile_from_tool: CPU-DATA')
      deallocate(empty_array_i1)
      deallocate(empty_array_i4)
      deallocate(empty_array_i8)
      deallocate(empty_array_rp)
      deallocate(empty_array2d_rp)
      call nvtxEndRange()
   end subroutine dummy_write_mshRank_data_in_hdf5_meshFile_from_tool

   subroutine load_hdf5_meshFile(mnnode,mnpbou)
      implicit none
      integer(4),intent(in) :: mnnode,mnpbou
      character(256) :: groupname,dsetname
      integer(hid_t) :: file_id,dset_id,fspace_id
      integer(4) :: h5err
      integer(hsize_t),dimension(1) :: fs_dims,fs_maxdims
      integer(1) :: aux_array_i1(1)

      if(mpi_rank.eq.0) write(*,*) '# Loading hdf5 mesh: ',trim(adjustl(meshFile_h5_name))

      call nvtxStartRange('load_hdf5_meshFile: open_hdf5_file')
      call open_hdf5_file(meshFile_h5_name,file_id)
      call nvtxEndRange()

      call nvtxStartRange('load_hdf5_meshFile: CPU-DATA')
      dsetname = '/globalIds/globalIdSrl'
      call h5dopen_f(file_id, dsetname, dset_id, h5err)
      call h5dget_space_f(dset_id, fspace_id, h5err) !get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err) !get dimensions of the filespace
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
      totalNumNodesPar = fs_dims(1)
      dsetname = '/globalIds/elemGid'
      call h5dopen_f(file_id, dsetname, dset_id, h5err)
      call h5dget_space_f(dset_id, fspace_id, h5err) !get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err) !get dimensions of the filespace
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
      totalNumElements = fs_dims(1)
      call nvtxEndRange()

      call nvtxStartRange('load_hdf5_meshFile: load_meshOutputInfo_hdf5')      
      call load_meshOutputInfo_hdf5(file_id)! meshOutputInfo
      call nvtxEndRange()

      call nvtxStartRange('load_hdf5_meshFile: load_porder_hdf5')
      call load_porder_hdf5(file_id,mnnode,mnpbou) !load porder
      call nvtxEndRange()
      
      !load the parallel data
      if(mpi_size.ge.2) then
         call nvtxStartRange('load_hdf5_meshFile: load_parallel_data_hdf5')
         call load_parallel_data_hdf5(file_id)
         call nvtxEndRange()
      else !only 1 rank
         call nvtxStartRange('load_hdf5_meshFile: CPU-DATA')
         numNodesRankPar = totalNumNodesPar
         numElemsRankPar = totalNumElements
         rankNodeStart = 1
         rankNodeEnd = numNodesRankPar
         rankElemStart = 1
         rankElemEnd = numElemsRankPar

         numNodesToComm=0
         numRanksWithComms=0
         allocate(nodesToComm(numNodesToComm))
         allocate(ranksToComm(numRanksWithComms))
         allocate(commsMemPosInLoc(numRanksWithComms))
         allocate(commsMemPosInNgb(numRanksWithComms))
         allocate(commsMemSize(numRanksWithComms))
         call nvtxEndRange()
      end if

      call nvtxStartRange('load_hdf5_meshFile: load_periodic_data_hdf5')
      call load_periodic_data_hdf5(file_id) !load periodic data
      call nvtxEndRange()

      call nvtxStartRange('load_hdf5_meshFile: load_boundary_data_hdf5')
      call load_boundary_data_hdf5(file_id,mnnode,mnpbou) !load boundary data
      call nvtxEndRange()

      if((isMeshBoundaries).and.(mpi_size.ge.2)) then
         call nvtxStartRange('load_hdf5_meshFile: load_parallel_data_boundary_hdf5')
         call load_parallel_data_boundary_hdf5(file_id)
         call nvtxEndRange()
      else
         call nvtxStartRange('load_hdf5_meshFile: CPU-DATA')
         bnd_numNodesToComm=0
         bnd_numRanksWithComms=0
         allocate(bnd_nodesToComm(bnd_numNodesToComm))
         allocate(bnd_ranksToComm(bnd_numRanksWithComms))
         allocate(bnd_commsMemPosInLoc(bnd_numRanksWithComms))
         allocate(bnd_commsMemPosInNgb(bnd_numRanksWithComms))
         allocate(bnd_commsMemSize(bnd_numRanksWithComms))
         call nvtxEndRange()
      end if

      call nvtxStartRange('load_hdf5_meshFile: load_coordinates_hdf5')
      call load_coordinates_hdf5(file_id) !load the coordinates
      call nvtxEndRange()

      call nvtxStartRange('load_hdf5_meshFile: load_connectivity_hdf5')
      call load_connectivity_hdf5(file_id,mnnode)!load connectivity
      call nvtxEndRange()

      call nvtxStartRange('load_hdf5_meshFile: load_globalIds_hdf5')
      call load_globalIds_hdf5(file_id)!load globalIds
      call nvtxEndRange()

      call nvtxStartRange('load_hdf5_meshFile: close_hdf5_file')
      call close_hdf5_file(file_id) !close h5 file
      call nvtxEndRange()

      if(mpi_rank.eq.0) write(*,*) '# Mesh ',trim(adjustl(meshFile_h5_name)),' succesfully loaded!'
      mesh_isLoaded = .true.

   end subroutine load_hdf5_meshFile

   subroutine create_group_hdf5(file_id,groupname)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*), intent(in) :: groupname
      integer(hid_t) :: group_id
      integer(4) :: h5err

      call nvtxStartRange('create_group_hdf5: CPU-DATA')
      call h5gcreate_f(file_id,groupname,group_id,h5err)
      call h5gclose_f(group_id, h5err)
      call nvtxEndRange()

   end subroutine create_group_hdf5

   subroutine create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ds_rank
      integer(hsize_t),intent(in) :: ds_dims(ds_rank)
      integer(hid_t),intent(in) :: dtype
      integer(hid_t) :: dset_id,dspace_id
      integer(4) :: h5err

      call nvtxStartRange('create_dataspace_hdf5: CPU-DATA')
      ! Create the data space for the  dataset.
      call h5screate_simple_f(ds_rank,ds_dims,dspace_id,h5err)
      ! Create the dataset with default properties.
      call h5dcreate_f(file_id,dsetname,dtype,dspace_id,dset_id,h5err)
      call h5sclose_f(dspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
      call nvtxEndRange()

   end subroutine create_dataspace_hdf5

   subroutine create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims, chunk_dims, dtype)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ds_rank
      integer(hsize_t),intent(in) :: ds_dims(ds_rank),max_dims(ds_rank), chunk_dims(ds_rank)
      integer(hid_t),intent(in) :: dtype
      integer(hid_t) :: dset_id,dspace_id,plist_id
      integer :: h5err

      call nvtxStartRange('create_dataspace_maxdims_hdf5: CPU-DATA')
      ! Create the data space for the  dataset.
      call h5screate_simple_f(ds_rank,ds_dims,dspace_id,h5err,max_dims)
      ! Create the dataset with default properties.
      call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,h5err)
      call h5pset_chunk_f(plist_id,ds_rank,chunk_dims,h5err)
      call h5dcreate_f(file_id, dsetname,dtype,dspace_id,dset_id, h5err,plist_id)
      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(dspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
      call nvtxEndRange()

   end subroutine create_dataspace_maxdims_hdf5

   subroutine extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ds_rank
      integer(hsize_t),intent(in) :: ds_dims(ds_rank)
      integer(hid_t) :: dset_id
      integer :: h5err

      call nvtxStartRange('extend_dataset_hdf5: CPU-DATA')
      ! Open dataset
      call h5dopen_f(file_id, dsetname,dset_id,h5err)
      ! Extend the dataset to ds_dims.
      call h5dextend_f(dset_id,ds_dims,h5err)
      call h5dclose_f(dset_id,h5err)
      call nvtxEndRange()
      
   end subroutine extend_dataset_hdf5

   subroutine create_chunked_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,chunk_dims,dtype)
      !BE CAREFUL: THIS SUBROUTINE MUST BE DOUBLECHECKED
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ds_rank
      integer(hsize_t),intent(in) :: ds_dims(ds_rank),chunk_dims(ds_rank)
      integer(hid_t),intent(in) :: dtype
      integer(hid_t) :: dset_id,dspace_id,plist_id
      integer(4) :: h5err

      call nvtxStartRange('create_chunked_dataspace_hdf5: CPU-DATA')
      ! Create the data space for the  dataset.
      call h5screate_simple_f(ds_rank,ds_dims,dspace_id,h5err)
      call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,h5err)
      call h5pset_chunk_f(plist_id,ds_rank,chunk_dims,h5err)
      call h5dcreate_f(file_id, dsetname,dtype,dspace_id,dset_id, h5err,plist_id)
      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(dspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
      call nvtxEndRange()

   end subroutine create_chunked_dataspace_hdf5

!-------------------------------------------------------------------------------------------------------------------
   subroutine open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
      dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset
      integer(hid_t),intent(out) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank),intent(out) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('open_create_dataspace_hyperslab_parallel: CPU-DATA')
      call h5dopen_f(file_id, dsetname, dset_id, h5err)
      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)
      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)
      ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err)
      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)
      call nvtxEndRange()

   end subroutine open_create_dataspace_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
   subroutine close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      implicit none
      integer(hid_t),intent(inout) :: dset_id,fspace_id,mspace_id,plist_id
      integer(4) :: h5err

      call nvtxStartRange('close_dataspace_hyperslab_parallel: CPU-DATA')
      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
      call nvtxEndRange()

   end subroutine close_dataspace_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!  FP(rp 4/8)

   subroutine select_dtype_rp(dtype)
      implicit none
      integer(hid_t),intent(inout) :: dtype

      call nvtxStartRange('select_dtype_rp: CPU-DATA')
      if(rp.eq.4) then
         dtype = h5_datatype_real4
      else if(rp.eq.8) then
         dtype = h5_datatype_real8
      else
         write(*,*) 'Fatal error in select_dtype_rp! rp is not 4 or 8 >> CRASH!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call nvtxEndRange()

   end subroutine select_dtype_rp

   subroutine select_dtype_rp_vtk(dtype)
      implicit none
      integer(hid_t),intent(inout) :: dtype

      call nvtxStartRange('select_dtype_rp_vtk: CPU-DATA')
      if(rp_vtk.eq.4) then
         dtype = h5_datatype_real4
      else if(rp_vtk.eq.8) then
         dtype = h5_datatype_real8
      else
         write(*,*) 'Fatal error in select_dtype_rp_vtk! rp_vtk is not 4 or 8 >> CRASH!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call nvtxEndRange()

   end subroutine select_dtype_rp_vtk

   subroutine write_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('write_dataspace_1d_real_rp_hyperslab_parallel: select_dtype_rp')
      call select_dtype_rp(dtype)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_real_rp_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_real_rp_hyperslab_parallel: CPU-DATA')
      call h5dwrite_f(dset_id,dtype,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_real_rp_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine write_dataspace_1d_real_rp_hyperslab_parallel

   subroutine read_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('read_dataspace_1d_real_rp_hyperslab_parallel: select_dtype_rp')
      call select_dtype_rp(dtype)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_real_rp_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_real_rp_hyperslab_parallel: CPU-DATA')
      call h5dread_f(dset_id,dtype,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_real_rp_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine read_dataspace_1d_real_rp_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------------------
   subroutine write_dataspace_2d_tr_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array2d)
      implicit none
      integer(4),parameter :: ms_rank = 2
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp),intent(in) :: array2d(ms_dims(2),ms_dims(1)) !FORTRAN is COLUMN-MAJOR & HDF5 is ROW-MAJOR
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: ii,jj,h5err
      real(rp) :: array2d_tr(ms_dims(1),ms_dims(2)) !FORTRAN is COLUMN-MAJOR & HDF5 is ROW-MAJOR

      call nvtxStartRange('write_dataspace_2d_tr_real_rp_hyperslab_parallel: CPU-DATA')
      !For the moment, this loop is done with data in the host (before saving in hdf5)
      !!!!$acc kernels
      do ii=1,ms_dims(1)
         do jj=1,ms_dims(2)
            array2d_tr(ii,jj)=array2d(jj,ii)
         end do
      end do
      !!!!$acc end kernels
      !!!!$acc update device(host(:,:))
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_2d_tr_real_rp_hyperslab_parallel: select_dtype_rp')
      call select_dtype_rp(dtype)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_2d_tr_real_rp_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_2d_tr_real_rp_hyperslab_parallel: CPU-DATA')
      call h5dwrite_f(dset_id,dtype,array2d_tr,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_2d_tr_real_rp_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine write_dataspace_2d_tr_real_rp_hyperslab_parallel

   subroutine read_dataspace_2d_tr_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array2d)
      implicit none
      integer(4),parameter :: ms_rank = 2
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp),intent(out) :: array2d(ms_dims(2),ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: ii,jj,h5err
      real(rp) :: array2d_tr(ms_dims(1),ms_dims(2)) !FORTRAN is COLUMN-MAJOR & HDF5 is ROW-MAJOR

      call nvtxStartRange('read_dataspace_2d_tr_real_rp_hyperslab_parallel: select_dtype_rp')
      call select_dtype_rp(dtype)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_2d_tr_real_rp_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_2d_tr_real_rp_hyperslab_parallel: CPU-DATA')
      call h5dread_f(dset_id,dtype,array2d_tr,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_2d_tr_real_rp_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_2d_tr_real_rp_hyperslab_parallel: CPU-DATA')
      !For the moment, this loop is done with data in the host
      !!!!$acc update device(array2d(:,:)) !put in device memory the value
      !!!!$acc kernels
      do ii=1,ms_dims(1)
         do jj=1,ms_dims(2)
            array2d(jj,ii) = array2d_tr(ii,jj)
         end do
      end do
      !!!!$acc end kernels
      !!!!$acc update host(array2d(:,:)) ! put in device memory the value
      call nvtxEndRange()

   end subroutine read_dataspace_2d_tr_real_rp_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------------------------
   subroutine write_dataspace_1d_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp_vtk),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('write_dataspace_1d_real_rp_vtk_hyperslab_parallel: select_dtype_rp_vtk')
      call select_dtype_rp_vtk(dtype)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_real_rp_vtk_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_real_rp_vtk_hyperslab_parallel: CPU-DATA')
      call h5dwrite_f(dset_id,dtype,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_real_rp_vtk_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine write_dataspace_1d_real_rp_vtk_hyperslab_parallel

   subroutine read_dataspace_1d_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter:: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp_vtk),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('read_dataspace_1d_real_rp_vtk_hyperslab_parallel: select_dtype_rp_vtk')
      call select_dtype_rp_vtk(dtype)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_real_rp_vtk_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_real_rp_vtk_hyperslab_parallel: CPU-DATA')
      call h5dread_f(dset_id,dtype,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_real_rp_vtk_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine read_dataspace_1d_real_rp_vtk_hyperslab_parallel
!--------------------------------------------------------------------------------------------------------------------------------------
   subroutine write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array2d)
      implicit none
      integer(4),parameter :: ms_rank = 2
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp_vtk),intent(in) :: array2d(ms_dims(2),ms_dims(1)) !FORTRAN is COLUMN-MAJOR & HDF5 is ROW-MAJOR
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: ii,jj,h5err
      real(rp_vtk) :: array2d_tr(ms_dims(1),ms_dims(2)) !FORTRAN is COLUMN-MAJOR & HDF5 is ROW-MAJOR

      call nvtxStartRange('write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel: CPU-DATA')
      !For the moment, this loop is done with data in the host
      !!!!!$acc kernels !!BE CAREFUL HtoD/DtoH requrired if to be used!!!
      do ii=1,ms_dims(1)
         do jj=1,ms_dims(2)
            array2d_tr(ii,jj)=array2d(jj,ii)
         end do
      end do
      !!!!!$acc end kernels
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel: select_dtype_rp_vtk')
      call select_dtype_rp_vtk(dtype)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel: CPU-DATA')
      call h5dwrite_f(dset_id,dtype,array2d_tr,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel

   subroutine read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array2d)
      implicit none
      integer(4),parameter :: ms_rank = 2
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp_vtk),intent(out) :: array2d(ms_dims(2),ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: ii,jj,h5err
      real(rp_vtk) :: array2d_tr(ms_dims(1),ms_dims(2)) !FORTRAN is COLUMN-MAJOR & HDF5 is ROW-MAJOR

      call nvtxStartRange('read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel: select_dtype_rp_vtk')
      call select_dtype_rp_vtk(dtype)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel: CPU-DATA')
      call h5dread_f(dset_id,dtype,array2d_tr,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel: CPU-DATA')
      !For the moment, this loop is done with data in the host
      !!!!$acc kernels !!BE CAREFUL HtoD/DtoH requrired if to be used!!!
      do ii=1,ms_dims(1)
         do jj=1,ms_dims(2)
            array2d(jj,ii) = array2d_tr(ii,jj)
         end do
      end do
      !!!$acc end kernels
      call nvtxEndRange()
   end subroutine read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel

!-------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------------------------
!  FP32
   subroutine write_dataspace_1d_fp32_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(4),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('write_dataspace_1d_fp32_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_fp32_hyperslab_parallel: CPU-DATA')
      call h5dwrite_f(dset_id,h5_datatype_real4,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_fp32_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine write_dataspace_1d_fp32_hyperslab_parallel

   subroutine read_dataspace_1d_fp32_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(4),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(4) :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims

      call nvtxStartRange('read_dataspace_1d_fp32_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_fp32_hyperslab_parallel: CPU-DATA')
      call h5dread_f(dset_id,h5_datatype_real4,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_fp32_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine read_dataspace_1d_fp32_hyperslab_parallel

   subroutine read_dataspace_2d_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset
      real(4),intent(out) :: data(ms_dims(1), ms_dims(2))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(hid_t) :: dtype
      dtype = H5T_NATIVE_REAL

      call nvtxStartRange('read_dataspace_2d_fp32_hyperslab_parallel: CPU-DATA')
      call h5dopen_f(file_id, dsetname, dset_id, h5err)
      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)
      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)
      ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err)
      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)
      call h5dread_f(dset_id,dtype,data,fs_dims,h5err,&
         file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
      call nvtxEndRange()

   end subroutine read_dataspace_2d_fp32_hyperslab_parallel

!  FP64
   subroutine write_dataspace_1d_fp64_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(8),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('write_dataspace_1d_fp64_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_fp64_hyperslab_parallel: CPU-DATA')
      call h5dwrite_f(dset_id,h5_datatype_real8,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_fp64_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine write_dataspace_1d_fp64_hyperslab_parallel

   subroutine read_dataspace_1d_fp64_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(8),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(4) :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims

      call nvtxStartRange('read_dataspace_1d_fp64_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_fp64_hyperslab_parallel: CPU-DATA')
      call h5dread_f(dset_id,h5_datatype_real8,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_fp64_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine read_dataspace_1d_fp64_hyperslab_parallel

!  INT1
   subroutine write_dataspace_1d_int1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter:: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(1),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('write_dataspace_1d_int1_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_int1_hyperslab_parallel: CPU-DATA')
      call h5dwrite_f(dset_id,h5_datatype_int1,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_int1_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine write_dataspace_1d_int1_hyperslab_parallel

   subroutine read_dataspace_1d_int1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(1),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('read_dataspace_1d_int1_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_int1_hyperslab_parallel: CPU-DATA')
      call h5dread_f(dset_id,h5_datatype_int1,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_int1_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine read_dataspace_1d_int1_hyperslab_parallel

!  UINT1
   subroutine write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(1),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('write_dataspace_1d_uint1_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_uint1_hyperslab_parallel: CPU-DATA')
      call h5dwrite_f(dset_id,h5_datatype_uint1,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_uint1_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine write_dataspace_1d_uint1_hyperslab_parallel

   subroutine read_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(1),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('read_dataspace_1d_uint1_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_uint1_hyperslab_parallel: CPU-DATA')
      call h5dread_f(dset_id,h5_datatype_uint1,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_uint1_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine read_dataspace_1d_uint1_hyperslab_parallel

!  INT4
   subroutine write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(4),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('write_dataspace_1d_int4_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_int4_hyperslab_parallel: CPU-DATA')
      call h5dwrite_f(dset_id,h5_datatype_int4,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_int4_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine write_dataspace_1d_int4_hyperslab_parallel

   subroutine read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(4),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('read_dataspace_1d_int4_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_int4_hyperslab_parallel: CPU-DATA')
      call h5dread_f(dset_id,h5_datatype_int4,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_int4_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine read_dataspace_1d_int4_hyperslab_parallel

!  INT8
   subroutine write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(8),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('write_dataspace_1d_int8_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_int8_hyperslab_parallel: CPU-DATA')
      call h5dwrite_f(dset_id,h5_datatype_int8,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('write_dataspace_1d_int8_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine write_dataspace_1d_int8_hyperslab_parallel


   subroutine read_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(8),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call nvtxStartRange('read_dataspace_1d_int8_hyperslab_parallel: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_int8_hyperslab_parallel: CPU-DATA')
      call h5dread_f(dset_id,h5_datatype_int8,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_dataspace_1d_int8_hyperslab_parallel: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

   end subroutine read_dataspace_1d_int8_hyperslab_parallel

   subroutine load_connectivity_hdf5(file_id,mnnode)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(4),intent(in) :: mnnode
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(4) :: iElemL,i,m
      integer(hssize_t), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      call nvtxStartRange('load_connectivity_hdf5: CPU-COMPUTE')
      ms_dims(1) = int(numElemsRankPar,hsize_t)*int(mnnode,hsize_t)
      ms_offset(1) = int((rankElemStart-1),hssize_t)*int(mnnode,hssize_t)
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: CPU-DATA')
      allocate( connecParOrig(numElemsRankPar,mnnode) )
      allocate( connecParWork(numElemsRankPar,mnnode) )
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: GPU-DATA')
      !$acc enter data create(connecParOrig(:,:))
      !$acc enter data create(connecParWork(:,:))
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: CPU-DATA')
      allocate(aux_array(numElemsRankPar*mnnode))
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Connectivity/connecParOrig'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: CPU-DATA')
      i=1
      do iElemL=1,numElemsRankPar
         do m=1,mnnode
            connecParOrig(iElemL,m)=aux_array(i) !it contains the iNodeL
            i=i+1
         end do
      end do
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: GPU-DATA')
      !$acc update device(connecParOrig(:,:))
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Connectivity/connecParWork'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: CPU-DATA')
      i=1
      do iElemL=1,numElemsRankPar
         do m=1,mnnode
            connecParWork(iElemL,m)=aux_array(i) !it contains the iNodeL
            i=i+1
         end do
      end do
      deallocate(aux_array)
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: GPU-DATA')
      !$acc update device(connecParWork(:,:))
      call nvtxEndRange()
      
      call nvtxStartRange('load_connectivity_hdf5: CPU-DATA')
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      allocate(aux_array(1))
      call nvtxEndRange()
      
      call nvtxStartRange('load_connectivity_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      numWorkingNodesRankPar=aux_array(1)
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: CPU-DATA')
      deallocate(aux_array)
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: CPU-DATA')
      allocate(workingNodesPar(numWorkingNodesRankPar))
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: GPU-DATA')
      !$acc enter data create(workingNodesPar(:))
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: CPU-DATA')
      allocate(aux_array(mpi_size))
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      !read data set numWorkingNodesRankPar of all ranks
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: CPU-COMPUTE')
      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numWorkingNodesRankPar,hsize_t)
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Connectivity/workingNodesPar'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,workingNodesPar)
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: CPU-DATA')
      deallocate(aux_array)
      call nvtxEndRange()

      call nvtxStartRange('load_connectivity_hdf5: GPU-DATA')
      !$acc update device(workingNodesPar(:))
      call nvtxEndRange()

   end subroutine load_connectivity_hdf5

   subroutine load_parallel_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(4) :: i,h5err
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      !write(*,*) 'Loading parallel data hdf5...'

      call nvtxStartRange('load_parallel_data_hdf5: CPU-DATA')
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      allocate(aux_array(1))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data/rankNodeStart'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      rankNodeStart=aux_array(1)

      dsetname = '/Parallel_data/rankNodeEnd'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      rankNodeEnd=aux_array(1)

      dsetname = '/Parallel_data/rankElemStart'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      rankElemStart=aux_array(1)

      dsetname = '/Parallel_data/rankElemEnd'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      rankElemEnd=aux_array(1)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: CPU-COMPUTE')
      numNodesRankPar = rankNodeEnd - rankNodeStart + 1
      numElemsRankPar  = rankElemEnd - rankElemStart + 1
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data/numRanksWithComms'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      numRanksWithComms=aux_array(1)

      dsetname = '/Parallel_data/numNodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      numNodesToComm=aux_array(1)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: CPU-DATA')
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      !read data set numRanksWithComms of all ranks
      dsetname = '/Parallel_data/numRanksWithComms'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: CPU-COMPUTE')
      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numRanksWithComms,hsize_t)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: CPU-DATA')
      allocate(ranksToComm(numRanksWithComms))
      allocate(commsMemPosInLoc(numRanksWithComms))
      allocate(commsMemPosInNgb(numRanksWithComms))
      allocate(commsMemSize(numRanksWithComms))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: GPU-DATA')
      !$acc enter data create(ranksToComm(:))
      !$acc enter data create(commsMemPosInLoc(:))
      !$acc enter data create(commsMemSize(:))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data/ranksToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,ranksToComm)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: GPU-DATA')
      !$acc update device(ranksToComm(:))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data/commsMemPosInLoc'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,commsMemPosInLoc)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: GPU-DATA')
      !$acc update device(commsMemPosInLoc(:))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data/commsMemSize'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,commsMemSize)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: GPU-DATA')
      !$acc update device(commsMemSize(:))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data/commsMemPosInNgb'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,commsMemPosInNgb)

      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      dsetname = '/Parallel_data/numNodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      numNodesToComm=aux_array(1)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: CPU-COMPUTE')
      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numNodesToComm,hsize_t)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: CPU-DATA')
      allocate(nodesToComm(numNodesToComm))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: GPU-DATA')
      !$acc enter data create(nodesToComm(:))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data/nodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,nodesToComm)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: GPU-DATA')
      !$acc update device(nodesToComm(:))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_hdf5: CPU-DATA')
      deallocate(aux_array)
      call nvtxEndRange()

   end subroutine load_parallel_data_hdf5

   subroutine load_parallel_data_boundary_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(4) :: i,h5err
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      
      call nvtxStartRange('load_parallel_data_boundary_hdf5: CPU-DATA')
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      allocate(aux_array(1))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      bnd_numRanksWithComms=aux_array(1)

      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      bnd_numNodesToComm=aux_array(1)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: CPU-DATA')
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      !read data set numRanksWithComms of all ranks
      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: CPU-COMPUTE')
      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(bnd_numRanksWithComms,hsize_t)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: CPU-DATA')
      allocate(bnd_ranksToComm(bnd_numRanksWithComms))
      allocate(bnd_commsMemPosInLoc(bnd_numRanksWithComms))
      allocate(bnd_commsMemPosInNgb(bnd_numRanksWithComms))
      allocate(bnd_commsMemSize(bnd_numRanksWithComms))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: GPU-DATA')
      !$acc enter data create(bnd_ranksToComm(:))
      !$acc enter data create(bnd_commsMemPosInLoc(:))
      !$acc enter data create(bnd_commsMemSize(:))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data_boundary/ranksToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bnd_ranksToComm)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: GPU-DATA')
      !$acc update device(bnd_ranksToComm(:))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bnd_commsMemPosInLoc)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: GPU-DATA')
      !$acc update device(bnd_commsMemPosInLoc(:))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data_boundary/commsMemSize'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bnd_commsMemSize)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: GPU-DATA')
      !$acc update device(bnd_commsMemSize(:))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bnd_commsMemPosInNgb)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: CPU-DATA')
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      bnd_numNodesToComm=aux_array(1)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: CPU-COMPUTE')
      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(bnd_numNodesToComm,hsize_t)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: CPU-DATA')
      allocate(bnd_nodesToComm(bnd_numNodesToComm))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: GPU-DATA')
      !$acc enter data create(bnd_nodesToComm(:))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/Parallel_data_boundary/nodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bnd_nodesToComm)
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: GPU-DATA')
      !$acc update device(bnd_nodesToComm(:))
      call nvtxEndRange()

      call nvtxStartRange('load_parallel_data_boundary_hdf5: CPU-DATA')
      deallocate(aux_array)
      call nvtxEndRange()

   end subroutine load_parallel_data_boundary_hdf5

   subroutine load_periodic_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer(4) :: ds_rank,h5err
      integer(4) :: i,accumVal,iBound,m
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)
      logical :: isPeriodicFolder

      call nvtxStartRange('load_periodic_data_hdf5: CPU-DATA')
      groupname = trim('/Periodic_data')
      call h5lexists_f(file_id,groupname,isPeriodicFolder,h5err)
      call nvtxEndRange()

      if(mpi_rank.eq.0) write(*,*) 'Loading Periodic data hdf5. -> isPeriodic:',isPeriodicFolder

      if(isPeriodicFolder) then
         call nvtxStartRange('load_periodic_data_hdf5: CPU-DATA')
         isMeshPeriodic = .true.
         dtype = h5_datatype_int4
         ds_rank = 1
         ds_dims(1) = mpi_size
         ms_dims(1) = 1
         ms_offset(1) = int(mpi_rank,hssize_t)
         allocate(aux_array(1))
         call nvtxEndRange()

         call nvtxStartRange('load_periodic_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Periodic_data/nPerRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         call nvtxEndRange()

         call nvtxStartRange('load_periodic_data_hdf5: CPU-DATA')
         nPerRankPar=aux_array(1)
         allocate(masSlaRankPar(nPerRankPar,2))
         deallocate(aux_array)
         allocate(aux_array(mpi_size))
         ms_dims(1) = int(mpi_size,hsize_t)
         ms_offset(1) = 0
         call nvtxEndRange()

         call nvtxStartRange('load_periodic_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
         !read data set numBoundsRankPar of all ranks
         dsetname = '/Periodic_data/nPerRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         call nvtxEndRange()

         call nvtxStartRange('load_periodic_data_hdf5: CPU-COMPUTE')
         ds_dims(1)=0
         do i=1,mpi_size
            ds_dims(1)=ds_dims(1)+aux_array(i)
         end do

         ms_offset(1)=0
         do i=1,(mpi_rank) !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
         end do
         ms_dims(1)=int(nPerRankPar,hsize_t)
         call nvtxEndRange()

         call nvtxStartRange('load_periodic_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Periodic_data/masSlaRankPar1'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,masSlaRankPar(:,1))

         dsetname = '/Periodic_data/masSlaRankPar2'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,masSlaRankPar(:,2))
         call nvtxEndRange()

         call nvtxStartRange('load_periodic_data_hdf5: CPU-DATA')
         deallocate(aux_array)
         call nvtxEndRange()
      else
         isMeshPeriodic = .false.
      end if

   end subroutine load_periodic_data_hdf5

   subroutine load_boundary_data_hdf5(file_id,mnnode,mnpbou)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(4),intent(in) :: mnnode,mnpbou
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer(4) :: ds_rank,h5err
      integer(4) :: i,accumVal,iBound,m,iNodeL
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)
      logical :: isBoundaryFolder

      call nvtxStartRange('load_boundary_data_hdf5: CPU-DATA')
      groupname = trim('/Boundary_data')
      call h5lexists_f(file_id,groupname,isBoundaryFolder,h5err)
      isMeshBoundaries = isBoundaryFolder
      call nvtxEndRange()

      if(mpi_rank.eq.0) write(*,*) 'Loading Boundary data hdf5. -> isBoundaryFolder:',isBoundaryFolder

      if(isMeshBoundaries) then
         call nvtxStartRange('load_boundary_data_hdf5: CPU-DATA')
         dtype = h5_datatype_int4
         ds_rank = 1
         ds_dims(1) = mpi_size
         ms_dims(1) = 1
         ms_offset(1) = 0
         allocate(aux_array(1))
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Boundary_data/numBoundCodes'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         numBoundCodes=aux_array(1)

         ms_offset(1) = int(mpi_rank,hssize_t)

         dsetname = '/Boundary_data/numBoundsRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         numBoundsRankPar=aux_array(1)

         dsetname = '/Boundary_data/ndofRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         ndofRankPar=aux_array(1)

         dsetname = '/Boundary_data/numBoundaryNodesRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         numBoundaryNodesRankPar=aux_array(1)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: CPU-DATA')
         allocate(boundPar(numBoundsRankPar,mnpbou))
         allocate(boundParOrig(numBoundsRankPar,mnpbou))
         allocate(bouCodesPar(numBoundsRankPar))
         allocate(ldofPar(ndofRankPar))
         allocate(lbnodesPar(numBoundaryNodesRankPar))
         deallocate(aux_array)
         allocate( aux_array(mpi_size) )
         ms_dims(1) = int(mpi_size,hsize_t)
         ms_offset(1) = 0
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
         !read data set numBoundsRankPar of all ranks
         dsetname = '/Boundary_data/numBoundsRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: CPU-COMPUTE')
         !taking advantadge and setting totalNumBoundsSrl
         totalNumBoundsSrl=0
         ds_dims(1)=0
         do i=1,mpi_size
            ds_dims(1)=ds_dims(1)+aux_array(i)
            totalNumBoundsSrl=totalNumBoundsSrl+aux_array(i)
         end do

         ms_offset(1)=0
         do i=1,(mpi_rank) !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
         end do
         ms_dims(1)=int(numBoundsRankPar,hsize_t)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Boundary_data/bouCodesPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bouCodesPar)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: CPU-DATA')
         deallocate(aux_array)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: CPU-COMPUTE')
         ds_dims(1) = int(ds_dims(1),hsize_t)*int(mnpbou,hsize_t)!totalNumBoundsSrl*mnpbou
         ms_dims(1) = int(ms_dims(1),hsize_t)*int(mnpbou,hsize_t)!numBoundsRankPar*mnpbou
         ms_offset(1) = ms_offset(1)*int(mnpbou,hssize_t)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: CPU-DATA')
         allocate(aux_array(numBoundsRankPar*mnnode))
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
         !boundPar
         dsetname = '/Boundary_data/boundPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: CPU-DATA')
         i=1
         do iBound=1,numBoundsRankPar
            do m=1,mnpbou
               boundPar(iBound,m)=aux_array(i)
               i=i+1
            end do
         end do
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
         !boundParOrig
         dsetname = '/Boundary_data/boundParOrig'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: CPU-DATA')
         i=1
         do iBound=1,numBoundsRankPar
            do m=1,mnpbou
               boundParOrig(iBound,m)=aux_array(i)
               i=i+1
            end do
         end do

         deallocate(aux_array)
         allocate( aux_array(mpi_size) )
         ms_dims(1) = int(mpi_size,hsize_t)
         ms_offset(1) = 0
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
         !read data set ndofRankPar of all ranks
         dsetname = '/Boundary_data/ndofRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: CPU-COMPUTE')
         ds_dims(1)=0
         do i=1,mpi_size
            ds_dims(1)=ds_dims(1)+aux_array(i)
         end do

         ms_offset(1)=0
         do i=1,(mpi_rank) !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
         end do
         ms_dims(1)=int(ndofRankPar,hsize_t)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Boundary_data/ldofPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,ldofPar)

         ms_dims(1) = int(mpi_size,hsize_t)
         ms_offset(1) = 0
         !read data set numBoundaryNodesRankPar of all ranks
         dsetname = '/Boundary_data/numBoundaryNodesRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: CPU-COMPUTE')
         ds_dims(1)=0
         do i=1,mpi_size
            ds_dims(1)=ds_dims(1)+aux_array(i)
         end do

         ms_offset(1)=0
         do i=1,(mpi_rank) !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
         end do
         ms_dims(1)=int(numBoundaryNodesRankPar,hsize_t)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
         dsetname = '/Boundary_data/lbnodesPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,lbnodesPar)
         call nvtxEndRange()

         call nvtxStartRange('load_boundary_data_hdf5: CPU-DATA')
         deallocate(aux_array)
         call nvtxEndRange()
      else
         call nvtxStartRange('load_boundary_data_hdf5: CPU-DATA')
         numBoundsRankPar=0
         allocate(boundPar(numBoundsRankPar,mnpbou))
         allocate(boundParOrig(numBoundsRankPar,mnpbou))
         allocate(bouCodesPar(numBoundsRankPar))

         numBoundaryNodesRankPar=0
         ndofRankPar = numNodesRankPar
         allocate(ldofPar(ndofRankPar))
         allocate(lbnodesPar(numBoundaryNodesRankPar))
         do iNodeL = 1,ndofRankPar
            ldofPar(iNodeL) = iNodeL
         end do
         call nvtxEndRange()
      end if
   end subroutine load_boundary_data_hdf5

   subroutine load_coordinates_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(2) :: ms_dims2d
      integer(hssize_t), dimension(2) :: ms_offset2d

      call nvtxStartRange('load_coordinates_hdf5: CPU-DATA')
      allocate(coordPar(numNodesRankPar,ndime))
      call nvtxEndRange()

      call nvtxStartRange('load_coordinates_hdf5: GPU-DATA')
      !$acc enter data create(coordPar(:,:))
      call nvtxEndRange()

      call nvtxStartRange('load_coordinates_hdf5: read_array2D_tr_rp_in_dataset_hdf5_file')
      ms_dims2d(1) = int(ndime,hsize_t)
      ms_dims2d(2) = int(numNodesRankPar,hsize_t)
      ms_offset2d(1) = 0
      ms_offset2d(2) = int(rankNodeStart,hssize_t)-1

      if(isMeshLinealOutput) then
         dsetname = '/VTKHDF/Points'
      else
         dsetname = '/Coords/Points'
      end if
      call read_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims2d,ms_offset2d,coordPar)
      call nvtxEndRange()

      call nvtxStartRange('load_coordinates_hdf5: GPU-DATA')
      !$acc update device(coordPar(:,:))
      call nvtxEndRange()

   end subroutine load_coordinates_hdf5

   subroutine load_globalIds_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4) :: iNodeL
      integer(8) :: iNodeGSrl,max_iNodeGSrl_l,max_iNodeGSrl_g

      call nvtxStartRange('load_globalIds_hdf5: CPU-DATA')
      allocate(globalIdSrl(numNodesRankPar))
      allocate(globalIdPar(numNodesRankPar))
      allocate(elemGid(numElemsRankPar))
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1
      call nvtxEndRange()

      call nvtxStartRange('load_globalIds_hdf5: read_dataspace_1d_int8_hyperslab_parallel')
      dsetname = '/globalIds/globalIdSrl'
      call read_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,globalIdSrl)

      dsetname = '/globalIds/globalIdPar'
      call read_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,globalIdPar)
      call nvtxEndRange()

      call nvtxStartRange('load_globalIds_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      ms_dims(1) = int(numElemsRankPar,hsize_t)
      ms_offset(1) = int(rankElemStart,hssize_t)-1

      dsetname = '/globalIds/elemGid'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,elemGid)
      call nvtxEndRange()

      call nvtxStartRange('load_globalIds_hdf5: CPU-DATA')
      !setting totalNumNodesSrl
      max_iNodeGSrl_l=0
      do iNodeL=1,numNodesRankPar
         iNodeGSrl = globalIdSrl(iNodeL)
         max_iNodeGSrl_l = max(iNodeGSrl,max_iNodeGSrl_l)
      end do
      call nvtxEndRange()

      call nvtxStartRange('load_globalIds_hdf5: MPI_Allreduce')
      call MPI_Allreduce(max_iNodeGSrl_l,max_iNodeGSrl_g,1,mpi_datatype_int8,MPI_MAX,app_comm,mpi_err)
      call nvtxEndRange()

      totalNumNodesSrl = max_iNodeGSrl_g
      
   end subroutine load_globalIds_hdf5

   subroutine load_meshOutputInfo_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(hssize_t), dimension(1) :: ms_offset
      integer(1) :: aux_array_i1(1)
      integer(4) :: aux_array_i4(1)

      ms_dims(1) = 1
      ms_offset(1) = 0

      call nvtxStartRange('load_meshOutputInfo_hdf5: read_dataspace_1d_uint1_hyperslab_parallel')
      dsetname = '/meshOutputInfo/isLinealOutput'
      call read_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i1)
      call nvtxEndRange()

      call nvtxStartRange('load_meshOutputInfo_hdf5: CPU-DATA')
      isMeshLinealOutput = .false.
      if(aux_array_i1(1).eq.1) isMeshLinealOutput = .true.
      if(mpi_rank.eq.0) write(*,*) 'Lineal Output:',isMeshLinealOutput
      call nvtxEndRange()

      call nvtxStartRange('load_meshOutputInfo_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/meshOutputInfo/nnodeVTK'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      mesh_VTKnnode = aux_array_i4(1)

      dsetname = '/meshOutputInfo/numVTKElemsPerMshElem'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      mesh_numVTKElemsPerMshElem = aux_array_i4(1)

      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)

      dsetname = '/meshOutputInfo/numElemsVTKMshRank'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      numElemsVTKRankPar = aux_array_i4(1)

      dsetname = '/meshOutputInfo/sizeConnecVTKMshRank'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      sizeConnecVTKRankPar = aux_array_i4(1)
      call nvtxEndRange()

   end subroutine load_meshOutputInfo_hdf5

   subroutine load_porder_hdf5(file_id,mnnode,mnpbou)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(4),intent(in) :: mnnode,mnpbou
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(hssize_t), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      call nvtxStartRange('load_porder_hdf5: CPU-DATA')
      allocate(aux_array(1))
      ms_dims(1) = 1
      ms_offset(1) = 0
      call nvtxEndRange()

      call nvtxStartRange('load_porder_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/order/porder'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      call nvtxEndRange()

      call nvtxStartRange('load_porder_hdf5: CPU-DATA')
      mesh_porder = aux_array(1)
      deallocate(aux_array)

      if(mesh_porder .ne. porder) then
         write(*,*) 'FATAL ERROR! mesh_porder',mesh_porder,' different to porder',porder
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      allocate(mesh_a2ijk(mnnode))
      allocate(mesh_gmsh2ijk(mnnode))
      allocate(mesh_vtk2ijk(mnnode))
      ms_dims(1) = mnnode
      ms_offset(1) = 0
      call nvtxEndRange()

      call nvtxStartRange('load_porder_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/order/a2ijk'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,mesh_a2ijk)

      dsetname = '/order/gmsh2ijk'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,mesh_gmsh2ijk)

      dsetname = '/order/vtk2ijk'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,mesh_vtk2ijk)
      call nvtxEndRange()

      call nvtxStartRange('load_porder_hdf5: CPU-DATA')
      allocate(mesh_a2ij(mnpbou))
      allocate(mesh_gmsh2ij(mnpbou))
      allocate(mesh_vtk2ij(mnpbou))
      ms_dims(1) = mnpbou
      ms_offset(1) = 0
      call nvtxEndRange()

      call nvtxStartRange('load_porder_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname = '/order/a2ij'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,mesh_a2ij)

      dsetname = '/order/gmsh2ij'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,mesh_gmsh2ij)

      dsetname = '/order/vtk2ij'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,mesh_vtk2ij)
      call nvtxEndRange()

   end subroutine load_porder_hdf5

!---------------------------------------------------------------------------------------------------------
!           RESULTS FILE
!---------------------------------------------------------------------------------------------------------

   subroutine set_hdf5_resultsFile_name(iStep,full_fileName)
      implicit none
      integer, intent(in) :: iStep
      character(len=*), intent(out) :: full_fileName
      character(len=12) :: aux_step

      call nvtxStartRange('set_hdf5_resultsFile_name: CPU-DATA')
      write(aux_step,'(I0)') iStep
      full_fileName = trim(adjustl(base_resultsFile_h5_name))//trim(aux_step)//'.hdf'
      call nvtxEndRange()

   end subroutine set_hdf5_resultsFile_name

   subroutine set_hdf5_restartFile_name(iStep,full_fileName)
      implicit none
      integer, intent(in) :: iStep
      character(len=*), intent(out) :: full_fileName
      character(len=12) :: aux_step

      call nvtxStartRange('set_hdf5_restartFile_name: CPU-DATA')
      write(aux_step,'(I0)') iStep
      full_fileName = trim(adjustl(base_restartFile_h5_name))//trim(aux_step)//'.h5'
      call nvtxEndRange()

   end subroutine set_hdf5_restartFile_name

   subroutine set_hdf5_avgResultsFile_name(iStep,full_fileName)
      implicit none
      integer, intent(in) :: iStep
      character(len=*), intent(out) :: full_fileName
      character(len=12) :: aux_step

      call nvtxStartRange('set_hdf5_avgResultsFile_name: CPU-DATA')
      write(aux_step,'(I0)') iStep
      full_fileName = trim(adjustl(base_avgResultsFile_h5_name))//trim(aux_step)//'.hdf'
      call nvtxEndRange()

   end subroutine set_hdf5_avgResultsFile_name

   subroutine set_hdf5_surface_resultsFile_name(surf_res_fileName,res_fileName)
      implicit none
      character(len=*), intent(out) :: surf_res_fileName
      character(len=*), intent(in)  :: res_fileName

      call nvtxStartRange('set_hdf5_surface_resultsFile_name: CPU-DATA')
      surf_res_fileName = 'surface_'//trim(adjustl(res_fileName))
      call nvtxEndRange()

   end subroutine set_hdf5_surface_resultsFile_name

   subroutine create_dataspace_for_rp_vtk_hdf5(file_id,dsetname,ds_rank,ds_dims)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(4),intent(in) :: ds_rank
      integer(hsize_t),dimension(ds_rank),intent(in) :: ds_dims
      integer(hid_t) :: dtype

      call nvtxStartRange('create_dataspace_for_rp_vtk_hdf5: CPU-DATA')
      if(rp_vtk.eq.4) then
         dtype = h5_datatype_real4
      else if(rp_vtk.eq.8) then
         dtype = h5_datatype_real8
      else
         write(*,*) 'Fatal error in create_dataspace_for_rp_vtk_hdf5! rp is not 4 or 8 >> CRASH!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call nvtxEndRange()

      call nvtxStartRange('create_dataspace_for_rp_vtk_hdf5: create_dataspace_hdf5')
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

   end subroutine create_dataspace_for_rp_vtk_hdf5

   subroutine save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,data_array_rp,isCreateDataspaceOpt)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(hsize_t),dimension(1),intent(in) :: ds_dims
      integer(hsize_t),dimension(1),intent(in) :: ms_dims
      integer(hssize_t),dimension(1),intent(in) :: ms_offset
      real(rp),intent(in) :: data_array_rp(ms_dims(1))
      logical, intent(in), optional :: isCreateDataspaceOpt
      integer(4) :: ds_rank = 1 !it is forced
      logical :: isCreateDataspace
      integer(4) :: h5err
      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:)

      call nvtxStartRange('save_array1D_rp_in_dataset_hdf5_file: CPU-DATA')
      isCreateDataspace = .true.
      if (present(isCreateDataspaceOpt)) then
         isCreateDataspace = isCreateDataspaceOpt
      end if
      call nvtxEndRange()

      if (isCreateDataspace) then
         call nvtxStartRange('save_array1D_rp_in_dataset_hdf5_file: create_dataspace_for_rp_vtk_hdf5')
         call create_dataspace_for_rp_vtk_hdf5(file_id,dsetname,ds_rank,ds_dims)
         call nvtxEndRange()
      end if

      if(rp .eq. rp_vtk) then
         call nvtxStartRange('save_array1D_rp_in_dataset_hdf5_file: write_dataspace_1d_real_rp_hyperslab_parallel')
         call write_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,data_array_rp)
         call nvtxEndRange()
      else
         !copying in the aux array
         call nvtxStartRange('save_array1D_rp_in_dataset_hdf5_file: CPU-DATA')
         allocate(aux_data_array_rp_vtk(ms_dims(1)))
         aux_data_array_rp_vtk(:) = real(data_array_rp(:),rp_vtk)
         call nvtxEndRange()

         call nvtxStartRange('save_array1D_rp_in_dataset_hdf5_file: write_dataspace_1d_real_rp_vtk_hyperslab_parallel')
         call write_dataspace_1d_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_data_array_rp_vtk)
         call nvtxEndRange()

         call nvtxStartRange('save_array1D_rp_in_dataset_hdf5_file: CPU-DATA')
         deallocate(aux_data_array_rp_vtk)
         call nvtxEndRange()
      end if

   end subroutine save_array1D_rp_in_dataset_hdf5_file

   subroutine save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,data_array_rp,isCreateDataspaceOpt)
      implicit none
      integer(4),parameter :: ds_rank = 2 !it is forced
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(hsize_t),dimension(ds_rank),intent(in) :: ds_dims2d,ms_dims2d
      integer(hssize_t),dimension(ds_rank),intent(inout) :: ms_offset2d
      logical, intent(in), optional :: isCreateDataspaceOpt
      real(rp),intent(in) :: data_array_rp(ms_dims2d(2),ms_dims2d(1)) !fortran is column-major & hdf5 writes in row-major
      logical :: isCreateDataspace
      integer(4) :: h5err
      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:,:)

      call nvtxStartRange('save_array2D_tr_rp_in_dataset_hdf5_file: CPU-DATA')
      isCreateDataspace = .true.
      if (present(isCreateDataspaceOpt)) then
         isCreateDataspace = isCreateDataspaceOpt
      end if
      call nvtxEndRange()

      if (isCreateDataspace) then
         call nvtxStartRange('save_array2D_tr_rp_in_dataset_hdf5_file: create_dataspace_for_rp_vtk_hdf5')
         call create_dataspace_for_rp_vtk_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
         call nvtxEndRange()
      end if

      if(rp .eq. rp_vtk) then
         call nvtxStartRange('save_array2D_tr_rp_in_dataset_hdf5_file: write_dataspace_2d_tr_real_rp_hyperslab_parallel')
         call write_dataspace_2d_tr_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,data_array_rp)
         call nvtxEndRange()
      else
         call nvtxStartRange('save_array2D_tr_rp_in_dataset_hdf5_file: CPU-DATA')
         allocate(aux_data_array_rp_vtk(ms_dims2d(2),ms_dims2d(1)))
         aux_data_array_rp_vtk(:,:) = real(data_array_rp(:,:),rp_vtk)
         call nvtxEndRange()

         call nvtxStartRange('save_array2D_tr_rp_in_dataset_hdf5_file: write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel')
         call write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,aux_data_array_rp_vtk)
         call nvtxEndRange()

         call nvtxStartRange('save_array2D_tr_rp_in_dataset_hdf5_file: CPU-DATA')
         deallocate(aux_data_array_rp_vtk)
         call nvtxEndRange()
      end if

   end subroutine save_array2D_tr_rp_in_dataset_hdf5_file

   subroutine save_array1D_rp_vtk_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,data_array_rp_vtk,isCreateDataspaceOpt)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(hsize_t),dimension(1),intent(in) :: ds_dims
      integer(hsize_t),dimension(1),intent(in) :: ms_dims
      integer(hssize_t),dimension(1),intent(in) :: ms_offset
      real(rp_vtk),intent(in) :: data_array_rp_vtk(ms_dims(1))
      logical, intent(in), optional :: isCreateDataspaceOpt
      integer(4) :: ds_rank = 1 !it is forced
      logical :: isCreateDataspace
      integer(4) :: h5err

      call nvtxStartRange('save_array1D_rp_vtk_in_dataset_hdf5_file: CPU-DATA')
      isCreateDataspace = .true.
      if (present(isCreateDataspaceOpt)) then
         isCreateDataspace = isCreateDataspaceOpt
      end if
      call nvtxEndRange()

      if (isCreateDataspace) then
         call nvtxStartRange('save_array1D_rp_vtk_in_dataset_hdf5_file: create_dataspace_for_rp_vtk_hdf5')
         call create_dataspace_for_rp_vtk_hdf5(file_id,dsetname,ds_rank,ds_dims)
         call nvtxEndRange()
      end if

      call nvtxStartRange('save_array1D_rp_vtk_in_dataset_hdf5_file: write_dataspace_1d_real_rp_vtk_hyperslab_parallel')
      call write_dataspace_1d_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,data_array_rp_vtk)
      call nvtxEndRange()

   end subroutine save_array1D_rp_vtk_in_dataset_hdf5_file

   subroutine save_array2D_tr_rp_vtk_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,data_array_rp_vtk,isCreateDataspaceOpt)
      implicit none
      integer(4),parameter :: ds_rank = 2 !it is forced
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(hsize_t),dimension(ds_rank),intent(in) :: ds_dims2d,ms_dims2d
      integer(hssize_t),dimension(ds_rank),intent(inout) :: ms_offset2d
      logical, intent(in), optional :: isCreateDataspaceOpt
      real(rp_vtk),intent(in) :: data_array_rp_vtk(ms_dims2d(2),ms_dims2d(1)) !fortran is column-major & hdf5 writes in row-major
      logical :: isCreateDataspace
      integer(4) :: h5err

      call nvtxStartRange('save_array2D_tr_rp_vtk_in_dataset_hdf5_file: CPU-DATA')
      isCreateDataspace = .true.
      if (present(isCreateDataspaceOpt)) then
         isCreateDataspace = isCreateDataspaceOpt
      end if
      call nvtxEndRange()

      if (isCreateDataspace) then
         call nvtxStartRange('save_array2D_tr_rp_vtk_in_dataset_hdf5_file: create_dataspace_for_rp_vtk_hdf5')
         call create_dataspace_for_rp_vtk_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
         call nvtxEndRange()
      end if

      call nvtxStartRange('save_array2D_tr_rp_vtk_in_dataset_hdf5_file: write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel')
      call write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,data_array_rp_vtk)
      call nvtxEndRange()

   end subroutine save_array2D_tr_rp_vtk_in_dataset_hdf5_file
!------------------------------------------------------------------------------------------------------------------------------

   subroutine read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,data_array_rp)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(HSIZE_T),dimension(1),intent(in) :: ms_dims
      integer(HSSIZE_T),dimension(1),intent(in) :: ms_offset
      real(rp),intent(out) :: data_array_rp(ms_dims(1))

      integer(4) :: h5err
      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:)

      if(rp .eq. rp_vtk) then
         call nvtxStartRange('read_array1D_rp_in_dataset_hdf5_file: read_dataspace_1d_real_rp_hyperslab_parallel')
         call read_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,data_array_rp)
         call nvtxEndRange()
      else
         call nvtxStartRange('read_array1D_rp_in_dataset_hdf5_file: CPU-DATA')
         allocate(aux_data_array_rp_vtk(ms_dims(1)))
         call nvtxEndRange()

         call nvtxStartRange('read_array1D_rp_in_dataset_hdf5_file: read_dataspace_1d_real_rp_vtk_hyperslab_parallel')
         call read_dataspace_1d_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_data_array_rp_vtk)
         call nvtxEndRange()

         call nvtxStartRange('read_array1D_rp_in_dataset_hdf5_file: CPU-DATA')
         data_array_rp(:) = real(aux_data_array_rp_vtk(:),rp)
         deallocate(aux_data_array_rp_vtk)
         call nvtxEndRange()
      end if

   end subroutine read_array1D_rp_in_dataset_hdf5_file

   subroutine read_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims2d,ms_offset2d,data_array_rp)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(hsize_t),dimension(2),intent(in) :: ms_dims2d
      integer(hssize_t),dimension(2),intent(inout) :: ms_offset2d
      real(rp),intent(out) :: data_array_rp(ms_dims2d(2),ms_dims2d(1)) !fortran is column-major & hdf5 writes in row-major

      integer(4) :: h5err
      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:,:)

      if(rp .eq. rp_vtk) then
         call nvtxStartRange('read_array2D_tr_rp_in_dataset_hdf5_file: read_dataspace_2d_tr_real_rp_hyperslab_parallel')
         call read_dataspace_2d_tr_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,data_array_rp)
         call nvtxEndRange()
      else
         call nvtxStartRange('read_array2D_tr_rp_in_dataset_hdf5_file: CPU-DATA')
         allocate(aux_data_array_rp_vtk(ms_dims2d(2),ms_dims2d(1)))
         call nvtxEndRange()

         call nvtxStartRange('read_array2D_tr_rp_in_dataset_hdf5_file: read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel')
         call read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,aux_data_array_rp_vtk)
         call nvtxEndRange()

         call nvtxStartRange('read_array2D_tr_rp_in_dataset_hdf5_file: CPU-DATA')
         data_array_rp(:,:) = real(aux_data_array_rp_vtk(:,:),rp)
         deallocate(aux_data_array_rp_vtk)
         call nvtxEndRange()
      end if

   end subroutine read_array2D_tr_rp_in_dataset_hdf5_file

   subroutine save_int4_in_dataset_hdf5_file(file_id,dsetname,int2save)
      implicit none
      integer(4),parameter :: ds_rank = 1, ms_rank = 1 !it is forced
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(4),intent(in) :: int2save
      integer(hsize_t),dimension(ms_rank) :: ds_dims,ms_dims,fs_dims,fs_maxdims
      integer(hssize_t),dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(4) :: h5err
      integer(4),allocatable :: aux_data_array_int4(:)

      call nvtxStartRange('save_int4_in_dataset_hdf5_file: CPU-DATA')
      ms_dims(1) = 0
      ds_dims = 1
      ms_offset(1) = 0
      if(mpi_rank.eq.0) then
         ms_dims(1) = 1
      endif

      allocate(aux_data_array_int4(ms_dims(1)))

      if(mpi_rank.eq.0) then
         aux_data_array_int4(1) = int2save
      end if
      call nvtxEndRange()

      call nvtxStartRange('save_int4_in_dataset_hdf5_file: create_dataspace_hdf5')
      dtype = h5_datatype_int4
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('save_int4_in_dataset_hdf5_file: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('save_int4_in_dataset_hdf5_file: CPU-DATA')
      call h5dwrite_f(dset_id,dtype,aux_data_array_int4,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('save_int4_in_dataset_hdf5_file: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

      call nvtxStartRange('save_int4_in_dataset_hdf5_file: CPU-DATA')
      deallocate(aux_data_array_int4)
      call nvtxEndRange()

   end subroutine save_int4_in_dataset_hdf5_file

   subroutine read_int4_in_dataset_hdf5_file(file_id,dsetname,int2read)
      implicit none
      integer(4),parameter :: ds_rank = 1, ms_rank = 1 !it is forced
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(4),intent(out) :: int2read
      integer(hsize_t),dimension(ms_rank) :: ds_dims,ms_dims,fs_dims,fs_maxdims
      integer(hssize_t),dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(4) :: h5err
      integer(4),allocatable :: aux_data_array_int4(:)

      call nvtxStartRange('read_int4_in_dataset_hdf5_file: CPU-DATA')
      ms_dims(1) = 1
      ds_dims = 1
      ms_offset(1) = 0
      allocate(aux_data_array_int4(ms_dims(1)))
      dtype = h5_datatype_int4
      call nvtxEndRange()

      call nvtxStartRange('read_int4_in_dataset_hdf5_file: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('read_int4_in_dataset_hdf5_file: CPU-DATA')
      call h5dread_f(dset_id,dtype,aux_data_array_int4,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_int4_in_dataset_hdf5_file: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_int4_in_dataset_hdf5_file: CPU-DATA')
      int2read = aux_data_array_int4(1)
      deallocate(aux_data_array_int4)
      call nvtxEndRange()

   end subroutine read_int4_in_dataset_hdf5_file

   subroutine save_real_rp_in_dataset_hdf5_file(file_id,dsetname,real2save)
      implicit none
      integer(4),parameter :: ds_rank = 1, ms_rank = 1 !it is forced
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      real(rp),intent(in) :: real2save
      integer(hsize_t),dimension(ms_rank) :: ds_dims,ms_dims,fs_dims,fs_maxdims
      integer(hssize_t),dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(4) :: h5err

      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:)

      call nvtxStartRange('save_real_rp_in_dataset_hdf5_file: CPU-DATA')
      ms_dims(1) = 0
      ds_dims = 1
      ms_offset(1) = 0
      if(mpi_rank.eq.0) then
         ms_dims(1) = 1
      endif

      allocate(aux_data_array_rp_vtk(ms_dims(1)))

      if(mpi_rank.eq.0) then
         aux_data_array_rp_vtk(1) = real(real2save,rp_vtk)
      end if
      call nvtxEndRange()

      call nvtxStartRange('save_real_rp_in_dataset_hdf5_file: create_dataspace_for_rp_vtk_hdf5')
      call create_dataspace_for_rp_vtk_hdf5(file_id,dsetname,ds_rank,ds_dims)
      call nvtxEndRange()

      call nvtxStartRange('save_real_rp_in_dataset_hdf5_file: select_dtype_rp_vtk')
      call select_dtype_rp_vtk(dtype)
      call nvtxEndRange()

      call nvtxStartRange('save_real_rp_in_dataset_hdf5_file: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('save_real_rp_in_dataset_hdf5_file: CPU-DATA')
      call h5dwrite_f(dset_id,dtype,aux_data_array_rp_vtk,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('save_real_rp_in_dataset_hdf5_file: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

      call nvtxStartRange('save_real_rp_in_dataset_hdf5_file: CPU-DATA')
      deallocate(aux_data_array_rp_vtk)
      call nvtxEndRange()

   end subroutine save_real_rp_in_dataset_hdf5_file

   subroutine read_real_rp_in_dataset_hdf5_file(file_id,dsetname,real2read)
      implicit none
      integer(4),parameter :: ds_rank = 1, ms_rank = 1 !it is forced
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      real(rp),intent(out) :: real2read
      integer(hsize_t),dimension(ms_rank) :: ds_dims,ms_dims,fs_dims,fs_maxdims
      integer(hssize_t),dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(4) :: h5err

      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:)

      call nvtxStartRange('read_real_rp_in_dataset_hdf5_file: CPU-DATA')
      ms_dims(1) = 1
      ds_dims = 1
      ms_offset(1) = 0

      allocate(aux_data_array_rp_vtk(ms_dims(1)))
      call nvtxEndRange()

      call nvtxStartRange('read_real_rp_in_dataset_hdf5_file: select_dtype_rp_vtk')
      call select_dtype_rp_vtk(dtype)
      call nvtxEndRange()

      call nvtxStartRange('read_real_rp_in_dataset_hdf5_file: open_create_dataspace_hyperslab_parallel')
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call nvtxEndRange()

      call nvtxStartRange('read_real_rp_in_dataset_hdf5_file: CPU-DATA')
      call h5dread_f(dset_id,dtype,aux_data_array_rp_vtk,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_real_rp_in_dataset_hdf5_file: close_dataspace_hyperslab_parallel')
      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      call nvtxEndRange()

      call nvtxStartRange('read_real_rp_in_dataset_hdf5_file: CPU-DATA')
      real2read = real(aux_data_array_rp_vtk(1),rp)
      deallocate(aux_data_array_rp_vtk)
      call nvtxEndRange()

   end subroutine read_real_rp_in_dataset_hdf5_file

!----------------------------------------------------------------------------------------------------------------------------------

   subroutine save_hdf5_restartFile(mnnode,mngaus,restartCnt,iStep,flag_walave_value,time,rho,u,pr,E,mu_e,mu_t,walave_u)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus
      integer(4),intent(in) :: restartCnt,iStep
      logical,intent(in) :: flag_walave_value
      real(rp),intent(in) :: time
      real(rp),intent(inout),dimension(numNodesRankPar)       :: rho,pr,E
      real(rp),intent(inout),dimension(numNodesRankPar,ndime) :: u,walave_u
      real(rp),intent(inout),dimension(numElemsRankPar,mngaus) :: mu_e,mu_t

      integer(hid_t) :: file_id,plist_id,dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4) :: ds_rank,ms_rank,h5err
      character(512) :: full_fileName,dsetname

      integer(4) :: aux_array_i4(1)

      real(rp),dimension(numNodesRankPar) :: aux_mu_e,aux_mu_t
      integer(4) :: iElem,iNode

      call nvtxStartRange('save_hdf5_restartFile: GPU-DATA')
      !$acc kernels
      do iElem = 1,numElemsRankPar
         do iNode = 1,mnnode
            aux_mu_e(connecParOrig(iElem,iNode)) = mu_e(iElem,iNode)
            aux_mu_t(connecParOrig(iElem,iNode)) = mu_t(iElem,iNode)
         end do
      end do
      !$acc end kernels

      !$acc update host(rho(:))
      !$acc update host(u(:,:))
      !$acc update host(pr(:))
      !$acc update host(E(:))
      if(flag_walave_value) then
         !$acc update host(walave_u(:,:))
      end if
      call nvtxEndRange()

      ! Writing HDF5 Files
      call nvtxStartRange('save_hdf5_restartFile: set_hdf5_restartFile_name')
      call set_hdf5_restartFile_name(restartCnt,full_fileName)
      call nvtxEndRange()

      call nvtxStartRange('save_hdf5_restartFile: create_hdf5_file')
      call create_hdf5_file(full_fileName,file_id)
      call nvtxEndRange()

      ds_dims(1) = int(totalNumNodesPar,hsize_t)
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1
      
      call nvtxStartRange('save_hdf5_restartFile: save_array1D_rp_in_dataset_hdf5_file')
      dsetname = 'rho'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,rho)

      dsetname = 'u_x'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u(:,1))

      dsetname = 'u_y'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u(:,2))

      dsetname = 'u_z'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u(:,3))

      dsetname = 'pr'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,pr)

      dsetname = 'E'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,E)

      dsetname = 'mue'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,aux_mu_e)

      dsetname = 'mut'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,aux_mu_t)

      if(flag_walave_value) then
         dsetname = 'walave_u_x'
         call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,walave_u(:,1))

         dsetname = 'walave_u_y'
         call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,walave_u(:,2))

         dsetname = 'walave_u_z'
         call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,walave_u(:,3))
      end if
      call nvtxEndRange()

      call nvtxStartRange('save_hdf5_restartFile: save_real_rp_in_dataset_hdf5_file')
      ! ----  time  -----
      dsetname = 'time'
      call save_real_rp_in_dataset_hdf5_file(file_id,dsetname,time)
      call nvtxEndRange()

      call nvtxStartRange('save_hdf5_restartFile: save_int4_in_dataset_hdf5_file')
      ! ---- istep -----
      dsetname = 'istep'
      call save_int4_in_dataset_hdf5_file(file_id,dsetname,iStep)
      call nvtxEndRange()

      call nvtxStartRange('save_hdf5_restartFile: close_hdf5_file')
      !close the file.
      call close_hdf5_file(file_id)
      call nvtxEndRange()

   end subroutine save_hdf5_restartFile

   subroutine load_hdf5_restartFile(mnnode,mngaus,restartCnt,load_step,flag_walave_value,time,rho,u,pr,E,mu_e,mu_t,walave_u)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,restartCnt
      logical,intent(in) :: flag_walave_value
      integer(4),intent(inout) :: load_step
      real(rp),intent(inout) :: time
      real(rp),intent(inout),dimension(numNodesRankPar)       :: rho,pr,E
      real(rp),intent(inout),dimension(numNodesRankPar,ndime) :: u,walave_u
      real(rp),intent(inout),dimension(numElemsRankPar,mngaus) :: mu_e,mu_t

      character(512) :: full_restartFileName
      integer(hid_t) :: file_id,plist_id
      integer(hsize_t),dimension(1) :: ms_dims
      integer(hssize_t),dimension(1) :: ms_offset
      integer(4) :: ms_rank,iPer,h5err
      character(128) :: dsetname

      real(rp) :: aux_array_rp(1)
      integer(4) :: aux_array_i4(1)

      real(rp),dimension(numNodesRankPar) :: aux_mu_e,aux_mu_t
      integer(4) :: iElem,iNode
      logical    :: link_exists

      if((restartCnt.ne.1).and.(restartCnt.ne.2)) then
         write(*,*) 'FATAL ERROR in load_hdf5_restartFile! restartFile to load must be 1 or 2 and is',restartCnt,'CRASHING!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      call nvtxStartRange('load_hdf5_restartFile: set_hdf5_restartFile_name')
      call set_hdf5_restartFile_name(restartCnt,full_restartFileName)
      call nvtxEndRange()

      if(mpi_rank.eq.0) write(*,*) '# Loading restart file: ',trim(adjustl(full_restartFileName))

      call nvtxStartRange('load_hdf5_restartFile: open_hdf5_file')
      call open_hdf5_file(full_restartFileName,file_id)
      call nvtxEndRange()

      call nvtxStartRange('load_hdf5_restartFile: read_real_rp_in_dataset_hdf5_file')
      ! ----  read time  --------------------------------------------------------------------------
      dsetname = 'time'
      call read_real_rp_in_dataset_hdf5_file(file_id,dsetname,time)
      call nvtxEndRange()

      call nvtxStartRange('load_hdf5_restartFile: read_int4_in_dataset_hdf5_file')
      ! ----  read istep --------------------------------------------------------------------------
      dsetname = 'istep'
      call read_int4_in_dataset_hdf5_file(file_id,dsetname,load_step)
      call nvtxEndRange()

      ! ----  read arrays  --------------------------------------------------------------------------
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      call nvtxStartRange('load_hdf5_restartFile: read_array1D_rp_in_dataset_hdf5_file')
      dsetname = 'rho'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,rho)

      dsetname = 'u_x'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u(:,1))

      dsetname = 'u_y'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u(:,2))

      dsetname = 'u_z'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u(:,3))

      dsetname = 'pr'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,pr)

      dsetname = 'E'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,E)

      dsetname = 'mue'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,aux_mu_e)

      dsetname = 'mut'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,aux_mu_t)
      call nvtxEndRange()

      if(flag_walave_value) then
         dsetname = 'walave_u_x'
         call nvtxStartRange('load_hdf5_restartFile: CPU-DATA')
         call h5lexists_f(file_id,dsetname, link_exists, h5err)
         if(h5err /= 0) then
            write(*,*) ' error checking if walave_u_x exists in restart file'
            call MPI_Abort(app_comm,-1,mpi_err)
         end if
         call nvtxEndRange()

         if(link_exists) then
            if(mpi_rank.eq.0) write(111,*) ' walave_u exists'

            call nvtxStartRange('load_hdf5_restartFile: read_array1D_rp_in_dataset_hdf5_file')
            call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,walave_u(:,1))
            ! I suposse that if walave_u_x exists then y and z also exist - no need to check again
            dsetname = 'walave_u_y'
            call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,walave_u(:,2))
            dsetname = 'walave_u_z'
            call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,walave_u(:,3))
            call nvtxEndRange()

         else
            call nvtxStartRange('load_hdf5_restartFile: CPU-DATA')
            if(mpi_rank.eq.0) write(111,*) ' walave_u does not exist - using u'
            walave_u(:,:) = u(:,:)
            call nvtxEndRange()
         end if

         call nvtxStartRange('load_hdf5_restartFile: GPU-DATA')
         !$acc update device(walave_u(:,:))
         call nvtxEndRange()
      end if

      call nvtxStartRange('load_hdf5_restartFile: GPU-DATA')
      !$acc update device(rho(:))
      !$acc update device(pr(:))
      !$acc update device(E(:))
      !$acc update device(u(:,:))

      !$acc kernels
      do iElem = 1,numElemsRankPar
         do iNode = 1, mnnode
            mu_e(iElem,iNode) = aux_mu_e(connecParOrig(iElem,iNode))
            mu_t(iElem,iNode) = aux_mu_t(connecParOrig(iElem,iNode))
         end do
      end do
      !$acc end kernels
      !$acc update device(mu_e(:,:))
      call nvtxEndRange()
      !-------------------------------------------
      !2. If case is periodic, adjust slave nodes

      if (isMeshPeriodic) then
         call nvtxStartRange('load_hdf5_restartFile: GPU-COMPUTE')
         !$acc parallel loop
         do iPer = 1,nPerRankPar
            rho(masSlaRankPar(iPer,2)) = rho(masSlaRankPar(iPer,1))
            u(masSlaRankPar(iPer,2),1) = u(masSlaRankPar(iPer,1),1)
            u(masSlaRankPar(iPer,2),2) = u(masSlaRankPar(iPer,1),2)
            u(masSlaRankPar(iPer,2),3) = u(masSlaRankPar(iPer,1),3)
            pr(masSlaRankPar(iPer,2))  = pr(masSlaRankPar(iPer,1))
            E(masSlaRankPar(iPer,2))   = E(masSlaRankPar(iPer,1))
         end do
         !$acc end parallel loop
         if(flag_walave_value) then
            !$acc parallel loop
            do iPer = 1,nPerRankPar
               walave_u(masSlaRankPar(iPer,2),1) = walave_u(masSlaRankPar(iPer,1),1)
               walave_u(masSlaRankPar(iPer,2),2) = walave_u(masSlaRankPar(iPer,1),2)
               walave_u(masSlaRankPar(iPer,2),3) = walave_u(masSlaRankPar(iPer,1),3)
            end do
            !$acc end parallel loop
         end if
         call nvtxEndRange()
      end if

      call nvtxStartRange('load_hdf5_restartFile: close_hdf5_file')
      call close_hdf5_file(file_id)
      call nvtxEndRange()

   end subroutine load_hdf5_restartFile

   subroutine interpolate_scalarField_in_nodes(mnnode,mngaus,Ngp,connecParW,connecParO,origNodeScalarField,interpNodeScalarField)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,connecParW(numElemsRankPar,mnnode),connecParO(numElemsRankPar,mnnode)
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      real(rp_vtk),intent(in) :: origNodeScalarField(numNodesRankPar)
      real(rp_vtk),intent(out) :: interpNodeScalarField(numNodesRankPar)
      integer(4) :: iElem,igp,inode
      real(rp_vtk) :: var_a,Ngp_vtk(mngaus,mnnode)

      call nvtxStartRange('interpolate_scalarField_in_nodes: GPU-DATA')
      !$acc kernels
      Ngp_vtk(:,:) = real(Ngp(:,:),rp_vtk)
      !$acc end kernels
      call nvtxEndRange()

      call nvtxStartRange('interpolate_scalarField_in_nodes: GPU-COMPUTE')
      !$acc parallel loop gang
      do iElem = 1,numElemsRankPar
         !$acc loop vector
         do igp = 1,mngaus
            var_a = 0.0_rp_vtk
            !$acc loop seq
            do inode = 1,mnnode
               var_a = var_a+Ngp_vtk(igp,inode)*origNodeScalarField(connecParW(iElem,inode))
            end do
            !$acc atomic write
            interpNodeScalarField(connecParO(iElem,igp)) = var_a
            !$acc end atomic
         end do
      end do
      !$acc end parallel loop
      call nvtxEndRange()

   end subroutine interpolate_scalarField_in_nodes

   subroutine interpolate_vectorField_in_nodes(mnnode,mngaus,Ngp,connecParW,connecParO,origNodeVectorField,interpNodeVectorField)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,connecParW(numElemsRankPar,mnnode),connecParO(numElemsRankPar,mnnode)
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      real(rp_vtk),intent(in) :: origNodeVectorField(numNodesRankPar,ndime)
      real(rp_vtk),intent(out) :: interpNodeVectorField(numNodesRankPar,ndime)
      integer(4) :: iElem,igp,inode,idime
      real(rp_vtk) :: var_a,Ngp_vtk(mngaus,mnnode)

      call nvtxStartRange('interpolate_vectorField_in_nodes: GPU-DATA')
      !$acc kernels
      Ngp_vtk(:,:) = real(Ngp(:,:),rp_vtk)
      !$acc end kernels
      call nvtxEndRange()

      call nvtxStartRange('interpolate_vectorField_in_nodes: GPU-COMPUTE')
      !$acc parallel loop gang
      do iElem = 1,numElemsRankPar
         !$acc loop vector collapse(2)
         do igp = 1,mngaus
            do idime = 1,ndime
               var_a = 0.0_rp_vtk
               !$acc loop seq
               do inode = 1,mnnode
                  var_a = var_a+Ngp_vtk(igp,inode)*origNodeVectorField(connecParW(iElem,inode),idime)
               end do
               !$acc atomic write
               interpNodeVectorField(connecParO(iElem,igp),idime) = var_a
               !$acc end atomic
            end do
         end do
      end do
      !$acc end parallel loop
      call nvtxEndRange()

   end subroutine interpolate_vectorField_in_nodes

   subroutine interpolate_elemGpScalarField_in_nodes_for_inst(mnnode,mngaus,Ngp,connecParW,connecParO,origElemGpScalarField,interpNodeScalarField)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,connecParW(numElemsRankPar,mnnode),connecParO(numElemsRankPar,mnnode)
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      real(rp),intent(in) :: origElemGpScalarField(numElemsRankPar,mngaus)
      real(rp_vtk),intent(out) :: interpNodeScalarField(numNodesRankPar)
      integer(4) :: iElem,igp,inode,iPer
      real(rp) :: var_a

      call nvtxStartRange('interpolate_elemGpScalarField_in_nodes_for_inst: GPU-COMPUTE')
      !$acc parallel loop gang
      do iElem = 1,numElemsRankPar

         !$acc loop vector
         do igp = 1,mngaus
            var_a = 0.0_rp
            !$acc loop seq
            do inode = 1,mnnode
               var_a = var_a+Ngp(igp,inode)*origElemGpScalarField(iElem,inode)
            end do
            !$acc atomic write
            interpNodeScalarField(connecParO(iElem,igp)) = real(var_a,rp_vtk)
            !$acc end atomic
         end do
      end do
      !$acc end parallel loop

      if(isMeshPeriodic) then
         !$acc parallel loop
         do iPer = 1,nPerRankPar
            interpNodeScalarField(masSlaRankPar(iPer,2)) = interpNodeScalarField(masSlaRankPar(iPer,1))
         end do
         !$acc end parallel loop
      end if
      call nvtxEndRange()

   end subroutine interpolate_elemGpScalarField_in_nodes_for_inst

   subroutine interpolate_elemGpScalarField_in_nodes_for_avg(mnnode,mngaus,Ngp,connecParW,connecParO,origElemGpScalarField,interpNodeScalarField)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,connecParW(numElemsRankPar,mnnode),connecParO(numElemsRankPar,mnnode)
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      real(rp_avg),intent(in) :: origElemGpScalarField(numElemsRankPar,mngaus)
      real(rp_vtk),intent(out) :: interpNodeScalarField(numNodesRankPar)
      integer(4) :: iElem,igp,inode,iPer
      real(rp_avg) :: var_a,Ngp_avg(mngaus,mnnode)

      call nvtxStartRange('interpolate_elemGpScalarField_in_nodes_for_avg: GPU-DATA')
      !$acc kernels
      Ngp_avg(:,:) = real(Ngp(:,:),rp_avg)
      !$acc end kernels
      call nvtxEndRange()

      call nvtxStartRange('interpolate_elemGpScalarField_in_nodes_for_avg: GPU-COMPUTE')
      !$acc parallel loop gang
      do iElem = 1,numElemsRankPar

         !$acc loop vector
         do igp = 1,mngaus
            var_a = 0.0_rp_avg
            !$acc loop seq
            do inode = 1,mnnode
               var_a = var_a+Ngp_avg(igp,inode)*origElemGpScalarField(iElem,inode)
            end do
            !$acc atomic write
            interpNodeScalarField(connecParO(iElem,igp)) = real(var_a,rp_vtk)
            !$acc end atomic
         end do
      end do
      !$acc end parallel loop

      if(isMeshPeriodic) then
         !$acc parallel loop
         do iPer = 1,nPerRankPar
            interpNodeScalarField(masSlaRankPar(iPer,2)) = interpNodeScalarField(masSlaRankPar(iPer,1))
         end do
         !$acc end parallel loop
      end if
      call nvtxEndRange()

   end subroutine interpolate_elemGpScalarField_in_nodes_for_avg

   subroutine interpolate_scalarField_in_elemGp(mnnode,mngaus,Ngp,connecParW,connecParO,origNodeScalarField,interpElemGpScalarField)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,connecParW(numElemsRankPar,mnnode),connecParO(numElemsRankPar,mnnode)
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      real(rp_vtk),intent(in) :: origNodeScalarField(numNodesRankPar)
      real(rp_vtk),intent(out) :: interpElemGpScalarField(numElemsRankPar,mngaus)
      integer(4) :: iElem,igp,inode,iPer
      real(rp_vtk) :: var_a,Ngp_vtk(mngaus,mnnode)

      call nvtxStartRange('interpolate_scalarField_in_elemGp: GPU-DATA')
      !$acc kernels
      Ngp_vtk(:,:) = real(Ngp(:,:),rp_vtk)
      !$acc end kernels
      call nvtxEndRange()

      call nvtxStartRange('interpolate_scalarField_in_elemGp: GPU-COMPUTE')
      !$acc parallel loop gang
      do iElem = 1,numElemsRankPar
         !$acc loop vector
         do igp = 1,mngaus
            var_a = 0.0_rp_vtk
            !$acc loop seq
            do inode = 1,mnnode
               var_a = var_a+Ngp_vtk(igp,inode)*origNodeScalarField(connecParW(iElem,inode))
            end do
            !$acc atomic write
            interpElemGpScalarField(iElem,igp) = var_a
            !$acc end atomic
         end do
      end do
      !$acc end parallel loop
      call nvtxEndRange()

   end subroutine interpolate_scalarField_in_elemGp

   subroutine copyPeriodicNodes_scalarField(nodeScalarField)
      implicit none
      real(rp_vtk),intent(inout) :: nodeScalarField(numNodesRankPar)
      integer(4) :: iPer

      call nvtxStartRange('copyPeriodicNodes_scalarField: GPU-COMPUTE')
      if(isMeshPeriodic) then
         !$acc parallel loop
         do iPer = 1,nPerRankPar
            nodeScalarField(masSlaRankPar(iPer,2)) = nodeScalarField(masSlaRankPar(iPer,1))
         end do
         !$acc end parallel loop
      end if
      call nvtxEndRange()

   end subroutine copyPeriodicNodes_scalarField

   subroutine copyPeriodicNodes_vectorField(nodeVectorField)
      implicit none
      real(rp_vtk),intent(inout) :: nodeVectorField(numNodesRankPar,ndime)
      integer(4) :: iPer

      call nvtxStartRange('copyPeriodicNodes_vectorField: GPU-COMPUTE')
      if(isMeshPeriodic) then
         !$acc parallel loop
         do iPer = 1,nPerRankPar
            nodeVectorField(masSlaRankPar(iPer,2),1) = nodeVectorField(masSlaRankPar(iPer,1),1)
            nodeVectorField(masSlaRankPar(iPer,2),2) = nodeVectorField(masSlaRankPar(iPer,1),2)
            nodeVectorField(masSlaRankPar(iPer,2),3) = nodeVectorField(masSlaRankPar(iPer,1),3)
         end do
         !$acc end parallel loop
      end if
      call nvtxEndRange()

   end subroutine copyPeriodicNodes_vectorField

   subroutine copy_elemGpScalarField_in_nodes_for_inst(mnnode,connecParW,connecParO,origElemGpScalarField,interpNodeScalarField)
      implicit none
      integer(4),intent(in) :: mnnode,connecParW(numElemsRankPar,mnnode),connecParO(numElemsRankPar,mnnode)
      real(rp),intent(in) :: origElemGpScalarField(numElemsRankPar,mnnode)
      real(rp_vtk),intent(out) :: interpNodeScalarField(numNodesRankPar)
      integer(4) :: iElem,inode,iPer

      call nvtxStartRange('copy_elemGpScalarField_in_nodes_for_inst: GPU-COMPUTE')
      !$acc parallel loop gang
      do iElem = 1,numElemsRankPar
         !$acc loop vector
         do inode = 1,mnnode
            interpNodeScalarField(connecParO(iElem,inode)) = real(origElemGpScalarField(iElem,inode),rp_vtk)
         end do
      end do
      !$acc end parallel loop

      if(isMeshPeriodic) then
         !$acc parallel loop
         do iPer = 1,nPerRankPar
            interpNodeScalarField(masSlaRankPar(iPer,2)) = interpNodeScalarField(masSlaRankPar(iPer,1))
         end do
         !$acc end parallel loop
      end if
      call nvtxEndRange()

   end subroutine copy_elemGpScalarField_in_nodes_for_inst

   subroutine copy_elemGpScalarField_in_nodes_for_avg(mnnode,connecParW,connecParO,origElemGpScalarField,interpNodeScalarField)
      implicit none
      integer(4),intent(in) :: mnnode,connecParW(numElemsRankPar,mnnode),connecParO(numElemsRankPar,mnnode)
      real(rp_avg),intent(in) :: origElemGpScalarField(numElemsRankPar,mnnode)
      real(rp_vtk),intent(out) :: interpNodeScalarField(numNodesRankPar)
      integer(4) :: iElem,inode,iPer

      call nvtxStartRange('copy_elemGpScalarField_in_nodes_for_avg: GPU-COMPUTE')
      !$acc parallel loop gang
      do iElem = 1,numElemsRankPar
         !$acc loop vector
         do inode = 1,mnnode
            interpNodeScalarField(connecParO(iElem,inode)) = real(origElemGpScalarField(iElem,inode),rp_vtk)
         end do
      end do
      !$acc end parallel loop

      if(isMeshPeriodic) then
         !$acc parallel loop
         do iPer = 1,nPerRankPar
            interpNodeScalarField(masSlaRankPar(iPer,2)) = interpNodeScalarField(masSlaRankPar(iPer,1))
         end do
         !$acc end parallel loop
      end if
      call nvtxEndRange()

   end subroutine copy_elemGpScalarField_in_nodes_for_avg

   subroutine copy_scalarField_in_elemGp(mnnode,connecParW,connecParO,origNodeScalarField,interpElemGpScalarField)
      implicit none
      integer(4),intent(in) :: mnnode,connecParW(numElemsRankPar,mnnode),connecParO(numElemsRankPar,mnnode)
      real(rp_vtk),intent(in) :: origNodeScalarField(numNodesRankPar)
      real(rp_vtk),intent(out) :: interpElemGpScalarField(numElemsRankPar,mnnode)
      integer(4) :: iElem,inode

      call nvtxStartRange('copy_scalarField_in_elemGp: GPU-COMPUTE')
      !$acc parallel loop gang
      do iElem = 1,numElemsRankPar
         !$acc loop vector
         do inode = 1,mnnode
            interpElemGpScalarField(iElem,inode) = origNodeScalarField(connecParO(iElem,inode))
         end do
      end do
      !$acc end parallel loop
      call nvtxEndRange()

   end subroutine copy_scalarField_in_elemGp

!----------------------------------------------------------------------------------------------------------------------------------
   subroutine copy_nodeScalarField2save_in_aux_for_inst(nodeScalarField)
      implicit none
      real(rp),intent(inout) :: nodeScalarField(numNodesRankPar)

      call nvtxStartRange('copy_nodeScalarField2save_in_aux_for_inst: GPU-DATA')
      !$acc kernels
      auxNodeScalarField_vtk(:) = real(nodeScalarField(:),rp_vtk)
      !$acc end kernels
      call nvtxEndRange()

   end subroutine copy_nodeScalarField2save_in_aux_for_inst

   subroutine copy_nodeScalarField2save_in_aux_for_avg(nodeScalarField)
      implicit none
      real(rp_avg),intent(inout) :: nodeScalarField(numNodesRankPar)

      call nvtxStartRange('copy_nodeScalarField2save_in_aux_for_avg: GPU-DATA')
      !$acc kernels
      auxNodeScalarField_vtk(:) = real(nodeScalarField(:),rp_vtk)
      !$acc end kernels
      call nvtxEndRange()

   end subroutine copy_nodeScalarField2save_in_aux_for_avg

   subroutine copy_nodeVectorField2save_in_aux_for_inst(nodeVectorField)
      implicit none
      real(rp),intent(inout) :: nodeVectorField(numNodesRankPar,ndime)

      call nvtxStartRange('copy_nodeVectorField2save_in_aux_for_inst: GPU-DATA')
      !$acc kernels
      auxNodeVectorField_vtk(:,:) = real(nodeVectorField(:,:),rp_vtk)
      !$acc end kernels
      call nvtxEndRange()

   end subroutine copy_nodeVectorField2save_in_aux_for_inst

   subroutine copy_nodeVectorField2save_in_aux_for_avg(nodeVectorField)
      implicit none
      real(rp_avg),intent(inout) :: nodeVectorField(numNodesRankPar,ndime)

      call nvtxStartRange('copy_nodeVectorField2save_in_aux_for_avg: GPU-DATA')
      !$acc kernels
      auxNodeVectorField_vtk(:,:) = real(nodeVectorField(:,:),rp_vtk)
      !$acc end kernels
      call nvtxEndRange()

   end subroutine copy_nodeVectorField2save_in_aux_for_avg

   subroutine save_hdf5_resultsFile_baseFunc(mnnode,mngaus,Ngp_equi,hdf5_fileId,save_type_inst,numNodeScalarFields2save,nodeScalarFields2save,&
      numNodeVectorFields2save,nodeVectorFields2save,&
      numElemGpScalarFields2save,elemGpScalarFields2save)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus
      real(rp),intent(in) :: Ngp_equi(mngaus,mnnode)
      integer(hid_t),intent(in) :: hdf5_fileId
      logical,intent(in) :: save_type_inst          ! Denote if instaneous (True) or average (False) data is to be saved
      integer(4),intent(in) :: numNodeScalarFields2save,numNodeVectorFields2save,numElemGpScalarFields2save
      type(ptr_array1d_rp_save),intent(in) :: nodeScalarFields2save(:)
      type(ptr_array2d_rp_save),intent(in) :: nodeVectorFields2save(:),elemGpScalarFields2save(:)

      integer(hsize_t) :: ds_dims(1),ms_dims(1),ds_dims2d(2),ms_dims2d(2)
      integer(hssize_t) :: ms_offset(1),ms_offset2d(2)
      integer(4) :: h5err
      character(512) :: groupname,dsetname
      integer(4) :: iElem,iGp,iPer,iField

      call nvtxStartRange('save_hdf5_resultsFile_baseFunc: create_vtkhdf_unstructuredGrid_struct_for_resultsFile')
      !   Creating the VTK-HDF structure
      call create_vtkhdf_unstructuredGrid_struct_for_resultsFile(mnnode,hdf5_fileId)
      call nvtxEndRange()

      ds_dims(1) = int(totalNumNodesPar,hsize_t)
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1
      groupname = '/VTKHDF/PointData/'

      ! Scalar Fields
      do iField=1,numNodeScalarFields2save
         dsetname = trim(adjustl(groupname))//trim(nodeScalarFields2save(iField)%nameField)
         !if(mpi_rank.eq.0) write(*,*) 'saving field',iField,'name',dsetname,'type_inst',save_type_inst

         if(save_type_inst) then    ! Instantaneous fields
            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: copy_nodeScalarField2save_in_aux_for_inst')
            call copy_nodeScalarField2save_in_aux_for_inst(nodeScalarFields2save(iField)%ptr_rp)
            call nvtxEndRange()
         else                       ! Averages fields
            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: copy_nodeScalarField2save_in_aux_for_avg')
            call copy_nodeScalarField2save_in_aux_for_avg(nodeScalarFields2save(iField)%ptr_avg)
            call nvtxEndRange()
         end if

         if(isMeshLinealOutput) then
            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: copyPeriodicNodes_scalarField')
            call copyPeriodicNodes_scalarField(auxNodeScalarField_vtk)
            call nvtxEndRange()

            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: GPU-DATA')
            !$acc update host(auxNodeScalarField_vtk(:))
            call nvtxEndRange()

            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: save_array1D_rp_vtk_in_dataset_hdf5_file')
            call save_array1D_rp_vtk_in_dataset_hdf5_file(hdf5_fileId,dsetname,ds_dims,ms_dims,ms_offset,auxNodeScalarField_vtk)
            call nvtxEndRange()
         else
            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: interpolate_scalarField_in_nodes')
            call interpolate_scalarField_in_nodes(mnnode,mngaus,Ngp_equi,connecParWork,connecParOrig,auxNodeScalarField_vtk,auxInterpNodeScalarField)
            call nvtxEndRange()

            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: GPU-DATA')
            !$acc update host(auxInterpNodeScalarField(:))
            call nvtxEndRange()

            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: save_array1D_rp_vtk_in_dataset_hdf5_file')
            call save_array1D_rp_vtk_in_dataset_hdf5_file(hdf5_fileId,dsetname,ds_dims,ms_dims,ms_offset,auxInterpNodeScalarField)
            call nvtxEndRange()
         end if


      end do

      ! Vector Fields
      ds_dims2d(1) = int(ndime,hsize_t)
      ds_dims2d(2) = int(totalNumNodesPar,hsize_t)
      ms_dims2d(1) = int(ndime,hsize_t)
      ms_dims2d(2) = int(numNodesRankPar,hsize_t)
      ms_offset2d(1) = 0
      ms_offset2d(2) = int(rankNodeStart,hssize_t)-1

      do iField=1,numNodeVectorFields2save
         dsetname = trim(adjustl(groupname))//trim(nodeVectorFields2save(iField)%nameField)

         if(save_type_inst) then
            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: copy_nodeVectorField2save_in_aux_for_inst')
            call copy_nodeVectorField2save_in_aux_for_inst(nodeVectorFields2save(iField)%ptr_rp)
            call nvtxEndRange()
         else
            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: copy_nodeVectorField2save_in_aux_for_avg')
            call copy_nodeVectorField2save_in_aux_for_avg(nodeVectorFields2save(iField)%ptr_avg)
            call nvtxEndRange()
         end if

         if(isMeshLinealOutput) then
            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: copyPeriodicNodes_vectorField')
            call copyPeriodicNodes_vectorField(auxNodeVectorField_vtk)
            call nvtxEndRange()

            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: GPU-DATA')
            !$acc update host(auxNodeVectorField_vtk(:,:))
            call nvtxEndRange()

            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: save_array2D_tr_rp_vtk_in_dataset_hdf5_file')
            call save_array2D_tr_rp_vtk_in_dataset_hdf5_file(hdf5_fileId,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxNodeVectorField_vtk)
            call nvtxEndRange()
         else
            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: interpolate_vectorField_in_nodes')
            call interpolate_vectorField_in_nodes(mnnode,mngaus,Ngp_equi,connecParWork,connecParOrig,auxNodeVectorField_vtk,auxInterpNodeVectorField)
            call nvtxEndRange()

            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: GPU-DATA')
            !$acc update host(auxInterpNodeVectorField(:,:))
            call nvtxEndRange()

            call nvtxStartRange('save_hdf5_resultsFile_baseFunc: save_array2D_tr_rp_vtk_in_dataset_hdf5_file')
            call save_array2D_tr_rp_vtk_in_dataset_hdf5_file(hdf5_fileId,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxInterpNodeVectorField)
            call nvtxEndRange()
         end if
      end do

      !--------------------------------------------------------------------------------------------------------------------------------------
      ! ElemGP Fields
      !--------------------------------------------------------------------------------------------------------------------------------------
      do iField=1,numElemGpScalarFields2save
         dsetname = trim(adjustl(groupname))//trim(elemGpScalarFields2save(iField)%nameField)

         if(isMeshLinealOutput) then
            if(save_type_inst) then
               call nvtxStartRange('save_hdf5_resultsFile_baseFunc: copy_elemGpScalarField_in_nodes_for_inst')
               call copy_elemGpScalarField_in_nodes_for_inst(mnnode,connecParWork,connecParOrig,elemGpScalarFields2save(iField)%ptr_rp,auxInterpNodeScalarField)
               call nvtxEndRange()
            else
               call nvtxStartRange('save_hdf5_resultsFile_baseFunc: copy_elemGpScalarField_in_nodes_for_avg')
               call copy_elemGpScalarField_in_nodes_for_avg(mnnode,connecParWork,connecParOrig,elemGpScalarFields2save(iField)%ptr_avg,auxInterpNodeScalarField)
               call nvtxEndRange()
            end if
         else
            if(save_type_inst) then
               call nvtxStartRange('save_hdf5_resultsFile_baseFunc: interpolate_elemGpScalarField_in_nodes_for_inst')
               call interpolate_elemGpScalarField_in_nodes_for_inst(mnnode,mngaus,Ngp_equi,connecParWork,connecParOrig,elemGpScalarFields2save(iField)%ptr_rp,auxInterpNodeScalarField)
               call nvtxEndRange()
            else
               call nvtxStartRange('save_hdf5_resultsFile_baseFunc: interpolate_elemGpScalarField_in_nodes_for_avg')
               call interpolate_elemGpScalarField_in_nodes_for_avg(mnnode,mngaus,Ngp_equi,connecParWork,connecParOrig,elemGpScalarFields2save(iField)%ptr_avg,auxInterpNodeScalarField)
               call nvtxEndRange()
            end if
         end if

         call nvtxStartRange('save_hdf5_resultsFile_baseFunc: GPU-DATA')
         !$acc update host(auxInterpNodeScalarField(:))
         call nvtxEndRange()

         call nvtxStartRange('save_hdf5_resultsFile_baseFunc: save_array1D_rp_vtk_in_dataset_hdf5_file')
         call save_array1D_rp_vtk_in_dataset_hdf5_file(hdf5_fileId,dsetname,ds_dims,ms_dims,ms_offset,auxInterpNodeScalarField)
         call nvtxEndRange()
      end do

   end subroutine save_hdf5_resultsFile_baseFunc

!----------------------------------------------------------------------------------------------------------------------------------
   subroutine copy_aux_in_nodeScalarField2save_for_inst(nodeScalarField)
      implicit none
      real(rp),intent(inout) :: nodeScalarField(numNodesRankPar)

      call nvtxStartRange('copy_aux_in_nodeScalarField2save_for_inst: GPU-DATA')
      !$acc kernels
      nodeScalarField(:) = real(auxNodeScalarField_vtk(:),rp)
      !$acc end kernels
      call nvtxEndRange()

   end subroutine copy_aux_in_nodeScalarField2save_for_inst

   subroutine copy_aux_in_nodeScalarField2save_for_avg(nodeScalarField)
      implicit none
      real(rp_avg),intent(inout) :: nodeScalarField(numNodesRankPar)

      call nvtxStartRange('copy_aux_in_nodeScalarField2save_for_avg: GPU-DATA')
      !$acc kernels
      nodeScalarField(:) = real(auxNodeScalarField_vtk(:),rp_avg)
      !$acc end kernels
      call nvtxEndRange()

   end subroutine copy_aux_in_nodeScalarField2save_for_avg

   subroutine copy_aux_in_nodeVectorField2save_for_inst(nodeVectorField)
      implicit none
      real(rp),intent(inout) :: nodeVectorField(numNodesRankPar,ndime)

      call nvtxStartRange('copy_aux_in_nodeVectorField2save_for_inst: GPU-DATA')
      !$acc kernels
      nodeVectorField(:,:) = real(auxNodeVectorField_vtk(:,:),rp)
      !$acc end kernels
      call nvtxEndRange()

   end subroutine copy_aux_in_nodeVectorField2save_for_inst

   subroutine copy_aux_in_nodeVectorField2save_for_avg(nodeVectorField)
      implicit none
      real(rp_avg),intent(inout) :: nodeVectorField(numNodesRankPar,ndime)

      call nvtxStartRange('copy_aux_in_nodeVectorField2save_for_avg: GPU-DATA')
      !$acc kernels
      nodeVectorField(:,:) = real(auxNodeVectorField_vtk(:,:),rp_avg)
      !$acc end kernels
      call nvtxEndRange()

   end subroutine copy_aux_in_nodeVectorField2save_for_avg
!----------------------------------------------------------------------------------------------------------------------------------

   subroutine load_hdf5_resultsFile_baseFunc(mnnode,mngaus,Ngp,hdf5_fileId,load_type_inst,numNodeScalarFields2load,nodeScalarFields2load,&
      numNodeVectorFields2load,nodeVectorFields2load,&
      numElemGpScalarFields2load,elemGpScalarFields2load)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      integer(hid_t),intent(in) :: hdf5_fileId
      logical,intent(in) :: load_type_inst          ! Denote if instaneous (True) or average (False) data is to be loaded
      integer(4),intent(in) :: numNodeScalarFields2load,numNodeVectorFields2load,numElemGpScalarFields2load
      type(ptr_array1d_rp_save),intent(inout) :: nodeScalarFields2load(:)
      type(ptr_array2d_rp_save),intent(inout) :: nodeVectorFields2load(:),elemGpScalarFields2load(:)

      integer(hsize_t) :: ds_dims(1),ms_dims(1),ds_dims2d(2),ms_dims2d(2)
      integer(hssize_t) :: ms_offset(1),ms_offset2d(2)
      integer(4) :: h5err
      character(512) :: groupname,dsetname
      integer(4) :: iElem,iGp,iPer,iField
      real(rp) :: aux_nodeScalarField(numNodesRankPar)
      real(rp) :: aux_array_time(1)

      ds_dims(1) = int(totalNumNodesPar,hsize_t)
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1
      groupname = '/VTKHDF/PointData/'

      do iField=1,numNodeScalarFields2load
         dsetname = trim(adjustl(groupname))//trim(nodeScalarFields2load(iField)%nameField)

         call nvtxStartRange('load_hdf5_resultsFile_baseFunc: read_dataspace_1d_real_rp_vtk_hyperslab_parallel')
         call read_dataspace_1d_real_rp_vtk_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,auxNodeScalarField_vtk)
         call nvtxEndRange()

         call nvtxStartRange('load_hdf5_resultsFile_baseFunc: GPU-DATA')
         !$acc update device(auxNodeScalarField_vtk(:))
         call nvtxEndRange()

         if(.not.(isMeshLinealOutput)) then
            call nvtxStartRange('load_hdf5_resultsFile_baseFunc: interpolate_scalarField_in_nodes')
            call interpolate_scalarField_in_nodes(mnnode,mngaus,Ngp,connecParWork,connecParOrig,auxNodeScalarField_vtk,auxInterpNodeScalarField)
            call nvtxEndRange()

            call nvtxStartRange('load_hdf5_resultsFile_baseFunc: CPU-DATA')
            auxNodeScalarField_vtk(:) = auxInterpNodeScalarField(:)
            call nvtxEndRange()
         end if

         if(load_type_inst) then    ! Instantaneous fields
            call nvtxStartRange('load_hdf5_resultsFile_baseFunc: copy_aux_in_nodeScalarField2save_for_inst')
            call copy_aux_in_nodeScalarField2save_for_inst(nodeScalarFields2load(iField)%ptr_rp)
            call nvtxEndRange()
         else                       ! Averages fields
            call nvtxStartRange('load_hdf5_resultsFile_baseFunc: copy_aux_in_nodeScalarField2save_for_avg')
            call copy_aux_in_nodeScalarField2save_for_avg(nodeScalarFields2load(iField)%ptr_avg)
            call nvtxEndRange()
         end if

      end do

      ds_dims2d(1) = int(ndime,hsize_t)
      ds_dims2d(2) = int(totalNumNodesPar,hsize_t)
      ms_dims2d(1) = int(ndime,hsize_t)
      ms_dims2d(2) = int(numNodesRankPar,hsize_t)
      ms_offset2d(1) = 0
      ms_offset2d(2) = int(rankNodeStart,hssize_t)-1

      do iField=1,numNodeVectorFields2load
         dsetname = trim(adjustl(groupname))//trim(nodeVectorFields2load(iField)%nameField)

         call nvtxStartRange('load_hdf5_resultsFile_baseFunc: read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel')
         call read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims2d,ms_offset2d,auxNodeVectorField_vtk)
         call nvtxEndRange()

         call nvtxStartRange('load_hdf5_resultsFile_baseFunc: GPU-DATA')
         !$acc update device(auxNodeVectorField_vtk(:,:))
         call nvtxEndRange()

         if(.not.(isMeshLinealOutput)) then
            call nvtxStartRange('load_hdf5_resultsFile_baseFunc: interpolate_vectorField_in_nodes')
            call interpolate_vectorField_in_nodes(mnnode,mngaus,Ngp,connecParWork,connecParOrig,auxNodeVectorField_vtk,auxInterpNodeVectorField)
            auxNodeVectorField_vtk(:,:) = auxInterpNodeVectorField(:,:)
            call nvtxEndRange()
         end if

         if(load_type_inst) then    ! Instantaneous fields
            call nvtxStartRange('load_hdf5_resultsFile_baseFunc: copy_aux_in_nodeVectorField2save_for_inst')
            call copy_aux_in_nodeVectorField2save_for_inst(nodeVectorFields2load(iField)%ptr_rp)
            call nvtxEndRange()
         else                       ! Averages fields
            call nvtxStartRange('load_hdf5_resultsFile_baseFunc: copy_aux_in_nodeVectorField2save_for_avg')
            call copy_aux_in_nodeVectorField2save_for_avg(nodeVectorFields2load(iField)%ptr_avg)
            call nvtxEndRange()
         end if
      end do

      !--------------------------------------------------------------------------------------------------------------------------------------
      do iField=1,numElemGpScalarFields2load
         dsetname = trim(adjustl(groupname))//trim(elemGpScalarFields2load(iField)%nameField)

         call nvtxStartRange('load_hdf5_resultsFile_baseFunc: read_dataspace_1d_real_rp_vtk_hyperslab_parallel')
         call read_dataspace_1d_real_rp_vtk_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,auxNodeScalarField_vtk)
         call nvtxEndRange()

         call nvtxStartRange('load_hdf5_resultsFile_baseFunc: GPU-DATA')
         !$acc update device(auxNodeScalarField_vtk(:))
         call nvtxEndRange()

         if(isMeshLinealOutput) then
            call nvtxStartRange('load_hdf5_resultsFile_baseFunc: copy_scalarField_in_elemGp')
            call copy_scalarField_in_elemGp(mnnode,connecParWork,connecParOrig,auxNodeScalarField_vtk,auxNodeVectorField_vtk)
            call nvtxEndRange()
         else
            call nvtxStartRange('load_hdf5_resultsFile_baseFunc: interpolate_scalarField_in_elemGp')
            call interpolate_scalarField_in_elemGp(mnnode,mngaus,Ngp,connecParWork,connecParOrig,auxNodeScalarField_vtk,auxNodeVectorField_vtk)
            call nvtxEndRange()
         end if

         if(load_type_inst) then    ! Instantaneous fields
            call nvtxStartRange('load_hdf5_resultsFile_baseFunc: copy_aux_in_nodeVectorField2save_for_inst')
            call copy_aux_in_nodeVectorField2save_for_inst(nodeVectorFields2load(iField)%ptr_rp)
            call nvtxEndRange()
         else                       ! Averages fields
            call nvtxStartRange('load_hdf5_resultsFile_baseFunc: copy_aux_in_nodeVectorField2save_for_avg')
            call copy_aux_in_nodeVectorField2save_for_avg(nodeVectorFields2load(iField)%ptr_avg)
            call nvtxEndRange()
         end if
      end do

   end subroutine load_hdf5_resultsFile_baseFunc

   subroutine save_instResults_hdf5_file(mnnode,mngaus,Ngp,iStep,time,numNodeScalarFields2save,nodeScalarFields2save,&
      numNodeVectorFields2save,nodeVectorFields2save,&
      numElemGpScalarFields2save,elemGpScalarFields2save)
      implicit none
      integer(4), intent(in) :: mnnode,mngaus,iStep
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      real(rp),intent(in) :: time
      integer(4),intent(in) :: numNodeScalarFields2save,numNodeVectorFields2save,numElemGpScalarFields2save
      type(ptr_array1d_rp_save),intent(in) :: nodeScalarFields2save(:)
      type(ptr_array2d_rp_save),intent(in) :: nodeVectorFields2save(:),elemGpScalarFields2save(:)

      integer(hid_t) :: hdf5_fileId
      character(512) :: full_hdf5_fileName,dsetname

      call nvtxStartRange('save_instResults_hdf5_file: set_hdf5_resultsFile_name')
      call set_hdf5_resultsFile_name(iStep,full_hdf5_fileName)
      call nvtxEndRange()

      call nvtxStartRange('save_instResults_hdf5_file: create_hdf5_file')
      call create_hdf5_file(full_hdf5_fileName,hdf5_fileId)
      call nvtxEndRange()

      call nvtxStartRange('save_instResults_hdf5_file: save_hdf5_resultsFile_baseFunc')
      call save_hdf5_resultsFile_baseFunc(mnnode,mngaus,Ngp,hdf5_fileId,.true.,numNodeScalarFields2save,nodeScalarFields2save,&
         numNodeVectorFields2save,nodeVectorFields2save,&
         numElemGpScalarFields2save,elemGpScalarFields2save)
      call nvtxEndRange()

      call nvtxStartRange('save_instResults_hdf5_file: save_real_rp_in_dataset_hdf5_file')
      dsetname = 'time'
      call save_real_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,time)
      call nvtxEndRange()

      call nvtxStartRange('save_instResults_hdf5_file: close_hdf5_file')
      call close_hdf5_file(hdf5_fileId)
      call nvtxEndRange()

   end subroutine save_instResults_hdf5_file

   subroutine save_avgResults_hdf5_file(mnnode,mngaus,Ngp,restartCnt,initial_avgTime,elapsed_avgTime,numAvgNodeScalarFields2save,avgNodeScalarFields2save,&
      numAvgNodeVectorFields2save,avgNodeVectorFields2save,&
      numAvgElemGpScalarFields2save,avgElemGpScalarFields2save)
      implicit none
      integer(4), intent(in) :: mnnode,mngaus,restartCnt
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      real(rp),intent(in) :: initial_avgTime,elapsed_avgTime
      integer(4),intent(in) :: numAvgNodeScalarFields2save,numAvgNodeVectorFields2save,numAvgElemGpScalarFields2save
      type(ptr_array1d_rp_save),intent(in) :: avgNodeScalarFields2save(:)
      type(ptr_array2d_rp_save),intent(in) :: avgNodeVectorFields2save(:),avgElemGpScalarFields2save(:)

      integer(hid_t) :: hdf5_fileId
      character(512) :: full_hdf5_fileName,dsetname

      call nvtxStartRange('save_avgResults_hdf5_file: set_hdf5_avgResultsFile_name')
      call set_hdf5_avgResultsFile_name(restartCnt,full_hdf5_fileName)
      call nvtxEndRange()

      call nvtxStartRange('save_avgResults_hdf5_file: create_hdf5_file')
      call create_hdf5_file(full_hdf5_fileName,hdf5_fileId)
      call nvtxEndRange()

      call nvtxStartRange('save_avgResults_hdf5_file: save_hdf5_resultsFile_baseFunc')
      call save_hdf5_resultsFile_baseFunc(mnnode,mngaus,Ngp,hdf5_fileId,.false.,numAvgNodeScalarFields2save,avgNodeScalarFields2save,&
         numAvgNodeVectorFields2save,avgNodeVectorFields2save,&
         numAvgElemGpScalarFields2save,avgElemGpScalarFields2save)
      call nvtxEndRange()

      call nvtxStartRange('save_avgResults_hdf5_file: save_real_rp_in_dataset_hdf5_file')
      dsetname = 'elapsed_avgTime'
      call save_real_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,elapsed_avgTime)
      call nvtxEndRange()

      call nvtxStartRange('save_avgResults_hdf5_file: save_real_rp_in_dataset_hdf5_file')
      dsetname = 'initial_avgTime'
      call save_real_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,initial_avgTime)
      call nvtxEndRange()
      
      call nvtxStartRange('save_avgResults_hdf5_file: close_hdf5_file')
      call close_hdf5_file(hdf5_fileId)
      call nvtxEndRange()

   end subroutine save_avgResults_hdf5_file

   subroutine load_avgResults_hdf5_file(mnnode,mngaus,Ngp,restartCnt,initial_avgTime,elapsed_avgTime,numAvgNodeScalarFields2load,avgNodeScalarFields2load,&
      numAvgNodeVectorFields2load,avgNodeVectorFields2load,&
      numAvgElemGpScalarFields2load,avgElemGpScalarFields2load)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,restartCnt
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      real(rp),intent(inout) :: initial_avgTime,elapsed_avgTime
      integer(4),intent(in) :: numAvgNodeScalarFields2load,numAvgNodeVectorFields2load,numAvgElemGpScalarFields2load
      type(ptr_array1d_rp_save),intent(inout) :: avgNodeScalarFields2load(:)
      type(ptr_array2d_rp_save),intent(inout) :: avgNodeVectorFields2load(:),avgElemGpScalarFields2load(:)

      integer(hid_t) :: hdf5_fileId
      character(512) :: full_hdf5_fileName,dsetname

      call nvtxStartRange('load_avgResults_hdf5_file: set_hdf5_avgResultsFile_name')
      call set_hdf5_avgResultsFile_name(restartCnt,full_hdf5_fileName)
      call nvtxEndRange()

      call nvtxStartRange('load_avgResults_hdf5_file: open_hdf5_file')
      call open_hdf5_file(full_hdf5_fileName,hdf5_fileId)
      call nvtxEndRange()

      call nvtxStartRange('load_avgResults_hdf5_file: load_hdf5_resultsFile_baseFunc')
      call load_hdf5_resultsFile_baseFunc(mnnode,mngaus,Ngp,hdf5_fileId,.False.,numAvgNodeScalarFields2load,avgNodeScalarFields2load,&
         numAvgNodeVectorFields2load,avgNodeVectorFields2load,&
         numAvgElemGpScalarFields2load,avgElemGpScalarFields2load)
      call nvtxEndRange()

      call nvtxStartRange('load_avgResults_hdf5_file: read_real_rp_in_dataset_hdf5_file')
      dsetname = 'elapsed_avgTime'
      call read_real_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,elapsed_avgTime)
      call nvtxEndRange()

      call nvtxStartRange('load_avgResults_hdf5_file: read_real_rp_in_dataset_hdf5_file')
      dsetname = 'initial_avgTime'
      call read_real_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,initial_avgTime)
      call nvtxEndRange()

      call nvtxStartRange('load_avgResults_hdf5_file: close_hdf5_file')
      call close_hdf5_file(hdf5_fileId)
      call nvtxEndRange()

   end subroutine load_avgResults_hdf5_file

   !-------------------------------------------------------------------------------------------------------------------------------

   subroutine save_surface_instResults_hdf5_file(iStep,numNodeScalarFields2save,nodeScalarFields2save,&
      numNodeVectorFields2save,nodeVectorFields2save,&
      numElemGpScalarFields2save,elemGpScalarFields2save)
      implicit none
      integer(4), intent(in) :: iStep
      integer(4),intent(in) :: numNodeScalarFields2save,numNodeVectorFields2save,numElemGpScalarFields2save
      character(512) :: res_hdf5_fileName
      type(ptr_array1d_rp_save),intent(in) :: nodeScalarFields2save(:)
      type(ptr_array2d_rp_save),intent(in) :: nodeVectorFields2save(:),elemGpScalarFields2save(:)

      call nvtxStartRange('save_surface_instResults_hdf5_file: set_hdf5_resultsFile_name')
      call set_hdf5_resultsFile_name(iStep,res_hdf5_fileName)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_instResults_hdf5_file: save_surface_results_hdf5_file')
      call save_surface_results_hdf5_file(res_hdf5_fileName,numNodeScalarFields2save,nodeScalarFields2save,&
         numNodeVectorFields2save,nodeVectorFields2save,&
         numElemGpScalarFields2save,elemGpScalarFields2save)
      call nvtxEndRange()

   end subroutine save_surface_instResults_hdf5_file

   subroutine save_surface_avgResults_hdf5_file(restartCnt,numAvgNodeScalarFields2save,avgNodeScalarFields2save,&
      numAvgNodeVectorFields2save,avgNodeVectorFields2save,&
      numAvgElemGpScalarFields2save,avgElemGpScalarFields2save)
      implicit none
      integer(4), intent(in) :: restartCnt
      integer(4),intent(in) :: numAvgNodeScalarFields2save,numAvgNodeVectorFields2save,numAvgElemGpScalarFields2save
      type(ptr_array1d_rp_save),intent(in) :: avgNodeScalarFields2save(:)
      type(ptr_array2d_rp_save),intent(in) :: avgNodeVectorFields2save(:),avgElemGpScalarFields2save(:)
      character(512) :: res_hdf5_fileName

      call nvtxStartRange('save_surface_avgResults_hdf5_file: set_hdf5_avgResultsFile_name')
      call set_hdf5_avgResultsFile_name(restartCnt,res_hdf5_fileName)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_avgResults_hdf5_file: save_surface_results_hdf5_file')
      call save_surface_results_hdf5_file(res_hdf5_fileName,numAvgNodeScalarFields2save,avgNodeScalarFields2save,&
         numAvgNodeVectorFields2save,avgNodeVectorFields2save,&
         numAvgElemGpScalarFields2save,avgElemGpScalarFields2save) 
      call nvtxEndRange()

   end subroutine save_surface_avgResults_hdf5_file

   subroutine save_surface_results_hdf5_file(res_hdf5_fileName,numNodeScalarFields2save,nodeScalarFields2save,&
      numNodeVectorFields2save,nodeVectorFields2save,&
      numElemGpScalarFields2save,elemGpScalarFields2save)
      implicit none
      character(512),intent(in) :: res_hdf5_fileName
      integer(4),intent(in) :: numNodeScalarFields2save,numNodeVectorFields2save,numElemGpScalarFields2save
      type(ptr_array1d_rp_save),intent(in) :: nodeScalarFields2save(:)
      type(ptr_array2d_rp_save),intent(in) :: nodeVectorFields2save(:),elemGpScalarFields2save(:)
      integer(4) :: h5err

      integer(hid_t) :: hdf5_fileId
      character(512) :: surf_res_hdf5_fileName,groupname,dsetname
      integer(4) :: iField

      call nvtxStartRange('save_surface_results_hdf5_file: set_hdf5_surface_resultsFile_name')
      call set_hdf5_surface_resultsFile_name(surf_res_hdf5_fileName,res_hdf5_fileName)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_results_hdf5_file: create_hdf5_file')
      call create_hdf5_file(surf_res_hdf5_fileName,hdf5_fileId)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_results_hdf5_file: set_vtkhdf_attributes_and_basic_groups')
      call set_vtkhdf_attributes_and_basic_groups(hdf5_fileId)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_results_hdf5_file: CPU-DATA')
      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/Points',hdf5_fileId,'/VTKHDF/Points',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/NumberOfPoints',hdf5_fileId,'/VTKHDF/NumberOfPoints',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/NumberOfCells',hdf5_fileId,'/VTKHDF/NumberOfCells',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/NumberOfConnectivityIds',hdf5_fileId,'/VTKHDF/NumberOfConnectivityIds',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/Offsets',hdf5_fileId,'/VTKHDF/Offsets',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/Connectivity',hdf5_fileId,'/VTKHDF/Connectivity',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/Types',hdf5_fileId,'/VTKHDF/Types',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/CellData/mpi_rank',hdf5_fileId,'/VTKHDF/CellData/mpi_rank',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/CellData/boundCode',hdf5_fileId,'/VTKHDF/CellData/boundCode',h5err)

      groupname = '/VTKHDF/PointData/'
      
      do iField=1,numNodeScalarFields2save
         dsetname = trim(adjustl(groupname))//trim(nodeScalarFields2save(iField)%nameField)
         call h5lcreate_external_f(res_hdf5_fileName,dsetname,hdf5_fileId,dsetname,h5err)
      end do
      
      do iField=1,numNodeVectorFields2save
         dsetname = trim(adjustl(groupname))//trim(nodeVectorFields2save(iField)%nameField)
         call h5lcreate_external_f(res_hdf5_fileName,dsetname,hdf5_fileId,dsetname,h5err)
      end do
      
      do iField=1,numElemGpScalarFields2save
         dsetname = trim(adjustl(groupname))//trim(elemGpScalarFields2save(iField)%nameField)
         call h5lcreate_external_f(res_hdf5_fileName,dsetname,hdf5_fileId,dsetname,h5err)
      end do
      call nvtxEndRange()

      call nvtxStartRange('save_surface_results_hdf5_file: close_hdf5_file')
      call close_hdf5_file(hdf5_fileId)
      call nvtxEndRange()

   end subroutine save_surface_results_hdf5_file
   !-------------------------------------------------------------------------------------------------------------------------------

   subroutine save_surface_mesh_hdf5_file(mnpbou,gmsh2ij,vtk2ij)
      implicit none
      integer(4),intent(in) :: mnpbou,gmsh2ij(mnpbou),vtk2ij(mnpbou)
      integer(4) :: ds_rank,h5err
      integer(hsize_t),dimension(1) :: ds_dims,ms_dims
      integer(hssize_t),dimension(1) :: ms_offset

      integer(hid_t) :: hdf5_fileId
      character(512) :: groupname,dsetname
      integer(hid_t) :: dtype
      integer(1),allocatable :: aux_array_i1(:)
      integer(8),allocatable :: aux_array_i8(:)
      integer(4), dimension(0:mpi_size-1) :: vecNumBoundsRankPar
      integer(4) :: mpiRankBoundStart
      integer(4) :: ii,jj,iBound,iRank,iNodeL

      call nvtxStartRange('save_surface_mesh_hdf5_file: set_hdf5_surface_meshFile_name')
      call set_hdf5_surface_meshFile_name()
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: create_hdf5_file')
      call create_hdf5_file(surface_meshFile_h5_name,hdf5_fileId)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: set_vtkhdf_attributes_and_basic_groups')
      call set_vtkhdf_attributes_and_basic_groups(hdf5_fileId)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-DATA')
      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/Points',hdf5_fileId,'/VTKHDF/Points',h5err)

      ds_rank = 1
      ds_dims(1) = int(mpi_size,hsize_t)
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      dtype = h5_datatype_int8
      allocate(aux_array_i8(1))
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: create_dataspace_hdf5')
      dsetname = '/VTKHDF/NumberOfPoints'
      aux_array_i8(1) = numNodesRankPar
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()  

      call nvtxStartRange('save_surface_mesh_hdf5_file: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: create_dataspace_hdf5')
      dsetname = '/VTKHDF/NumberOfCells'
      aux_array_i8(1) = numBoundsRankPar
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: create_dataspace_hdf5')
      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      aux_array_i8(1) = numBoundsRankPar*mnpbou
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-DATA')
      deallocate(aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: MPI_Allgather')
      call MPI_Allgather(numBoundsRankPar,1,mpi_datatype_int4,vecNumBoundsRankPar,1,mpi_datatype_int4,app_comm,mpi_err)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-COMPUTE')
      mpiRankBoundStart = 1
      do iRank=1,mpi_rank
         mpiRankBoundStart = mpiRankBoundStart +vecNumBoundsRankPar(iRank-1)
      end do
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-DATA')
      allocate(aux_array_i8(numBoundsRankPar+1))
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-COMPUTE')
      ms_dims(1)   = int(numBoundsRankPar,hsize_t)  + 1
      ms_offset(1) = int(mpiRankBoundStart,hssize_t) - 1 + int(mpi_rank,hssize_t)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: create_dataspace_hdf5')
      dsetname = '/VTKHDF/Offsets'
      ds_rank    = 1
      ds_dims(1) = int(totalNumBoundsSrl,hsize_t) + int(mpi_size,hsize_t)
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-COMPUTE')
      aux_array_i8(1) = 0
      do iBound = 2,(numBoundsRankPar+1)
         aux_array_i8(iBound) = aux_array_i8(iBound-1)+mnpbou
      end do
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-DATA')
      deallocate(aux_array_i8)
      allocate(aux_array_i8(numBoundsRankPar*mnpbou))
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: create_dataspace_hdf5')
      dsetname   = '/VTKHDF/Connectivity'
      ds_rank    = 1
      ds_dims(1) = int(totalNumBoundsSrl,hsize_t)  * int(mnpbou,hsize_t)
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-COMPUTE')
      ms_dims(1)   = int(numBoundsRankPar,hsize_t)   * int(mnpbou,hsize_t)
      ms_offset(1) = int((mpiRankBoundStart-1),hssize_t)* int(mnpbou,hssize_t)

      do iBound=1,numBoundsRankPar
         do ii=1,mnpbou
            iNodeL = boundParOrig(iBound,gmsh2ij(ii))
            jj = (iBound-1)*mnpbou + vtk2ij(ii)
            aux_array_i8(jj) = iNodeL - 1
         end do
      end do
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-DATA')
      deallocate(aux_array_i8)
      allocate(aux_array_i1(numBoundsRankPar))
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: create_dataspace_hdf5')
      dsetname = '/VTKHDF/Types'
      ds_rank    = 1
      ds_dims(1) = int(totalNumBoundsSrl,hsize_t)
      dtype      = h5_datatype_uint1
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-DATA')
      aux_array_i1(:) = 70
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: write_dataspace_1d_uint1_hyperslab_parallel')
      ms_dims(1)   = int(numBoundsRankPar ,hsize_t)
      ms_offset(1) = int(mpiRankBoundStart,hssize_t) - 1
      call write_dataspace_1d_uint1_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i1)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: create_dataspace_hdf5')
      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-DATA')
      do iBound = 1,numBoundsRankPar
         aux_array_i1(iBound) = mpi_rank
      end do
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: write_dataspace_1d_uint1_hyperslab_parallel')
      call write_dataspace_1d_uint1_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i1)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: create_dataspace_hdf5')
      ! ## bou_codes ##
      dsetname = '/VTKHDF/CellData/boundCode'
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-DATA')
      do iBound = 1,numBoundsRankPar
         aux_array_i1(iBound) = bouCodesPar(iBound)
      end do
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: write_dataspace_1d_uint1_hyperslab_parallel')
      call write_dataspace_1d_uint1_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i1)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: CPU-DATA')
      deallocate(aux_array_i1)
      call nvtxEndRange()

      call nvtxStartRange('save_surface_mesh_hdf5_file: close_hdf5_file')
      call close_hdf5_file(hdf5_fileId)
      call nvtxEndRange()

   end subroutine save_surface_mesh_hdf5_file

!---------------------------------------------------------------------------------------------------------
!        VTKHDF5 FILE
!---------------------------------------------------------------------------------------------------------

   subroutine set_vtkhdf_resultsFile_name(iStep,full_fileName)
      implicit none
      integer, intent(in) :: iStep
      character(len=*), intent(out) :: full_fileName
      character(len=12) :: aux_step

      call nvtxStartRange('set_vtkhdf_resultsFile_name: CPU-DATA')
      write(aux_step,'(I0)') iStep
      full_fileName = trim(adjustl(base_resultsFile_h5_name))//trim(aux_step)//'-vtk.hdf'
      call nvtxEndRange()
   end subroutine set_vtkhdf_resultsFile_name

   subroutine set_vtkhdf_avgResultsFile_name(iStep,full_fileName)
      implicit none
      integer, intent(in) :: iStep
      character(len=*), intent(out) :: full_fileName
      character(len=12) :: aux_step

      call nvtxStartRange('set_vtkhdf_avgResultsFile_name: CPU-DATA')
      write(aux_step,'(I0)') iStep
      full_fileName = trim(adjustl(base_avgResultsFile_h5_name))//trim(aux_step)//'-vtk.hdf'
      call nvtxEndRange()
   end subroutine set_vtkhdf_avgResultsFile_name

   subroutine set_vtkhdf_finalAvgResultsFile_name(full_fileName)
      implicit none
      character(len=*), intent(out) :: full_fileName

      call nvtxStartRange('set_vtkhdf_finalAvgResultsFile_name: CPU-DATA')
      full_fileName = trim(adjustl(base_avgResultsFile_h5_name))//'final-vtk.hdf'
      call nvtxEndRange()

   end subroutine set_vtkhdf_finalAvgResultsFile_name

   subroutine set_vtkhdf_attributes_and_basic_groups(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(512) :: groupname
      integer(hid_t) :: group_id,aspace_id,attr_id,atype
      integer(size_t) :: attr_length
      integer(hsize_t) :: a_dims(1)
      integer(4) :: a_rank,h5err
      character(16) :: attr_value
      integer(4),allocatable :: aux_array_i4(:)

      call nvtxStartRange('set_vtkhdf_attributes_and_basic_groups: CPU-DATA')
      groupname = '/VTKHDF'
      call h5gcreate_f(file_id,groupname,group_id,h5err)

      ! Attribute 'Type'
      call h5screate_f(H5S_SCALAR_F,aspace_id,h5err)

      attr_value = "UnstructuredGrid"
      attr_length = len_trim(attr_value)

      a_dims(1) = 1
      call h5tcopy_f(H5T_C_S1,atype,h5err)
      call h5tset_size_f(atype, attr_length,h5err)
      call h5tset_strpad_f(atype,H5T_STR_NULLPAD_F,h5err)

      call h5acreate_f(group_id,'Type',atype,aspace_id,attr_id,h5err)

      call h5awrite_f(attr_id,atype,attr_value,a_dims,h5err)
      call h5tclose_f(atype,h5err)

      call h5aclose_f(attr_id, h5err)

      call h5sclose_f(aspace_id, h5err)

      ! Attribute 'Version'
      a_rank = 1
      a_dims(1) = 2
      call h5screate_simple_f(a_rank,a_dims,aspace_id,h5err)

      call h5acreate_f(group_id,'Version',h5_datatype_int4,aspace_id,attr_id,h5err)

      allocate(aux_array_i4(2))
      aux_array_i4(1) = 1
      aux_array_i4(2) = 0
      call h5awrite_f(attr_id,h5_datatype_int4,aux_array_i4,a_dims,h5err)
      deallocate(aux_array_i4)

      call h5aclose_f(attr_id,h5err)

      call h5sclose_f(aspace_id,h5err)

      call h5gclose_f(group_id, h5err)
      call nvtxEndRange()

      call nvtxStartRange('set_vtkhdf_attributes_and_basic_groups: create_group_hdf5')
      groupname = '/VTKHDF/CellData'
      call create_group_hdf5(file_id,groupname)

      groupname = '/VTKHDF/FieldData'
      call create_group_hdf5(file_id,groupname)

      groupname = '/VTKHDF/PointData'
      call create_group_hdf5(file_id,groupname)
      call nvtxEndRange()

   end subroutine

   subroutine create_vtkhdf_unstructuredGrid_meshFile(mnnode,file_id)
      implicit none
      integer(4),intent(in) :: mnnode
      integer(hid_t),intent(in) :: file_id
      integer(hid_t) :: dtype
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hsize_t), dimension(2) :: ds_dims2d,ms_dims2d
      integer(hssize_t), dimension(1) :: ms_offset
      integer(hssize_t), dimension(2) :: ms_offset2d
      integer(4) :: ds_rank,h5err
      character(512) :: dsetname
      integer(4) :: ii,iElemL
      integer(1),allocatable :: aux_array_i1(:)
      integer(8),allocatable :: aux_array_i8(:)

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: set_vtkhdf_attributes_and_basic_groups')
      call set_vtkhdf_attributes_and_basic_groups(file_id)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      ds_dims2d(1) = int(ndime,hsize_t)
      ds_dims2d(2) = int(totalNumNodesPar,hsize_t)
      ms_dims2d(1) = int(ndime,hsize_t)
      ms_dims2d(2) = int(numNodesRankPar,hsize_t)
      ms_offset2d(1) = 0
      ms_offset2d(2) = int(rankNodeStart,hssize_t)-1
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: save_array2D_tr_rp_in_dataset_hdf5_file')
      dsetname = '/VTKHDF/Points'
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,coordPar)
      call nvtxEndRange()
      
      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      dtype = h5_datatype_int8
      allocate(aux_array_i8(1))
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: create_dataspace_hdf5')
      dsetname = '/VTKHDF/NumberOfPoints'
      aux_array_i8(1) = numNodesRankPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: create_dataspace_hdf5')
      dsetname = '/VTKHDF/NumberOfCells'
      aux_array_i8(1) = numElemsRankPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: create_dataspace_hdf5')
      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      aux_array_i8(1) = numElemsRankPar*mnnode
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      deallocate(aux_array_i8)
      allocate(aux_array_i8(numElemsRankPar+1))
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: CPU-COMPUTE')
      ds_dims(1)   = int(totalNumElements,hsize_t) + int(mpi_size,hsize_t)
      ms_dims(1)   = int(numElemsRankPar,hsize_t)  + 1
      ms_offset(1) = int(rankElemStart,hssize_t)   - 1 + int(mpi_rank,hssize_t)

      dsetname = '/VTKHDF/Offsets'
      aux_array_i8(1) = 0
      do iElemL = 2,(numElemsRankPar+1)
         aux_array_i8(iElemL) = aux_array_i8(iElemL-1)+mnnode
      end do
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: create_dataspace_hdf5')
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      deallocate(aux_array_i8)
      allocate(aux_array_i8(numElemsRankPar*mnnode))
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: CPU-COMPUTE')
      ds_dims(1)   = int(totalNumElements,hsize_t)  * int(mnnode,hsize_t)
      ms_dims(1)   = int(numElemsRankPar,hsize_t)   * int(mnnode,hsize_t)
      ms_offset(1) = int((rankElemStart-1),hssize_t)* int(mnnode,hssize_t)

      do ii = 1,numElemsRankPar*mnnode
         aux_array_i8(ii) = connecVTK(ii)-1
      end do
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: create_dataspace_hdf5')
      dsetname = '/VTKHDF/Connectivity'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      deallocate(aux_array_i8)
      allocate(aux_array_i1(numElemsRankPar))
      dtype = h5_datatype_uint1
      dsetname = '/VTKHDF/Types'
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: GPU-DATA')
      !$acc kernels
      aux_array_i1(:) = 72
      !$acc end kernels
      call nvtxEndRange()

      ds_dims(1)   = int(totalNumElements,hsize_t)
      ms_dims(1)   = int(numElemsRankPar ,hsize_t)
      ms_offset(1) = int(rankElemStart,hssize_t) - 1

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: create_dataspace_hdf5')
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_uint1_hyperslab_parallel')
      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i1)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      deallocate(aux_array_i1)
      call nvtxEndRange()

   end subroutine create_vtkhdf_unstructuredGrid_meshFile

   subroutine create_groups_datasets_vtkhdf_unstructuredGrid_meshFile(mporder,mnnode,file_id,isLinealOutput,evalMeshQuality,numMshRanks2Part,numElemsGmsh,numNodesParTotal_i8,mnnodeVTK,numVTKElemsPerMshElem)
      implicit none
      integer(hid_t),intent(in) :: file_id
      logical,intent(in) :: isLinealOutput,evalMeshQuality
      integer(4),intent(in) :: mporder,mnnode,numMshRanks2Part,numElemsGmsh,mnnodeVTK,numVTKElemsPerMshElem
      integer(8),intent(in) :: numNodesParTotal_i8
      integer(hid_t) :: dtype
      integer(hsize_t) :: ds_dims(1),ds_dims2d(2),aux_ds_dims
      integer(4) :: ds_rank,h5err
      character(512) :: dsetname

      CALL nvtxStartRange('create_groups_datasets_vtkhdf_unstructuredGrid_meshFile: set_vtkhdf_attributes_and_basic_groups')
      call set_vtkhdf_attributes_and_basic_groups(file_id)
      CALL nvtxEndRange()

      CALL nvtxStartRange('create_groups_datasets_vtkhdf_unstructuredGrid_meshFile: select_dtype_rp')
      call select_dtype_rp(dtype)
      CALL nvtxEndRange()

      ds_rank = 2
      ds_dims2d(1) = ndime
      ds_dims2d(2) = numNodesParTotal_i8

      call nvtxStartRange('create_groups_datasets_vtkhdf_unstructuredGrid_meshFile: create_dataspace_hdf5')
      dsetname = '/VTKHDF/Points'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims2d,dtype)

      ds_rank = 1
      ds_dims(1) = numMshRanks2Part
      dtype = h5_datatype_int8

      dsetname = '/VTKHDF/NumberOfPoints'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/VTKHDF/NumberOfCells'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = int(numElemsGmsh*numVTKElemsPerMshElem+numMshRanks2Part,hsize_t)

      dsetname = '/VTKHDF/Offsets'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = int(numElemsGmsh*numVTKElemsPerMshElem,hsize_t)*int(mnnodeVTK,hsize_t)

      dsetname = '/VTKHDF/Connectivity'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = int(numElemsGmsh*numVTKElemsPerMshElem,hsize_t)

      dtype = h5_datatype_uint1

      dsetname = '/VTKHDF/Types'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      if (evalMeshQuality) then
         call nvtxStartRange('create_groups_datasets_vtkhdf_unstructuredGrid_meshFile: select_dtype_rp')
         call select_dtype_rp(dtype)
         call nvtxEndRange()

         call nvtxStartRange('create_groups_datasets_vtkhdf_unstructuredGrid_meshFile: create_dataspace_hdf5')
         dsetname = '/VTKHDF/CellData/mesh_quality'
         call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
         call nvtxEndRange()
      end if

   end subroutine create_groups_datasets_vtkhdf_unstructuredGrid_meshFile

   subroutine write_mshRank_data_vtkhdf_unstructuredGrid_meshFile(mporder,mnnode,file_id,eval_mesh_quality,mshRank,numMshRanks2Part,numElemsMshRank,numElemsVTKMshRank,sizeConnecVTKMshRank,mnnodeVTK,numVTKElemsPerMshElem,mshRankElemStart,mshRankElemEnd,mshRankNodeStart_i8,mshRankNodeEnd_i8,numNodesMshRank,coordVTKMshRank,connecVTKMshRank,quality)
      implicit none
      integer(hid_t),intent(in) :: file_id
      logical, intent(in) :: eval_mesh_quality
      integer(4),intent(in) :: mporder,mnnode,mshRank,numMshRanks2Part
      integer(4),intent(in) :: numElemsMshRank,numElemsVTKMshRank,sizeConnecVTKMshRank,mnnodeVTK,numVTKElemsPerMshElem,mshRankElemStart,mshRankElemEnd
      integer(8),intent(in) :: mshRankNodeStart_i8,mshRankNodeEnd_i8
      integer(4),intent(in) :: numNodesMshRank
      real(rp),intent(in)   :: coordVTKMshRank(numNodesMshRank,3), quality(numVTKElemsPerMshElem)
      integer(4),intent(in) :: connecVTKMshRank(sizeConnecVTKMshRank)

      integer(hsize_t) :: ms_dims(1),ms_dims2d(2),aux_ms_dims
      integer(hssize_t) :: ms_offset(1),ms_offset2d(2),aux_ms_offset
      character(512) :: dsetname
      integer(4) :: ii,iElemL
      integer(1),allocatable :: aux_array_i1(:)
      integer(8),allocatable :: aux_array_i8(:)

      ms_dims2d(1) = ndime
      ms_dims2d(2) = numNodesMshRank
      ms_offset2d(1) = 0
      ms_offset2d(2) = mshRankNodeStart_i8-1

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: write_dataspace_2d_tr_real_rp_hyperslab_parallel')
      dsetname = '/VTKHDF/Points'
      call write_dataspace_2d_tr_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,coordVTKMshRank)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      ms_dims(1) = 1
      ms_offset(1) = int(mshRank,hssize_t)
      allocate(aux_array_i8(1))
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_int8_hyperslab_parallel')
      dsetname = '/VTKHDF/NumberOfPoints'
      aux_array_i8(1) = int(numNodesMshRank,8)
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

      dsetname = '/VTKHDF/NumberOfCells'
      aux_array_i8(1) = int(numElemsVTKMshRank,8)
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      aux_array_i8(1) = sizeConnecVTKMshRank
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      deallocate(aux_array_i8)
      dsetname = '/VTKHDF/Offsets'
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: CPU-COMPUTE')
      aux_ms_dims   = int(numElemsVTKMshRank + 1,hsize_t)
      aux_ms_offset = int((mshRankElemStart-1)*(numVTKElemsPerMshElem) + mshRank,hssize_t)

      allocate(aux_array_i8(aux_ms_dims))
      aux_array_i8(1) = 0
      do iElemL = 2,aux_ms_dims
         aux_array_i8(iElemL) = aux_array_i8(iElemL-1)+int(mnnodeVTK,8)
      end do

      ms_dims(1)   = aux_ms_dims
      ms_offset(1) = aux_ms_offset
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      deallocate(aux_array_i8)
      allocate(aux_array_i8(sizeConnecVTKMshRank))
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: CPU-COMPUTE')
      dsetname = '/VTKHDF/Connectivity'

      aux_ms_dims   = int(sizeConnecVTKMshRank,hsize_t)
      aux_ms_offset = int((mshRankElemStart-1)*numVTKElemsPerMshElem,hssize_t)*int(mnnodeVTK,hssize_t)

      ms_dims(1)   = aux_ms_dims
      ms_offset(1) = aux_ms_offset

      do ii = 1,sizeConnecVTKMshRank
         aux_array_i8(ii) = connecVTKMshRank(ii)-1
      end do
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      deallocate(aux_array_i8)

      dsetname = '/VTKHDF/Types'

      aux_ms_dims   = int(numElemsMshRank*numVTKElemsPerMshElem,hsize_t)
      aux_ms_offset = int((mshRankElemStart-1)*numVTKElemsPerMshElem,hssize_t)

      allocate(aux_array_i1(aux_ms_dims))
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: GPU-DATA')
      !$acc kernels
      aux_array_i1(:) = 72
      !$acc end kernels
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_uint1_hyperslab_parallel')
      ms_dims(1)   = aux_ms_dims
      ms_offset(1) = aux_ms_offset
      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i1)
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      do iElemL = 1,numElemsVTKMshRank
         aux_array_i1(iElemL) = mshRank
      end do
      call nvtxEndRange()

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_uint1_hyperslab_parallel')
      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i1)
      call nvtxEndRange()

      !!! Save mesh quality
      if (eval_mesh_quality) then
         call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_real_rp_hyperslab_parallel')
         dsetname = '/VTKHDF/CellData/mesh_quality'
         call write_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,quality)
         call nvtxEndRange()
      end if

      call nvtxStartRange('write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      deallocate(aux_array_i1)
      call nvtxEndRange()

   end subroutine write_mshRank_data_vtkhdf_unstructuredGrid_meshFile

   subroutine dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile(file_id,eval_mesh_quality,numMshRanks2Part)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(4),intent(in) :: numMshRanks2Part
      logical,intent(in) :: eval_mesh_quality
      integer(hsize_t), dimension(1) :: ms_dims
      integer(hsize_t), dimension(2) :: ms_dims2d
      integer(hssize_t), dimension(1) :: ms_offset
      integer(hssize_t), dimension(2) :: ms_offset2d
      character(512) :: dsetname
      integer(1),allocatable :: empty_array_i1(:)
      integer(8),allocatable :: empty_array_i8(:)
      real(rp),allocatable :: empty_array2d_rp(:,:)
      real(rp),allocatable :: empty_array1d_rp(:)

      call nvtxStartRange('dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      allocate(empty_array_i1(0))
      allocate(empty_array_i8(0))
      allocate(empty_array2d_rp(0,0))
      allocate(empty_array1d_rp(0))

      ms_dims2d(1) = 0
      ms_dims2d(2) = 0
      ms_offset2d(1) = 0
      ms_offset2d(2) = 0
      call nvtxEndRange()

      call nvtxStartRange('dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: write_dataspace_2d_tr_real_rp_hyperslab_parallel')
      dsetname = '/VTKHDF/Points'
      call write_dataspace_2d_tr_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,empty_array2d_rp)
      call nvtxEndRange()

      call nvtxStartRange('dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_int8_hyperslab_parallel')
      ms_dims(1) = 0
      ms_offset(1) = 0
      dsetname = '/VTKHDF/NumberOfPoints'
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

      dsetname = '/VTKHDF/NumberOfCells'
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

      dsetname = '/VTKHDF/Offsets'
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

      dsetname = '/VTKHDF/Connectivity'
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_uint1_hyperslab_parallel')
      dsetname = '/VTKHDF/Types'
      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i1)

      dsetname = '/VTKHDF/CellData/mpi_rank'
      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i1)
      call nvtxEndRange()

      if (eval_mesh_quality) then
         call nvtxStartRange('dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: write_dataspace_1d_real_rp_hyperslab_parallel')
         dsetname = '/VTKHDF/CellData/mesh_quality'
         call write_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array1d_rp)
         call nvtxEndRange()
      end if

      call nvtxStartRange('dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile: CPU-DATA')
      deallocate(empty_array_i1)
      deallocate(empty_array_i8)
      deallocate(empty_array2d_rp)
      deallocate(empty_array1d_rp)
      call nvtxEndRange()

   end subroutine dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile

   subroutine create_vtkhdf_unstructuredGrid_struct_for_resultsFile(mnnode,file_id)
      implicit none
      integer(4),intent(in) :: mnnode
      integer(hid_t),intent(in) :: file_id
      integer(hid_t) :: dtype
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hssize_t), dimension(1) :: ms_offset
      integer(4) :: ds_rank,h5err
      character(512) :: dsetname
      integer(8),allocatable :: aux_array_i8(:)

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_struct_for_resultsFile: set_vtkhdf_attributes_and_basic_groups')
      call set_vtkhdf_attributes_and_basic_groups(file_id)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_struct_for_resultsFile: CPU-DATA')
      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/Points',file_id,'/VTKHDF/Points',h5err)

      ds_rank = 1
      ds_dims(1) = int(mpi_size,hsize_t)
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      dtype = h5_datatype_int8
      allocate(aux_array_i8(1))
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_struct_for_resultsFile: create_dataspace_hdf5')
      dsetname = '/VTKHDF/NumberOfPoints'
      aux_array_i8(1) = int(numNodesRankPar,8)
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_struct_for_resultsFile: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_struct_for_resultsFile: create_dataspace_hdf5')
      dsetname = '/VTKHDF/NumberOfCells'
      aux_array_i8(1) = int(numElemsVTKRankPar,8)!numElemsRankPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_struct_for_resultsFile: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_struct_for_resultsFile: create_dataspace_hdf5')
      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      aux_array_i8(1) = int(sizeConnecVTKRankPar,8)!numElemsRankPar*mnnode
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_struct_for_resultsFile: write_dataspace_1d_int8_hyperslab_parallel')
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      call nvtxEndRange()

      call nvtxStartRange('create_vtkhdf_unstructuredGrid_struct_for_resultsFile: CPU-DATA')
      deallocate(aux_array_i8)
      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/Offsets',file_id,'/VTKHDF/Offsets',h5err)
      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/Connectivity',file_id,'/VTKHDF/Connectivity',h5err)
      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/Types',file_id,'/VTKHDF/Types',h5err)
      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/CellData/mpi_rank',file_id,'/VTKHDF/CellData/mpi_rank',h5err)
      call nvtxEndRange()

   end subroutine create_vtkhdf_unstructuredGrid_struct_for_resultsFile

   subroutine write_vtkhdf_realFieldFile(mnnode,full_fileName,realField)
      implicit none
      integer(4),intent(in) :: mnnode
      character(512),intent(in) :: full_fileName
      real(rp),dimension(numNodesRankPar),intent(in) :: realField

      character(512) :: dsetname
      integer(hid_t) :: file_id,plist_id,dset_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer :: ds_rank,h5err

      integer(4) :: iElemL
      integer(1),allocatable :: aux_array_i1(:)

      call nvtxStartRange('write_vtkhdf_realFieldFile: CPU-DATA')
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create VTKHDF ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)
      call nvtxEndRange()

      call nvtxStartRange('write_vtkhdf_realFieldFile: create_vtkhdf_unstructuredGrid_meshFile')
      call create_vtkhdf_unstructuredGrid_meshFile(mnnode,file_id)
      call nvtxEndRange()

      ds_rank = 1
      ds_dims(1)   = int(totalNumNodesPar,hsize_t)
      ms_dims(1)   = int(numNodesRankPar ,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      call nvtxStartRange('write_vtkhdf_realFieldFile: save_array1D_rp_in_dataset_hdf5_file')
      ! ## realField ##
      dsetname = '/VTKHDF/PointData/realField'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,realField)
      call nvtxEndRange()

      call nvtxStartRange('write_vtkhdf_realFieldFile: CPU-DATA')
      allocate(aux_array_i1(numNodesRankPar))
      dtype =  h5_datatype_uint1
      ds_rank = 1
      ds_dims(1)   = int(totalNumElements,hsize_t)
      ms_dims(1)   = int(numElemsRankPar ,hsize_t)
      ms_offset(1) = int(rankElemStart,hssize_t)-1
      call nvtxEndRange()
      
      call nvtxStartRange('write_vtkhdf_realFieldFile: create_dataspace_hdf5')
      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('write_vtkhdf_realFieldFile: CPU-DATA')
      do iElemL = 1,numElemsRankPar
         aux_array_i1(iElemL) = mpi_rank
      end do
      call nvtxEndRange()

      call nvtxStartRange('write_vtkhdf_realFieldFile: write_dataspace_1d_uint1_hyperslab_parallel')
      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i1)
      call nvtxEndRange()

      call nvtxStartRange('write_vtkhdf_realFieldFile: CPU-DATA')
      deallocate(aux_array_i1)
      call h5fclose_f(file_id,h5err)
      call nvtxEndRange()

   end subroutine write_vtkhdf_realFieldFile


   subroutine save_vtkhdf_realFieldFile(mnnode,realField)
      implicit none
      integer(4),intent(in) :: mnnode
      real(rp),dimension(numNodesRankPar),intent(in) :: realField
      character(512) :: full_fileName

      call nvtxStartRange('save_vtkhdf_realFieldFile: set_vtkhdf_resultsFile_name')
      call set_vtkhdf_resultsFile_name(0,full_fileName)
      call nvtxEndRange()
      
      if(mpi_rank.eq.0) write(*,*) '# Saving VTKHDF realField file: ',trim(adjustl(full_fileName))

      call nvtxStartRange('save_vtkhdf_realFieldFile: write_vtkhdf_realFieldFile')
      call write_vtkhdf_realFieldFile(mnnode,full_fileName,realField)
      call nvtxEndRange()

   end subroutine save_vtkhdf_realFieldFile

   subroutine get_dims_hdf5(file_id,dsetname, ms_rank, fs_dims, fs_maxdims)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer, intent(in) :: ms_rank
      integer(hid_t) :: dset_id,fspace_id
      integer :: h5err
      integer(hsize_t), intent(out) :: fs_dims(ms_rank),fs_maxdims(ms_rank)
      integer(hid_t) :: dtype
      dtype = H5T_NATIVE_INTEGER

      call nvtxStartRange('get_dims_hdf5: CPU-DATA')
      call h5dopen_f(file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
      call nvtxEndRange()

   end subroutine get_dims_hdf5

!-------------------------------WITNESS POINTS-------------------------------!
   subroutine create_witness_hdf5(full_fileName, mnnode, xyz, witel, witxi, shapewit, nwit, nwitPar, witGlob, save_u_i, save_pr, save_rho)
      implicit none
      character(512),intent(in)  :: full_fileName
      integer(4),intent(in)      :: mnnode, nwit, nwitPar
      integer(4),intent(in)      :: witel(nwit), witGlob(nwit)
      real(rp),intent(in)        :: witxi(nwit, ndime), shapewit(nwit,mnnode)
      real(rp),intent(in)        :: xyz(nwit,ndime)
      logical,intent(in)         :: save_u_i, save_pr, save_rho
      integer(4)                 :: aux(1), nwitParAllRanks(mpi_size), nwitOffset=0, inode
      integer(hid_t)             :: file_id,plist_id,dset_id,dspace_id,group_id, dtype
      integer(hsize_t)           :: ds_dims(1),ms_dims(1),max_dims(1),chunk_dims(1)
      integer(hssize_t)          :: ms_offset(1)
      integer(hsize_t)           :: ds_dims2d(2),ms_dims2d(2),max_dims2d(2),chunk_dims2d(2)
      integer(hssize_t)          :: ms_offset2d(2)
      integer(4)                 :: ds_rank,h5err,irank,iwit
      character(256)             :: dsetname
      real(rp)                   :: auxwitxyz(nwitPar, ndime), auxwitxi(nwitPar,ndime), auxshapefunc(nwitPar,mnnode)

      call nvtxStartRange('create_witness_hdf5: CPU-DATA')
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create results file ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      !Create dataspece for nwitPar and save it!
      dtype        = h5_datatype_int4
      ds_rank      = 1
      dsetname     = 'nwitPar'
      ds_dims(1)   = mpi_size
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      aux(1)       = nwitPar
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: create_dataspace_hdf5')
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux)
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: MPI_Barrier')
      call MPI_Barrier(app_comm, mpi_err)
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      !Compute sum of nwitPar until that rank and save it on its dataspace!
      dtype        = h5_datatype_int4
      ms_dims(1)   = mpi_size
      ms_offset(1) = 0
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,nwitParAllRanks)
      call nvtxEndRange()
      
      call nvtxStartRange('create_witness_hdf5: CPU-COMPUTE')
      if (mpi_rank > 0) then
         do irank = 1, mpi_rank
            nwitOffset = nwitOffset + nwitParAllRanks(irank)
         end do
      end if
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: CPU-DATA')
      ds_rank      = 1
      dsetname     = 'nwitOffset'
      ds_dims(1)   = mpi_size
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      aux(1)       = nwitOffset
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: create_dataspace_hdf5')
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux)
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: MPI_Barrier')
      call MPI_Barrier(app_comm, mpi_err)
      call nvtxEndRange()

      !Create dataspece for global numeration and save it!
      call nvtxStartRange('create_witness_hdf5: CPU-DATA')
      dtype        = h5_datatype_int4
      dsetname     = 'global'
      ds_rank      = 1
      ds_dims(1)   = nwit
      ms_dims(1)   = nwitPar
      ms_offset(1) = nwitOffset
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: create_dataspace_hdf5')
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,witGlob)
      call nvtxEndRange()

      !Create dataspece for element containing the witness and save it!
      call nvtxStartRange('create_witness_hdf5: CPU-DATA')
      dtype        = h5_datatype_int4
      dsetname     = 'element'
      ds_rank      = 1
      ds_dims(1)   = nwit
      ms_dims(1)   = nwitPar
      ms_offset(1) = nwitOffset
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: create_dataspace_hdf5')
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: write_dataspace_1d_int4_hyperslab_parallel')
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,witel)
      call nvtxEndRange()

      !Create dataspace for witness coordinates and save them!
      call nvtxStartRange('create_witness_hdf5: CPU-DATA')
      do iwit = 1, nwitPar
         auxwitxyz(iwit,:) = xyz(iwit,:)
         auxwitxi(iwit,:)  = witxi(iwit,:)
         auxshapefunc(iwit,:) = shapewit(iwit,:)
      end do
      call nvtxEndRange()
      
      call nvtxStartRange('create_witness_hdf5: save_array2D_tr_rp_in_dataset_hdf5_file')
      dsetname       = 'xyz'
      ds_dims2d(1)   = ndime
      ds_dims2d(2)   = nwit
      ms_dims2d(1)   = ndime
      ms_dims2d(2)   = nwitPar
      ms_offset2d(1) = 0
      ms_offset2d(2) = nwitOffset
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwitxyz)

      !Create dataspace for witness isoparametric coordinates and save them!
      dsetname       = 'witxi'
      ds_dims2d(1)   = ndime
      ds_dims2d(2)   = nwit
      ms_dims2d(1)   = ndime
      ms_dims2d(2)   = nwitPar
      ms_offset2d(1) = 0
      ms_offset2d(2) = nwitOffset
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwitxi)

      !Create dataspace for the shape functions evaluated on the witness points and save them!
      dsetname       = 'shape_functions'
      ds_dims2d(1)   = mnnode
      ds_dims2d(2)   = nwit
      ms_dims2d(1)   = mnnode
      ms_dims2d(2)   = nwitPar
      ms_offset2d(1) = 0
      ms_offset2d(2) = nwitOffset
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxshapefunc)
      call nvtxEndRange()

      !Create time dataset!
      call nvtxStartRange('create_witness_hdf5: select_dtype_rp_vtk')
      call select_dtype_rp_vtk(dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: create_dataspace_maxdims_hdf5')
      ds_rank       = 1
      dsetname      = 'time'
      ds_dims(1)    = 1
      max_dims(1)   = H5S_UNLIMITED_F
      chunk_dims(1) = 1
      call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims,chunk_dims,dtype)

      !Create istep dataset!
      dtype         = h5_datatype_int4
      ds_rank       = 1
      dsetname      = 'istep'
      ds_dims(1)    = 1
      max_dims(1)   = H5S_UNLIMITED_F
      chunk_dims(1) = 1
      call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims,chunk_dims,dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: CPU-DATA')
      call h5fclose_f(file_id,h5err)
      call nvtxEndRange()

      !Create dataspaces for the magnitudes to save!
      call nvtxStartRange('create_witness_hdf5: select_dtype_rp_vtk')
      call select_dtype_rp_vtk(dtype)
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: create_dataspace_maxdims_hdf5')
      ds_rank         = 2
      ds_dims2d(1)    = 1
      ds_dims2d(2)    = nwit
      max_dims2d(1)   = H5S_UNLIMITED_F
      max_dims2d(2)   = nwit
      chunk_dims2d(1) = 1
      chunk_dims2d(2) = nwit
      if (save_u_i) then
         dsetname = 'u_x'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims2d,max_dims2d,chunk_dims2d,dtype)

         dsetname = 'u_y'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims2d,max_dims2d,chunk_dims2d,dtype)

         dsetname = 'u_z'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims2d,max_dims2d,chunk_dims2d,dtype)
      end if

      if (save_pr) then
         dsetname = 'pr'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims2d,max_dims2d,chunk_dims2d,dtype)
      end if

      if (save_rho) then
         dsetname = 'rho'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims2d,max_dims2d,chunk_dims2d,dtype)
      end if
      call nvtxEndRange()

      call nvtxStartRange('create_witness_hdf5: CPU-DATA')
      call h5fclose_f(file_id,h5err)
      call nvtxEndRange()

   end subroutine create_witness_hdf5

   subroutine load_witness_hdf5(full_fileName, mnnode, nwit, loadstep, load_stepwit, nwitPar, witel, witxi, shapefunc)
      implicit none
      character(512), intent(in)  :: full_fileName
      integer(4),     intent(in)  :: mnnode, nwit, loadstep
      integer(4),     intent(out) :: witel(nwit)
      real(rp),       intent(out) :: witxi(nwit,ndime), shapefunc(nwit,mnnode)
      integer(4),     intent(out) :: nwitPar, load_stepwit
      integer(hid_t)              :: file_id,plist_id,dset_id,dspace_id,group_id
      integer(hsize_t)            :: ds_dims(1),ms_dims(1),max_dims(1),chunk_dims(1)
      integer(hssize_t)           :: ms_offset(1)
      integer(hsize_t)            :: ds_dims2d(2),ms_dims2d(2),max_dims2d(2),chunk_dims2d(2)
      integer(hssize_t)           :: ms_offset2d(2)
      integer                     :: ms_rank, h5err, iwit, istep
      integer(hsize_t)            :: nsteps(2), maxnsteps(2)
      character(256)              :: dsetname
      real(rp), allocatable       :: auxwitxi(:,:), auxshapefunc(:,:)
      integer(4)                  :: nwitOffset, auxread(1)
      integer(4), allocatable     :: steps(:)

      call nvtxStartRange('load_witness_hdf5: CPU-DATA')
      witel(:)   = 0
      witxi(:,:) = 0.0_rp
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_fileName, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load results file ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)
      call nvtxEndRange()


      !Read step to continue!
      call nvtxStartRange('load_witness_hdf5: get_dims_hdf5')
      dsetname     = 'istep'
      ms_rank      = 1
      call get_dims_hdf5(file_id, dsetname, ms_rank, nsteps, maxnsteps)
      call nvtxEndRange()

      call nvtxStartRange('load_witness_hdf5: CPU-DATA')
      ms_dims(1)   = nsteps(1)
      ms_offset(1) = 0
      allocate(steps(nsteps(1)))
      call nvtxEndRange()

      call nvtxStartRange('load_witness_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,steps)
      call nvtxEndRange()

      call nvtxStartRange('load_witness_hdf5: CPU-COMPUTE')
      do istep = nsteps(1), 1, -1
         if (steps(istep) < loadstep) then
            load_stepwit = istep+1
            exit
         end if
      end do
      call nvtxEndRange()

      call nvtxStartRange('load_witness_hdf5: CPU-DATA')
      deallocate(steps)
      call nvtxEndRange()

      !Read nwitPar!
      call nvtxStartRange('load_witness_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname     = 'nwitPar'
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,auxread)
      call nvtxEndRange()

      call nvtxStartRange('load_witness_hdf5: CPU-DATA')
      nwitPar = auxread(1)
      allocate(auxwitxi(nwitPar,ndime))
      allocate(auxshapefunc(nwitPar,mnnode))
      call nvtxEndRange()

      !Read nwitOffset!
      call nvtxStartRange('load_witness_hdf5: read_dataspace_1d_int4_hyperslab_parallel')
      dsetname     = 'nwitOffset'
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,auxread)
      nwitOffset = auxread(1)

      !Read elements containing the witness points
      dsetname     = 'element'
      ms_dims(1)   = nwitPar
      ms_offset(1) = nwitOffset
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,witel)
      call nvtxEndRange()

      !Read witness isoparametric coordinates!
      call nvtxStartRange('load_witness_hdf5: read_array2D_tr_rp_in_dataset_hdf5_file')
      dsetname       = 'witxi'
      ms_dims2d(1)   = ndime
      ms_dims2d(2)   = nwitPar
      ms_offset2d(1) = 0
      ms_offset2d(2) = nwitOffset
      call read_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims2d,ms_offset2d,auxwitxi)

      !Read the shape functions coordinates!
      dsetname       = 'shape_functions'
      ms_dims2d(1)   = mnnode
      ms_dims2d(2)   = nwitPar
      ms_offset2d(1) = 0
      ms_offset2d(2) = nwitOffset
      call read_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims2d,ms_offset2d,auxshapefunc)
      call nvtxEndRange()

      call nvtxStartRange('load_witness_hdf5: CPU-DATA')
      do iwit = 1,nwitPar
         witxi(iwit,:) = auxwitxi(iwit,:)
         shapefunc(iwit,:) = auxshapefunc(iwit,:)
      end do

      deallocate(auxwitxi)
      deallocate(auxshapefunc)

      call h5fclose_f(file_id,h5err)
      call nvtxEndRange()
   end subroutine load_witness_hdf5

   subroutine update_witness_hdf5(itewit, leapwitsave, witval, nwit, nwitPar, nvarwit, full_fileName, t, steps, save_u_i, save_pr, save_rho)
      integer(4),     intent(in) :: itewit, nwit, nwitPar, nvarwit, leapwitsave, steps(leapwitsave)
      real(rp),       intent(in) :: witval(leapwitsave, nwitPar, nvarwit), t(leapwitsave)
      logical,        intent(in) :: save_u_i, save_pr, save_rho
      character(512), intent(in) :: full_fileName
      character(256)             :: dsetname
      real(rp)                   :: auxwrite(leapwitsave,nwitPar)
      integer(4)                 :: h5err, iwit, ds_rank, ileap
      integer(4)                 :: nwitOffset, auxread(1)
      integer(hsize_t)           :: ds_dims(1),ms_dims(1),max_dims(1),chunk_dims(1)
      integer(hssize_t)          :: ms_offset(1)
      integer(hsize_t)           :: ds_dims2d(2),ms_dims2d(2),max_dims2d(2),chunk_dims2d(2)
      integer(hssize_t)          :: ms_offset2d(2)
      integer(hid_t)             :: file_id,plist_id, dtype

      ! Setup file access property list with parallel I/O access.
      call nvtxStartRange('update_witness_hdf5: CPU-DATA')
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_fileName, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load results file ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)
      call nvtxEndRange()

      !Read nwitOffset!
      call nvtxStartRange('update_witness_hdf5: read_dataspace_1d_int4_hyperslab_parallel') 
      dsetname     = 'nwitOffset'
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,auxread)
      nwitOffset = auxread(1)
      call nvtxEndRange()

      !Save variables!
      ms_dims        = nwitPar
      ds_rank        = 2
      ds_dims2d(1)   = itewit
      ds_dims2d(2)   = nwit
      ms_dims2d(1)   = leapwitsave
      ms_dims2d(2)   = nwitPar
      ms_offset2d(1) = itewit-leapwitsave
      ms_offset2d(2) = nwitOffset
      if (save_u_i) then
         dsetname = 'u_x'
         call nvtxStartRange('update_witness_hdf5: GPU-DATA')
         !$acc kernels
         auxwrite(:,:) = witval(:,:,1)
         !$acc end kernels
         call nvtxEndRange()

         call nvtxStartRange('update_witness_hdf5: extend_dataset_hdf5')
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
         call nvtxEndRange()

         call nvtxStartRange('update_witness_hdf5: save_array2D_tr_rp_in_dataset_hdf5_file')
         call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwrite,isCreateDataspaceOpt=.False.)
         call nvtxEndRange()

         dsetname = 'u_y'
         call nvtxStartRange('update_witness_hdf5: GPU-DATA')
         !$acc kernels
         auxwrite(:,:) = witval(:,:,2)
         !$acc end kernels
         call nvtxEndRange()

         call nvtxStartRange('update_witness_hdf5: extend_dataset_hdf5')
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
         call nvtxEndRange()

         call nvtxStartRange('update_witness_hdf5: save_array2D_tr_rp_in_dataset_hdf5_file')
         call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwrite,isCreateDataspaceOpt=.False.)
         call nvtxEndRange()

         dsetname = 'u_z'
         call nvtxStartRange('update_witness_hdf5: GPU-DATA')
         !$acc kernels
         auxwrite(:,:) = witval(:,:,3)
         !$acc end kernels
         call nvtxEndRange()

         call nvtxStartRange('update_witness_hdf5: extend_dataset_hdf5')
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
         call nvtxEndRange()

         call nvtxStartRange('update_witness_hdf5: save_array2D_tr_rp_in_dataset_hdf5_file')
         call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwrite,isCreateDataspaceOpt=.False.)
         call nvtxEndRange()
      end if
      if (save_pr) then
         dsetname = 'pr'
         call nvtxStartRange('update_witness_hdf5: GPU-DATA')
         !$acc kernels
         auxwrite(:,:) = witval(:,:,4)
         !$acc end kernels
         call nvtxEndRange()

         call nvtxStartRange('update_witness_hdf5: extend_dataset_hdf5')
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
         call nvtxEndRange()

         call nvtxStartRange('update_witness_hdf5: save_array2D_tr_rp_in_dataset_hdf5_file')
         call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwrite,isCreateDataspaceOpt=.False.)
         call nvtxEndRange()
      end if
      if (save_rho) then
         dsetname = 'rho'
         call nvtxStartRange('update_witness_hdf5: GPU-DATA')
         !$acc kernels
         auxwrite(:,:) = witval(:,:,5)
         !$acc end kernels
         call nvtxEndRange()

         call nvtxStartRange('update_witness_hdf5: extend_dataset_hdf5')
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
         call nvtxEndRange()

         call nvtxStartRange('update_witness_hdf5: save_array2D_tr_rp_in_dataset_hdf5_file')
         call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwrite,isCreateDataspaceOpt=.False.)
         call nvtxEndRange()
      end if

      !Save time!
      call nvtxStartRange('update_witness_hdf5: extend_dataset_hdf5')
      dsetname     = 'time'
      ms_dims(1)   = leapwitsave
      ms_offset(1) = itewit - leapwitsave
      ds_rank      = 1
      ds_dims(1)   = itewit
      call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
      call nvtxEndRange()

      call nvtxStartRange('update_witness_hdf5: save_array1D_rp_in_dataset_hdf5_file')
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,t,isCreateDataspaceOpt=.False.)
      call nvtxEndRange()

      !Save istep!
      call nvtxStartRange('update_witness_hdf5: extend_dataset_hdf5')
      dsetname     = 'istep'
      ms_dims(1)   = leapwitsave
      ms_offset(1) = itewit - leapwitsave
      ds_rank      = 1
      ds_dims(1)   = itewit
      call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
      call nvtxEndRange()

      call nvtxStartRange('update_witness_hdf5: save_array1D_int4_in_dataset_hdf5_file')
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,steps)
      call nvtxEndRange()

      call nvtxStartRange('update_witness_hdf5: CPU-DATA')
      call h5fclose_f(file_id,h5err)
      call nvtxEndRange()

   end subroutine update_witness_hdf5

end module mod_hdf5
