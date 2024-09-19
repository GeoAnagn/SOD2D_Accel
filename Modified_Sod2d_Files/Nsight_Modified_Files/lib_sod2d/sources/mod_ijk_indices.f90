module mod_ijk_indices
   use mod_constants
   use mod_mpi
   use mod_nvtx
   use elem_qua ! Using only the ordering table for edges
   use elem_hex ! Using only the ordering tables for edges and faces
   implicit none

   !--------------------------------------------------------------------------------------------
   ! GMSH Indices
   integer(4),allocatable :: gmshHexahedraHO_ijkTable(:,:,:),gmshQuadrilateralHO_ijTable(:,:)
   integer(4) :: gmsh_porder=0
   !--------------------------------------------------------------------------------------------

contains

   subroutine get_porder_values(mporder,mnnode,mngaus,mnpbou)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),intent(out) :: mnnode,mngaus,mnpbou

      call nvtxStartRange('get_porder_values: CPU-COMPUTE')
      mnnode = (mporder+1)**3
      call nvtxEndRange()

      call nvtxStartRange('get_porder_values: CPU-DATA')
      mngaus = mnnode
      call nvtxEndRange()

      call nvtxStartRange('get_porder_values: CPU-COMPUTE')
      mnpbou = (mporder+1)**2
      call nvtxEndRange()

   end subroutine get_porder_values

   subroutine set_allocate_array_ijk_sod2d_criteria(mporder,ijk_sod2d_to_gmsh,ijk_gmsh_to_sod2d)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),intent(inout),dimension(:),allocatable :: ijk_sod2d_to_gmsh,ijk_gmsh_to_sod2d
      integer(4) :: i,j

      !--------------------------------------------------------------------
      !defining criteria for ijk
      call nvtxStartRange('set_allocate_array_ijk_sod2d_criteria: CPU-DATA')
      allocate(ijk_sod2d_to_gmsh(mporder+1))
      allocate(ijk_gmsh_to_sod2d(0:mporder))

      i=1
      ijk_sod2d_to_gmsh(i) = 0
      ijk_gmsh_to_sod2d(ijk_sod2d_to_gmsh(i)) = i

      i=i+1
      ijk_sod2d_to_gmsh(i) = mporder
      ijk_gmsh_to_sod2d(ijk_sod2d_to_gmsh(i)) = i
      call nvtxEndRange()

      call nvtxStartRange('set_allocate_array_ijk_sod2d_criteria: CPU-COMPUTE')
      do j=1,(mporder-1)
         i=i+1
         ijk_sod2d_to_gmsh(i)=j
         ijk_gmsh_to_sod2d(ijk_sod2d_to_gmsh(i)) = i
      end do
      call nvtxEndRange()
   end subroutine set_allocate_array_ijk_sod2d_criteria

   function get_indexIJK_sod2d(mporder,i,j,k) result(indexIJK)
      implicit none
      integer(4),intent(in) :: mporder,i,j,k
      integer :: indexIJK

      call nvtxStartRange('get_indexIJK_sod2d: CPU-COMPUTE')
      indexIJK = ((mporder+1)**2)*k+(mporder+1)*i+j+1
      call nvtxEndRange()
   end function get_indexIJK_sod2d

   subroutine set_allocate_hexahedronHO_ijk_indices(mporder,gmsh2ijk,vtk2ijk)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),allocatable,dimension(:),intent(inout) :: gmsh2ijk,vtk2ijk
      !integer(4),allocatable,intent(inout),optional :: ijk2gmsh(:,:,:),ijk2vtk(:,:,:)
      integer(4) :: mnnode,mngaus,mnpbou
      integer(4) :: i,j,k,ip,jp,kp,gmshCnt,vtkCnt,gmshIndex,vtkIndex,pIndex
      integer(4),allocatable :: ijk_sod2d_to_gmsh(:),ijk_gmsh_to_sod2d(:)

      !-----------------------------------------------------------------------------------------
      call nvtxStartRange('set_allocate_hexahedronHO_ijk_indices: get_porder_values')
      call get_porder_values(mporder,mnnode,mngaus,mnpbou)
      call nvtxEndRange()
      !if(mpi_rank.eq.0) write(*,*) 'mporder',mporder,'mnnode',mnnode,'mngaus',mngaus,'mnpbou',mnpbou

      call nvtxStartRange('set_allocate_hexahedronHO_ijk_indices: set_allocate_array_ijk_sod2d_criteria')
      call set_allocate_array_ijk_sod2d_criteria(mporder,ijk_sod2d_to_gmsh,ijk_gmsh_to_sod2d)
      call nvtxEndRange()
      !-----------------------------------------------------------------------------------------

      call nvtxStartRange('set_allocate_hexahedronHO_ijk_indices: CPU-DATA')
      allocate(gmsh2ijk(mnnode))
      allocate(vtk2ijk(mnnode))
      call nvtxEndRange()

      !if(present(ijk2gmsh)) allocate(ijk2gmsh(0:porder,0:porder,0:porder))
      !if(present(ijk2vtk))  allocate(ijk2vtk(0:porder,0:porder,0:porder))

      if(mporder<2) then
         write(*,*) 'SOD2D is not ready to work for mporder < 2... You know, #gobigorgohome and set mporder >= 2'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      !--------------------------------------------------------------------
      !Filling high order hexahedra
      pIndex=0
      do kp=1,mporder+1
         k=ijk_sod2d_to_gmsh(kp)
         do ip=1,mporder+1
            i=ijk_sod2d_to_gmsh(ip)
            do jp=1,mporder+1
               j=ijk_sod2d_to_gmsh(jp)
               
               call nvtxStartRange('set_allocate_hexahedronHO_ijk_indices: vtkHigherOrderHexahedron_pointIndexFromIJK')
               call vtkHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,vtkIndex)
               call nvtxEndRange()

               call nvtxStartRange('set_allocate_hexahedronHO_ijk_indices: gmshHigherOrderHexahedron_pointIndexFromIJK')
               call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
               call nvtxEndRange()

               pIndex               = pIndex + 1
               gmsh2ijk(pIndex)     = gmshIndex
               vtk2ijk(pIndex)      = vtkIndex

            end do
         end do
      end do

      !--------------------------------------------------------------------

      call nvtxStartRange('set_allocate_hexahedronHO_ijk_indices: CPU-DATA')
      deallocate(ijk_sod2d_to_gmsh)
      deallocate(ijk_gmsh_to_sod2d)
      call nvtxEndRange()
   end subroutine set_allocate_hexahedronHO_ijk_indices

   subroutine set_allocate_quadrilateralHO_ij_indices(mporder,gmsh2ij,vtk2ij)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),allocatable,dimension(:),intent(inout) :: gmsh2ij,vtk2ij
      integer(4) :: mnnode,mngaus,mnpbou
      integer(4) :: i,j,ip,jp,gmshCnt,vtkCnt,gmshIndex,vtkIndex,pIndex
      integer(4),allocatable :: ijk_sod2d_to_gmsh(:),ijk_gmsh_to_sod2d(:)

      !-----------------------------------------------------------------------------------------
      call nvtxStartRange('set_allocate_quadrilateralHO_ij_indices: get_porder_values')
      call get_porder_values(mporder,mnnode,mngaus,mnpbou)
      call nvtxEndRange()
      !if(mpi_rank.eq.0) write(*,*) 'mporder',mporder,'mnnode',mnnode,'mngaus',mngaus,'mnpbou',mnpbou

      call nvtxStartRange('set_allocate_quadrilateralHO_ij_indices: set_allocate_array_ijk_sod2d_criteria')
      call set_allocate_array_ijk_sod2d_criteria(mporder,ijk_sod2d_to_gmsh,ijk_gmsh_to_sod2d)
      call nvtxEndRange()
      !-----------------------------------------------------------------------------------------

      call nvtxStartRange('set_allocate_quadrilateralHO_ij_indices: CPU-DATA')
      allocate(gmsh2ij(mnpbou))
      allocate(vtk2ij(mnpbou))
      call nvtxEndRange()

      !--------------------------------------------------------------------
      !filling high order quads
      pIndex=0
      do ip=1,mporder+1
         i=ijk_sod2d_to_gmsh(ip)
         do jp=1,mporder+1
            j=ijk_sod2d_to_gmsh(jp)
            
            call nvtxStartRange('set_allocate_quadrilateralHO_ij_indices: vtkHigherOrderQuadrilateral_pointIndexFromIJ')
            call vtkHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,vtkIndex)
            call nvtxEndRange()

            call nvtxStartRange('set_allocate_quadrilateralHO_ij_indices: gmshHigherOrderQuadrilateral_pointIndexFromIJ')
            call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
            call nvtxEndRange()

            pIndex = pIndex + 1
            gmsh2ij(pIndex) = gmshIndex
            vtk2ij(pIndex)  = vtkIndex

         end do
      end do

      call nvtxStartRange('set_allocate_quadrilateralHO_ij_indices: CPU-DATA')
      deallocate(ijk_sod2d_to_gmsh)
      deallocate(ijk_gmsh_to_sod2d)
      call nvtxEndRange()
   end subroutine set_allocate_quadrilateralHO_ij_indices

   subroutine get_gmshHexHOVertIndex(mporder,gmshHexVertInd)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),intent(inout) :: gmshHexVertInd(8)
      integer(4) :: i,j,k,gmshIndex

      !--------------------------------------------------------------------
      !Filling gmshHexVertInd(:)
      !--------------------------------------------------------------------
      i=0
      j=0
      k=0
      call nvtxStartRange('get_gmshHexHOVertIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      call nvtxEndRange()
      gmshHexVertInd(1) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder
      j=0
      k=0
      call nvtxStartRange('get_gmshHexHOVertIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      call nvtxEndRange()
      gmshHexVertInd(2) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder
      j=mporder
      k=0
      call nvtxStartRange('get_gmshHexHOVertIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      call nvtxEndRange()
      gmshHexVertInd(3) = gmshIndex
      !--------------------------------------------------------------------
      i=0
      j=mporder
      k=0
      call nvtxStartRange('get_gmshHexHOVertIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      call nvtxEndRange()
      gmshHexVertInd(4) = gmshIndex
      !--------------------------------------------------------------------
      i=0
      j=0
      k=mporder
      call nvtxStartRange('get_gmshHexHOVertIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      call nvtxEndRange()
      gmshHexVertInd(5) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder
      j=0
      k=mporder
      call nvtxStartRange('get_gmshHexHOVertIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      call nvtxEndRange()
      gmshHexVertInd(6) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder
      j=mporder
      k=mporder
      call nvtxStartRange('get_gmshHexHOVertIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      call nvtxEndRange()
      gmshHexVertInd(7) = gmshIndex
      !--------------------------------------------------------------------
      i=0
      j=mporder
      k=mporder
      call nvtxStartRange('get_gmshHexHOVertIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      call nvtxEndRange()
      gmshHexVertInd(8) = gmshIndex
      !--------------------------------------------------------------------

   end subroutine get_gmshHexHOVertIndex

   subroutine get_gmshHexHOFacesIndex(mporder,mnpbou,gmshHexFaceFrontInd,gmshHexFaceBackInd,gmshHexFaceBottomInd,gmshHexFaceTopInd,gmshHexFaceLeftInd,gmshHexFaceRightInd)
      implicit none
      integer(4),intent(in) :: mporder,mnpbou
      integer(4),intent(inout),dimension(mnpbou) :: gmshHexFaceFrontInd,gmshHexFaceBackInd,gmshHexFaceBottomInd,gmshHexFaceTopInd,gmshHexFaceLeftInd,gmshHexFaceRightInd
      integer(4) :: i,j,k,gmshIndex,gmshCnt

      !--------------------------------------------------------------------
      !Filling faceFront2ijk(:)
      gmshCnt=0
      j=0
      do i=0,mporder
         do k=0,mporder
            call nvtxStartRange('get_gmshHexHOFacesIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
            call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
            call nvtxEndRange()

            gmshCnt = gmshCnt+1
            gmshHexFaceFrontInd(gmshCnt) = gmshIndex
         end do
      end do
      !--------------------------------------------------------------------
      !Filling faceBack2ijk(:)
      gmshCnt=0
      j=mporder
      do i=0,mporder
         do k=0,mporder
            call nvtxStartRange('get_gmshHexHOFacesIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
            call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
            call nvtxEndRange()

            gmshCnt = gmshCnt+1
            gmshHexFaceBackInd(gmshCnt) = gmshIndex
         end do
      end do
      !--------------------------------------------------------------------
      !Filling faceBottom2ijk(:)
      gmshCnt=0
      k=0
      do i=0,mporder
         do j=0,mporder
            call nvtxStartRange('get_gmshHexHOFacesIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
            call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
            call nvtxEndRange()
            gmshCnt = gmshCnt+1
            gmshHexFaceBottomInd(gmshCnt) = gmshIndex
         end do
      end do
      !--------------------------------------------------------------------
      !Filling faceTop2ijk(:)
      gmshCnt=0
      k=mporder
      do i=0,mporder
         do j=0,mporder
            call nvtxStartRange('get_gmshHexHOFacesIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
            call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
            call nvtxEndRange()

            gmshCnt = gmshCnt+1
            gmshHexFaceTopInd(gmshCnt) = gmshIndex
         end do
      end do
      !--------------------------------------------------------------------
      !Filling faceLeft2ijk(:)
      gmshCnt=0
      i=0
      do j=0,mporder
         do k=0,mporder
            call nvtxStartRange('get_gmshHexHOFacesIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
            call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
            call nvtxEndRange()

            gmshCnt = gmshCnt+1
            gmshHexFaceLeftInd(gmshCnt) = gmshIndex
         end do
      end do
      !--------------------------------------------------------------------
      !Filling faceRight2ijk(:)
      gmshCnt=0
      i=mporder
      do j=0,mporder
         do k=0,mporder
            call nvtxStartRange('get_gmshHexHOFacesIndex: gmshHigherOrderHexahedron_pointIndexFromIJK')
            call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
            call nvtxEndRange()

            gmshCnt = gmshCnt+1
            gmshHexFaceRightInd(gmshCnt) = gmshIndex
         end do
      end do
      !--------------------------------------------------------------------

   end subroutine get_gmshHexHOFacesIndex

   subroutine get_gmshQuadHOVertIndex(mporder,gmshQuadVertInd)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),intent(inout) :: gmshQuadVertInd(4)
      integer(4) :: i,j,gmshIndex

      i=0
      j=0
      call nvtxStartRange('get_gmshQuadHOVertIndex: gmshHigherOrderQuadrilateral_pointIndexFromIJ')
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      call nvtxEndRange()
      gmshQuadVertInd(1) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder
      j=0
      call nvtxStartRange('get_gmshQuadHOVertIndex: gmshHigherOrderQuadrilateral_pointIndexFromIJ')
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      call nvtxEndRange()
      gmshQuadVertInd(2) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder
      j=mporder
      call nvtxStartRange('get_gmshQuadHOVertIndex: gmshHigherOrderQuadrilateral_pointIndexFromIJ')
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      call nvtxEndRange()
      gmshQuadVertInd(3) = gmshIndex
      !--------------------------------------------------------------------
      i=0
      j=mporder
      call nvtxStartRange('get_gmshQuadHOVertIndex: gmshHigherOrderQuadrilateral_pointIndexFromIJ')
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      call nvtxEndRange()
      gmshQuadVertInd(4) = gmshIndex
      !--------------------------------------------------------------------

   end subroutine get_gmshQuadHOVertIndex

   subroutine get_gmshQuadHOInnerVertIndex(mporder,gmshQuadInnerVertInd)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),intent(inout) :: gmshQuadInnerVertInd(4)
      integer(4) :: i,j,gmshIndex


      !--------------------------------------------------------------------
      !Filling gmshQuadInnerVertInd(:)
      !--------------------------------------------------------------------
      i=1
      j=1
      call nvtxStartRange('get_gmshQuadHOInnerVertIndex: gmshHigherOrderQuadrilateral_pointIndexFromIJ')
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      call nvtxEndRange()
      gmshQuadInnerVertInd(1) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder-1
      j=1
      call nvtxStartRange('get_gmshQuadHOInnerVertIndex: gmshHigherOrderQuadrilateral_pointIndexFromIJ')
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      call nvtxEndRange()
      gmshQuadInnerVertInd(2) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder-1
      j=mporder-1
      call nvtxStartRange('get_gmshQuadHOInnerVertIndex: gmshHigherOrderQuadrilateral_pointIndexFromIJ')
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      call nvtxEndRange()
      gmshQuadInnerVertInd(3) = gmshIndex
      !--------------------------------------------------------------------
      i=1
      j=mporder-1
      call nvtxStartRange('get_gmshQuadHOInnerVertIndex: gmshHigherOrderQuadrilateral_pointIndexFromIJ')
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      call nvtxEndRange()
      gmshQuadInnerVertInd(4) = gmshIndex
      !--------------------------------------------------------------------

   end subroutine get_gmshQuadHOInnerVertIndex

   subroutine vtkHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,pointIndex)
      implicit none
      integer(4),intent(in) :: mporder,i,j,k
      integer(4),intent(out) :: pointIndex
      integer(4) :: ibdy,jbdy,kbdy,nbdy
      integer(4) :: aux_pi

      !----------------------------------------------------------------------------
      call nvtxStartRange('vtkHigherOrderHexahedron_pointIndexFromIJK: CPU-DATA')
      if((i.eq.0).or.(i.eq.mporder)) then
         ibdy = 1
      else
         ibdy = 0
      endif

      if((j.eq.0).or.(j.eq.mporder)) then
         jbdy = 1
      else
         jbdy = 0
      endif

      if((k.eq.0).or.(k.eq.mporder)) then
         kbdy = 1
      else
         kbdy = 0
      endif
      call nvtxEndRange()

      call nvtxStartRange('vtkHigherOrderHexahedron_pointIndexFromIJK: CPU-COMPUTE')
      nbdy = ibdy + jbdy + kbdy

      !----------------------------------------------------------------------------
      pointIndex = 1

      if(nbdy .eq. 3) then !Vertex
         !return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);

         if(i.ne.0) then
            if(j.ne.0) then
               pointIndex = pointIndex + 2
            else
               pointIndex = pointIndex + 1
            end if
         else
            if(j.ne.0) then
               pointIndex = pointIndex + 3
            end if
         end if

         if(k.ne.0) then
            pointIndex = pointIndex + 4
         end if

         return
      end if

      pointIndex = pointIndex + 8
      if(nbdy .eq. 2) then !Edge

         if(ibdy.eq.0) then !iaxis
            !return (i - 1) + (j ? order[0] + order[1] - 2 : 0) + (k ? 2 * (order[0] + order[1] - 2) : 0) + offset;
            pointIndex = pointIndex + (i-1)
            if(j.ne.0) then
               pointIndex = pointIndex + (mporder*2 - 2)
            end if
            if(k.ne.0) then
               pointIndex = pointIndex + 2*(mporder*2 - 2)
            end if
         else if(jbdy.eq.0) then !jaxis
            !return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) + (k ? 2 * (order[0] + order[1] - 2) : 0) + offset;
            pointIndex = pointIndex + (j-1)
            if(i.ne.0) then
               pointIndex = pointIndex + (mporder - 1)
            else
               pointIndex = pointIndex + (2*(mporder - 1) + mporder - 1)
            end if
            if(k.ne.0) then
               pointIndex = pointIndex + 2*(2*mporder - 2)
            end if
         else !kaxis
            !offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
            !return (k - 1) + (order[2] - 1) * (i ? (j ? 2 : 1) : (j ? 3 : 0)) + offset;
            pointIndex = pointIndex + 4*(mporder-1)+ 4*(mporder-1)

            aux_pi = 0
            if(i.ne.0) then
               if(j.ne.0) then
                  aux_pi = 2
               else
                  aux_pi = 1
               end if
            else
               if(j.ne.0) then
                  aux_pi = 3
               end if
            end if

            pointIndex = pointIndex + (k-1) + (mporder - 1)*aux_pi
         end if

         return
      end if

      pointIndex = pointIndex + 4*(3*mporder - 3)
      if(nbdy .eq. 1) then !Face
         if(ibdy.ne.0) then ! on i-normal face
            pointIndex = pointIndex + (j-1) + ((mporder-1)*(k-1))
            if(i.ne.0) then
               pointIndex = pointIndex + (mporder-1)*(mporder-1)
            end if
            return
         end if

         pointIndex = pointIndex + 2*(mporder-1)*(mporder-1)
         if(jbdy.ne.0) then! on j-normal face
            pointIndex = pointIndex + (i-1) + ((mporder-1)*(k-1))
            if(j.ne.0) then
               pointIndex = pointIndex + (mporder-1)*(mporder-1)
            end if
            return
         end if

         ! on k-normal face
         pointIndex = pointIndex + 2*(mporder-1)*(mporder-1)
         pointIndex = pointIndex + (i-1) + ((mporder-1)*(j-1))
         if(k.ne.0) then
            pointIndex = pointIndex + (mporder-1)*(mporder-1)
         end if
         return

      end if

      ! nbdy == 0: Body DOF
      pointIndex = pointIndex + 2* ((mporder-1)*(mporder-1)+(mporder-1)*(mporder-1)+(mporder-1)*(mporder-1))
      pointIndex = pointIndex + (i-1)+(mporder-1)*((j-1)+(mporder-1)*((k - 1)))
      call nvtxEndRange()
   end subroutine vtkHigherOrderHexahedron_pointIndexFromIJK

   subroutine vtkHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,pointIndex)
      implicit none
      integer(4),intent(in) :: mporder,i,j
      integer(4),intent(out) :: pointIndex
      integer(4) :: ibdy,jbdy,nbdy
      integer(4) :: aux_pi

      !----------------------------------------------------------------------------
      call nvtxStartRange('vtkHigherOrderQuadrilateral_pointIndexFromIJ: CPU-DATA')
      if((i.eq.0).or.(i.eq.mporder)) then
         ibdy = 1
      else
         ibdy = 0
      endif

      if((j.eq.0).or.(j.eq.mporder)) then
         jbdy = 1
      else
         jbdy = 0
      endif
      call nvtxEndRange()

      call nvtxStartRange('vtkHigherOrderQuadrilateral_pointIndexFromIJ: CPU-COMPUTE')
      nbdy = ibdy + jbdy

      !----------------------------------------------------------------------------
      pointIndex = 1

      if(nbdy .eq. 2) then !Vertex
         !return (i ? (j ? 2 : 1) : (j ? 3 : 0));

         if(i.ne.0) then
            if(j.ne.0) then
               pointIndex = pointIndex + 2
            else
               pointIndex = pointIndex + 1
            end if
         else
            if(j.ne.0) then
               pointIndex = pointIndex + 3
            end if
         end if

         return
      end if

      pointIndex = pointIndex + 4
      if(nbdy .eq. 1) then !Edge

         if(ibdy.eq.0) then !iaxis
            !return (i - 1) + (j ? order[0] - 1 + order[1] - 1 : 0) + offset;

            pointIndex = pointIndex + (i-1)
            if(j.ne.0) then
               pointIndex = pointIndex + (mporder*2 - 2)
            end if

         else !jaxis
            !return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) + offset;

            pointIndex = pointIndex + (j-1)
            if(i.ne.0) then
               pointIndex = pointIndex + (mporder - 1)
            else
               pointIndex = pointIndex + (2*(mporder - 1) + mporder - 1)
            end if

         end if

         return
      end if

      ! nbdy == 0: Body DOF
      pointIndex = pointIndex + (4*mporder-4)
      pointIndex = pointIndex + (i-1)+(mporder-1)*(j-1)
      call nvtxEndRange()
   end subroutine vtkHigherOrderQuadrilateral_pointIndexFromIJ

!--------------------------------------------------------------------------------------------------

   recursive subroutine genHighOrderHexGmsh(p,indexTable)
      implicit none
      integer(4), intent(in)  :: p
      integer(4), intent(out) :: indexTable((p+1)**3,3) ! Set as I, J, K
      integer(4)              :: tableFace((p-1)**2,2), tableVolume((p-1)**3,3)
      integer(4)              :: inode, iedge, iface, i0, i1, i3, u(3), v(3), i

      ! Initialize corner node to 0th position, or generate element of order 0
      indexTable(1,:) = [0,0,0]

      ! Generate element of order 1 (corner nodes)
      if (p .gt. 0) then
         !                  i, j, k
         call nvtxStartRange('genHighOrderHexGmsh: CPU-DATA')
         indexTable(2,:) = [p, 0, 0]
         indexTable(3,:) = [p, p, 0]
         indexTable(4,:) = [0, p, 0]
         indexTable(5,:) = [0, 0, p]
         indexTable(6,:) = [p, 0, p]
         indexTable(7,:) = [p, p, p]
         indexTable(8,:) = [0, p, p]
         call nvtxEndRange()
         if (p .gt. 1) then
            ! Generate high-order edges
            call nvtxStartRange('genHighOrderHexGmsh: CPU-COMPUTE')
            inode = 9
            do iedge = 1,12
               i0 = hex_order_edges(iedge,1)
               i1 = hex_order_edges(iedge,2)
               u(:) = (indexTable(i1,:) - indexTable(i0,:))/p
               do i = 1,p-1
                  indexTable(inode,:) = indexTable(i0,:) + i*u(:)
                  inode = inode + 1
               end do
            end do
            call nvtxEndRange()

            ! Generate a generic high-order face with p` = p-2
            call nvtxStartRange('genHighOrderHexGmsh: genHighOrderQuadGmsh')
            call genHighOrderQuadGmsh(p-2,tableFace)
            call nvtxEndRange()

            call nvtxStartRange('genHighOrderHexGmsh: CPU-COMPUTE')
            tableFace = tableFace + 1
            ! Generate faces interior nodes
            do iface = 1,6
               i0 = hex_order_faces(iface,1)
               i1 = hex_order_faces(iface,2)
               i3 = hex_order_faces(iface,4)
               u(:) = (indexTable(i1,:) - indexTable(i0,:))/p
               v(:) = (indexTable(i3,:) - indexTable(i0,:))/p
               do i = 1,((p-1)**2)
                  indexTable(inode,:) = indexTable(i0,:) + u(:)*tableFace(i,1) + v(:)*tableFace(i,2)
                  inode = inode + 1
               end do
            end do
            call nvtxEndRange()

            ! Generate volume nodes
            call nvtxStartRange('genHighOrderHexGmsh: genHighOrderHexGmsh')
            call genHighOrderHexGmsh(p-2,tableVolume)
            call nvtxEndRange()
            tableVolume = tableVolume + 1

            call nvtxStartRange('genHighOrderHexGmsh: joinTablesGmsh')
            call joinTablesGmsh([(p-1)**3,3],tableVolume,inode,[(p+1)**3,3],indexTable)
            call nvtxEndRange()
         end if
      end if
   end subroutine genHighOrderHexGmsh

   recursive subroutine genHighOrderQuadGmsh(p,indexTable)
      implicit none
      integer(4), intent(in)  :: p
      integer(4), intent(out) :: indexTable((p+1)**2,2) ! Set as I, J
      integer(4)              :: tableFace((p-1)**2,2)
      integer(4)              :: inode, iedge, iface, i0, i1, u(2), i

      indexTable(1,:) = [0,0]
      if (p .gt. 0) then
         call nvtxStartRange('genHighOrderQuadGmsh: CPU-DATA')
         indexTable(2,:) = [p,0]
         indexTable(3,:) = [p,p]
         indexTable(4,:) = [0,p]
         call nvtxEndRange()
         if (p .gt. 1) then
            call nvtxStartRange('genHighOrderQuadGmsh: CPU-COMPUTE')
            inode = 5
            do iedge = 1,4
               i0 = quad_order_edges(iedge,1)
               i1 = quad_order_edges(iedge,2)
               u(:) = (indexTable(i1,:) - indexTable(i0,:))/p
               do i = 1,p-1
                  indexTable(inode,:) = indexTable(i0,:) + i*u(:)
                  inode = inode + 1
               end do
            end do
            call nvtxEndRange()

            call nvtxStartRange('genHighOrderQuadGmsh: genHighOrderQuadGmsh')
            call genHighOrderQuadGmsh(p-2,tableFace)
            call nvtxEndRange()

            tableFace = tableFace + 1
            
            call nvtxStartRange('genHighOrderQuadGmsh: joinTablesGmsh')
            call joinTablesGmsh([(p-1)**2,2],tableFace,inode,[(p+1)**2,2],indexTable)
            call nvtxEndRange()
         end if
      end if
   end subroutine genHighOrderQuadGmsh

   subroutine joinTablesGmsh(size1,table1,indexDesti,size2,table2)
      implicit none
      integer(4), intent(in)    :: indexDesti, size1(2), size2(2)
      integer(4), intent(in)    :: table1(size1(1),size1(2))
      integer(4), intent(inout) :: table2(size2(1),size2(2))
      integer(4)                :: i, j

      call nvtxStartRange('joinTablesGmsh: CPU-DATA')
      j = indexDesti
      do i = 1,size1(1)
         table2(j,:) = table1(i,:)
         j = j + 1
      end do
      call nvtxEndRange()
   end subroutine joinTablesGmsh

!--------------------------------------------------------------------------------------------------

   subroutine initGmshIJKTables(mporder)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4) :: mnnode,mngaus,mnpbou,pIndex,i,j,k
      integer(4),allocatable :: auxHexHOtable(:,:),auxQuadHOtable(:,:)

      call nvtxStartRange('initGmshIJKTables: get_porder_values')
      call get_porder_values(mporder,mnnode,mngaus,mnpbou)
      call nvtxEndRange()
      !if(mpi_rank.eq.0) write(*,*) 'mporder',mporder,'mnnode',mnnode,'mngaus',mngaus,'mnpbou',mnpbou

      if(gmsh_porder.eq.0) then !arrays not initialized
         if(mpi_rank.eq.0) write(*,*) 'Initialising GMSH IJK Tables to order',mporder

      else if(gmsh_porder.eq.mporder) then !arrays already initalized to current order, do nothing and exit!
         if(mpi_rank.eq.0) write(*,*) 'GMSH IJK Tables already initialised to order',mporder,'doing nothing! :)'
         return
      else !arrays initalized to another other
         if(mpi_rank.eq.0) write(*,*) 'GMSH IJK Tables initialised to order',gmsh_porder,'! -> changing to order',mporder

         call nvtxStartRange('initGmshIJKTables: CPU-DATA')
         deallocate(gmshHexahedraHO_ijkTable)
         deallocate(gmshQuadrilateralHO_ijTable)
         call nvtxEndRange()
      end if

      call nvtxStartRange('initGmshIJKTables: CPU-DATA')
      allocate(auxHexHOtable(mnnode,3))
      allocate(auxQuadHOtable(mnpbou,2))
      allocate(gmshHexahedraHO_ijkTable(0:mporder,0:mporder,0:mporder))
      allocate(gmshQuadrilateralHO_ijTable(0:mporder,0:mporder) )
      call nvtxEndRange()

      call nvtxStartRange('initGmshIJKTables: init_quad_info')
      call init_quad_info()
      call nvtxEndRange()

      call nvtxStartRange('initGmshIJKTables: init_hex_info')
      call init_hex_info()
      call nvtxEndRange()

      call nvtxStartRange('initGmshIJKTables: genHighOrderHexGmsh')
      call genHighOrderHexGmsh(mporder,auxHexHOtable)
      call nvtxEndRange()

      call nvtxStartRange('initGmshIJKTables: genHighOrderQuadGmsh')
      call genHighOrderQuadGmsh(mporder,auxQuadHOtable)
      call nvtxEndRange()

      gmsh_porder = mporder

      !write(*,*) 'generating HexHO_ijktable'

      call nvtxStartRange('initGmshIJKTables: CPU-DATA')
      do pIndex=1,mnnode

         i = auxHexHOtable(pIndex,1)
         j = auxHexHOtable(pIndex,2)
         k = auxHexHOtable(pIndex,3)

         !write(*,*) 'ijk',i,j,k,'pI',pIndex

         gmshHexahedraHO_ijkTable(i,j,k) = pIndex
      end do

      !write(*,*) 'generating QuadHO_ijtable'

      do pIndex=1,mnpbou
         i = auxQuadHOtable(pIndex,1)
         j = auxQuadHOtable(pIndex,2)

         !write(*,*) 'ij',i,j,'pI',pIndex
         gmshQuadrilateralHO_ijTable(i,j) = pIndex
      end do

      deallocate(auxHexHOtable)
      deallocate(auxQuadHOtable)
      call nvtxEndRange()
   end subroutine initGmshIJKTables

   subroutine gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,pointIndex)
      implicit none
      integer(4),intent(in) :: mporder,i,j,k
      integer(4),intent(out) :: pointIndex

      call nvtxStartRange('gmshHigherOrderHexahedron_pointIndexFromIJK: CPU-DATA')
      if(gmsh_porder.ne.mporder) then !arrays not initialized
         if(mpi_rank.eq.0) write(*,*) 'ERROR! GMSH IJK TABLES NOT PROPERLY INITALISED!! gmsh_porder',gmsh_porder,'!= mporder',mporder
         pointIndex = 0
         return
      end if

      pointIndex = gmshHexahedraHO_ijkTable(i,j,k)
      call nvtxEndRange()

   end subroutine gmshHigherOrderHexahedron_pointIndexFromIJK

   subroutine gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,pointIndex)
      implicit none
      integer(4),intent(in) :: mporder,i,j
      integer(4),intent(out) :: pointIndex

      call nvtxStartRange('gmshHigherOrderQuadrilateral_pointIndexFromIJ: CPU-DATA')
      if(gmsh_porder.ne.mporder) then !arrays not initialized
         if(mpi_rank.eq.0) write(*,*) 'ERROR! GMSH IJK TABLES NOT PROPERLY INITALISED!! gmsh_porder',gmsh_porder,'!= mporder',mporder
         pointIndex = 0
         return
      end if

      pointIndex = gmshQuadrilateralHO_ijTable(i,j)
      call nvtxEndRange()
   end subroutine gmshHigherOrderQuadrilateral_pointIndexFromIJ

end module mod_ijk_indices
