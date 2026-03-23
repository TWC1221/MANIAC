module m_wall_current
  use mesh_types
  use BC_types
  implicit none
  private
  public :: compute_wall_currents

contains

subroutine compute_wall_currents(phi, mesh, Dg, BCs, Jx_min, Jx_max, Jy_min, Jy_max, Jz_min, Jz_max)
  implicit none

  ! Inputs
  type(BoundarySet), intent(in)  :: BCs   
  type(MeshGrid),    intent(in)  :: mesh
  real(8),           intent(in)  :: phi(:,:)    ! (G, nx*ny*nz)
  real(8),           intent(in)  :: Dg(:)         ! diffusion for this group

  ! Outputs (for all groups, but this call only fills slice gg,:)
  real(8), intent(inout) :: Jx_min(:,:), Jx_max(:,:)  ! (G, ny*nz)
  real(8), intent(inout) :: Jy_min(:,:), Jy_max(:,:)  ! (G, nx*nz)
  real(8), intent(inout) :: Jz_min(:,:), Jz_max(:,:)  ! (G, nx*ny)

  integer :: nx, ny, nz, G, gg
  integer :: j, k, i
  integer :: rowP, rowN
  integer :: idx_face
  real(8) :: dx, dy, dz, h, phi_bc, gamma

  nx = mesh%nx
  ny = mesh%ny
  nz = mesh%nz
  G  = size(phi,1)

  dx = mesh%dx
  dy = mesh%dy
  dz = mesh%dz

  ! Initialize group currents
  Jx_min(:,:) = 0.0d0
  Jx_max(:,:) = 0.0d0
  Jy_min(:,:) = 0.0d0
  Jy_max(:,:) = 0.0d0
  Jz_min(:,:) = 0.0d0
  Jz_max(:,:) = 0.0d0

  !===========================================================
  ! XMIN face: outward normal points -x, outflow is positive
  !===========================================================
  do gg = 1,G
    do k = 1, nz
      do j = 1, ny
        rowP     = 1 + (j-1)*nx + (k-1)*nx*ny
        idx_face = j + (k-1)*ny

        select type (bc => BCs%face(FACE_XMIN)%p)
        type is (BC_Dirichlet)
          h      = dx
          phi_bc = bc%value
          Jx_min(gg,idx_face) = (2.0d0*Dg(gg)/h) * (phi(gg,rowP) - phi_bc) * dy * dz

        type is (BC_Neumann)
          Jx_min(gg,idx_face) = 0.0d0

        type is (BC_Vacuum)
          h = dx
          Jx_min(gg,idx_face) = (2.0d0 * Dg(gg)) / (4.0d0 * Dg(gg) + h) * phi(gg,rowP) * dy * dz

        type is (BC_Albedo)
          h = dx
          gamma = (1.0d0/(2.0d0*Dg(gg)))*(1.0d0 - bc%alpha)/(1.0d0 + bc%alpha)
          Jx_min(gg,idx_face) = (2.0d0*Dg(gg)*gamma)/(2.0d0 + gamma*h) * phi(gg,rowP) * dy * dz * (1.0d0 - bc%alpha)/(1.0d0 + bc%alpha)

        class default
          if (nx > 1) then
            rowN = 2 + (j-1)*nx + (k-1)*nx*ny
            Jx_min(gg,idx_face) = Dg(gg) * (phi(gg,rowP) - phi(gg,rowN)) * dy * dz / dx
          else
            Jx_min(gg,idx_face) = 0.0d0
          end if
        end select
      end do
    end do

    !===========================================================
    ! XMAX face: outward normal points +x, outflow is positive
    !===========================================================
    do k = 1, nz
      do j = 1, ny
        rowP     = nx + (j-1)*nx + (k-1)*nx*ny
        idx_face = j + (k-1)*ny

        select type (bc => BCs%face(FACE_XMAX)%p)
        type is (BC_Dirichlet)
          h      = dx
          phi_bc = bc%value
          Jx_max(gg,idx_face) = (2.0d0*Dg(gg)/h) * (phi(gg,rowP) - phi_bc) * dy * dz

        type is (BC_Neumann)
          Jx_max(gg,idx_face) = 0.0d0

        type is (BC_Vacuum)
          h = dx
          Jx_max(gg,idx_face) = (2.0d0 * Dg(gg)) / (4.0d0 * Dg(gg) + h) * phi(gg,rowP) * dy * dz

        type is (BC_Albedo)
          h = dx
          gamma = (1.0d0/(2.0d0*Dg(gg)))*(1.0d0 - bc%alpha)/(1.0d0 + bc%alpha)
          Jx_max(gg,idx_face) = (2.0d0*Dg(gg)*gamma)/(2.0d0 + gamma*h) * phi(gg,rowP) * dy * dz * (1.0d0 - bc%alpha)/(1.0d0 + bc%alpha)

        class default
          if (nx > 1) then
            rowN = (nx-1) + (j-1)*nx + (k-1)*nx*ny
            Jx_max(gg,idx_face) = Dg(gg) * (phi(gg,rowP) - phi(gg,rowN)) * dy * dz / dx
          else
            Jx_max(gg,idx_face) = 0.0d0
          end if
        end select
      end do
    end do

    !===========================================================
    ! YMIN face: outward normal -y, outflow positive
    !===========================================================
    do k = 1, nz
      do i = 1, nx
        rowP     = i + (1-1)*nx + (k-1)*nx*ny
        idx_face = i + (k-1)*nx

        select type (bc => BCs%face(FACE_YMIN)%p)
        type is (BC_Dirichlet)
          h      = dy
          phi_bc = bc%value
          Jy_min(gg,idx_face) = (2.0d0*Dg(gg)/h) * (phi(gg,rowP) - phi_bc) * dx * dz

        type is (BC_Neumann)
          Jy_min(gg,idx_face) = 0.0d0

        type is (BC_Vacuum)
          h = dy
          Jy_min(gg,idx_face) = (2.0d0 * Dg(gg)) / (4.0d0 * Dg(gg) + h) * phi(gg,rowP) * dx * dz

        type is (BC_Albedo)
          h = dy
          gamma = (1.0d0/(2.0d0*Dg(gg)))*(1.0d0 - bc%alpha)/(1.0d0 + bc%alpha)
          Jy_min(gg,idx_face) = (2.0d0*Dg(gg)*gamma)/(2.0d0 + gamma*h) * phi(gg,rowP) * dx * dz * (1.0d0 - bc%alpha)/(1.0d0 + bc%alpha)

        class default
          if (ny > 1) then
            rowN = i + (2-1)*nx + (k-1)*nx*ny
            Jy_min(gg,idx_face) = Dg(gg) * (phi(gg,rowP) - phi(gg,rowN)) * dx * dz / dy
          else
            Jy_min(gg,idx_face) = 0.0d0
          end if
        end select
      end do
    end do

    !===========================================================
    ! YMAX face: outward normal +y, outflow positive
    !===========================================================
    do k = 1, nz
      do i = 1, nx
        rowP     = i + (ny-1)*nx + (k-1)*nx*ny
        idx_face = i + (k-1)*nx

        select type (bc => BCs%face(FACE_YMAX)%p)
        type is (BC_Dirichlet)
          h      = dy
          phi_bc = bc%value
          Jy_max(gg,idx_face) = (2.0d0*Dg(gg)/h) * (phi(gg,rowP) - phi_bc) * dx * dz

        type is (BC_Neumann)
          Jy_max(gg,idx_face) = 0.0d0

        type is (BC_Vacuum)
          h = dy
          Jy_max(gg,idx_face) = (2.0d0 * Dg(gg)) / (4.0d0 * Dg(gg) + h) * phi(gg,rowP) * dx * dz

        type is (BC_Albedo)
          h = dy
          gamma = (1.0d0/(2.0d0*Dg(gg)))*(1.0d0 - bc%alpha)/(1.0d0 + bc%alpha)
          Jy_max(gg,idx_face) = (2.0d0*Dg(gg)*gamma)/(2.0d0 + gamma*h) * phi(gg,rowP) * dx * dz * (1.0d0 - bc%alpha)/(1.0d0 + bc%alpha)

        class default
          if (ny > 1) then
            rowN = i + (ny-2)*nx + (k-1)*nx*ny
            Jy_max(gg,idx_face) = Dg(gg) * (phi(gg,rowP) - phi(gg,rowN)) * dx * dz / dy
          else
            Jy_max(gg,idx_face) = 0.0d0
          end if
        end select
      end do
    end do

    !===========================================================
    ! ZMIN face: outward normal -z, outflow positive
    !===========================================================
    do j = 1, ny
      do i = 1, nx
        rowP     = i + (j-1)*nx + (1-1)*nx*ny
        idx_face = i + (j-1)*nx

        select type (bc => BCs%face(FACE_ZMIN)%p)
        type is (BC_Dirichlet)
          h      = dz
          phi_bc = bc%value
          Jz_min(gg,idx_face) = (2.0d0*Dg(gg)/h) * (phi(gg,rowP) - phi_bc) * dx * dy

        type is (BC_Neumann)
          Jz_min(gg,idx_face) = 0.0d0

        type is (BC_Vacuum)
          h = dz
          Jz_min(gg,idx_face) = (2.0d0 * Dg(gg)) / (4.0d0 * Dg(gg) + h) * phi(gg,rowP) * dx * dy

        type is (BC_Albedo)
          h = dz
          gamma = (1.0d0/(2.0d0*Dg(gg)))*(1.0d0 - bc%alpha)/(1.0d0 + bc%alpha)
          Jz_min(gg,idx_face) = (2.0d0*Dg(gg)*gamma)/(2.0d0 + gamma*h) * phi(gg,rowP) * dx * dy * (1.0d0 - bc%alpha)/(1.0d0 + bc%alpha)

        class default
          if (nz > 1) then
            rowN = i + (j-1)*nx + (2-1)*nx*ny
            Jz_min(gg,idx_face) = Dg(gg) * (phi(gg,rowP) - phi(gg,rowN)) * dx * dy / dz
          else
            Jz_min(gg,idx_face) = 0.0d0
          end if
        end select
      end do
    end do

    !===========================================================
    ! ZMAX face: outward normal +z, outflow positive
    !===========================================================
    do j = 1, ny
      do i = 1, nx
        rowP     = i + (j-1)*nx + (nz-1)*nx*ny
        idx_face = i + (j-1)*nx

        select type (bc => BCs%face(FACE_ZMAX)%p)
        type is (BC_Dirichlet)
          h      = dz
          phi_bc = bc%value
          Jz_max(gg,idx_face) = (2.0d0*Dg(gg)/h) * (phi(gg,rowP) - phi_bc) * dx * dy

        type is (BC_Neumann)
          Jz_max(gg,idx_face) = 0.0d0

        type is (BC_Vacuum)
          h = dz
          Jz_max(gg,idx_face) = (2.0d0 * Dg(gg)) / (4.0d0 * Dg(gg) + h) * phi(gg,rowP) * dx * dy

        type is (BC_Albedo)
          h = dz
          gamma = (1.0d0/(2.0d0*Dg(gg)))*(1.0d0 - bc%alpha)/(1.0d0 + bc%alpha)
          Jz_max(gg,idx_face) = (2.0d0*Dg(gg)*gamma)/(2.0d0 + gamma*h) * phi(gg,rowP) * dx * dy * (1.0d0 - bc%alpha)/(1.0d0 + bc%alpha)
        class default
          if (nz > 1) then
            rowN = i + (j-1)*nx + (nz-2)*nx*ny
            Jz_max(gg,idx_face) = Dg(gg) * (phi(gg,rowP) - phi(gg,rowN)) * dx * dy / dz
          else
            Jz_max(gg,idx_face) = 0.0d0
          end if
        end select
      end do
    end do
  end do 

end subroutine compute_wall_currents


end module m_wall_current