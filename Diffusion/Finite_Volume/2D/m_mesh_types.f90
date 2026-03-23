module mesh_types
  implicit none
  private
  public :: MeshGrid
  type :: MeshGrid
    integer :: nx, ny, nz
    real(8) :: X_domain, Y_domain, Z_domain
    integer :: N
    real(8) :: dx, dy, dz
    real(8) :: dV
  contains
    procedure :: init
  end type MeshGrid

contains
  subroutine init(mesh, x_domain, y_domain, z_domain, nx, ny, nz)
    class(MeshGrid), intent(inout) :: mesh
    real(8), intent(in) :: x_domain, y_domain, z_domain
    integer, intent(in) :: nx, ny, nz

    mesh%nx = nx
    mesh%ny = ny
    mesh%nz = nz
    mesh%x_domain = x_domain
    mesh%y_domain = y_domain
    mesh%z_domain = z_domain
    mesh%dx = x_domain/nx
    mesh%dy = y_domain/ny
    mesh%dz = z_domain/nz
    mesh%N = nx*ny*nz
    mesh%dV = x_domain/nx*y_domain/ny*z_domain/nz
  end subroutine

end module 
