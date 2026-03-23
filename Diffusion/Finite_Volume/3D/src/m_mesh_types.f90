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

  subroutine init(mesh, filename)
    class(MeshGrid), intent(inout) :: mesh
    character(len=*), intent(in) :: filename

    integer :: unit, ios
    character(len=256) :: line
    real(8) :: x_dom, y_dom, z_dom
    integer :: nx_in, ny_in, nz_in
    logical :: found

    ! Open the file
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print*, 'Error opening mesh file: ', filename
        stop
    end if

    ! --- Read domain sizes ---
    found = .false.
    do
        read(unit,'(A)', iostat=ios) line
        if (ios /= 0) then
            print*, 'Unexpected EOF while reading domain sizes!'
            stop
        end if
        line = adjustl(line)
        if (len_trim(line) == 0) cycle    ! skip blank lines
        if (line(1:1) == '#') cycle      ! skip comment lines
        read(line,*, iostat=ios) x_dom, y_dom, z_dom
        if (ios == 0) then
            found = .true.
            exit
        end if
    end do
    if (.not. found) then
        print*, 'Could not read domain sizes from file.'
        stop
    end if

    ! --- Read number of nodes ---
    found = .false.
    do
        read(unit,'(A)', iostat=ios) line
        if (ios /= 0) then
            print*, 'Unexpected EOF while reading discretization!'
            stop
        end if
        line = adjustl(line)
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#') cycle
        read(line,*, iostat=ios) nx_in, ny_in, nz_in
        if (ios == 0) then
            found = .true.
            exit
        end if
    end do
    if (.not. found) then
        print*, 'Could not read discretization (nx,ny,nz) from file.'
        stop
    end if

    close(unit)

    ! --- Assign to mesh object ---
    mesh%X_domain = x_dom
    mesh%Y_domain = y_dom
    mesh%Z_domain = z_dom
    mesh%nx = nx_in
    mesh%ny = ny_in
    mesh%nz = nz_in
    mesh%dx = x_dom / nx_in
    mesh%dy = y_dom / ny_in
    mesh%dz = z_dom / nz_in
    mesh%N = nx_in * ny_in * nz_in
    mesh%dV = mesh%dx * mesh%dy * mesh%dz

  end subroutine init

end module mesh_types
