module m_assembly
    implicit none
    contains
subroutine assemble_3D_diffusion_matrix_g(D, Sigma_r, A, mesh, BCs, use_adjoint, periodic_x, periodic_y, periodic_z)
    use CSR_types, only: CSRMatrix
    use Mesh_types, only: MeshGrid
    use BC_types, only: BoundarySet, BC_Base, &
                        FACE_XMIN, FACE_XMAX, FACE_YMIN, FACE_YMAX, FACE_ZMIN, FACE_ZMAX
    implicit none

    type(MeshGrid), intent(in) :: mesh
    type(CSRMatrix), intent(out) :: A
    type(BoundarySet), intent(in) :: BCs
    
    real(8), intent(in) :: D, Sigma_r
    logical, intent(in) :: use_adjoint
    logical, intent(in), optional :: periodic_x, periodic_y, periodic_z 
  
    integer :: ii, jj, kk, row, ptr
    real(8) :: aL, aR, aB, aT, aBack, aFront, aC, Dloc
    logical :: on_xmin, on_xmax, on_ymin, on_ymax, on_zmin, on_zmax, px, py, pz

    ! allocate CSR arrays
    allocate(A%ia(mesh%N+1))
    allocate(A%ja(7*mesh%N))
    allocate(A%aa(7*mesh%N))

    ! Periodic flags (default: false)
    px =.false.; if (present(periodic_x)) px = periodic_x
    py =.false.; if (present(periodic_y)) py = periodic_y
    pz =.false.; if (present(periodic_z)) pz = periodic_z

    Dloc = D
    if (use_adjoint) Dloc = D

    ! initialize
    A%ia = 0
    A%ja = 0
    A%aa = 0

    ptr = 1
    A%ia(1) = 1  ! 1-based CSR convention: IA(1) = index of first element = 1

    ! loop rows in k-j-i order (consistent with your flattening)
    do kk = 1, mesh%nz
      do jj = 1, mesh%ny
        do ii = 1, mesh%nx

          ! flatten to row index (1-based)
          row = ii + (jj-1)*mesh%nx + (kk-1)*mesh%nx*mesh%ny

          ! stencil coefficients (off-diagonals)
          aL = -Dloc * (mesh%dy*mesh%dz)/mesh%dx
          aR = aL
          aB = -Dloc * (mesh%dx*mesh%dz)/mesh%dy
          aT = aB
          aBack  = -Dloc * (mesh%dx*mesh%dy)/mesh%dz
          aFront = aBack

          ! diagonal contribution (positive) - note sign convention: off-diagonals are negative
          aC = Sigma_r*mesh%dV - (aL + aR + aB + aT + aBack + aFront)

          !--------------------------------------------------
          ! Identify if this cell is on geometric faces
          !--------------------------------------------------
          on_xmin = (ii == 1)
          on_xmax = (ii == mesh%nx)
          on_ymin = (jj == 1)
          on_ymax = (jj == mesh%ny)
          on_zmin = (kk == 1)
          on_zmax = (kk == mesh%nz)

          !--------------------------------------------------
          ! Apply BCs ONLY on non-periodic directions
          ! For periodic directions, faces are connected, no BC.
          !--------------------------------------------------
          if (.not. px) then
            call apply_face_bc_oop(BCs%face(FACE_XMIN)%p, on_xmin, D, mesh%dx, mesh%dy, mesh%dz, aC, aL)
            call apply_face_bc_oop(BCs%face(FACE_XMAX)%p, on_xmax, D, mesh%dx, mesh%dy, mesh%dz, aC, aR)
          end if

          if (.not. py) then
            call apply_face_bc_oop(BCs%face(FACE_YMIN)%p, on_ymin, D, mesh%dx, mesh%dy, mesh%dz, aC, aB)
            call apply_face_bc_oop(BCs%face(FACE_YMAX)%p, on_ymax, D, mesh%dx, mesh%dy, mesh%dz, aC, aT)
          end if

          if (.not. pz) then
            call apply_face_bc_oop(BCs%face(FACE_ZMIN)%p, on_zmin, D, mesh%dx, mesh%dy, mesh%dz, aC, aBack)
            call apply_face_bc_oop(BCs%face(FACE_ZMAX)%p, on_zmax, D, mesh%dx, mesh%dy, mesh%dz, aC, aFront)
          end if

          ! --- center (column = row) ---
          if (ptr > 7*mesh%N) then
            print *, "ERROR: ptr exceeds est_nnz in assemble_3D_diffusion_matrix_g, increase est_nnz"
            stop
          end if
          A%ja(ptr) = row       ! column index of diagonal
          A%aa(ptr) = aC
          ptr = ptr + 1

          ! LEFT neighbour (col = row-1)
          if (ii > 1) then
            A%ja(ptr) = row - 1
            A%aa(ptr) = aL
            ptr = ptr + 1
          else if (px) then
            ! periodic: ii=1 left neighbour is ii=nx
            A%ja(ptr) = row + (mesh%nx - 1)
            A%aa(ptr) = aL
            ptr = ptr + 1
          end if

          ! RIGHT neighbour
          if (ii < mesh%nx) then
            A%ja(ptr) = row + 1
            A%aa(ptr) = aR
            ptr = ptr + 1
          else if (px) then
            ! periodic: ii=1 left neighbour is ii=nx
            A%ja(ptr) = row - (mesh%nx - 1)
            A%aa(ptr) = aL
            ptr = ptr + 1
          end if

          ! BOTTOM neighbour (j-1)
          if (jj > 1) then
            A%ja(ptr) = row - mesh%nx
            A%aa(ptr) = aB
            ptr = ptr + 1
          else if (py) then
            ! periodic: jj=1 bottom neighbour is jj=ny
            A%ja(ptr) = row + mesh%nx * (mesh%ny - 1)
            A%aa(ptr) = aB
            ptr = ptr + 1
          end if

          ! TOP neighbour (j+1)
          if (jj < mesh%ny) then
            A%ja(ptr) = row + mesh%nx
            A%aa(ptr) = aT
            ptr = ptr + 1
          else if (py) then
            ! periodic: jj=1 bottom neighbour is jj=ny
            A%ja(ptr) = row - mesh%nx * (mesh%ny - 1)
            A%aa(ptr) = aB
            ptr = ptr + 1
          end if

          ! BACK neighbour (k-1)
          if (kk > 1) then
            A%ja(ptr) = row - mesh%nx*mesh%ny
            A%aa(ptr) = aBack
            ptr = ptr + 1
          else if (pz) then
            ! periodic: kk=1 back neighbour is kk=nz
            A%ja(ptr) = row + mesh%nx*mesh%ny*(mesh%nz - 1)
            A%aa(ptr) = aBack
            ptr = ptr + 1
          end if

          ! FRONT neighbour (k+1)
          if (kk < mesh%nz) then
            A%ja(ptr) = row + mesh%nx*mesh%ny
            A%aa(ptr) = aFront
            ptr = ptr + 1
          else if (pz) then
            ! periodic: kk=nz front neighbour is kk=1
            A%ja(ptr) = row - mesh%nx*mesh%ny*(mesh%nz - 1)
            A%aa(ptr) = aFront
            ptr = ptr + 1
          end if

          ! set IA for next row (first index after this row's entries)
          A%ia(row+1) = ptr

        end do
      end do
    end do

    ! final NNZ index pointer: ptr is (nnz + 1)
    A%ia(mesh%N+1) = ptr

    ! trim JA and AA to actual nnz
    if (ptr - 1 < size(A%ja)) A%ja = A%ja(1:ptr-1)
    if (ptr - 1 < size(A%aa)) A%aa = A%aa(1:ptr-1)

  contains

    subroutine apply_face_bc_oop(bc, is_on_face, D_in, h_in, h1_in, h2_in, aC_inout, aN_inout)
      class(BC_Base), intent(in) :: bc
      logical,        intent(in) :: is_on_face
      real(8),        intent(in) :: D_in, h_in, h1_in, h2_in
      real(8),        intent(inout) :: aC_inout, aN_inout
      call bc%apply_matrix(is_on_face, D_in, h_in, h1_in, h2_in, aC_inout, aN_inout)
    end subroutine apply_face_bc_oop

  end subroutine assemble_3D_diffusion_matrix_g

end module m_assembly