module m_types
    use m_constants
    implicit none
    
    type t_patch
        integer, allocatable  :: cp_ids(:)
        integer               :: material_id
        real(dp), allocatable :: knots_xi(:), knots_eta(:), knots_zeta(:)
        integer               :: face_to_surface(6) ! Links faces to surface indices
    end type t_patch

    type t_surface
        integer, allocatable  :: cp_ids(:)
        integer               :: bc_id
        real(dp), allocatable :: knots_xi(:), knots_eta(:)
    end type t_surface

    type t_mesh
        integer               :: dim
        integer               :: n_groups
        integer               :: n_nodes          ! Total number of unique control points
        real(dp), allocatable :: nodes(:,:)       ! Global coordinates of control points
        real(dp), allocatable :: weights(:)       ! Weights of control points (NURBS)

        ! --- Original NURBS Patch Data ---
        type(t_patch),   allocatable :: patches(:)
        type(t_surface), allocatable :: surfaces(:)

        ! --- Solver Elements (Knot Spans) ---
        integer :: n_elems          
        integer :: order            ! Global polynomial order

        integer, allocatable  :: elems(:,:)        ! Global CP IDs per knot span element
        integer, allocatable  :: material_ids(:)
        integer, allocatable  :: elem_patch_id(:)  ! Parent patch reference
        integer, allocatable  :: elem_span_indices(:,:) ! (dim, n_elems) -> span_u, span_v, span_w
        integer, allocatable  :: elem_map_to_id(:,:,:,:) ! (patch_id, span_u, span_v, span_w) -> elem_id
        real(dp), allocatable :: elem_u_min(:), elem_u_max(:) 
        real(dp), allocatable :: elem_v_min(:), elem_v_max(:)
        real(dp), allocatable :: elem_w_min(:), elem_w_max(:)

        ! --- Transport Physics & Connectivity ---
        integer               :: n_faces_per_elem
        integer, allocatable  :: face_connectivity(:,:,:)
        real(dp), allocatable :: face_normals(:,:,:)
        integer, allocatable  :: upwind_idx(:,:,:)
        integer, allocatable  :: reflect_map(:,:,:)

        real(dp), allocatable :: elem_mass_matrix(:,:,:)
        real(dp), allocatable :: elem_stiffness_x(:,:,:)
        real(dp), allocatable :: elem_stiffness_y(:,:,:)
        real(dp), allocatable :: elem_stiffness_z(:,:,:)
        real(dp), allocatable :: face_mass_x(:,:,:,:)
        real(dp), allocatable :: face_mass_y(:,:,:,:)
        real(dp), allocatable :: face_mass_z(:,:,:,:)
        real(dp), allocatable :: basis_integrals_vol(:,:)
        real(dp), allocatable :: local_lu(:,:,:,:,:)
        integer,  allocatable :: local_pivots(:,:,:,:)
    end type t_mesh

    type t_finite
        integer :: order
        integer :: p_order, q_order, r_order
        integer :: n_basis
        integer :: n_nodes_per_face
        integer, allocatable :: face_node_map(:,:)
    end type t_finite

    type t_material
        character(len=32)     :: name
        real(dp), allocatable :: D(:)       ! Diffusion coefficient [cm]
        real(dp), allocatable :: SigmaT(:)  ! Total Cross Section [cm^-1]
        real(dp), allocatable :: SigmaR(:)  ! Removal Cross Section [cm^-1]
        real(dp), allocatable :: SigA(:)    ! Absorption Cross Section [cm^-1]
        real(dp), allocatable :: Nu(:)
        real(dp), allocatable :: SigF(:)
        real(dp), allocatable :: NuSigF(:)  ! Fission Production [cm^-1]
        real(dp), allocatable :: Chi(:)     ! Fission Spectrum [-]
        real(dp), allocatable :: Src(:)     ! Fixed Source [cm^-2]
        real(dp), allocatable :: SigmaS(:,:)! Scattering Matrix [from, to]
    end type t_material

    type :: t_quadrature
        integer :: n_points
        real(dp), dimension(:), allocatable :: xi, eta, zeta, weights
    end type t_quadrature

    type t_sn_quadrature
        integer :: n_angles, order
        real(dp), dimension(:,:), allocatable :: dirs
        real(dp), dimension(:), allocatable   :: weights
    end type t_sn_quadrature

end module m_types