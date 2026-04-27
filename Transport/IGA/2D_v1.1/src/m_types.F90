module m_types
    use m_constants
    implicit none
    
    type t_mesh
        integer               :: n_nodes
        integer               :: n_elems
        integer               :: n_edges
        integer               :: n_faces_per_elem
        integer               :: order
        integer               :: nloc               ! Max control points per patch
        integer               :: dim

        real(dp), allocatable :: nodes(:,:)
        real(dp), allocatable :: weights(:)         ! Control point weights
        integer,  allocatable :: elems(:,:)
        integer,  allocatable :: edges(:,:)
        integer,  allocatable :: material_ids(:)      
        integer,  allocatable :: boundary_ids(:)  

        real(dp), allocatable :: span_range(:,:)            ! (4, n_elems) -> u_min, u_max, v_min, v_max
        integer,  allocatable :: n_knots_xi_patch(:), n_knots_eta_patch(:)
        integer,  allocatable :: n_cp_edge(:), n_knots_edge(:)
        real(dp), allocatable :: edge_knots(:,:)
        integer,  allocatable :: n_cp_xi(:), n_cp_eta(:)
        real(dp), allocatable :: knot_vectors_xi(:,:), knot_vectors_eta(:,:)
        
        integer, allocatable  :: face_connectivity(:,:,:)       
        real(dp), allocatable :: face_normals(:,:,:)        ! (dim, n_faces_per_elem, n_elems)
        real(dp), allocatable :: face_mass_x(:,:,:,:)           ! (n_basis, n_basis, n_faces, n_elems)
        real(dp), allocatable :: face_mass_y(:,:,:,:)           ! Normal-weighted mass matrices
        real(dp), allocatable :: elem_mass_matrix(:,:,:)    ! (n_basis, n_basis, n_elems)
        real(dp), allocatable :: basis_integrals_vol(:,:)       ! (n_basis, n_elems)
        real(dp), allocatable :: elem_stiffness_x(:,:,:)        ! (n_basis, n_basis, n_elems) gradient x
        real(dp), allocatable :: elem_stiffness_y(:,:,:)        ! (n_basis, n_basis, n_elems) gradient y
        
        real(dp), allocatable :: local_lu(:,:,:,:,:)        ! local_lu: (Basis, Basis, Element, Angle, Group)
        integer,  allocatable :: local_pivots(:,:,:,:)      ! local_pivots: (Basis, Element, Angle, Group)
        integer,  allocatable :: reflect_map(:,:,:)         ! reflect_map: (Angle, Face, Element) -> Index of reflected angle
        integer, allocatable  :: upwind_idx(:,:,:)          ! upwind_idx(face_node, face, element) -> global angular flux index
    end type t_mesh

    type :: t_finite
        integer               :: order, p_order, q_order
        integer               :: n_basis
        integer               :: n_nodes_per_face                   ! Max nodes on any face across mesh
        integer, allocatable  :: face_node_map(:,:,:)       ! (n_nodes_per_face, 4, n_elems)
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
        real(dp), allocatable :: dir(:,:)
    end type t_quadrature

    type t_sn_quadrature
        integer :: n_angles, order
        real(dp), dimension(:), allocatable   :: mu, eta, zeta
        real(dp), dimension(:,:), allocatable :: dirs
        real(dp), dimension(:), allocatable   :: weights
    end type t_sn_quadrature

    type :: t_edge_sort
        integer :: n1, n2
    end type t_edge_sort

end module m_types