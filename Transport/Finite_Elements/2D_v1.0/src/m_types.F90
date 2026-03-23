module m_types
    use m_constants
    implicit none
    
    type t_mesh
        integer :: n_nodes, n_elems, n_edges, dim, n_faces_per_elem
        real(dp), allocatable :: nodes(:,:)
        integer,  allocatable :: elems(:,:)
        integer,  allocatable :: edges(:,:)
        integer,  allocatable :: mats(:)       ! For Quads
        integer,  allocatable :: edge_mats(:)  ! For Lines/Boundaries

    ! Consolidated per-face metadata: dims = (4, n_faces_per_elem, n_elems)
    ! face_info(1, f, e) = neighbor element id (or -1)
    ! face_info(2, f, e) = neighbor face index
    ! face_info(3, f, e) = orientation (+1/-1/0)
    ! face_info(4, f, e) = boundary id (0 = internal)
        integer, allocatable  :: face_info(:,:,:)
        real(dp), allocatable :: face_normals(:,:,:) ! (dim, n_faces_per_elem, n_elems)
        real(dp), allocatable :: face_mass_matrices(:,:,:,:)

        ! Dimensions: (Basis, Element, Angle, Group)
        real(dp), allocatable :: local_LU(:,:,:,:,:) 
        integer,  allocatable :: local_pivots(:,:,:,:)

        real(dp), allocatable :: elem_mass_matrix(:,:,:) ! (n_basis, n_basis, n_elems)
        real(dp), allocatable :: elem_integral_v(:,:)   ! (n_basis, n_elems)
        real(dp), allocatable :: elem_stiff_x(:,:,:)    ! (n_basis, n_basis, n_elems) gradient x
        real(dp), allocatable :: elem_stiff_y(:,:,:)    ! (n_basis, n_basis, n_elems) gradient y
    end type t_mesh

    type :: t_finite
        integer :: order
        integer :: n_basis
        integer, allocatable :: p(:)           ! Nodal Mappings
        real(dp), allocatable :: Xi(:)         ! 1D Node positions
        integer :: n_nodes_per_face
        integer, allocatable :: face_node_map(:,:) ! (n_nodes_per_face, 4)

        real(dp), allocatable :: N(:,:)        ! Shape functions
        real(dp), allocatable :: dN_dxi(:,:)   ! Derivatives w.r.t xi
        real(dp), allocatable :: dN_deta(:,:)  ! Derivatives w.r.t eta
        real(dp), allocatable :: N_B(:,:)      
        real(dp), allocatable :: dN_B(:,:)     
    end type t_finite

    type t_material
        character(len=32) :: name
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

    type :: t_PCG_CSR
        real(dp), allocatable :: val(:)
        integer,  allocatable :: col(:)
        integer,  allocatable :: row_ptr(:)
    end type t_PCG_CSR

    type :: t_vec
        real(dp), allocatable :: vec(:)
    end type t_vec

    type :: t_bc_config
        integer  :: mat_id 
        integer  :: bc_type  
        real(dp) :: value    
    end type t_bc_config

    type t_sweep_plan
        integer, allocatable :: sweep_order(:,:) ! (n_elems, n_angles)
        integer, allocatable :: colors(:,:)      ! (n_elems, n_angles)
        integer, allocatable :: num_colors(:)    ! (n_angles)
    end type t_sweep_plan

    type :: t_quadrature
    real(dp), dimension(:), allocatable :: Xi, Eta, Zeta, W
    real(dp), allocatable :: dir(:,:)
    integer :: NoPoints
end type t_quadrature

type t_sn_quadrature
    integer :: NoAngles, P, Q, order
    real(dp), dimension(:), allocatable   :: mu, eta, zeta
    real(dp), dimension(:,:), allocatable :: Angles
    real(dp), dimension(:), allocatable   :: w
    real(dp), dimension(:), allocatable   :: alpha1D
    real(dp), dimension(:,:), allocatable :: alpha2D

    integer, allocatable :: pp(:), qq(:)            ! maps: angle nn -> level p, order q
    integer, allocatable :: zeta_pp(:), qq_len(:)   ! per-level metadata

    real(dp), allocatable :: mu_pq(:,:)    ! (P,Q): μ at (pp,qq)
    real(dp), allocatable :: zeta_p(:)  ! (P,Q): ζ at (pp,qq)
    real(dp), allocatable :: w_pq(:,:)     ! (P,Q): weight at (pp,qq)
    integer, allocatable :: n_pq(:,:)
    integer, allocatable :: reflect_x(:), reflect_y(:)
end type t_sn_quadrature

end module m_types