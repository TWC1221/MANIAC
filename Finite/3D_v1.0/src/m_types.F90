module m_types
    use m_constants
    implicit none

    type t_timer
        character(len=32) :: name
        real(dp) :: start_t
        real(dp) :: accumulated_t = 0.0_dp
        integer  :: calls = 0
        logical  :: active = .false.
    end type t_timer
    
    type t_mesh
        integer :: n_nodes, n_elems, n_edges, nloc, dim
        real(dp), allocatable :: nodes(:,:)
        integer,  allocatable :: elems(:,:)
        integer,  allocatable :: edges(:,:)
        integer,  allocatable :: mats(:)       ! For Quads
        integer,  allocatable :: edge_mats(:)  ! For Lines/Boundaries
    end type t_mesh

    type :: t_finite
        integer :: order
        integer :: n_basis
        integer, allocatable :: p(:)           ! Nodal Mappings
        real(dp), allocatable :: Xi(:)         ! 1D Node positions

        real(dp), allocatable :: N(:,:)        ! Shape functions
        real(dp), allocatable :: dN_dxi(:,:)   ! Derivatives w.r.t xi
        real(dp), allocatable :: dN_deta(:,:)  ! Derivatives w.r.t eta
        real(dp), allocatable :: dN_dzeta(:,:)  ! Derivatives w.r.t zeta
        real(dp), allocatable :: N_mat(:,:,:)  ! Pre-computed outer products

        real(dp), allocatable :: N_B(:,:)      ! Boundary Shape functions
        real(dp), allocatable :: dN_B_xi(:,:)  ! Boundary derivatives xi
        real(dp), allocatable :: dN_B_eta(:,:) ! Boundary derivatives eta
    end type t_finite

    type :: t_mat
        real(dp), allocatable :: val(:)
        integer,  allocatable :: col(:)
        integer,  allocatable :: row_ptr(:)
    end type t_mat

    type :: t_vec
        real(dp), allocatable :: vec(:)
    end type t_vec

    type :: t_bc_config
        integer  :: mat_id 
        integer  :: bc_type  
        real(dp) :: value    
    end type t_bc_config


end module m_types