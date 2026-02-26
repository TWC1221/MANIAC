module m_types
    use m_constants
    implicit none
    
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
        real(dp), allocatable :: N_mat(:,:,:)  ! Pre-computed outer products (Mass matrix template)

        real(dp), allocatable :: N_B(:,:)      
        real(dp), allocatable :: dN_B(:,:)     
    end type t_finite

    type :: t_mat
        real(dp), allocatable :: val(:)
        integer,  allocatable :: col(:)
        integer,  allocatable :: row_ptr(:)
    end type t_mat

    type :: t_vec
        real(dp), allocatable :: vec(:)
    end type t_vec


end module m_types