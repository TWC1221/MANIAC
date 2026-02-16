module m_types
    use m_constants
    implicit none
    
    type :: t_mesh
        integer :: dim=2, n_nodes=0, n_elems=0, n_edges=0, nloc = 9
        real(dp),    allocatable :: nodes(:,:)         ! (2, n_nodes)
        integer,     allocatable :: elems(:,:)         ! (n_elems, 4)  (vtk=9)
        integer,     allocatable :: edges(:,:)         ! (n_edges, 2)  (vtk=3)
        integer,     allocatable :: mats(:)            ! (n_elems)
    end type t_mesh

    type :: t_finite
        integer :: order
        integer :: n_basis
        integer, allocatable :: p(:) ! Nodal Mappings
        real(dp), allocatable :: Xi(:) ! 1D Node positions

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