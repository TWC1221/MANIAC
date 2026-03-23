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



end module m_types