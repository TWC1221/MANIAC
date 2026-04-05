module m_types
    use omp_lib
    use m_constants
    implicit none
    private
    public :: t_mesh, t_finite, t_material, t_quadrature, t_sn_quadrature, &
              t_edge_sort, t_timer, timers, num_timers, &
              int_to_str, timer_start, timer_stop, dp

    type t_mesh
        integer               :: n_nodes
        integer               :: n_elems
        integer               :: n_edges
        integer               :: n_faces_per_elem
        integer               :: dim

        real(dp), allocatable :: nodes(:,:)
        integer,  allocatable :: elems(:,:)
        integer,  allocatable :: edges(:,:)
        integer,  allocatable :: material_ids(:)      
        integer,  allocatable :: boundary_ids(:)  
        integer,  allocatable :: pin_ids(:)

        integer, allocatable  :: face_connectivity(:,:,:)       
        real(dp), allocatable :: face_normals(:,:,:)        ! (dim, n_faces_per_elem, n_elems)
        real(dp), allocatable :: face_mass_x(:,:,:,:)           ! (n_basis, n_basis, n_faces, n_elems)
        real(dp), allocatable :: face_mass_y(:,:,:,:)           ! Normal-weighted mass matrices
        real(dp), allocatable :: face_mass_z(:,:,:,:)
        real(dp), allocatable :: elem_mass_matrix(:,:,:)    ! (n_basis, n_basis, n_elems)
        real(dp), allocatable :: basis_integrals_vol(:,:)       ! (n_basis, n_elems)
        real(dp), allocatable :: elem_stiffness_x(:,:,:)        ! (n_basis, n_basis, n_elems) gradient x
        real(dp), allocatable :: elem_stiffness_y(:,:,:)        ! (n_basis, n_basis, n_elems) gradient y
        real(dp), allocatable :: elem_stiffness_z(:,:,:)
        
        real(dp), allocatable :: local_lu(:,:,:,:,:)        ! local_lu: (Basis, Basis, Element, Angle, Group)
        integer,  allocatable :: local_pivots(:,:,:,:)      ! local_pivots: (Basis, Element, Angle, Group)
        integer,  allocatable :: reflect_map(:,:,:)         ! reflect_map: (Angle, Face, Element) -> Index of reflected angle
        integer, allocatable  :: upwind_idx(:,:,:)          ! upwind_idx(face_node, face, element) -> global angular flux index
    end type t_mesh

    type :: t_finite
        integer               :: order
        integer               :: n_basis
        integer               :: n_nodes_per_face
        real(dp), allocatable :: node_roots(:)                      ! 1D Node positions
        integer, allocatable  :: face_node_map(:,:)         ! (n_nodes_per_face, 4)
        real(dp), allocatable :: basis_at_quad(:,:)                     ! Shape functions
        real(dp), allocatable :: dbasis_dxi(:,:)                ! Derivatives w.r.t xi
        real(dp), allocatable :: dbasis_deta(:,:)               ! Derivatives w.r.t eta
        real(dp), allocatable :: dbasis_dzeta(:,:)              ! Derivatives w.r.t zeta
        real(dp), allocatable :: basis_at_bound_quad(:,:)      
        real(dp), allocatable :: dbasis_bound_dxi(:,:)
        real(dp), allocatable :: dbasis_bound_deta(:,:)
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

    type t_timer
        character(len=32) :: name
        real(dp) :: start_t
        real(dp) :: accumulated_t = 0.0_dp
        integer  :: calls = 0
        logical  :: active = .false.
    end type t_timer
    
    ! --- Timing Parameters ---
    integer, parameter :: MAX_TIMERS = 50
    type(t_timer), save :: timers(MAX_TIMERS)
    integer, save       :: num_timers = 0

contains

    pure function int_to_str(i) result(res)
        integer, intent(in) :: i
        character(len=12) :: res
        write(res, '(I0)') i
        res = adjustl(res)
    end function int_to_str

    subroutine timer_start(label)
        character(len=*), intent(in) :: label
        integer :: i
        do i = 1, num_timers
            if (trim(timers(i)%name) == trim(label)) then
                if (.not. timers(i)%active) then
                    timers(i)%start_t = omp_get_wtime()
                    timers(i)%active = .true.
                    timers(i)%calls = timers(i)%calls + 1
                end if
                return
            end if
        end do
        if (num_timers < MAX_TIMERS) then
            num_timers = num_timers + 1
            timers(num_timers)%name = label
            timers(num_timers)%start_t = omp_get_wtime()
            timers(num_timers)%active = .true.
            timers(num_timers)%calls = 1
        end if
    end subroutine timer_start
    
    subroutine timer_stop(label)
        character(len=*), intent(in) :: label
        integer :: i
        real(dp) :: t_now
        t_now = omp_get_wtime()
        do i = 1, num_timers
            if (trim(timers(i)%name) == trim(label)) then
                if (timers(i)%active) then
                    timers(i)%accumulated_t = timers(i)%accumulated_t + (t_now - timers(i)%start_t)
                    timers(i)%active = .false.
                end if
                return
            end if
        end do
    end subroutine timer_stop

end module m_types