module BC_types
    implicit none
    private

    !===========================================================
    ! Public face identifiers
    !===========================================================
    integer, public, parameter :: FACE_XMIN = 1
    integer, public, parameter :: FACE_XMAX = 2
    integer, public, parameter :: FACE_YMIN = 3
    integer, public, parameter :: FACE_YMAX = 4
    integer, public, parameter :: FACE_ZMIN = 5
    integer, public, parameter :: FACE_ZMAX = 6
    integer, public, parameter :: N_FACES  = 6

    !===========================================================
    ! Abstract base type for all boundary conditions
    !===========================================================
    type, abstract, public :: BC_Base
    contains
        procedure(apply_bc_if), deferred :: apply_matrix
    end type BC_Base

    !-----------------------------------------------------------
    ! Polymorphic BC interface
    !-----------------------------------------------------------
    abstract interface
        subroutine apply_bc_if(self, is_on_face, D, h, h1, h2, aC, aN)
            import :: BC_Base
            class(BC_Base), intent(in)    :: self
            logical,        intent(in)    :: is_on_face
            real(8),        intent(in)    :: D, h, h1, h2
            real(8),        intent(inout) :: aC, aN
        end subroutine apply_bc_if
    end interface

    !===========================================================
    ! Concrete BC types
    !===========================================================

    !---------------------------
    ! "None" BC: does nothing
    !---------------------------
    type, extends(BC_Base), public :: BC_None
        character(len=30) :: description = "N/A"
    contains
        procedure :: apply_matrix => apply_matrix_none
    end type BC_None

    !---------------------------
    ! Neumann BC:
    ! prescribed outward normal flux = self%flux
    ! If flux = 0 => homogeneous Neumann
    ! Matrix effect: aC = aC + aN (for homogeneous case)
    !---------------------------
    type, extends(BC_Base), public :: BC_Neumann
        real(8)        :: flux = 0.0d0        ! outward normal flux
        character(len=30) :: description = "Neumann"
    contains
        procedure :: apply_matrix => apply_neumann_matrix
    end type BC_Neumann

    !---------------------------
    ! Dirichlet BC:
    ! Enforced value in RHS, matrix modification:
    ! aC = aC - aN
    !---------------------------
    type, extends(BC_Base), public :: BC_Dirichlet
        real(8)        :: value = 0.0d0       ! prescribed field value
        character(len=30) :: description = "Dirichlet"
    contains
        procedure :: apply_matrix => apply_dirichlet_matrix
    end type BC_Dirichlet

    !---------------------------
    ! Vacuum BC:
    ! aC = aC + aN + (2*D*h1*h2)/(4*D + h)
    !---------------------------
    type, extends(BC_Base), public :: BC_Vacuum
        character(len=30) :: description = "Vacuum"
    contains
        procedure :: apply_matrix => apply_vacuum_matrix
    end type BC_Vacuum

    !---------------------------
    ! Albedo BC:
    ! gamma = (1/(2D))*(1-alpha)/(1+alpha)
    ! aC = aC + aN + (2*D*gamma*h1*h2)/(2 + gamma*h)
    !---------------------------
    type, extends(BC_Base), public :: BC_Albedo
        real(8)        :: alpha = 0.0d0
        character(len=30) :: description = "Albedo"
    contains
        procedure :: apply_matrix => apply_albedo_matrix
    end type BC_Albedo

    !===========================================================
    ! Container types
    !===========================================================
    type :: BC_ptr
        class(BC_Base), pointer :: p => null()
    end type BC_ptr

    type, public :: BoundarySet
        type(BC_ptr), allocatable :: face(:)   ! dimension N_FACES
    end type BoundarySet

    !===========================================================
    ! Public helper procedures (optional but convenient)
    !===========================================================
    public :: init_boundary_set, set_bc_dirichlet, set_bc_neumann, &
              set_bc_vacuum, set_bc_albedo

contains

    !===========================================================
    ! Matrix-level BC applications
    !===========================================================

    !---------------- Neumann ----------------
    subroutine apply_neumann_matrix(self, is_on_face, D, h, h1, h2, aC, aN)
        class(BC_Neumann), intent(in)    :: self
        logical,           intent(in)    :: is_on_face
        real(8),           intent(in)    :: D, h, h1, h2
        real(8),           intent(inout) :: aC, aN

        if (.not. is_on_face) return

        ! Homogeneous Neumann (zero-normal-flux): aC += aN
        ! Inhomogeneous part (self%flux) should be handled in RHS assembly.
        aC = aC + aN
    end subroutine apply_neumann_matrix

    !---------------- Dirichlet ----------------
    subroutine apply_dirichlet_matrix(self, is_on_face, D, h, h1, h2, aC, aN)
        class(BC_Dirichlet), intent(in)    :: self
        logical,             intent(in)    :: is_on_face
        real(8),             intent(in)    :: D, h, h1, h2
        real(8),             intent(inout) :: aC, aN

        if (.not. is_on_face) return

        ! Matrix: remove neighbour coupling and push it into RHS
        aC = aC - aN !+ 10**6
        ! self%value used separately in RHS assembly
    end subroutine apply_dirichlet_matrix

    !---------------- Vacuum --------------------
    subroutine apply_vacuum_matrix(self, is_on_face, D, h, h1, h2, aC, aN)
        class(BC_Vacuum), intent(in)    :: self
        logical,          intent(in)    :: is_on_face
        real(8),          intent(in)    :: D, h, h1, h2
        real(8),          intent(inout) :: aC, aN
        real(8) :: term

        if (.not. is_on_face) return

        term = (2.0d0 * D * h1 * h2) / (4.0d0 * D + h)
        aC   = aC + aN + term
    end subroutine apply_vacuum_matrix

    !---------------- Albedo --------------------
    subroutine apply_albedo_matrix(self, is_on_face, D, h, h1, h2, aC, aN)
        class(BC_Albedo), intent(in)    :: self
        logical,          intent(in)    :: is_on_face
        real(8),          intent(in)    :: D, h, h1, h2
        real(8),          intent(inout) :: aC, aN
        real(8) :: gamma, term

        if (.not. is_on_face) return

        gamma = (1.0d0 / (2.0d0 * D)) * (1.0d0 - self%alpha) / (1.0d0 + self%alpha)
        term  = (2.0d0 * D * gamma * h1 * h2) / (2.0d0 + gamma * h)

        aC = aC + aN + term
    end subroutine apply_albedo_matrix

    !---------------- None: do nothing ----------
    subroutine apply_matrix_none(self, is_on_face, D, h, h1, h2, aC, aN)
        class(BC_None), intent(in)    :: self
        logical,        intent(in)    :: is_on_face
        real(8),        intent(in)    :: D, h, h1, h2
        real(8),        intent(inout) :: aC, aN

        if (.not. is_on_face) return
    end subroutine apply_matrix_none

    !===========================================================
    ! Helper routines for BoundarySet management
    !===========================================================

    !---------------- Initialize BoundarySet ----------------
    subroutine init_boundary_set(BCs, default_type)
        ! Initialize BCs%face(:) and assign a default BC type for all faces.
        ! default_type: 1=None, 2=Neumann, 3=Dirichlet, 4=Vacuum, 5=Albedo
        type(BoundarySet), intent(inout) :: BCs
        integer,           intent(in), optional :: default_type

        integer :: f, dt

        dt = 4   ! default to Vacuum if not supplied
        if (present(default_type)) dt = default_type

        if (.not. allocated(BCs%face)) then
            allocate(BCs%face(N_FACES))
        end if

        do f = FACE_XMIN, FACE_ZMAX
            call allocate_face_bc(BCs%face(f)%p, dt)
        end do
    end subroutine init_boundary_set

    !---------------- Set Dirichlet on one face --------------
    subroutine set_bc_dirichlet(BCs, face_id, value)
        type(BoundarySet), intent(inout) :: BCs
        integer,           intent(in)    :: face_id
        real(8),           intent(in)    :: value

        call allocate_face_bc(BCs%face(face_id)%p, 3)
        select type(p => BCs%face(face_id)%p)
        type is (BC_Dirichlet)
            p%value = value
        end select
    end subroutine set_bc_dirichlet

    !---------------- Set Neumann on one face ----------------
    subroutine set_bc_neumann(BCs, face_id, flux)
        type(BoundarySet), intent(inout) :: BCs
        integer,           intent(in)    :: face_id
        real(8),           intent(in)    :: flux

        call allocate_face_bc(BCs%face(face_id)%p, 2)
        select type(p => BCs%face(face_id)%p)
        type is (BC_Neumann)
            p%flux = flux
        end select
    end subroutine set_bc_neumann

    !---------------- Set Vacuum on one face -----------------
    subroutine set_bc_vacuum(BCs, face_id)
        type(BoundarySet), intent(inout) :: BCs
        integer,           intent(in)    :: face_id

        call allocate_face_bc(BCs%face(face_id)%p, 4)
    end subroutine set_bc_vacuum

    !---------------- Set Albedo on one face -----------------
    subroutine set_bc_albedo(BCs, face_id, alpha)
        type(BoundarySet), intent(inout) :: BCs
        integer,           intent(in)    :: face_id
        real(8),           intent(in)    :: alpha

        call allocate_face_bc(BCs%face(face_id)%p, 5)
        select type(p => BCs%face(face_id)%p)
        type is (BC_Albedo)
            p%alpha = alpha
        end select
    end subroutine set_bc_albedo

    !---------------- Internal allocation helper -------------
    subroutine allocate_face_bc(ptr, bc_type)
        class(BC_Base), pointer :: ptr
        integer,        intent(in) :: bc_type
        ! bc_type: 1=None, 2=Neumann, 3=Dirichlet, 4=Vacuum, 5=Albedo

        if (associated(ptr)) then
            deallocate(ptr)
        end if

        select case (bc_type)
        case (1)
            allocate(BC_None     :: ptr)
        case (2)
            allocate(BC_Neumann  :: ptr)
        case (3)
            allocate(BC_Dirichlet:: ptr)
        case (4)
            allocate(BC_Vacuum   :: ptr)
        case (5)
            allocate(BC_Albedo   :: ptr)
        case default
            allocate(BC_None     :: ptr)
        end select
    end subroutine allocate_face_bc

end module BC_types