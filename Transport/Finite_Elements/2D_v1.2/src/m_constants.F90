module m_constants
    implicit none
  
    !! Private definition
  
    !! Precision Definitions
    Integer, Parameter        :: sp=selected_Real_Kind(4)   !! single precision
    Integer, Parameter        :: dp=selected_Real_Kind(8)  !! double precision
    Integer, Parameter        :: qp=selected_Real_Kind(16)  !! quad precision
  
    !! Constant Values
    Real(Kind=dp),  Parameter :: PI=4.D0*DATAN(1.0D0)
    Real(Kind=dp),  Parameter :: dp_EPSILON = 1.0E-12_dp
    Real(Kind=dp),  Parameter :: VSMALL_NUMBER = 1.0E-9_dp
    Real(Kind=dp),  Parameter :: SMALL_NUMBER = 1.0E-6_dp
    Real(Kind=dp),  Parameter :: LARGE_NUMBER = 1.0E+6_dp
    ! Real(Kind=dp),  Parameter :: VLARGE_NUMBER = 1.0E+9_dp
    real(Kind=dp),  parameter :: VLARGE_NUMBER = 1.0E+9_dp
    INTEGER,PARAMETER         :: MAX_ITERATIONS = 10000
    real(dp),parameter        :: ADJUSTED_NUETRON_MASS = 1.04625e-8_dp
    
    !! Characters
    Character, Parameter          :: COMMENT_CHAR = '!'
    Character(len=2), Parameter   :: tab = "  "
    Character, Parameter          :: space = " "

    ! Solver Choices
    integer, parameter, public :: SOLVER_PCG = 1
    integer, parameter, public :: SOLVER_KSP_CG = 2
    integer, parameter, public :: SOLVER_KSP_GMRES = 3
    integer, parameter, public :: SOLVER_KSP_BCGS = 4

    ! Preconditioner Choices
    integer, parameter, public :: PRECON_NONE = 0
    integer, parameter, public :: PRECON_JACOBI = 1
    integer, parameter, public :: PRECON_ILU = 2
    integer, parameter, public :: PRECON_CHOLESKY = 3
    integer, parameter, public :: PRECON_GAMG = 4

    ! Boundary Params
    integer, parameter  :: BC_VACUUM     = 1  
    integer, parameter  :: BC_REFLECTIVE = 2  
    integer, parameter  :: BC_DIRICHLET  = 3  
    integer, parameter  :: BC_ALBEDO     = 4  
    real(dp), parameter :: PENALTY = 1.0e10_dp
  
contains

    ! Helper subroutines for NaN checking
    subroutine check_nan_scalar(value, var_name, location)
        real(dp), intent(in) :: value
        character(len=*), intent(in) :: var_name, location
        if (isnan(value)) then
            write(*,*) "FATAL: NaN detected in scalar ", trim(var_name), " at ", trim(location)
            stop "NaN detected."
        end if
    end subroutine check_nan_scalar

    subroutine check_nan_array(arr, var_name, location)
        real(dp), intent(in) :: arr(:)
        character(len=*), intent(in) :: var_name, location
        if (any(isnan(arr))) then
            write(*,*) "FATAL: NaN detected in array ", trim(var_name), " at ", trim(location)
        end if
    end subroutine check_nan_array

    subroutine check_nan_matrix(mat, var_name, location)
        real(dp), intent(in) :: mat(:,:)
        character(len=*), intent(in) :: var_name, location
        if (any(isnan(mat))) then
            write(*,*) "FATAL: NaN detected in matrix ", trim(var_name), " at ", trim(location)
            stop "NaN detected."
        end if
    end subroutine check_nan_matrix

    pure function int_to_str(i) result(res)
        integer, intent(in) :: i
        character(len=12) :: res
        write(res, '(I0)') i
        res = adjustl(res)
    end function int_to_str


end module