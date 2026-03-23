module m_constants

    Implicit None
  
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
    real(Kind=dp),  dimension(1), parameter :: VLARGE_NUMBER = 1.0E+9_dp
    INTEGER,PARAMETER         :: MAX_ITERATIONS = 10000
    real(dp),parameter        :: ADJUSTED_NUETRON_MASS = 1.04625e-8_dp
    
    !! Characters
    Character, Parameter          :: COMMENT_CHAR = '!'
    Character(len=2), Parameter   :: tab = "  "
    Character, Parameter          :: space = " "


    integer, parameter :: GEOM_SLAB=0, GEOM_CYL=1, GEOM_SPH=2
    integer, parameter :: BC_NONE=0, BC_DIRICHLET=1, BC_NEUMANN=2, BC_VACUUM=3, BC_ALBEDO=4, BC_PERIODIC=5
  
end module