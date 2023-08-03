MODULE machine_constants
  ! ... machine dependent parameters
  ! ... $Revision: 1.4 $//$Date: 2007/07/25 23:32:29 $
  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: kdp = SELECTED_REAL_KIND(14,60)
  ! ... BGREAL: A large real number representable in single precision
  ! ... BGINT:  A large integer number representable in 4 bytes
  INTEGER, PARAMETER :: BGINT=9999
  REAL(KIND=kdp), PARAMETER :: BGREAL=9.999E35_kdp
  REAL(KIND=kdp), PARAMETER :: epsmac=EPSILON(1._kdp)
END MODULE machine_constants
