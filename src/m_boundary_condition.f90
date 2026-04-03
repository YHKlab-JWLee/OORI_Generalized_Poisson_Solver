MODULE m_boundary_condition
  USE precision, ONLY: dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BC_DIRICHLET_ZERO, BC_PERIODIC
  PUBLIC :: parse_boundary_condition, apply_boundary_condition_3d
  PUBLIC :: build_periodic_halo

  INTEGER, PARAMETER :: BC_DIRICHLET_ZERO = 1
  INTEGER, PARAMETER :: BC_PERIODIC       = 2

CONTAINS

!======================================================================================
  SUBROUTINE parse_boundary_condition(bc_name, bc_type, is_known)
!======================================================================================
    CHARACTER(len=*), INTENT(IN)  :: bc_name
    INTEGER,          INTENT(OUT) :: bc_type
    LOGICAL, OPTIONAL, INTENT(OUT) :: is_known
    CHARACTER(len=32)             :: key

    key = adjustl(trim(bc_name))
    CALL lower_string(key)

    SELECT CASE (trim(key))
    CASE ('periodic','pbc')
      bc_type = BC_PERIODIC
      IF (present(is_known)) is_known = .true.
    CASE ('dirichlet_zero','dirichlet','zero','open')
      bc_type = BC_DIRICHLET_ZERO
      IF (present(is_known)) is_known = .true.
    CASE DEFAULT
      bc_type = BC_DIRICHLET_ZERO
      IF (present(is_known)) is_known = .false.
    END SELECT
  END SUBROUTINE parse_boundary_condition

!======================================================================================
  SUBROUTINE apply_boundary_condition_3d(bc_type, f)
!======================================================================================
    INTEGER,  INTENT(IN)    :: bc_type
    REAL(dp), INTENT(INOUT) :: f(:,:,:)

    IF (bc_type == BC_DIRICHLET_ZERO) THEN
      f(1,:,:) = 0.0_dp; f(size(f,1),:,:) = 0.0_dp
      f(:,1,:) = 0.0_dp; f(:,size(f,2),:) = 0.0_dp
      f(:,:,1) = 0.0_dp; f(:,:,size(f,3)) = 0.0_dp
    END IF
  END SUBROUTINE apply_boundary_condition_3d

!======================================================================================
  SUBROUTINE build_periodic_halo(f, nfd, halo)
!======================================================================================
    INTEGER,  INTENT(IN)  :: nfd
    REAL(dp), INTENT(IN)  :: f(:,:,:)
    REAL(dp), INTENT(OUT) :: halo(-nfd+1:,-nfd+1:,-nfd+1:)
    INTEGER               :: i, j, k
    INTEGER               :: ii, jj, kk
    INTEGER               :: n1, n2, n3

    n1 = size(f,1); n2 = size(f,2); n3 = size(f,3)

    DO k = lbound(halo,3), ubound(halo,3)
      kk = modulo(k-1, n3) + 1
      DO j = lbound(halo,2), ubound(halo,2)
        jj = modulo(j-1, n2) + 1
        DO i = lbound(halo,1), ubound(halo,1)
          ii = modulo(i-1, n1) + 1
          halo(i,j,k) = f(ii,jj,kk)
        END DO
      END DO
    END DO
  END SUBROUTINE build_periodic_halo

!======================================================================================
  SUBROUTINE lower_string(str)
!======================================================================================
    CHARACTER(len=*), INTENT(INOUT) :: str
    INTEGER                         :: i, ich
    DO i = 1, len_trim(str)
      ich = iachar(str(i:i))
      IF (ich >= iachar('A') .AND. ich <= iachar('Z')) str(i:i) = achar(ich + 32)
    END DO
  END SUBROUTINE lower_string

END MODULE m_boundary_condition
