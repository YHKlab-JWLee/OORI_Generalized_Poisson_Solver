MODULE m_poisson_bc_ops
  USE precision, ONLY: dp
  USE m_boundary_condition, ONLY: BC_PERIODIC, build_periodic_halo
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: is_periodic_bc, precondition_fd_3d

CONTAINS

!======================================================================================
  LOGICAL FUNCTION is_periodic_bc(bc_type)
!======================================================================================
    INTEGER, INTENT(IN) :: bc_type
    is_periodic_bc = (bc_type == BC_PERIODIC)
  END FUNCTION is_periodic_bc

!======================================================================================
  SUBROUTINE precondition_fd_3d(N1,N2,N3,bc_type,a,c)
!======================================================================================
    INTEGER,  INTENT(IN)  :: N1, N2, N3
    INTEGER,  INTENT(IN)  :: bc_type
    REAL(dp), INTENT(IN)  :: a(1:N1,1:N2,1:N3)
    REAL(dp), INTENT(OUT) :: c(1:N1,1:N2,1:N3)
    REAL(dp), ALLOCATABLE :: aux(:,:,:)
    INTEGER               :: i,j,k

    IF (is_periodic_bc(bc_type)) THEN
      ALLOCATE(aux(0:N1+1,0:N2+1,0:N3+1))
      CALL build_periodic_halo(a, 1, aux)
    ELSE
      ALLOCATE(aux(0:N1+1,0:N2+1,0:N3+1))
      aux = 0.0_dp
      aux(1:N1,1:N2,1:N3) = a(1:N1,1:N2,1:N3)
    END IF

    DO i=1,N1
      DO j=1,N2
        DO k=1,N3
          c(i,j,k)=(aux(i-1,j,k) + aux(i+1,j,k) + &
                    aux(i,j-1,k) + aux(i,j+1,k) + &
                    aux(i,j,k-1) + aux(i,j,k+1))/12.0_dp + &
                    aux(i,j,k)/2.0_dp
        END DO
      END DO
    END DO

    DEALLOCATE(aux)
  END SUBROUTINE precondition_fd_3d

END MODULE m_poisson_bc_ops
