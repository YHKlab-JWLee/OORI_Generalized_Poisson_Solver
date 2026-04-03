MODULE m_poison_mg
!**************************************************************************************
! Yong Hoon Kim, Department of Physics, University of Illinois
! Juho Lee, Graduate School of EEWS, KAIST
!
!-- DESCRIPTIONS --
! * Type definition for Poisson solver
!    MG FFT solver
!
!-- REVISION HISTORY --
! 990810 Written by YHK.
! 171024 Modified by JL.
!**************************************************************************************
  PRIVATE 
  PUBLIC :: generalized_poison_solver, cg_poisson_solver

CONTAINS

!======================================================================================
  SUBROUTINE generalized_poison_solver(CELL,N1,N2,N3,nfd,rho,esp,vhar,ehar,stress,bc_type)
!======================================================================================
    USE precision, ONLY: dp,grid_p
    USE m_boundary_condition, ONLY: BC_DIRICHLET_ZERO
    IMPLICIT NONE
    LOGICAL                  :: hbc
    INTEGER                  :: N1, N2, N3, ierr, nfd
    INTEGER, optional        :: bc_type
    INTEGER                  :: bc_type_eff
    integer                  :: IX, JX, mpc
    REAL(grid_p) :: rho(N1*N2*N3), vhar(N1*N2*N3), esp(N1*N2*N3)
    REAL(dp)                 :: CELL(3,3)
    REAL(dp), optional       :: stress(3,3)
    REAL(dp)     :: ehar
    REAL(dp)                 :: dx, dy, dz, dvol
!**************************************************************************************
! Yong-Hoon Kim, Department of Physics, University of Illinois
!
!-- DESCRIPTION --
! Calls Poisson equation solver(s) to generate Hartree potential and 
! calculates the Hartree energy if 'ehar' is present.
! It uses 1/ FFT, and 2/ minimize the error with preconditioned CG.
!
!-- REVISION HISTORY --
! 9710   Rewritten based on Inho's code.
! 980623 Written by YHK based on the previous code.
! 980717 `rho' and `vhar' is now on the subgrid on input.  
!        Copied to total grid and calculations will be performed.
!        (subgrid copy is outside of this routine now)
! 990411 Charge normalization checked with 'qtot'.
! 990615 Type 'grid_obj' used.
!**************************************************************************************
!     Initialize stress contribution
      do IX = 1,3
        do JX = 1,3
          STRESS(JX,IX) = 0.0_dp
        enddo
      enddo
! 1/ Initial `vhar2' by 3-point lower order solution of Poisson equation.
!    CALL poison_mg(hbc,CELL,N1,N2,N3,rho,vhar)

! 2/ Polish up by higher order, more accurate solution using CG.
    bc_type_eff = BC_DIRICHLET_ZERO
    IF (present(bc_type)) bc_type_eff = bc_type
    CALL cg_poisson_solver(CELL,N1,N2,N3,nfd,rho,esp,vhar,bc_type_eff)
  
!Grid Spacings
    dx = CELL(1,1)/real(N1-1, dp)
    dy = CELL(2,2)/real(N2-1, dp)
    dz = CELL(3,3)/real(N3-1, dp)

!Grid cell volume.
    dvol = dx*dy*dz   

    ehar = DOT_PRODUCT(vhar,rho)*dvol/2.0_dp 
  END SUBROUTINE generalized_poison_solver

!!======================================================================================
!  SUBROUTINE poison_mg(hbc,CELL,N1,N2,N3,rho,vhar)
!======================================================================================
!    USE precision, ONLY: dp, grid_p
!    IMPLICIT NONE
!    INTEGER                  :: N1, N2, N3, N1h, N2h, N3h, ierr, mpc
!    REAL(grid_p)   :: rho(1:N1,1:N2,1:N3), vhar(1:N1,1:N2,1:N3) ! Hartree potential
!    REAL(dp)                 :: CELL(3,3)
!**************************************************************************************
!Yong-Hoon Kim, Department of Physics, University of Illinois, Oct. 1997
!
!-- DESCRIPTION --
! Dummy interface for the use of subroutine `hw3crt'.
! `hw3crt' is the library routine for the solution of Helmholtz equation
! using FACR (combination of Fourier analysis and cyclic reduction) algorithm.
!
!-- NOTE --
! o Boundary values of 'vhar' should have been assigned
! o Presently only Dirichlet boundary condition 
!
!-- REVISION HISTORY --
! 970428 Original version written by IHL.
! 971023 Rewritten by YHK.
! 990615 Modified to use type 'grid'.
!**************************************************************************************
!    LOGICAL :: hbc
!    INTEGER  :: i,j,k,inx,iswitch,ipoi_s,ipoi_f, &
!                lbdcnd,mbdcnd,nbdcnd,ldimf,mdimf,&
!                ix, iy, iz
!    REAL(dp) :: pi,pi8,elmbda,pertrb,t1,t2
!!    REAL(dp), ALLOCATABLE, DIMENSION(:) :: w
!    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: bdxi,bdxf,bdyi,bdyf,bdzi,bdzf
!    external   hw3crt
!--------------------------------------------------------------------------------------
!    lbdcnd=0
!    mbdcnd=0
!    nbdcnd=0 
! 
!
!    elmbda=0.0_dp
!    ldimf=N1
!    mdimf=N2
!    inx = 30 + N1 + N2 + 5*N3 + MAX(N1,N2,N3) &
!         + 7*(INT((N1+1)/2) + INT((N2+1)/2))
!    ALLOCATE(bdxi(1:N2,1:N3),bdxf(1:N2,1:N3), &
!             bdyi(1:N1,1:N3),bdyf(1:N1,1:N3), &
!             bdzi(1:N1,1:N2),bdzf(1:N1,1:N2))!, stat=ierr)
!    ALLOCATE(w(inx))     !, stat=ierr)
!    pi = 4.0_dp*atan(1.0_dp)
!    pi8 = pi*8.0_dp

! 20171024 JLee START
! Here, we need a script to read the boundary potential at z-axis from previous
!    IF(hbc) THEN
!       N1h = nint(N1/2.0_dp)
!       N2h = nint(N2/2.0_dp)
!       N3h = nint(N3/2.0_dp)
!       rho(1:N1h,1:N2h,1:N3h) = 0.0_dp 
!    ENDIF
!    vhar(1:N1,1:N2,1:N3) = -pi8*rho(1:N1,1:N2,1:N3)
! 20171024 END
!    CALL hw3crt (0.0_dp, CELL(1,1), N1-1, lbdcnd, bdxi, bdxf, &
!                 0.0_dp, CELL(2,2), N2-1, mbdcnd, bdyi, bdyf, &
!                 0.0_dp, CELL(3,3), N3-1, nbdcnd, bdzi, bdzf, &
!                 elmbda, ldimf, mdimf, vhar, pertrb, ierr, w)
!    IF(ierr /= 0) WRITE(6,'(a,i4)') 'hw3crt_poisson_solver: %%% ERROR! -- ierr =',ierr
!
!    DEALLOCATE(w)
!    DEALLOCATE(bdxi,bdxf,bdyi,bdyf,bdzi,bdzf)
!
!    CALL CPU_TIME(t2)
!    IF(lpr) THEN
!       WRITE(6,'(1x,a,1x,f7.2)') '@ hw3crt_poisson_solver: COMPUTING TIME =>',(t2-t1)
!    ENDIF
!
!  END SUBROUTINE poison_mg
!
!======================================================================================
  SUBROUTINE cg_poisson_solver(CELL,N1,N2,N3, nfd,rho,esp,vhar,bc_type)
!======================================================================================
    USE precision,   ONLY: dp
!    USE m_laplacian, ONLY: laplacian_fd
    USE m_grid_obj,   only: gradient_fd
    IMPLICIT NONE
    logical   :: mg
    INTEGER                  :: N1, N2, N3
    REAL(dp)  :: rho(N1*N2*N3), vhar(N1*N2*N3), esp(N1*N2*N3)
    real(dp)  :: CELL(3,3)
    integer,  intent(in) :: nfd
    integer,  intent(in) :: bc_type
!    INTEGER,  INTENT(OUT) :: iter  ! Final iteration steps taken.
!    REAL(dp), INTENT(OUT) :: err   ! Final error.
!**************************************************************************************
!-- DESCRIPTIONS --
! Solves Poisson eq.: 
!    \nabla^2 Vhar = - 4 \pi \rho
! using preconditioned conjugate gradient method.
! So,
!    A = \nabla^2
!    x = V_{har}
!    b = -4 \pi \rho    
!
!-- NOTE --
! o Designed to be used together with subroutine 'hw3crt_poisson_solver', such that 'vhar'
!   should be initially good solution.
! o Preconditioning operator:
!
!-- CONTAINED SUBROUTINES --
! o <prcnd> 
! o <filbdz>
!
!--REFERENCES--
!   [0] "Numerical Recipes", Ch. 2.7.
!
!-- REVISION HISTORY --
! 981002 Written by YHK.
! 990615 Modified to use 'grid_obj' type.
! 990712 Modified to use 'poisson_pccg_prm' type.
!**************************************************************************************
    INTEGER  :: i,ng,ierr, iter
    REAL(dp) :: ak,akden,bk,bkden,bknum,bnrm
    real(dp)  :: pi, pi8, err
    real(dp)  :: dHtol=1.e-7_dp
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: p,r,z
    REAL(dp), allocatable, dimension(:) :: del1f, del2f, del3f
    REAL(dp), allocatable, dimension(:) :: ddel1fx,ddel1fy,ddel1fz
    REAL(dp), allocatable, dimension(:) :: ddel2fx,ddel2fy,ddel2fz
    REAL(dp), allocatable, dimension(:) :: ddel3fx,ddel3fy,ddel3fz
!--------------------------------------------------------------------------------------

    ng = N1*N2*N3
    
    ALLOCATE(p(ng),r(ng),z(ng),stat=ierr)
    IF(ierr/=0) STOP 'cg_poisson_solver: ALLOCATAION FAIL! -- p,r,z'   

    allocate(del1f(ng), del2f(ng), del3f(ng))
    allocate(ddel1fx(ng), ddel1fy(ng), ddel1fz(ng), &
             ddel2fx(ng), ddel2fy(ng), ddel2fz(ng), &
             ddel3fx(ng), ddel3fy(ng), ddel3fz(ng))

    pi = 4.0_dp*atan(1.0_dp)
    pi8 = pi*8.0_dp

    ! Initial CG step
    iter=1
     
    ! Initial residual vector 'r_0'
    call gradient_fd(.true.,CELL, N1, N2, N3, nfd, vhar,del1f, del2f, del3f) 

    do i = 1, ng
      del1f(i)=esp(i)*del1f(i)
      del2f(i)=esp(i)*del2f(i)
      del3f(i)=esp(i)*del3f(i)
    enddo
 
    call gradient_fd(.true., CELL, N1, N2, N3, nfd, del1f, ddel1fx, ddel1fy, ddel1fz) 
    call gradient_fd(.true., CELL, N1, N2, N3, nfd, del2f, ddel2fx, ddel2fy, ddel2fz) 
    call gradient_fd(.true., CELL, N1, N2, N3, nfd, del3f, ddel3fx, ddel3fy, ddel3fz) 

    do i = 1, ng
      r(i) = ddel1fx(i)+ddel2fy(i)+ddel3fz(i)
    enddo

    r = -pi8*(rho) - r

    CALL apply_bc_vector(N1,N2,N3,bc_type,r)

    ! Norm of vector '-4 \pi * rho'
    call prcnd(N1,N2,N3,r,z)
    CALL apply_bc_vector(N1,N2,N3,bc_type,z)
    bnrm = pi8* MAXVAL(ABS(rho))

    ! Main loop 
    DO

       bknum = DOT_PRODUCT(r,z)

       ! Calculate conjugate gradient direction.
       ! If first iteration, just steepest descent,
       ! else, calculate '\beta_{k}' and calculate conjugate gradient direction 'p_k'.
       IF (iter == 1) THEN
          p = z
       ELSE
          bk = bknum/bkden
          p = z + bk*p
       END IF
      !  CALL filbdzero(N1,N2,N3,p)
       bkden=bknum

       call gradient_fd(.true., CELL, N1, N2, N3, nfd, p,del1f, del2f, del3f) 

       do i = 1, ng
        del1f(i)=esp(i)*del1f(i)
        del2f(i)=esp(i)*del2f(i)
        del3f(i)=esp(i)*del3f(i)
       enddo
 
       call gradient_fd(.true., CELL, N1, N2, N3, nfd, del1f, ddel1fx, ddel1fy, ddel1fz) 
       call gradient_fd(.true., CELL, N1, N2, N3, nfd, del2f, ddel2fx, ddel2fy, ddel2fz) 
       call gradient_fd(.true., CELL, N1, N2, N3, nfd, del3f, ddel3fx, ddel3fy, ddel3fz) 

       do i = 1, ng
         z(i) = ddel1fx(i)+ddel2fy(i)+ddel3fz(i)
       enddo
      !  CALL filbdzero(N1,N2,N3,z)
       akden = DOT_PRODUCT(p,z)
       ak = bknum/akden

       ! Calculate new iterate 'x_{k+1}'
       vhar = vhar + ak*p

       ! Calculate new residual 'r_{k+1} = r_{k} - \alpha_{k} '
       r = r - ak*z

      call prcnd(N1,N2,N3,r,z)
      !  CALL filbdzero(N1,N2,N3,z)
       ! Check stopping criterion
       err = MAXVAL(ABS(r))/bnrm

       IF(err <= dHtol) EXIT

       ! Check iteration number 
       iter=iter+1
!       WRITE (6,*) 'cg_poisson_solver: step=',iter,' err=',err
       
       IF (iter > 2000) THEN
          iter=iter-1
          EXIT
       ENDIF
    END DO
    
    WRITE (6,*) 'cg_poisson_solver: iter=',iter

    DEALLOCATE(p,r,z,del1f, del2f, del3f, ddel1fx, ddel1fy, ddel1fz, &
               ddel2fx, ddel2fy, ddel2fz, ddel3fx, ddel3fy, ddel3fz)

  CONTAINS

!======================================================================================
  SUBROUTINE prcnd(N1,N2,N3,a,c)
!======================================================================================
    USE precision, ONLY: dp
    IMPLICIT NONE
    INTEGER                  :: N1, N2, N3
    REAL(dp),   INTENT(IN)  :: a(1:N1,1:N2,1:N3)  ! Input vector.
!    REAL(dp),   INTENT(IN)  :: esp(1:N1,1:N2,1:N3)  ! Input vector.
    REAL(dp),   INTENT(OUT) :: c(1:N1,1:N2,1:N3)  ! Preconditioned vector.
! JLee 200529
    REAL(dp), allocatable, dimension(:,:,:) :: aux
!    REAL(dp), allocatable, dimension(:,:,:) :: b, bux
!**************************************************************************************
!-- DESCRIPTIONS --
! Precondition vector 'a' to generate preconditoned vector 'b'
!
!--REFERENCES--
! [1] Hoshi et al. Phys. Rev. B {bf 52}, R5459 (1995).
! [2] Saas et al. BIT {\bf 36:3}, 563 (1996).
!**************************************************************************************

    INTEGER :: i,j,k,ipc
    allocate(aux(0:N1+1, 0:N2+1, 0:N3+1))!, bux(0:N1+1,0:N2+1,0:N3+1))
    aux(1:N1, 1:N2, 1:N3) = a(1:N1, 1:N2, 1:N3)
    aux(0, 1:N2, 1:N3) = a(N1, 1:N2, 1:N3)
    aux(N1+1, 1:N2, 1:N3) = a(1, 1:N2, 1:N3)
    aux(1:N1, 0, 1:N3) = a(1:N1, N2, 1:N3)
    aux(1:N1, N2+1, 1:N3) = a(1:N1, 1, 1:N3)
    aux(1:N1, 1:N2, 0) = a(1:N1, 1:N2, N3)
    aux(1:N1, 1:N2, N3+1) = a(1:N1, 1:N2, 1)

    DO i=1,N1
    DO j=1,N2
    DO k=1,N3
          c(i,j,k)=(aux(i-1,j,k) + aux(i+1,j,k) + &
                    aux(i,j-1,k) + aux(i,j+1,k) + &
                    aux(i,j,k-1) + aux(i,j,k+1))/12.0_dp + &
                    aux(i,j,k)/2.0_dp
    ENDDO
    ENDDO
    ENDDO

    deallocate(aux)!,bux,b)
  END SUBROUTINE prcnd

!======================================================================================
  SUBROUTINE apply_bc_vector(N1,N2,N3,bc_type,f)
!======================================================================================
    USE precision, ONLY: dp
    USE m_boundary_condition, ONLY: BC_PERIODIC, apply_boundary_condition_3d
    INTEGER, INTENT(IN)    :: N1,N2,N3
    INTEGER, INTENT(IN)    :: bc_type
    REAL(dp),   INTENT(INOUT) :: f(1:N1,1:N2,1:N3)
    REAL(dp), allocatable :: ftmp(:,:,:)
!**************************************************************************************
!-- DESCRIPTIONS --
! Set boundary values of array 'f' to zero.
!**************************************************************************************
    IF (bc_type == BC_PERIODIC) RETURN
    allocate(ftmp(N1,N2,N3))
    ftmp = f
    call apply_boundary_condition_3d(bc_type, ftmp)
    f = ftmp
    deallocate(ftmp)
  END SUBROUTINE apply_bc_vector

  END SUBROUTINE cg_poisson_solver


END MODULE m_poison_mg
