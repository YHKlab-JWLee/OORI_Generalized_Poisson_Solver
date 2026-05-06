MODULE m_grid_obj

  implicit none
  public :: laplacian_fd, gradient_fd, weights 

CONTAINS
!======================================================================================
  SUBROUTINE gradient_fd(pbc,CELL,N1,N2,N3,nfd,f,del1f,del2f,del3f)
!======================================================================================
    USE precision, ONLY: dp 
    IMPLICIT NONE
! JLee 191114
    integer :: N1, N2, N3
    real(dp) :: CELL(3,3)    
    logical  :: pbc
    INTEGER, INTENT(IN) :: nfd ! Higher-order finite difference order
    REAL(dp), INTENT(IN)  :: f(N1*N2*N3)
    REAL(dp), INTENT(OUT) :: del1f(N1*N2*N3), &
                             del2f(N1*N2*N3), &
                             del3f(N1*N2*N3)
!**************************************************************************************
! Yong Hoon Kim, Department of Physics, University of Illinois
!
!--DESCRIPTION--
! Calculates gradient of function f(x,y,z) using higher-order finite difference method:
!   del1f := grad_1 f(x,y,z) = d/dx f(x,y,z) 
!   del2f := grad_2 f(x,y,z) = d/dy f(x,y,z) 
!   del3f := grad_3 f(x,y,z) = d/dz f(x,y,z) 
!
!--NOTE--
! o 1D version: calls <grid_obj3d_gradient_fd>
! o This subroutine assumes that the grids are uniform along x,y, & z direction
!
!-- REFERENCES --
! o See <laplacian>
! o NOTE FDP-XC980420
!
!-- REVISION HISTORY --
! 990830 Written by YHK.
! 191114 Converted to SIESTA ver. by JLee.
!**************************************************************************************
    CALL grid_obj3d_gradient_fd(pbc,CELL,N1,N2,N3,nfd,f,del1f,del2f,del3f)
  END SUBROUTINE gradient_fd

!======================================================================================
  SUBROUTINE grid_obj3d_gradient_fd(pbc,CELL,N1,N2,N3,nfd,f,del1f,del2f,del3f)
!======================================================================================
    USE precision, ONLY: dp !,zero
    USE m_boundary_condition, ONLY: build_periodic_halo
!    USE m_nr, ONLY: weights
    IMPLICIT NONE
    integer :: N1, N2, N3
    real(dp) :: CELL(3,3), zero    
    INTEGER, INTENT(IN) :: nfd ! Higher-order finite difference order
    REAL(dp), INTENT(IN)  :: f(1:N1,1:N2,1:N3)
    REAL(dp), INTENT(OUT) :: del1f(1:N1,1:N2,1:N3), &
                             del2f(1:N1,1:N2,1:N3), &
                             del3f(1:N1,1:N2,1:N3)
!**************************************************************************************
! Yong Hoon Kim, Department of Physics, University of Illinois
!
!--DESCRIPTION--
! SUBROUTINE gradient CALCULATES GRADIENT OF FUNCTION f(x,y,z):
!   del1f := grad_1 f(x,y,z) = d/dx f(x,y,z) 
!   del2f := grad_2 f(x,y,z) = d/dy f(x,y,z) 
!   del3f := grad_3 f(x,y,z) = d/dz f(x,y,z) 
!
!--NOTE--
! o This subroutine assumes that the grids are uniform along x,y, & z direction
!
!-- REFERENCES --
! o See <laplacian>
! o NOTE FDP-XC980420
!
!-- REVISION HISTORY --
! 980420 Written & tested by YHK.
! 990830 Revised with new type definitions.
!**************************************************************************************
    INTEGER, PARAMETER :: mdrv=1 ! HIGHEST ORDER OF DERIVATIVE TO BE APPROXIMATED 
                                 ! (FOR SUBROUTINE weights)

    INTEGER :: nsten,ierr,i,il,ir,j,jl,jr,k,kl,kr,l
    REAL(dp) :: xi,yj,zk

    ! Allocatable arrays for subroutine <weights>
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: xsten, xcrd, ycrd, zcrd
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: cof,xcof,ycof,zcof
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: f_temp ! array for sliced f values
    real(dp) :: dx, dy, dz
    logical  :: pbc
!--------------------------------------------------------------------------------------

    nsten = nfd*2
    zero = 0.0_dp

    if (pbc) then
        allocate(f_temp(-nfd+1:N1+nfd,-nfd+1:N2+nfd,-nfd+1:N3+nfd))
        call build_periodic_halo(f, nfd, f_temp)
    endif

    ALLOCATE(xsten(0:nsten),zcrd(1:N3), &
             xcrd(1:N1),ycrd(1:N2), &
             cof(0:nsten,0:nsten,0:mdrv), &
             xcof(0:nsten,0:nsten,0:mdrv), &
             ycof(0:nsten,0:nsten,0:mdrv), &
             zcof(0:nsten,0:nsten,0:mdrv), &
             STAT=ierr)

!Grid Spacings
! SIESTA mesh files store N points over N mesh divisions:
! x_i = i * L/N, i = 0, ..., N-1.  Do not reinterpret them as endpoint grids.
    dx = CELL(1,1)/real(N1, dp)
    dy = CELL(2,2)/real(N2, dp)
    dz = CELL(3,3)/real(N3, dp)

    xcrd = arth_d(0.0_dp, dx, (N1))
    ycrd = arth_d(0.0_dp, dy, (N2))
    zcrd = arth_d(0.0_dp, dz, (N3))
 

! (1) ALONG THE X-DIR.  del1f := d/dx f(x,y,z) 

    ! Since we are using uniform grids, we need to calculate x,y, and z 
    ! coefficients only once. here we choose the calculation point 'xcrd(nfd)'
    ! but it can be any point between npnt and 'nx(y,z)-nfd'
    xsten(0:nsten) = xcrd(1:nsten+1)
    CALL weights(xcrd(nfd+1),xsten,nsten,mdrv,xcof) 

    x_dir: DO i = 1,N1
       xi = xcrd(i)
       il = i-nfd
       ir = i+nfd
       if(pbc) then
             do k = 1,N3
             do j = 1,N2
                del1f(i,j,k) = DOT_PRODUCT( xcof(0:nsten,nsten,mdrv),f_temp(il:ir,j,k) )
             enddo
             enddo
       else
       ! CASE 1: WHEN ALL THE STNECILS ARE IN THE GRID RANGE
          IF(il >= 1 .AND. ir <= N1) THEN
             DO k = 1,N3
             DO j = 1,N2
                del1f(i,j,k) = DOT_PRODUCT( xcof(0:nsten,nsten,mdrv),f(il:il+nsten,j,k) )
             ENDDO
             ENDDO

       ! CASE 2: WHEN THE LEFT STENCILS ARE OUTSIDE OF GRID RANGE
          ELSEIF(il < 1) THEN
             xsten(0:nsten) = xcrd(1:nsten+1)
             CALL weights(xi,xsten,nsten,mdrv,cof) 
             DO k = 1,N3
             DO j = 1,N2
                del1f(i,j,k) = DOT_PRODUCT( cof(0:nsten,nsten,mdrv),f(1:nsten+1,j,k) )
             ENDDO
             ENDDO

       ! CASE 3: WHEN RIGHT STNECILS ARE OUTSIDE OF THE GRID RANGE
          ELSE
             xsten(0:nsten) = xcrd(N1-nsten:N1)
             CALL weights(xi,xsten,nsten,mdrv,cof) 
             DO k = 1,N3
             DO j = 1,N2
                del1f(i,j,k) = DOT_PRODUCT( cof(0:nsten,nsten,mdrv),f(N1-nsten:N1,j,k) )
             ENDDO
             ENDDO
          ENDIF
       endif  ! pbc
    ENDDO x_dir

    ! (2) ALONG THE Y-DIR. del2f := d/dy  f(x,y,z) 
    xsten(0:nsten) = ycrd(1:nsten+1)
    CALL weights(ycrd(nfd+1),xsten,nsten,mdrv,ycof) 

    y_dir: DO j = 1,N2
       yj=ycrd(j)
       jl = j - nfd
       jr = j + nfd
       if(pbc) then
             DO i=1,N1
             DO k=1,N3
                del2f(i,j,k) = SUM( ycof(0:nsten,nsten,mdrv)*f_temp(i,jl:jr,k) )
             ENDDO
             ENDDO
       else
          IF(jl >= 1 .AND. jr <= N2)THEN
             DO i=1,N1
             DO k=1,N3
                del2f(i,j,k) = SUM( ycof(0:nsten,nsten,mdrv)*f(i,jl:jl+nsten,k) )
             ENDDO
             ENDDO
          ELSEIF(jl < 1) THEN
             xsten(0:nsten) = ycrd(1:nsten+1)
             CALL weights(yj,xsten,nsten,mdrv,cof) 
             DO k=1,N3
             DO i=1,N1
                del2f(i,j,k) = SUM( cof(0:nsten,nsten,mdrv)*f(i,1:nsten+1,k) )
             ENDDO
             ENDDO
          ELSE
             xsten(0:nsten) = ycrd(N2-nsten:N2)
             CALL weights(yj,xsten,nsten,mdrv,cof) 
             DO k=1,N3
             DO i=1,N1
                del2f(i,j,k) = SUM( cof(0:nsten,nsten,mdrv)*f(i,N2-nsten:N2,k) )
             ENDDO
             ENDDO
          ENDIF
       endif !pbc
    ENDDO y_dir

    ! (3) ALONG THE Z-DIR. : del3f := d/dz f(x,y,z)
    xsten(0:nsten) = zcrd(1:nsten+1)
    CALL weights(zcrd(nfd+1),xsten,nsten,mdrv,zcof)
    z_dir: DO k=1,N3
       zk = zcrd(k)
       kl = k - nfd
       kr = k + nfd
       if(pbc) then
             DO j=1,N2
             DO i=1,N1
                del3f(i,j,k) = SUM( zcof(0:nsten,nsten,mdrv)*f_temp(i,j,kl:kr) )
             ENDDO
             ENDDO
       else
          IF(kl >= 1 .AND. kr <= N3) THEN
             DO j=1,N2
             DO i=1,N1
                del3f(i,j,k) = SUM( zcof(0:nsten,nsten,mdrv)*f(i,j,kl:kl+nsten) )
             ENDDO
             ENDDO
          ELSEIF(kl < 1) THEN
             xsten(0:nsten) = zcrd(1:nsten+1)
             CALL weights(zk,xsten,nsten,mdrv,cof) 
             DO j=1,N2
             DO i=1,N1
                del3f(i,j,k) = SUM( cof(0:nsten,nsten,mdrv)*f(i,j,1:nsten+1) )
             ENDDO
             ENDDO
          ELSE
             xsten(0:nsten) = zcrd(N3-nsten:N3)
             CALL weights(zk,xsten,nsten,mdrv,cof) 
             DO j=1,N2
             DO i=1,N1
                del3f(i,j,k) = SUM( cof(0:nsten,nsten,mdrv)*f(i,j,N3-nsten:N3) )
             ENDDO
             ENDDO
          ENDIF
       endif !pbc
    ENDDO z_dir

  CONTAINS
!======================================================================================
  FUNCTION arth_d(first,increment,n)
!======================================================================================
    INTEGER, PARAMETER :: i4b = selected_int_KIND(9)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(i4b), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER(i4b) :: k,k2
    REAL(DP) :: temp
    integer(i4b), parameter  :: npar=16, npar2=8
!======================================================================================
    if (n > 0) arth_d(1)=first
    if (n <= npar) then
       do k=2,n
          arth_d(k)=arth_d(k-1)+increment
       end do
    else
       do k=2,npar2
          arth_d(k)=arth_d(k-1)+increment
       end do
       temp=increment*npar2
       k=npar2
       do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_d
 
  END SUBROUTINE grid_obj3d_gradient_fd


!======================================================================================
  SUBROUTINE laplacian_fd(mg,CELL,N1,N2,N3,nfd,f,lap_f)
!======================================================================================
    USE precision, ONLY: dp
    INTEGER, INTENT(in) :: nfd
    integer :: N1, N2, N3
    REAL(dp), INTENT(in) :: f(N1*N2*N3)
    REAL(dp), INTENT(out) :: lap_f(N1*N2*N3)
    logical :: mg
    real(dp) :: CELL(3,3)
!**************************************************************************************
!--DESCRIPTION--
! Calculates Laplacian of function f(x,y,z) on the grid 
!   lap_f(x,y,z) = (d^2/dx^2 + d^2/dy^2 + d^2d/dz^2) f(x,y,z) 
! using higher order finite difference formula.
!
!--NOTE--
! * SUBROUTINE USED : weights
! * This subroutine assumes that 
!   1/ 3 dim.
!   2/ the grids are uniform along x,y, & z direction.
!
! * 2 functions with 'lbdcnd' options:
!
!   1/ lbdcnd=1 => If the function has finite values at the boundary of box. (laplacian1)
!   2/ lbdcnd=2 => If the function is localized inside of the box, or f==0 outside 
!             of box boundary. (laplacian2)
!
!-- REFERENCES --
! NOTE 970927
! TEST 971007,FDP980204,11,20
!
!-- REVISION HISTORY --
! 970520 delsqp IN pcgpoi.f90 WRITTEN BY IHL
! 970927 Written by YHK
! 980306 Used Fortran 90 feature of subarray allocation and function 'sum'
! 980703 'sum'->'dot_product'
! 990203 Previous <laplacian> is moved to <laplacian1>. <laplacian2> written:
!        assumming ftn. is localized in the box.
! 990204 Now <laplacian> is interface to <laplacian1> and <laplacian2>.
! 990204 Algorithmic bug found! -- In <laplacian1>, z & y loops should be inside of 
!        the x loop, because of shifting of FD stenils in case 2 & 3. 
! 990819 Interface revised -- type 'grid', etc.
!**************************************************************************************

       CALL laplacian3(mg,CELL,N1,N2,N3,nfd,f,lap_f)

  END SUBROUTINE laplacian_fd


!======================================================================================
  SUBROUTINE laplacian3(mg,CELL,N1,N2,N3,nfd,f,lap_f)
!======================================================================================
    USE precision, ONLY: dp
    implicit none
    INTEGER, INTENT(in) :: nfd
    integer :: N1, N2, N3
    REAL(dp), intent(in) :: f(1:N1,1:N2,1:N3)
    real(dp)    :: lap_f(1:N1,1:N2,1:N3) 
    logical :: mg
    real(dp)  :: CELL(3,3)
    INTEGER, PARAMETER :: mdrv=2 ! HIGHEST ORDER OF DERIVATIVE TO BE APPROXIMATED 
                                 ! (FOR SUBROUTINE weights)

    INTEGER :: nsten,ierr,i,il,ir,j,jl,jr,k,kl,kr,l
    REAL(dp) :: xi,yj,zk, zero
 
    ! Allocatable arrays for subroutine <weights>
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: xsten, xcrd, ycrd, zcrd
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: cof, xcof,ycof,zcof
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: f_temp ! array for sliced f values
    real(dp) :: dx, dy, dz
!--------------------------------------------------------------------------------------

    nsten = nfd * 2
    zero = 0.0_dp
    ALLOCATE(f_temp(-nfd+1:N1+nfd,-nfd+1:N2+nfd,-nfd+1:N3+nfd), &
             STAT=ierr)
    IF(ierr/=0) STOP 'laplacian_node: Allocation fail! -- f_temp,ilap_f'
    f_temp = zero

! Later need to be modified for hybrid boundary condition 2017.10.28
    ! Make boundary values same for periodicity
    !IF(mg) THEN
    f_temp(1:N1,1:N2,1:N3)=f(1:N1,1:N2,1:N3)
    f_temp(-nfd+1:0,1:N2,1:N3)=f(N1-nfd+1:N1,1:N2,1:N3)
    f_temp(N1+1:N1+nfd,1:N2,1:N3)=f(1:nfd,1:N2,1:N3)
    f_temp(1:N1,-nfd+1:0,1:N3)=f(1:N1,N2-nfd+1:N2,1:N3)
    f_temp(1:N1,N2+1:N2+nfd,1:N3)=f(1:N1,1:nfd,1:N3)
    f_temp(1:N1,1:N2,-nfd+1:0)=f(1:N1,1:N2,N3-nfd+1:N3)
    f_temp(1:N1,1:N2,N3+1:N3+nfd)=f(1:N1,1:N2,1:nfd)
    !ENDIF

    ALLOCATE(xsten(0:nsten),zcrd(1:N3), &
             xcrd(1:N1),ycrd(1:N2), &
             xcof(0:nsten,0:nsten,0:mdrv), &
             ycof(0:nsten,0:nsten,0:mdrv), &
             zcof(0:nsten,0:nsten,0:mdrv), &
             STAT=ierr)
    !IF(ierr/=0) STOP 'laplacian: Allocation fail! -- xsten,cof,xcof,ycof,zcof'

    lap_f = zero

    !Grid Spacings
    ! Keep the Laplacian helper on the same SIESTA-compatible mesh convention.
    dx = CELL(1,1)/real(N1, dp)
    dy = CELL(2,2)/real(N2, dp)
    dz = CELL(3,3)/real(N3, dp)


    xcrd = arth_d(0.0_dp, dx, (N1))
    ycrd = arth_d(0.0_dp, dy, (N2))
    zcrd = arth_d(0.0_dp, dz, (N3))
     
    !write(6,*) 'test lapl-in 2 by J.Lee '
    xsten(0:nsten)=xcrd(1:nsten+1)
    CALL weights(xcrd(nfd+1), xsten,nsten,mdrv,xcof) 
    !write(6,*) 'test lapl: xsten by J.Lee ', xsten
    !write(6,*) 'test lapl:final xcrd', xcrd(N1), 'final box', CELL(1,1)
    !write(6,*) 'test lapl: xcof by J.Lee ', xcof
    
    !write(6,*) 'test lapl-in 3 by J.Lee '
    xxlp: DO i=1,N1
       il = i-nfd
       ir = i+nfd
          xzlp1: DO k=1,N3
          xylp1: DO j=1,N2
             lap_f(i,j,k) = lap_f(i,j,k) + &
                  DOT_PRODUCT( xcof(0:nsten,nsten,mdrv),f_temp(il:ir,j,k) )
          ENDDO xylp1
          ENDDO xzlp1
    ENDDO xxlp
    
   ! write(6,*) 'test lapl-in 4 by J.Lee '
    ! (2) ALONG THE Y-DIR. : d^2/dy^2  f(x,y,z) 
    
    xsten(0:nsten) = ycrd(1:nsten+1)
    CALL weights(ycrd(nfd+1),xsten,nsten,mdrv,ycof) 
    
    yylp: DO j=1,N2
       jl = j-nfd
       jr = j+nfd
          yzlp1: DO k=1,N3
          yxlp1: DO i=1,N1
             lap_f(i,j,k) = lap_f(i,j,k) + &
                  DOT_PRODUCT( ycof(0:nsten,nsten,mdrv),f_temp(i,jl:jr,k) )
          ENDDO yxlp1
          ENDDO yzlp1
    ENDDO yylp
    
    ! (3) ALONG THE Z-DIR. : d^2/dz^2  f(x,y,z) 
    
    xsten(0:nsten) = zcrd(1:nsten+1)
    CALL weights(zcrd(nfd+1),xsten,nsten,mdrv,zcof)
    
    zzlp: DO k=1,N3
       kl = k-nfd
       kr = k+nfd       
          zylp1: DO j=1,N2
          zxlp1: DO i=1,N1
             lap_f(i,j,k) = lap_f(i,j,k) + &
                  DOT_PRODUCT( zcof(0:nsten,nsten,mdrv),f_temp(i,j,kl:kr) )
          ENDDO zxlp1
          ENDDO zylp1
    ENDDO zzlp

    deallocate(f_temp) 
  CONTAINS
!======================================================================================
  FUNCTION arth_d(first,increment,n)
!======================================================================================
    INTEGER, PARAMETER :: i4b = selected_int_KIND(9)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(i4b), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER(i4b) :: k,k2
    REAL(DP) :: temp
    integer(i4b), parameter  :: npar=16, npar2=8
!======================================================================================
    if (n > 0) arth_d(1)=first
    if (n <= npar) then
       do k=2,n
          arth_d(k)=arth_d(k-1)+increment
       end do
    else
       do k=2,npar2
          arth_d(k)=arth_d(k-1)+increment
       end do
       temp=increment*npar2
       k=npar2
       do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_d
  
  END SUBROUTINE laplacian3

!======================================================================================
  SUBROUTINE weights(xi,x,n,m,c)
!======================================================================================
    USE precision, ONLY: dp
    IMPLICIT NONE
    INTEGER,  INTENT(in)  :: m,n
    REAL(dp), INTENT(in)  :: x(0:n),xi
    REAL(dp), INTENT(out) :: c(0:n,0:n,0:m)
!**************************************************************************************
!-- PURPOSE --
! COMPUTES WEIGHTS FOR THE MTH DERIVATIVE ON ARBITRARILY SPACED POINTS
!
!-- REFERENCE --
! B. FONBERG & D.M. SLOAN, "ACTA NUMERICA 94" EDITED BY A. ISERLES 
! (CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE)
!
!-- INPUT PARAMETERS --
!   XI   POINT AT WHICH THE APPROXIMATIONS ARE TO BE ACCURATE
!    X   X-COORDINATES FOR THE GRID POINTS, ARRAY DIMENSIONED X(0:N)
!    N   THE GRID POINTS ARE AT X(0),X(1),...X(N) (I.E. N+1 IN ALL)
!    M   HIGHEST ORDER OF DERIVATIVE TO BE APPROXIMATED 
!  
!-- OUTPUT PARAMETERS --
!    C   WEIGHTS, ARRAY DIMENSIONED C(0:N,0:N,0:M)
!        ON RETURN, THE ELEMENT C(J,K,I) CONTAINS THE WEIGHT TO BE 
!        APPLIED AT X(K) WHEN THE I:TH DERIVATIVE IS APPROXIMATED 
!        BY A STENCIL EXTENDING OVER X(0),X(1),...,X(J)
!**************************************************************************************
    INTEGER  :: i,j,k,mn
    REAL(dp) :: c1,c2,c3,c4,c5
    real(dp) :: one, zero
!--------------------------------------------------------------------------------------
    one = 1.0_dp
    zero = 0.0_dp

    c(0,0,0)= one
    c1      = one 
    c4      = x(0)-xi
    DO j=1,n
       mn = MIN(j,m) 
       c2 = one
       c5 = c4
       c4 = x(j)-xi
       DO k=0,j-1
          c3=x(j)-x(k)
          c2=c2*c3
          IF (j <= m) c(k,j-1,j) = zero
          c(k,j,0) =c4*c(k,j-1,0)/c3
          DO  i=1,mn
             c(k,j,i)=(c4*c(k,j-1,i)-i*c(k,j-1,i-1))/c3
          ENDDO
       ENDDO
       c(j,j,0)=-c1*c5*c(j-1,j-1,0)/c2
       DO i=1,mn
          c(j,j,i)=c1*(i*c(j-1,j-1,i-1)-c5*c(j-1,j-1,i))/c2
       ENDDO
       c1=c2
    ENDDO
    RETURN
    
  END SUBROUTINE weights
 

END MODULE m_grid_obj
