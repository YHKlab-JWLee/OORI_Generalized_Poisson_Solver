module m_iorho_serial
  use precision, only : dp, sp, grid_p
  implicit none
  private
  public :: read_rho_serial, write_rho_serial

contains

  subroutine read_rho_serial(fname, cell, mesh, nspin, rho)
    character(len=*), intent(in) :: fname
    real(dp), intent(out) :: cell(3,3)
    integer, intent(out) :: mesh(3)
    integer, intent(out) :: nspin
    real(grid_p), allocatable, intent(out) :: rho(:,:)

    integer :: iu, is, iz, iy, ind, np
    real(sp), allocatable :: temp(:)

    open(newunit=iu, file=fname, form='unformatted', status='old', action='read')
    read(iu) cell
    read(iu) mesh, nspin

    np = mesh(1)*mesh(2)*mesh(3)
    allocate(rho(np,nspin), temp(mesh(1)))

    do is = 1, nspin
      ind = 0
      do iz = 1, mesh(3)
        do iy = 1, mesh(2)
          read(iu) temp
          rho(ind+1:ind+mesh(1),is) = real(temp,kind=grid_p)
          ind = ind + mesh(1)
        end do
      end do
    end do

    close(iu)
    deallocate(temp)
  end subroutine read_rho_serial

  subroutine write_rho_serial(fname, cell, mesh, nspin, rho)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: mesh(3)
    integer, intent(in) :: nspin
    real(grid_p), intent(in) :: rho(:,:)

    integer :: iu, is, iz, iy, ind
    real(sp), allocatable :: temp(:)

    open(newunit=iu, file=fname, form='unformatted', status='replace', action='write')
    write(iu) cell
    write(iu) mesh, nspin

    allocate(temp(mesh(1)))
    do is = 1, nspin
      ind = 0
      do iz = 1, mesh(3)
        do iy = 1, mesh(2)
          temp = real(rho(ind+1:ind+mesh(1),is),kind=sp)
          write(iu) temp
          ind = ind + mesh(1)
        end do
      end do
    end do

    deallocate(temp)
    close(iu)
  end subroutine write_rho_serial

end module m_iorho_serial
