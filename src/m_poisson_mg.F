program m_poisson_mg_main
  use precision,      only : dp, grid_p
  use m_poison_mg,    only : generalized_poison_solver
  use m_boundary_condition, only : BC_PERIODIC, parse_boundary_condition
  use m_iorho_serial, only : read_rho_serial, write_rho_serial
  implicit none

  character(len=256) :: fin, fout, arg, bc_name
  logical :: bc_known
  integer :: narg, nfd, bc_type
  integer :: mesh(3), nspin, ng
  real(dp) :: cell(3,3), ehar, stress(3,3)
  real(grid_p), allocatable :: rho2d(:,:), rho(:), vhar(:), esp(:), vhar2d(:,:)

  narg = command_argument_count()
  if (narg < 2) then
    write(*,'(a)') 'Usage: m_poisson_mg.x <rho_input.bin> <vhar_output.bin> [nfd] [bc]'
    write(*,'(a)') '  bc: periodic|pbc|dirichlet_zero|dirichlet|open (default: periodic)'
    write(*,'(a)') '  optional env fallback: OORI_POISSON_BC=<bc>'
    stop 1
  end if

  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  bc_name = ''

  nfd = 2
  if (narg >= 3) then
    call get_command_argument(3, arg)
    read(arg,*) nfd
  end if

  bc_type = BC_PERIODIC
  if (narg >= 4) then
    call get_command_argument(4, bc_name)
    call parse_boundary_condition(trim(bc_name), bc_type, bc_known)
    if (.not. bc_known) then
      write(*,'(a,1x,a)') 'Warning: unknown bc from argv, fallback to dirichlet_zero =>', trim(bc_name)
    end if
  else
    call get_environment_variable('OORI_POISSON_BC', bc_name)
    if (len_trim(bc_name) > 0) then
      call parse_boundary_condition(trim(bc_name), bc_type, bc_known)
      if (.not. bc_known) then
        write(*,'(a,1x,a)') 'Warning: unknown bc from OORI_POISSON_BC, fallback to dirichlet_zero =>', trim(bc_name)
      end if
    end if
  end if

  call read_rho_serial(trim(fin), cell, mesh, nspin, rho2d)
  ng = mesh(1)*mesh(2)*mesh(3)

  allocate(rho(ng), vhar(ng), esp(ng), vhar2d(ng,nspin))
  rho = rho2d(:,1)
  esp = 1.0_grid_p
  vhar = 0.0_grid_p

  call generalized_poison_solver(cell, mesh(1), mesh(2), mesh(3), nfd, rho, esp, vhar, ehar, stress, bc_type)

  vhar2d = 0.0_grid_p
  vhar2d(:,1) = vhar
  call write_rho_serial(trim(fout), cell, mesh, nspin, vhar2d)

  write(*,'(a,1x,es16.8)') 'Hartree energy:', ehar

  deallocate(rho2d, rho, vhar, esp, vhar2d)
end program m_poisson_mg_main
