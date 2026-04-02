module precision
  implicit none
  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: grid_p = dp
end module precision
