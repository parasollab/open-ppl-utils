!-----------------------------------------------------------------------
program main
!-----------------------------------------------------------------------
  use energy; use utils
!-----------------------------------------------------------------------
  implicit none
  real(dp) :: f, norm, time1, time2
  real(dp), dimension(:), allocatable :: phi0, psi0, omega0, phi, psi, omega
  integer :: min, dim, mode
  character(len=len_fname) :: out_pdb_name

  ! setup potential and native 
  call setup_potential(mode)

  dim = protein%num_res - 1
  allocate(phi(dim), psi(dim), omega(dim), phi0(dim), psi0(dim), omega0(dim))

  ! get current (native) angles
  call get_phipsi(dim, phi, psi, omega)
  ! initial (native) energy
  call potential(dim, phi, psi, omega, f)
  write(*,997) ' Native Energy =', f

  ! set all angles to pi and calc potential
  phi0 = pi; psi0 = pi; omega0 = omega
  ! get potential for all-trans conformation
  call potential(dim, phi0, psi0, omega0, f)
  write(*,997) ' Energy for all-trans before min =', f

  call my_timer(time1)
  ! minimize in energy starting from all-trans
  call potential_after_min(dim, phi0, psi0, omega0, f, mode)
  call my_timer(time2)
  write(*,997) ' Energy for all-trans after  min =', f
  write(*,998) ' Time for min =', time2-time1, ' sec'

  ! get angles after minimization
  call get_phipsi(dim, phi, psi, omega)

  ! diff bet before and after min: shold not change if only side chain angles are varied
  call ang_norm(dim, phi0, psi0, omega0, phi, psi, omega, norm)
  write(*,997) ' Norm =', norm

!  out_pdb_name = 'out.pdb'
!  call write_pdb(out_pdb_name)

  deallocate(phi, psi, omega, phi0, psi0, omega0)

  call finalize_potential(mode)
997 format(a,e15.5)
998 format(a,e15.5,a)
end program main
!-----------------------------------------------------------------------
