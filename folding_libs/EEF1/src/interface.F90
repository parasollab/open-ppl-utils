!! file name: interface.F90
!! Guang Song, 03/21/02
!! This is an interface file which makes some subroutines inside modules
!! easily usable by F77/C/C++ by defining some wrapper non-module 
!! subroutines/functions.

	subroutine EEF1_setup_potential(mode)
	use utils
	implicit none
	integer, intent(out) :: mode
	call setup_potential(mode)
	return
	end subroutine EEF1_setup_potential


	
	subroutine EEF1_potential(dim, phi, psi, omega, f)
	use utils
	implicit none
	integer, intent(in) :: dim
	real(dp), dimension(dim), intent(in) :: phi, psi, omega
	real(dp), intent(out) :: f
	call potential(dim, phi, psi, omega, f)
	return
	end subroutine EEF1_potential

