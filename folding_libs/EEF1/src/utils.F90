!-----------------------------------------------------------------------
MODULE utils
!-----------------------------------------------------------------------
  use energy; use rmsd
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
subroutine setup_potential(mode)
!-----------------------------------------------------------------------
! setup potential and minimize native. 
!-----------------------------------------------------------------------
  implicit none
  real(dp) :: f
  integer, intent(out) :: mode

  call read_input_files    ! read input files

  call initialize_protein  ! initialize protein

  call read_pdb            ! read pdb file

  call cartesian2internal(0) ! get native angles
  call internal2cartesian  ! get coordinates for atoms not defined in pdb
    
  call initialize_energy(0)   ! initialize energy for cartesian min

  if (minimize_at_setup == 'YES') then
     ! minimize energy with all degrees of freedom starting from the current
     ! structure
     call minimize_energy(f, 0) ! cartesian min
     call cartesian2internal(1)    ! save minimized conf in internal coord
  end if
  protein_native = protein

  if (t_coherent == 1) then
     N_update = 0
     allocate(R_gr_prev(3,tot_num_gr))
  end if

  call finalize_energy(0) ! finalize energy for cartesian min

  if (include_ang_type == 0) then
     mode = 0 ! cartesian min
     include_ang_type = 1
  else
     mode = 1 ! internal min
  end if

  call initialize_energy(mode)   ! initialize energy for internal min

  ! one more minimization in the new minimization mode
  call minimize_energy(f, mode)

  if (mode == 0) then
     call cartesian2internal(1)    ! save minimized conf in internal coord
  end if
  protein_native = protein

end subroutine setup_potential
!-----------------------------------------------------------------------
subroutine finalize_potential(mode)
  integer, intent(in) :: mode

  call finalize_energy(mode)
  call finalize_protein
  if (t_coherent == 1) then
     deallocate(R_gr_prev)
  end if
end subroutine finalize_potential
!-----------------------------------------------------------------------
subroutine read_input_files
!-----------------------------------------------------------------------
! read input files
!-----------------------------------------------------------------------
  implicit none
  call read_user_input    ! read user input
  call read_topology_geo  ! read topology file for geometry setup
  call read_parameters    ! read parameter file
  call read_topology_eng  ! read topology file for energy setup
  call initialize_ss_bnd
end subroutine read_input_files
!-----------------------------------------------------------------------
subroutine initialize_protein
!-----------------------------------------------------------------------
! read coords from a pdb file, convert the cartesian coordinates into 
! internal coordinates.
!-----------------------------------------------------------------------
  implicit none
  integer :: num_res, i, res_type, n, j, res_no
  integer :: k, ang_type, imax, bond, ang_no, ang1
  character(len=len_fname) :: pdb_out
  type(reference_residue_type) :: ref
  
!-----------------------------------------------------------------------
  call get_sequence
! check sequence given in the input here!
  num_res = protein%num_res
!-----------------------------------------------------------------------
  ! setup reference residues 
  call setup_ter_ref('N')
  call setup_ter_ref('C')
  call setup_ter_ref_eng('N')
  call setup_ter_ref_eng('C')

  allocate(ref_res(num_res), ref_res_eng(num_res))
  do i = 1, num_res
     if (i==1) then  
        ref_res(i) = Nter_ref
        ref_res_eng(i) = Nter_ref_eng
     else if (i==num_res) then  
        ref_res(i) = Cter_ref
        ref_res_eng(i) = Cter_ref_eng
     else  
        res_type = protein%residue(i)%res_type
        ref_res(i) = residue_ref(res_type)
        ref_res_eng(i) = res_ref_eng(res_type)
     end if
  end do

  ! setup parameter links
  call setup_eng_para
  call setup_ss_bnd
  call fix_angles

  allocate(i_ww(2,3*protein%num_res-3))
  ang_no = 0; ang1 = 1
  do i = 1, protein%num_res
     ref = ref_res(i)
     ref_res(i)%ww_fixed = 1
     do j = 1, ref%num_ang
        ang_type = ref%ang_type(j); imax = ref%na_ang(j)-1
        bond = ref%a_idx(imax,j)
        if (ang_type == 1 .and. &
             (bond==1.or.bond==2.or.bond==3) .and. &
             (i/=1.or.ref%a_idx(1,j)/=0) .and. &
             (i/=num_res.or.bond/=3) ) then
           if (ang_no == 0 .and. ang1 == 1) then ! skip first ang
              ang1 = 0
           else
              ang_no = ang_no + 1
              i_ww(1:2,ang_no)=(/ i, j /)
              ref_res(i)%ww_fixed(j) = 0
           end if
        end if
     end do
  end do

  do i = 1, protein%num_res
     ref = ref_res(i)
     ref_res(i)%ww_change = 0
     do j = 1, ref%num_ang
        if (ref%ww_fixed(j)==0) then
           ang_type = ref%ang_type(j); bond = ref%a_idx(2,j)
           do k = 1, ref%num_ang
              if (ref%ww_fixed(k) == 1 .and. &
                   ref%ang_type(k) == ang_type .and. &
                   ref%a_idx(2,k) == bond) then
                 ref_res(i)%ww_change(k) = j
              end if
           end do
        end if
     end do
  end do

end subroutine initialize_protein
!-----------------------------------------------------------------------
subroutine finalize_protein
  deallocate(ref_res, ref_res_eng, i_ww)
end subroutine finalize_protein
!-----------------------------------------------------------------------
subroutine evaluate_energy(f,i_print)
  implicit none
  integer, intent(in) :: i_print
  real(dp), intent(out) :: f
  real(dp), dimension(:), allocatable :: g
  allocate(g(tot_num_dof))
  call energy_and_gradient(f,g,i_print,0,0,0)
  write(*,900)
  deallocate(g)
900 format(' ---------------------------------------------------')
end subroutine evaluate_energy
!-----------------------------------------------------------------------
subroutine minimize_energy(f, mode)
  implicit none
  integer, intent(in) :: mode
  real(dp), intent(out) :: f
  type(protein_type) :: protein0
  integer :: iprint_Emin, status
  real(dp) :: error
  real(dp), dimension(3,3) :: U

! local minimization
  protein0 = protein
  if (print_level > 0) then
     iprint_Emin = 0
     write(*,"(a)") ' Before Minimization'
     call evaluate_energy(f,1)
     write(*,"(a)") ' Minimizing Energy...'
  else
     iprint_Emin = -1
  end if

  call localmin_energy(f, iprint_Emin, status, mode)
  if (status == 0) then
     write(*,"(a)") ' Error occurred during local minimization of Energy.'
     !stop
  end if
  if (print_level > 0) then
     call calc_rmsd(protein0,protein,'aa',1,U,error)
     write(*,"(A, F10.3, A)") ' All atom rmsd after Energy minimization =', error, ' A'
     write(*,900)
     write(*,"(a)") ' After  Minimization'
     call evaluate_energy(f,1)
  end if

900 format(' ---------------------------------------------------')
end subroutine minimize_energy
!-----------------------------------------------------------------------
subroutine potential(dim, phi, psi, omega, f)
!-----------------------------------------------------------------------
! calculates energy
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: dim
  real(dp), dimension(dim), intent(in) :: phi, psi, omega 
  real(dp), intent(out) :: f
  real(dp) :: Evdw, Ees, Eslv, Eb
  real(dp), dimension(3) :: Eint
  integer ::  i, j
 
  ! set all angles to native angles, and change phi, psi, and omg
  protein = protein_native

  i = 1
  j = 1
  protein%residue(i_ww(1,j))%w(i_ww(2,j)) = psi(i)
  j = 2
  !protein%residue(i_ww(1,j))%w(i_ww(2,j)) = omega(i)
  do i = 2, protein%num_res - 1 
     j = j + 1
     protein%residue(i_ww(1,j))%w(i_ww(2,j)) = phi(i-1)
     j = j + 1
     protein%residue(i_ww(1,j))%w(i_ww(2,j)) = psi(i)
     j = j + 1
     !protein%residue(i_ww(1,j))%w(i_ww(2,j)) = omega(i)
  end do
  i = protein%num_res - 1
  j = j + 1
  protein%residue(i_ww(1,j))%w(i_ww(2,j)) = phi(i)

! need other internal2cartesian for phi,psi,omg change
  call internal2cartesian_ww

  call pairlist_update(t_coherent)

  call bond_energy_only(Eb)

  call angle_energy_only(Eint)
 
  call non_bonded_energy_only(Evdw,Ees,Eslv)

  f = Evdw + Eb + sum(Eint) + Ees + Eslv

end subroutine potential
!-----------------------------------------------------------------------
subroutine potential_after_min(dim, phi, psi, omega, f, mode)
!-----------------------------------------------------------------------
! calculates energy
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: dim, mode
  real(dp), dimension(dim), intent(inout) :: phi, psi, omega 
  real(dp), intent(out) :: f
  real(dp) :: Evdw, Ees, Eslv, Eb
  real(dp), dimension(3) :: Eint
  integer ::  i, j
 
  ! set all angles to native angles, and change phi, psi, and omg
  protein = protein_native

  i = 1
  j = 1
  protein%residue(i_ww(1,j))%w(i_ww(2,j)) = psi(i)
  j = 2
  !protein%residue(i_ww(1,j))%w(i_ww(2,j)) = omega(i)
  do i = 2, protein%num_res - 1 
     j = j + 1
     protein%residue(i_ww(1,j))%w(i_ww(2,j)) = phi(i-1)
     j = j + 1
     protein%residue(i_ww(1,j))%w(i_ww(2,j)) = psi(i)
     j = j + 1
     !protein%residue(i_ww(1,j))%w(i_ww(2,j)) = omega(i)
  end do
  i = protein%num_res - 1
  j = j + 1
  protein%residue(i_ww(1,j))%w(i_ww(2,j)) = phi(i)

  call internal2cartesian_ww
  call minimize_energy(f, mode)
  if (mode == 0) then
     call cartesian2internal(1)    ! save minimized conf in internal coord
  end if
end subroutine potential_after_min
!-----------------------------------------------------------------------
subroutine get_phipsi(dim, phi, psi, omega)
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: dim
  real(dp), dimension(dim), intent(out) :: phi, psi, omega 
  integer ::  i, j
 
  i = 1
  j = 1
  psi(i) = protein%residue(i_ww(1,j))%w(i_ww(2,j))
  j = 2
  omega(i) = protein%residue(i_ww(1,j))%w(i_ww(2,j))
  do i = 2, protein%num_res - 1 
     j = j + 1
     phi(i-1) = protein%residue(i_ww(1,j))%w(i_ww(2,j))
     j = j + 1
     psi(i) = protein%residue(i_ww(1,j))%w(i_ww(2,j))
     j = j + 1
     omega(i) = protein%residue(i_ww(1,j))%w(i_ww(2,j))
  end do
  i = protein%num_res - 1
  j = j + 1
  phi(i) = protein%residue(i_ww(1,j))%w(i_ww(2,j))


end subroutine get_phipsi
!-----------------------------------------------------------------------
subroutine ang_norm(dim, phi0, psi0, omega0, phi, psi, omega, norm)
  implicit none
  integer, intent(in) :: dim
  real(dp), dimension(dim), intent(in) :: phi0, psi0, omega0, phi, psi, omega 
  real(dp), intent(out) :: norm
  integer :: i, j
  real(dp) :: diff, sum

  sum = 0
  do i = 1, protein%num_res - 1
     do j = 1, 3
        if (j==1) then
           diff = phi(i) - phi0(i)
        else if (j==2) then
           diff = psi(i) - psi0(i)
        else if (j==3) then
           diff = omega(i) - omega0(i)
        end if
        diff = mod(diff,two_pi) 
        if (diff > pi) then
           diff = diff - two_pi
        else if (diff < -pi) then
           diff = diff + two_pi
        end if
        sum = sum + diff*diff
     end do
  end do
  norm = sqrt(sum)

end subroutine ang_norm
!-----------------------------------------------------------------------
subroutine test_timing
  implicit none
  real(dp) :: f
  real(dp) :: time1, time2
  integer :: i
  real(dp), dimension(:), allocatable :: phi, psi, omega
  integer :: dim

  call my_timer(time1)     
  dim = protein%num_res - 1
  allocate(phi(dim), psi(dim), omega(dim))
  call get_phipsi(dim, phi, psi, omega)

  do i = 1, 100
     call  potential(dim, phi, psi, omega, f)
  end do
  call my_timer(time2)
  write(*,*)'Total User Time = ', (time2-time1)/100.0*1000.0,'msec'

  deallocate(phi, psi, omega)

end subroutine test_timing
!-----------------------------------------------------------------------
subroutine my_timer(ttime)
!-----------------------------------------------------------------------
  use globals
  implicit none
  real(dp) :: ttime
  character(len=8) :: date
  character(len=10) :: time
  character(len=5) :: zone
  integer, dimension(8) :: values
  integer, dimension(12) :: days_in_month
  integer :: days

  days_in_month = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  call date_and_time(date,time,zone,values)
  if (mod(values(1),4) == 0) days_in_month(2) = 29
  days = sum(days_in_month(1:(values(2)-1)))+values(3)
  ttime = 1.0d-3*(values(8)+1.0d3*(values(7)+ &
       60*(values(6)+60*(values(5)+24*days))))
end subroutine my_timer
!-----------------------------------------------------------------------
end MODULE utils
!-----------------------------------------------------------------------
