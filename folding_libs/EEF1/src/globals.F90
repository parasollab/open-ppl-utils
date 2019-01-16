!----------------------------------------------------------------------------
MODULE GLOBALS
!----------------------------------------------------------------------------
implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: max_res = 1000  ! max no of residues in a protein
  integer, parameter :: max_atm = 19   ! max no of atoms in a residue
  integer, parameter :: max_gr  = 6    ! max no of atom groups in a residue
  integer, parameter :: max_bnd = 20   ! max no of bonds in a residue
  integer, parameter :: max_ang = 38   ! max no of angles in a residue
  integer, parameter :: max_ang_eng=55 ! max no of ang in a res for energy calc
  integer, parameter :: max_atm_typ=300! max no of atm typ with different charges
  integer, parameter :: max_atm_cls=19 ! max no of extended atm class
  integer, parameter :: max_ang_prm=140! max no of ang in parameters
  integer, parameter :: max_bnd_prm=40 ! max no of bnd in parameters
  integer, parameter :: max_crd = 4    ! max no of coordination number
  integer, parameter :: max_prm = 639  ! max no of parameters
  integer, parameter :: max_ss = 20    ! max no of disulfide bnds
  integer, parameter :: len_fname = 50 ! max length of a file name
  integer, parameter :: len_s = 20     ! max length of a short string
  integer, parameter :: len_l = 80     ! max length of a long string
  integer :: num_prm                   ! no of parameters
  integer :: tot_num_atm, tot_num_ang, tot_num_gr
  integer :: tot_num_ang_e, tot_num_bnd, tot_num_dof ! no of atm, ang, gr, ang_e in a protein
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197d0
  real(dp), parameter :: deg2rad=pi/180.0d0,rad2deg=180.0d0/pi,two_pi=2.0d0*pi
  real(dp), dimension(3), parameter :: origin=(/0.0d0,0.0d0,0.0d0/)
  real(dp), dimension(4), parameter :: unit_q=(/1.0d0,0.0d0,0.0d0,0.0d0/)
  real(dp), dimension(3,3), parameter :: unit_U=reshape( &
       (/1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))
  ! input file names
  character(len=len_fname), parameter :: &
       infile_topo_geo='topology_geo.in',infile_topo_eng='topology_eng.in', &
       infile_parameters='parameters.in',infile_user_input='user_input.in'
!----------------------------------------------------------------------------
  character(len=len_fname) :: pdb_file_name
  character(len=10) :: time_coherent, minimize_at_setup
  integer :: print_level, t_coherent
  integer :: N_update
!----------------------------------------------------------------------------
  integer :: include_ang_type
  ! index for backbone dihedral angles
  integer, dimension(:,:), allocatable :: i_ww
!----------------------------------------------------------------------------
  ! tolerances in energy minimization
  real(dp) :: factr_eng, pgtol_eng
  integer :: max_iter_lbfgsb
!----------------------------------------------------------------------------
type residue_type
  integer :: res_type                ! numeric representation of residue type
  real(dp), dimension(3, max_atm+1) :: R      ! absolute coordinate of atoms
  real(dp), dimension(3, -1:max_bnd) :: b     ! bond vector
  real(dp), dimension(-1:max_bnd) :: b_len    ! bond length
  real(dp), dimension(max_ang) :: w           ! angles
  real(dp), dimension(max_ang) :: sin_w       ! sin(w)
  integer, dimension(-1:max_atm+1) :: atom_pdb   ! 1, if in orig pdb, 0, if not
end type residue_type
!----------------------------------------------------------------------------
type ss_bnd_type
  real(dp) :: b_len
  real(dp),  dimension(3) :: b
end type ss_bnd_type
!----------------------------------------------------------------------------
type protein_type
  integer :: num_res, num_ss_bnd
  type(residue_type), dimension(max_res) :: residue
  type(ss_bnd_type), dimension(max_ss) :: ss_bnd
end type protein_type
!----------------------------------------------------------------------------
type(protein_type) :: protein, protein_native
!----------------------------------------------------------------------------
type reference_residue_type ! reference residue for geometry calc
  character(len=4) :: residue_type            ! 3 letter residue type
  integer :: num_atm, num_bnd, num_ang        ! no of atm, bnd, & ang 
  character(len=4), dimension(-1:max_atm+1) :: atom_type !name in pdb
  integer, dimension(2,-1:max_bnd) :: b_idx   ! head and tail of bnd
  real(dp), dimension(max_bnd) :: b_len       ! bond length
  integer, dimension(3,0:max_ang) :: a_idx    ! indices for bonds forming ang
  real(dp), dimension(max_ang) :: w0          ! reference angles
  integer, dimension(0:max_ang) :: ang_type   ! dihedral(1) or bond angle(2) 
  integer, dimension(0:max_ang) :: na_ang ! no of atms consisting of ang
  integer, dimension(max_ang) :: num_dpn_bnd  ! no of dependent bonds
  integer, dimension(max_crd-1,max_ang) :: i_dpn_bnd ! index for dpn bnds
  integer, dimension(0:max_ang) :: num_dpn_ang      ! no of dependent angles
  integer, dimension(max_crd-1,0:max_ang) :: i_dpn_ang ! index for dpn angs
  integer, dimension(max_ang) :: ang_fixed    ! 1 (fixed) or 0 (not fixed) 
  integer, dimension(max_ang) :: ang_change   ! dih ang to change together
  integer, dimension(max_ang) :: fixed_branch ! fixed ang sharing axis
  integer, dimension(max_ang) :: ww_fixed    ! 1 (fixed) or 0 (not fixed) 
  integer, dimension(max_ang) :: ww_change   ! dih ang to change together
end type reference_residue_type
!----------------------------------------------------------------------------
type ref_res_eng_type ! reference residue for energy calc
  integer :: num_atm, num_ang, num_gr         ! no of atm, ang, gr
  integer, dimension(max_atm) :: atm_cls      ! atom class no
  integer, dimension(max_atm) :: q_idx        ! indx for charge
  integer, dimension(max_atm) :: gr_no        ! group number
  integer, dimension(max_gr) :: na_gr         ! no of atoms in group
  integer, dimension(max_ang_eng) :: ang_type ! dih(1), bnd(2), imp(3)
  integer, dimension(max_ang_eng) :: na_ang   ! no of atms consisting of ang
  integer, dimension(4,max_ang_eng) :: aa_idx ! indx for atm forming ang
  integer, dimension(2,max_ang_eng) :: cor_ang! indx for corresponding ang in ref_res
  integer, dimension(max_ang_eng) :: prm      ! indx for prm
  integer, dimension(max_bnd) :: bnd_prm      ! indx for prm
  integer, dimension(4,max_ang_eng) :: atm_idx! indx for new atm no in ang
  integer, dimension(2,max_bnd) :: b_idx      ! indx for new atm no in bnd
end type ref_res_eng_type
!----------------------------------------------------------------------------
type eng_para_type ! energy parameters
  integer :: num_a_cls, num_q_typ, num_bnd, num_ang! no of atm class, charge type, ang
  integer, dimension(max_ang_prm) :: ang_type 
  character(len=4), dimension(2,max_atm_typ) :: q_type  ! residue and atom
  character(len=4), dimension(max_atm_cls) :: atom_cls  ! atom type
  character(len=4), dimension(2,max_bnd_prm) :: atm_in_bnd ! atm in bnd
  character(len=4), dimension(4,max_ang_prm) :: atm_in_ang ! atm in ang
  real(dp), dimension(max_atm_typ) :: charge        ! \ fix
  real(dp), dimension(4,max_atm_cls) :: LJ_para     !  \ fix
  real(dp), dimension(5,max_atm_cls) :: slv_para    !   \ fix
  real(dp), dimension(3,max_ang_prm) :: ang_para    !   / parameters  
  real(dp), dimension(2,max_bnd_prm) :: bnd_para    !  / 
  real(dp) :: E14fac                                ! / 
  integer, dimension(max_atm_typ) :: charge_on      ! \ fix 
  integer, dimension(4,max_atm_cls) :: LJ_on        !  \ whether parameters 
  integer, dimension(5,max_atm_cls) :: slv_on       !  / are fixed or not
  integer, dimension(3,max_ang_prm) :: ang_on       ! /
  integer, dimension(2,max_ang_prm) :: bnd_on       !/
  integer :: E14fac_on
  real(dp), dimension(5,max_atm_cls,max_atm_cls) :: aux_para
  ! aux_para(1:5,i,j): eps_LJ, sigma_sqr, eps14_LJ, sigma14_sqr, slv_factor
end type eng_para_type
!----------------------------------------------------------------------------
type ref_ss_bnd_type
  integer, dimension(2) :: res_no, s_atm_no, b_idx
end type ref_ss_bnd_type
!----------------------------------------------------------------------------
integer, parameter :: num_ref_res = 20
type(reference_residue_type), dimension(num_ref_res) :: residue_ref
type(reference_residue_type) :: Nter_ref, Cter_ref
type(reference_residue_type), dimension(:), allocatable :: ref_res
type(reference_residue_type) :: ss_ref
type(ref_res_eng_type), dimension(num_ref_res) :: res_ref_eng
type(ref_res_eng_type) :: Nter_ref_eng, Cter_ref_eng
type(ref_res_eng_type), dimension(:), allocatable :: ref_res_eng
type(ref_res_eng_type) :: ss_ref_eng
type(ref_ss_bnd_type), dimension(max_ss) :: ref_ss
type(eng_para_type) :: eng_para
!----------------------------------------------------------------------------
! numeric assignment based on Thomas-Dill classification
character(len=4), dimension(num_ref_res), parameter :: residue_name = (/ &
       'CYS ', 'MET ', 'PHE ', 'ILE ', 'LEU ', & !  1-5
       'VAL ', 'TRP ', 'TYR ', 'ALA ', 'GLY ', & !  6-10 
       'THR ', 'SER ', 'GLN ', 'ASN ', 'GLU ', & ! 11-15
       'ASP ', 'HIS ', 'ARG ', 'LYS ', 'PRO ' /) ! 16-20
!----------------------------------------------------------------------------
END MODULE GLOBALS
!----------------------------------------------------------------------------
