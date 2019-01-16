!-----------------------------------------------------------------------
MODULE in_out
  use globals
  type new_numbering_type
     ! atom numbered in a different way and atom at the end of the branch
     character(len=4), dimension(max_atm) :: atom, endc 
     integer, dimension(max_atm) :: ring ! 1, if member of a ring, and 0, if not
  end type new_numbering_type
  type(new_numbering_type), dimension(num_ref_res) :: new_numbering
!----------------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
subroutine get_sequence()
!-----------------------------------------------------------------------
! gets sequence from a pdb file
!-----------------------------------------------------------------------
  implicit none
  character(len=4) :: resname
  character(len=len_s) :: string, resnum, resnum0
  character(len=len_l) :: lstring
  character(len=len_s), dimension(10) :: sstring
  integer :: i, res_no, unit = 10, num_string, ioerror, ss_bnd_no
!-----------------------------------------------------------------------
  open(unit, file = trim(pdb_file_name))
! get ss bond
  ss_bnd_no = 0
  do    !endless loop, exit at ATOM
     read(unit, 50) string
     if (string == 'ATOM') exit
     if (string == 'SSBOND') then
        ss_bnd_no = ss_bnd_no + 1 
        backspace(unit)
        read(unit,*) string, i, &
             string, ref_ss(ss_bnd_no)%res_no(1), &
             string, ref_ss(ss_bnd_no)%res_no(2)
     end if
  end do
  protein%num_ss_bnd = ss_bnd_no
  backspace(unit)
  res_no = 0; resnum0 = '-100' ! start with something unlikely
  do    !endless loop, exit if not ATOM
     read(unit,65,iostat=ioerror)lstring
     if (ioerror < 0) exit ! end of file
     call find_num_string(lstring, num_string, sstring)
     string=sstring(1); resname=sstring(4); resnum=sstring(num_string)
     if (string /= 'ATOM') exit
     if (resnum /= resnum0) then
        res_no = res_no + 1; resnum0 = resnum
        call find_res(resname, protein%residue(res_no)%res_type)
     end if
  end do
  protein%num_res = res_no
  close(unit)
50 format (a6)
65 format (a28)
end subroutine get_sequence
!-----------------------------------------------------------------------
subroutine read_pdb()
!-----------------------------------------------------------------------
! reads a pdb file and gets coordinates and bond vectors.
!-----------------------------------------------------------------------
  implicit none
  type(residue_type) :: res
  type(reference_residue_type) :: ref
  character(len=4) :: atom_type, atom_saved, resname
  character(len=1) :: h_no
  character(len=len_s) :: string, resnum, resnum0, atmname
  character(len=len_l) :: lstring
  character(len=len_s), dimension(10) :: sstring
  integer :: i, res_no, unit = 10, j, idx1, idx2, num_res, num_string
  integer :: atom_tried, ioerror
  integer, dimension(2) :: ss_res, ss_atm
  real(dp) :: x, y, z
  real(dp), dimension(3) :: b
!-----------------------------------------------------------------------
  open(unit, file = trim(pdb_file_name))
  do    !endless loop, exit at ATOM
     read(unit, 50) string
     if (string == 'ATOM') exit
  end do
  backspace(unit)
!-----------------------------------------------------------------------
! get number of residues
  res_no = 0; resnum0 = '-100' ! start with something unlikely.
  do    !endless loop, exit if not ATOM
     read(unit,65,iostat=ioerror)lstring
     if (ioerror < 0) exit ! end of file
     call find_num_string(lstring, num_string, sstring)
     string = sstring(1); if (string /= 'ATOM') exit
     resnum = sstring(num_string)
     if(resnum /= resnum0) then
        res_no = res_no + 1; resnum0 = resnum
     end if
  end do
  num_res = res_no
  protein%num_res = num_res
  rewind(unit) ! goto the beginning of the file
!-----------------------------------------------------------------------
  ! initialize
  do i = 1, num_res 
     protein%residue(i)%atom_pdb = 0 ! atom_pdb: not read from pdb yet
     protein%residue(i)%b = 0.0d0
     protein%residue(i)%b_len(1:) = ref_res(i)%b_len
  end do
!-----------------------------------------------------------------------
! get atom coordinates
  do    !endless loop, exit at ATOM
     read(unit, 50) string; if(string == 'ATOM') exit
  end do
  backspace(unit)
  res_no = 0; resnum0 = '-100'
  do    !endless loop, exit if not ATOM
     read(unit,65,iostat=ioerror)lstring
     if (ioerror < 0) exit ! end of file
     call find_num_string(lstring, num_string, sstring)
     string = sstring(1); if (string /= 'ATOM') exit
     atmname= sstring(3); resname=sstring(4)
     resnum = sstring(num_string)
     backspace(unit); read(unit,*)(sstring(i),i=1,num_string),x,y,z
     if(resnum /= resnum0) then
        res_no = res_no + 1; resnum0 = resnum
        ref = ref_res(res_no)
     end if
     h_no = atmname(1:1)
     if (h_no=='1'.or.h_no=='2'.or.h_no=='3') then
        atmname = atmname(2:len(atmname))
     else
        h_no = '0'
     end if
     call fill_atom_type3(atmname)
     atom_type = atmname(1:3)//h_no
     atom_tried = 0
999  continue
     ! consider different representations for atom types
     if (atom_type(1:1) == 'D') atom_type = 'H'//atom_type(2:)
     if (atom_type(1:2) == 'HN') atom_type = 'H0'//atom_type(3:)
     if (res_no == 1 .and. atom_type == 'H001') atom_type = 'HT10'
     if (res_no == 1 .and. atom_type == 'H002') atom_type = 'HT20'
     if (res_no == 1 .and. atom_type == 'H003') atom_type = 'HT30'
     if (res_no == protein%num_res .and. (atom_type == 'O000' .or. &
          atom_type == 'O100')) atom_type = 'OT10'
     if (res_no == protein%num_res .and. (atom_type == 'OXT0' .or. &
          atom_type == 'O200')) atom_type = 'OT20'
     if (resname == 'ILE' .and. atom_type == 'CD00') atom_type = 'CD10'
     if (atom_tried == 0) atom_saved = atom_type
     do i = 1, ref%num_atm ! save coordinates
        if (atom_type == ref%atom_type(i) .and. &
             protein%residue(res_no)%atom_pdb(i) /= 1) then
           protein%residue(res_no)%R(:,i) = (/ x, y, z /)
           protein%residue(res_no)%atom_pdb(i) = 1
           atom_tried = 3
           exit
        end if
     end do
     if (atom_tried < 3 .and. atom_saved(1:1) == 'H' .and. &
          atom_saved(4:4) == '0' .and. (atom_saved(3:3) == '1' .or. &
          atom_saved(3:3) == '2' .or. atom_saved(3:3) == '3')) then
        if (atom_tried == 0) then
           atom_type = atom_saved(1:2)//'0'//atom_saved(3:3)
        else if (atom_tried == 1) then
           atom_type = atom_saved(1:3)//'1'
        else if (atom_tried == 2) then
           atom_type = atom_saved(1:3)//'2'
        end if
        atom_tried = atom_tried + 1
        goto 999
     end if
  end do
  close(unit)

  call find_atm(ref_res(num_res),'OT10',idx1)
  call find_atm(ref_res(num_res),'OT20',idx2)
  if (protein%residue(num_res)%atom_pdb(idx2) /= 1) then
     protein%residue(num_res)%atom_pdb(idx1) = 0 
  end if

50 format (a6)
55 format (a7, i4, 1x, a1, a3, 1x, a4, 1x, i4, 4x, 3f8.3)
65 format (a28)
!-----------------------------------------------------------------------
! more processes with coord, such as calc bnd vectors  
  do i = 1, num_res
     res = protein%residue(i); ref = ref_res(i)
     if (i < num_res) then ! get coord for +N 
        res%R(:,ref%num_atm+1) = protein%residue(i+1)%R(:,1)
        res%atom_pdb(ref%num_atm+1) = 1
     end if
     do j = 1, ref%num_bnd ! calc bond vectors
        idx1 = ref%b_idx(1,j); idx2 = ref%b_idx(2,j)
        if (res%atom_pdb(idx1) == 1 .and. res%atom_pdb(idx2) == 1) then
           b = res%R(:,idx2) - res%R(:,idx1)
           res%b(:,j) = b
           res%b_len(j) = sqrt(dot_product(b,b))
        end if
     end do
     if (i>1) then
        res%atom_pdb(-1:0) = 1
        do j = -1, 0 ! (-C,N) & (-CA,-C) bnd
           res%b(:,j) = protein%residue(i-1)%b(:,3+j)
           res%b_len(j) = protein%residue(i-1)%b_len(3+j)
        end do
     end if
     protein%residue(i) = res ! save res
  end do
  do i = 1, protein%num_ss_bnd
     ss_res(1:2) = ref_ss(i)%res_no(1:2)
     ss_atm(1:2) = ref_ss(i)%s_atm_no(1:2)
     b = protein%residue(ss_res(2))%R(:,ss_atm(2)) &
          - protein%residue(ss_res(1))%R(:,ss_atm(1))
     protein%ss_bnd(i)%b = b
     protein%ss_bnd(i)%b_len = sqrt(dot_product(b,b))
  end do
end subroutine read_pdb
!-----------------------------------------------------------------------
subroutine write_pdb(file_name)
!-----------------------------------------------------------------------
! This subroutine writes a pdb file.
! ATOM is the first string, followed by the atom serial number, atom name,
! residue name, chain identifier, residue sequence number,
! code for insertions of residues, X,Y,Z coords, occupancy,
! and temperature factors.
! write TER at the end.
!-----------------------------------------------------------------------
  implicit none
  character(len=len_fname) :: file_name
  character(len=7) :: string = 'ATOM', ter_string = 'TER'
  character(len=3) :: resname
  character(len=1) :: chain = 'A', s2, s3
  character(len=4) :: atom_type
  character(len=1), dimension(4) :: s
  integer :: i, j, k, m, unit = 20
  type(reference_residue_type) :: ref
  real(dp), dimension(3) :: R, center
!-----------------------------------------------------------------------
  open(unit, file = trim(file_name))
!-----------------------------------------------------------------------
  ! write ss bnd
  do i = 1, protein%num_ss_bnd 
     write(unit,69)i,ref_ss(i)%res_no(1),ref_ss(i)%res_no(2)
  end do
!-----------------------------------------------------------------------
  ! center of geometry
  center = 0; k = 0
  do i = 1, protein%num_res ! loop over residues
     ref = ref_res(i)      
     do j = 1, ref%num_atm ! loop over atoms
        k = k+1
        center = center + protein%residue(i)%R(:,j)
     end do
  end do
  center = center / k
!-----------------------------------------------------------------------
  ! write atoms
  k = 0
  do i = 1, protein%num_res ! loop over residues
     ref = ref_res(i)      
     resname = ref%residue_type
     do j = 1, ref%num_atm ! loop over atoms
        k = k+1
        R = protein%residue(i)%R(:,j)
        R = R - center ! translate
        atom_type = ref%atom_type(j)
        ! change atom_type in terminal residues
        if (atom_type(1:2)=='HT') atom_type = 'H00'//atom_type(3:3)
        if (atom_type=='OT10') atom_type = 'O000'
        if (atom_type=='OT20') atom_type = 'OXT0'
!-----------------------------------------------------------------------
        do m = 1, 4
           s(m) = atom_type(m:m)
           if (s(m) == '0') s(m) = ' '
        end do
        s2 = s(2); s3 = s(3)
        atom_type = s(1)//s2(1:len(s(2)))//s3(1:len(s(3)))
        write(unit,72) string, k, s(4), atom_type, &
             resname, chain, i, R(1), R(2), R(3)
     end do
  end do
  write(unit, 70) ter_string
  close (unit)
69 format('SSBOND',2x,i2,1x,'CYS',4x,i3,4x,'CYS',4x,i3)
72 format(a7, i4, 1x, a1, a3, 1x, a3, 1x, a1, i4, 4x, 3f8.3)
70 format(a7)
end subroutine write_pdb
!-----------------------------------------------------------------------
subroutine read_user_input()
  implicit none
  integer :: unit = 20, ioerror
  character(len=len_fname) :: string

  open(unit, file = infile_user_input)

  do
     read(unit,"(a20)",iostat=ioerror) string
     if (ioerror < 0 .or. string == 'END') exit 
     if (string(1:11) == 'IN_PDB_NAME') then
        backspace(unit)
        read(unit,*) string, pdb_file_name
     else if (string(1:11) == 'PRINT_LEVEL') then
        backspace(unit)
        read(unit,*) string, print_level
     else if (string(1:17) == 'MINIMIZE_AT_SETUP') then
        backspace(unit)
        read(unit,*) string, minimize_at_setup
     else if (string(1:13) == 'TIME_COHERENT') then
        backspace(unit)
        read(unit,*) string, time_coherent
        if (time_coherent == 'YES') then
           t_coherent = 1
        else
           t_coherent = 0
        end if
     else if (string(1:17) == 'MINIMIZATION_MODE') then
        backspace(unit)
        read(unit,*) string, include_ang_type
!	include_ang_type = 5
     else if (string(1:4) == 'FTOL') then
        backspace(unit)
        read(unit,*) string, factr_eng
     else if (string(1:4) == 'GTOL') then
        backspace(unit)
        read(unit,*) string, pgtol_eng
     else if (string(1:8) == 'MAX_ITER') then
        backspace(unit)
        read(unit,*) string, max_iter_lbfgsb
     end if
  end do

  close(unit)

  if (print_level > 0) then
     write(*,900) 
     write(*,"(2a)")   ' IN_PDB_NAME         = ', trim(pdb_file_name) 
     write(*,"(a,i4)") ' PRINT_LEVEL         = ', print_level
     write(*,"(2a)")   ' MINIMIZE_AT_SETUP   = ', trim(minimize_at_setup)
     write(*,"(2a)")   ' TIME_COHERENT       = ', time_coherent
     write(*,"(a,i4)") ' MINIMIZATION_MODE   = ', include_ang_type
     if (include_ang_type == 0) then
        write(*,*) 'Cartesian minimization'
     else if (include_ang_type == 1) then
        write(*,*) 'All dihedral and bond angles are varied.'
     else if (include_ang_type == 2) then
        write(*,*) 'Only dihedral angles are varied.'
     else if (include_ang_type == 3) then
        write(*,*) 'Only backbone dihedral angles are varied.'
     else if (include_ang_type == 4) then
        write(*,*) 'Only phi and psi angles are varied.'
     else if (include_ang_type == 5) then
        write(*,*) 'Only side chain dihedral angles are varied.'
     else
        write(*,*) 'Unknown angle type.'
        stop
     end if
     write(*,900) 
  end if

900 format(' ---------------------------------------------------')
end subroutine read_user_input
!-----------------------------------------------------------------------
subroutine read_parameters()
!-----------------------------------------------------------------------
! read parameters, except charges
!-----------------------------------------------------------------------
  implicit none
  character(len=16) :: lstring
  character(len=4) :: atom
  character(len=6) :: string
  integer :: unit = 23, i, j, atm_no, ang_no, bnd_no
  real(dp) :: para, b0
  real(dp) :: Emin,Rmin,Emin14,Rmin14,V,Gref,Gfree,Lambda,R0,k,n,phi0
!-----------------------------------------------------------------------
  num_prm = 0 ! initialize num of parameters
  open(unit, file = trim(infile_parameters))
! read ATOM CLASS
  do  
     read(unit,200) lstring
     if (lstring == 'ATOM CLASS')exit
  end do
  atm_no = 0
  do 
     read(unit,200) lstring
     if (lstring == 'END') then
        exit
     else
        atm_no = atm_no + 1; backspace(unit)
        read(unit,*) eng_para%atom_cls(atm_no)
     end if
  end do
  eng_para%num_a_cls = atm_no

! read bond parameters
  do  
     read(unit,200) lstring
     if (lstring == 'BOND')exit
  end do
  bnd_no = 0
  do 
     read(unit,200) lstring
     if (lstring == 'END') then
        exit
     else
        bnd_no = bnd_no + 1; backspace(unit)
        read(unit,*) (eng_para%atm_in_bnd(i,bnd_no),i=1,2),k,b0
        eng_para%bnd_para(1:2,bnd_no)=(/ k, b0 /)
        num_prm = num_prm + 2
     end if
  end do
  eng_para%num_bnd = bnd_no

! read bond angle parameters
  do  
     read(unit,200) lstring
     if (lstring == 'BOND ANGLE')exit
  end do
  ang_no = 0
  do 
     read(unit,200) lstring
     if (lstring == 'END') then
        exit
     else
        ang_no = ang_no + 1; backspace(unit)
        eng_para%ang_type(ang_no) = 2
        read(unit,*) (eng_para%atm_in_ang(i,ang_no),i=1,3),k,phi0
        eng_para%ang_para(1:2,ang_no)=(/ k, deg2rad*phi0 /)
        num_prm = num_prm + 2
     end if
  end do

! read dihedral angle parameters
  do  
     read(unit,200) lstring
     if (lstring == 'DIHEDRAL ANGLE')exit
  end do
  do 
     read(unit,200) lstring
     if (lstring == 'END') then
        exit
     else
        ang_no = ang_no + 1; backspace(unit)
        eng_para%ang_type(ang_no) = 1
        read(unit,*) (eng_para%atm_in_ang(i,ang_no),i=1,4),k,n,phi0
        eng_para%ang_para(1:3,ang_no)=(/ k, deg2rad*phi0, n /)
        num_prm = num_prm + 2
     end if
  end do

! read improper torsion angle parameters
  do  
     read(unit,200) lstring
     if (lstring == 'IMPROPER TROSION')exit
  end do
  do 
     read(unit,200) lstring
     if (lstring == 'END') then
        exit
     else
        ang_no = ang_no + 1; backspace(unit)
        eng_para%ang_type(ang_no) = 3
        read(unit,*) (eng_para%atm_in_ang(i,ang_no),i=1,4), &
             k,para,phi0
        eng_para%ang_para(1:2,ang_no)=(/ k, deg2rad*phi0 /)
        num_prm = num_prm + 2
     end if
  end do
  eng_para%num_ang = ang_no

! read LJ parameters
  do  
     read(unit,200) lstring
     if (lstring == 'VDW') exit
  end do
  atm_no = 0
  do 
     read(unit,200) lstring
     if (lstring == 'END') then
        exit
     else if (lstring(1:6)=='E14FAC') then
        backspace(unit); read(unit,*)string, eng_para%E14fac
        num_prm = num_prm + 1
     else
        backspace(unit); read(unit,*) atom 
        if (atom /= '!') then
           do i = 1, eng_para%num_a_cls
              if (atom == eng_para%atom_cls(i)) then
                 atm_no = i; exit
              end if
           end do
           if (atm_no == 0) then
              write(*,*) 'Error in finding atm class of ', atom 
           end if
           backspace(unit)
           if (atom(1:1) == 'C') then
              read(unit,*)atom,para,Emin,Rmin,para,Emin14,Rmin14
           else
              read(unit,*)atom,para,Emin,Rmin
              Emin14 = Emin; Rmin14 = Rmin
           end if
           eng_para%LJ_para(1:4,atm_no)=(/ Emin,Rmin,Emin14,Rmin14 /)
           num_prm = num_prm + 4
        end if
     end if
  end do

! calc sigma_ij and eps_ij for LJ intrxn using the combination rule
  do i = 1, eng_para%num_a_cls
     do j = 1, 4
        call get_aux_para('LJ',i,j)
     end do
  end do
  
! read eef1 solvation parameters
  do  
     read(unit,200) lstring
     if (lstring == 'SOLVATION') exit
  end do
  atm_no = 0
  do 
     read(unit,200) lstring
     if (lstring == 'END') then
        exit
     else
        backspace(unit); read(unit,*) atom 
        if (atom /= '!') then
           do i = 1, eng_para%num_a_cls
              if (atom == eng_para%atom_cls(i)) then
                 atm_no = i; exit
              end if
           end do
           if (atm_no == 0) then
              write(*,*) 'Error in finding atm class of ', atom 
           end if
           backspace(unit)
           read(unit,*)atom,V,Gref,Gfree,para,para,Lambda
           R0 = eng_para%LJ_para(2,atm_no)
           eng_para%slv_para(1:5,atm_no)=(/ V,Gref,Gfree,1.0/Lambda,R0 /)
           num_prm = num_prm + 5
        end if
     end if
  end do

! calc slv_factor = tmp*Gfree_i/Lambda_i*V_j
  do i = 1, eng_para%num_a_cls
     do j = 1, 5
        call get_aux_para('SL',i,j)
     end do
  end do
  close(unit)

200 format(a16)
!-----------------------------------------------------------------------
end subroutine read_parameters
!-----------------------------------------------------------------------
subroutine get_aux_para(para_type,i,j)
!-----------------------------------------------------------------------
  implicit none
  real(dp) :: para, sigma, tmp, tmp1, Gfree, Lambda, V
  integer :: i, ii, j
  character(len=2) :: para_type

  if (para_type == 'LJ') then
     para = eng_para%LJ_para(j,i)
     if (j==1.or.j==3) then
        eng_para%aux_para(j,i,i) = para
     else if (j==2.or.j==4) then
        eng_para%aux_para(j,i,i) = 4.0d0*para*para
     end if
     do ii = 1, eng_para%num_a_cls
        if (j==1.or.j==3) then
           eng_para%aux_para(j,i,ii)=-sqrt(para*eng_para%LJ_para(j,ii))
        else if (j==2.or.j==4) then
           sigma = para + eng_para%LJ_para(j,ii)
           eng_para%aux_para(j,i,ii) = sigma*sigma 
        end if
        eng_para%aux_para(j,ii,i) = eng_para%aux_para(j,i,ii)
     end do
  else if (para_type == 'SL') then
     tmp = 0.5d0/(pi*sqrt(pi))
     if (j==3.or.j==4) then
        Gfree=eng_para%slv_para(3,i)
        Lambda=eng_para%slv_para(4,i)
        do ii = 1, eng_para%num_a_cls
           V=eng_para%slv_para(1,ii)
           eng_para%aux_para(5,i,ii) = tmp*Gfree*Lambda*V
        end do
     else if (j==1) then
        V=eng_para%slv_para(1,i)
        do ii = 1, eng_para%num_a_cls
           Gfree=eng_para%slv_para(3,ii)
           Lambda=eng_para%slv_para(4,ii)
           eng_para%aux_para(5,ii,i) = tmp*Gfree*Lambda*V
        end do
     end if
  end if
end subroutine get_aux_para
!-----------------------------------------------------------------------
subroutine read_topology_geo()
!-----------------------------------------------------------------------
! builds the reference residue structure.  
!-----------------------------------------------------------------------
  implicit none
  type(reference_residue_type) :: ref
  character(len=4) :: string, resname
  integer :: unit = 19, res_no, atm_no, bnd_no, ang_no, dpn_no, i
  integer :: index, imax, ang_type
  integer, dimension(3) :: idx
  character(len=4), dimension(4) :: atm
  integer, dimension(0:max_ang) :: dpn_ang_no
!-----------------------------------------------------------------------
  open(unit, file = trim(infile_topo_geo))
  res_no = 0
  do  ! exit at "END" or "RES NTER"
     read(unit,*) string; backspace(unit)
!-----------------------------------------------------------------------
     if (string == 'END') then   ! complete structure of last residue
        ref%num_ang = ang_no; ref%num_dpn_ang = dpn_ang_no
        residue_ref(res_no) = ref
        exit
     end if
!-----------------------------------------------------------------------
! main: actions to set up various residues                       RES CYS  
     if (string == 'RES') then
        if (res_no > 0) then  ! previous residue done, complete structure
           ref%num_ang = ang_no; ref%num_dpn_ang = dpn_ang_no
           residue_ref(res_no) = ref
        end if
        read(unit,*) string, resname
        if (resname == 'NTER') exit ! terminal res taken care of elsewhere
        call find_res(resname, res_no)
        ref%residue_type = resname
        atm_no = 0; bnd_no = 0; ang_no = 0
        ref%atom_type = ' '; ref%b_idx = 0
        ref%atom_type(-1:0) = (/ '-CA0', '-C00' /)
        ref%b_idx(1:2,-1)= (/-1, 0 /) ! (-CA,-C) and 
        ref%b_idx(1:2,0) = (/ 0, 1 /) ! (-C,N) bnd
        ref%na_ang(0)=3;ref%ang_type(0)=2;ref%a_idx(1:2,0)=(/-1,0/) !(-CA,-C,N)
!-----------------------------------------------------------------------
! read atom: name                 ATM N 
     else if (string == 'ATM') then
        atm_no = atm_no + 1
        read(unit,*) string, ref%atom_type(atm_no)
        call fill_atom_type4(ref%residue_type, ref%atom_type(atm_no)) !replace blanks by 0's
!-----------------------------------------------------------------------
! Bonds: 1st, 2nd, length:  BND N    CA               1.4500000000000000
     else if (string == 'BND') then
        if (bnd_no == 0) then
           ref%num_atm = atm_no - 1
        end if
        bnd_no = bnd_no + 1
        read(unit,*) string, (atm(i),i=1,2), ref%b_len(bnd_no)
        do i = 1, 2
           call find_atm(ref,atm(i),idx(i))
        end do
        ref%b_idx(1:2,bnd_no) = (/ idx(1), idx(2) /)
!-----------------------------------------------------------------------
     else if (string == 'DIH' .or. string == 'ANG') then
        if (ang_no == 0) then  !just done reading bonds, start w. angles
           ref%num_bnd = bnd_no
           dpn_ang_no = 0
        end if
        ang_no = ang_no + 1                ! increment angle counter
        if (string == 'DIH') then
           imax = 4; ref%ang_type(ang_no) = 1 
        else if (string == 'ANG') then
           imax = 3; ref%ang_type(ang_no) = 2
        end if
        ref%na_ang(ang_no)=imax
        read(unit, *) string, (atm(i), i = 1, imax), ref%w0(ang_no)
        ref%w0(ang_no) = ref%w0(ang_no)*deg2rad 
        do i = 1, imax-1 ! find ref%a_idx
           call find_bnd(ref,atm(i),atm(i+1),idx(i),0)
        end do
        ref%a_idx(1:imax-1,ang_no) = idx(1:imax-1)
!-----------------------------------------------------------------------
        dpn_no = 0 ! get dpn bonds
        do 
           read(unit,*) string; backspace(unit)
           if (string /= 'DPN') exit
           dpn_no = dpn_no + 1
           read(unit,*) string, atm(1), atm(2)
           call find_bnd(ref,atm(1),atm(2),index,0)
           ref%i_dpn_bnd(dpn_no,ang_no)= index
        end do
        ref%num_dpn_bnd(ang_no) = dpn_no
!-----------------------------------------------------------------------
! get prn angles and assign dpn angles reversely
        read(unit,*) string; backspace(unit)
        if (string /= 'PRN') then
           write(*,*)'Error in ', trim(infile_topo_geo), ' expecting PRN'
           stop
        end if
        read(unit,*) string, string; backspace(unit)
        if (string == 'DIH') then
           imax = 4; ang_type = 1
        else if (string == 'ANG') then
           imax = 3; ang_type = 2
        end if
        read(unit,*) string, string, (atm(i), i = 1, imax)
        call find_ang(ref,ang_type,ang_no,atm(1:4),index,0)
        dpn_ang_no(index) = dpn_ang_no(index) + 1
        ref%i_dpn_ang(dpn_ang_no(index),index) = ang_no
     else
        read(unit,*)
     end if
  end do
  close(unit)
!-----------------------------------------------------------------------
end subroutine read_topology_geo
!-----------------------------------------------------------------------
subroutine setup_ter_ref(ter_type)
!-----------------------------------------------------------------------
! set up Nter_ref or Cter_ref depending on the sequence
!-----------------------------------------------------------------------
  implicit none
  character(len=1), intent(in) :: ter_type
  type(reference_residue_type) :: ref
  type(ref_res_eng_type) :: eref
  character(len=4) :: resname,del_atm,ter_name,string,atom_type
  character(len=4), dimension(3) :: add_atm
  character(len=4), dimension(4) :: atm, atom
  character(len=4), dimension(2,-1:max_bnd) :: atm_in_bnd, atom_in_bnd
  character(len=4), dimension(2,max_ang) :: atm_in_dpn, atom_in_dpn
  character(len=4), dimension(4,0:max_ang) :: atm_in_ang, atom_in_ang
  character(len=4), dimension(4,0:max_ang) :: atm_in_prn, atom_in_prn
  character(len=4), dimension(max_atm+1) :: atom_cls, atm_cls
  character(len=4), dimension(-1:max_atm+1) :: atm_type
  real(dp), dimension(max_bnd) :: b_len  
  real(dp), dimension(max_ang) :: w0
  integer, dimension(max_atm+1) :: delete_atm
  integer, dimension(max_bnd) :: delete_bnd
  integer, dimension(0:max_ang) :: angle_type, delete_ang, na_ang
  integer :: num_del_atm,num_add_atm,num_del_bnd,num_add_bnd
  integer :: num_del_ang,num_add_ang,del_atm_idx
  integer :: i,j,k,idx,imax,ang_type,unit=23,del,prm_no
  integer :: res_type,res_no,num_atm,num_bnd,num_ang,del_ang1
  integer, dimension(3) :: indx
  integer, dimension(0:max_ang) :: dpn_ang_no
!-----------------------------------------------------------------------
! info needed for terminal residue
  if (ter_type == 'N') res_no = 1
  if (ter_type == 'C') res_no = protein%num_res
  res_type = protein%residue(res_no)%res_type
  ref = residue_ref(res_type); eref = res_ref_eng(res_type)
  resname = ref%residue_type
  if (ter_type == 'N') then
     if (resname == 'GLY') then
        ter_name='GLYP'; del_atm='H000'; num_del_atm=1; num_add_atm=3
     else if (resname == 'PRO') then ! CD is not actually  deleted
        ter_name='PROP'; del_atm='CD00'; num_del_atm=0; num_add_atm=2
     else 
        ter_name='NTER'; del_atm='H000'; num_del_atm=1; num_add_atm=3
     end if
     add_atm = (/ 'HT10','HT20','HT30' /)
  else if (ter_type == 'C') then
     ter_name='CTER'; del_atm='O000'; num_del_atm=1; num_add_atm=2
     add_atm(1:2) = (/ 'OT10','OT20' /)
  else
     write(*,*)'Unknown terminal type'; stop
  end if
!-----------------------------------------------------------------------
! save data before changing indices
  num_atm = ref%num_atm; num_bnd = ref%num_bnd; num_ang = ref%num_ang
  atm_type=ref%atom_type; b_len=ref%b_len; angle_type=ref%ang_type
  na_ang = ref%na_ang; w0 = ref%w0
  atm_cls(1:num_atm) = eng_para%atom_cls(eref%atm_cls(1:num_atm))
  do i = -1, num_bnd
     atm_in_bnd(1:2,i) = ref%atom_type(ref%b_idx(1:2,i))
  end do
  do i = 0, num_ang
     imax = ref%na_ang(i)
     atm_in_ang(1,i) = atm_in_bnd(1,ref%a_idx(1,i))
     atm_in_ang(2:imax,i)=atm_in_bnd(2,ref%a_idx(1:(imax-1),i))
     if (i>0) then
        if (ref%num_dpn_bnd(i)==1) then
           atm_in_dpn(1:2,i)=atm_in_bnd(1:2,ref%i_dpn_bnd(1,i))
        end if
     end if
     do j = 1, ref%num_dpn_ang(i)
        atm_in_prn(1:imax,ref%i_dpn_ang(j,i)) = atm_in_ang(1:imax,i)
     end do
  end do
!-----------------------------------------------------------------------
  open(unit, file = infile_topo_geo)
  do ! find the start of terminal residue type
     read(unit,*) string
     if (string == 'RES') then
        backspace(unit); read(unit,*) string, string
        if (string == ter_name) exit
     end if
  end do
!-----------------------------------------------------------------------
! change atom data
  call find_atm(ref,del_atm,del_atm_idx) ! find del_atm no
  delete_atm = 0 ; del = 0
  if (num_del_atm/=0) delete_atm(del_atm_idx) = 1
  idx = 0
  do i = 1, num_atm+1 ! add and delete atom
     if (delete_atm(i)==1.and.del==0.or. &
          (num_del_atm==0.and.i==del_atm_idx)) then
        del = 1
        do j = 1, num_add_atm
           idx = idx + 1
           ref%atom_type(idx)=add_atm(j)
        end do
     end if
     if (delete_atm(i)==0) then
        idx = idx + 1
        ref%atom_type(idx)=atm_type(i)
        atom_cls(idx)=atm_cls(i)
     end if
  end do
  ref%num_atm = num_atm + num_add_atm - num_del_atm ! increase num_atm
  do ! change atom_cls
     read(unit,*) string, atom_type; backspace(unit)
     if (string /= 'ATM') exit
     call find_atm(ref,atom_type,idx)
     read(unit,*) string,atom_type,atom_cls(idx)
  end do
!  print*,'done with atoms in ', ref%residue_type,' ', ter_type, 'ter'
!  do i = -1, ref%num_atm+1
!     print*,ref%atom_type(i)
!  end do
!-----------------------------------------------------------------------
! change bond data
  delete_bnd = 0; num_del_bnd = 0
  do i = 1, num_bnd ! find bnd to be deleted 
     if ((num_del_atm/=0.and. (ref%b_idx(1,i)==del_atm_idx.or. &
          ref%b_idx(2,i)==del_atm_idx)).or. &
          (ter_type == 'C'.and.atm_in_bnd(2,i)=='+N00'))then
        delete_bnd(i) = 1; num_del_bnd = num_del_bnd + 1
     end if
  end do
  idx = 0; del = 0; num_add_bnd = 0
  do i = 1, num_bnd ! add and delete bond
     if((delete_bnd(i)==1.and.del==0).or.(num_del_atm==0.and.i==4))then
        ! the order of the first three bnd should not be disturbed 
         del = 1
        do ! add bnds
           read(unit,*) string; backspace(unit)
           if (string /= 'BND') exit
           idx = idx + 1; num_add_bnd = num_add_bnd + 1
           read(unit,*)string,(atom_in_bnd(j,idx),j=1,2),ref%b_len(idx)
        end do
     end if
     if (delete_bnd(i)==0) then
        idx = idx + 1
        atom_in_bnd(:,idx) = atm_in_bnd(:,i)
        ref%b_len(idx) = b_len(i)
     end if
  end do
  if (ter_type == 'N') then
     call find_atm(ref,'HT10',idx)
     ref%b_idx(1:2,0) = (/ 0, idx /) ! (-C,HT1) bnd
  end if
  ref%num_bnd = num_bnd + num_add_bnd - num_del_bnd ! increase num_bnd
  do i = 1, ref%num_bnd
     call find_atm(ref,atom_in_bnd(1,i),indx(1))
     call find_atm(ref,atom_in_bnd(2,i),indx(2))
     ref%b_idx(1:2,i) = indx(1:2)
  end do
!  print*,'done with bonds in ', ref%residue_type,' ', ter_type, 'ter'
!  do i = -1, ref%num_bnd
!     print*,ref%atom_type(ref%b_idx(1,i)),ref%atom_type(ref%b_idx(2,i))
!  end do
!-----------------------------------------------------------------------
! change angle data
  delete_ang = 0; num_del_ang = 0
  do i = 1, num_ang ! find ang to be deleted
     imax = na_ang(i)
     if((num_del_atm/=0.and.atm_in_ang(imax,i)==del_atm) .or. &
          (ter_type == 'C'.and.atm_in_ang(imax,i)=='+N00'))then
        delete_ang(i) = 1; num_del_ang = num_del_ang + 1
!        print*,num_del_ang,(atm_in_ang(j,i),j=1,imax)
     end if
  end do
  idx = 0; del = 0; num_add_ang = 0
  do i = 1, num_ang  ! add and delete angle
     if((delete_ang(i)==1.and.del==0).or.(num_del_atm==0.and.i==1))then
        del = 1; del_ang1 = i
        do ! add angles
           read(unit,*) string 
           if (string == 'END' .or. string == 'RES') exit
           if (string == 'DIH' .or. string == 'ANG') then
              backspace(unit)
              idx = idx + 1; num_add_ang = num_add_ang + 1
              if (string == 'DIH') then
                 imax = 4; ref%ang_type(idx) = 1
              else if (string == 'ANG') then
                 imax = 3; ref%ang_type(idx) = 2
              end if
              ref%na_ang(idx)=imax
              read(unit, *)string,(atom_in_ang(j,idx),j=1,imax),ref%w0(idx)
              ref%w0(idx) = ref%w0(idx)*deg2rad 
              call fill_atom_type4(resname, atom_in_ang(imax,idx))
              do j = i, num_ang ! delete ang to be redefined
                 if((angle_type(j)==1.and.imax==4).or.&
                      (angle_type(j)==2.and.imax==3)) then
                    if (atm_in_ang(imax,j) == atom_in_ang(imax,idx)) then
                       delete_ang(j) = 1; num_del_ang = num_del_ang + 1
!                       print*,num_del_ang,(atm_in_ang(k,j),k=1,imax)
                    end if
                 end if
              end do
               do 
                 read(unit,*) string; backspace(unit)
                 if (string /= 'DPN') exit
                 read(unit,*) string, (atom_in_dpn(j,idx),j=1,2)
              end do
              read(unit,*) string; backspace(unit)
              if (string /= 'PRN') then
                 write(*,*)'Error in topology_geo.in: expecting PRN'
                 stop
              end if
              read(unit,*) string, string; backspace(unit)
              if (string=='DIH')imax=4; if(string=='ANG')imax=3
              read(unit,*) string, string, (atom_in_prn(j,idx),j=1,imax)
           end if
        end do
        close(unit)
     end if
     if (delete_ang(i)==0) then
        idx = idx + 1
        ref%ang_type(idx)=angle_type(i)
        ref%na_ang(idx)=na_ang(i)
        atom_in_ang(:,idx)=atm_in_ang(:,i)
        atom_in_dpn(:,idx)=atm_in_dpn(:,i)
        atom_in_prn(:,idx)=atm_in_prn(:,i)
        ref%w0(idx)=w0(i)
     end if
  end do
  ref%num_ang = num_ang+num_add_ang-num_del_ang ! increase num_ang
  dpn_ang_no = 0
  do i = 1, ref%num_ang ! find indices: a_idx, i_dpn_bnd, i_dpn_ang
     ang_type = ref%ang_type(i)
     imax=ref%na_ang(i)
     do j = 1, imax-1 
        call find_bnd(ref,atom_in_ang(j,i),atom_in_ang(j+1,i),ref%a_idx(j,i),0)
     end do
     if (ang_type==2) then
        call find_bnd(ref,atom_in_dpn(1,i),atom_in_dpn(2,i),idx,0)
        ref%i_dpn_bnd(1,i)=idx; ref%num_dpn_bnd(i)=1
     else
        ref%num_dpn_bnd(i)=0
     end if
     if (ang_type==1) then
        ang_type=2
     else 
        ang_type=1
     end if
     call find_ang(ref,ang_type,i,atom_in_prn(1:4,i),idx,0)
     dpn_ang_no(idx) = dpn_ang_no(idx) + 1
     ref%i_dpn_ang(dpn_ang_no(idx),idx) = i
  end do
  ref%num_dpn_ang = dpn_ang_no
  ! get equil values for new angs 
  do i = del_ang1, num_add_ang+del_ang1-1
     atom_type = atom_in_ang(1,i)
     if (ref%ang_type(i)==2.and.atom_type(1:1)/='-'.and. &
          abs(ref%w0(i))<1.0d-10)then
        atm(1:3) = atom_in_ang(1:3,i)
        do j = 1, 3
           call find_atm(ref, atm(j), idx)
           atm(j) = atom_cls(idx)
        end do
        call find_ang_prm(atm,2,prm_no)
        ref%w0(i) = eng_para%ang_para(2,prm_no)
     end if
  end do
!  print*,'done with angles in ', ref%residue_type,' ', ter_type, 'ter'
!  do i = 1, ref%num_ang
!     if(ref%ang_type(i)==1)imax=4; if(ref%ang_type(i)==2)imax=3
!     print*,i,(atom_in_ang(j,i),j=1,imax),real(rad2deg*ref%w0(i)), &
!          ref%i_dpn_ang(1:ref%num_dpn_ang(i),i)
!  end do
!-----------------------------------------------------------------------
  if (ter_type == 'N') Nter_ref = ref ! save as Nter_ref
  if (ter_type == 'C') Cter_ref = ref ! save as Cter_ref
end subroutine setup_ter_ref
!-----------------------------------------------------------------------
subroutine read_topology_eng()
!-----------------------------------------------------------------------
! build the reference residue structure for energy calculation. 
!-----------------------------------------------------------------------
  implicit none
  type(reference_residue_type) :: ref
  type(ref_res_eng_type) :: eref
  character(len=4) :: string, resname, atom_type
  integer :: unit=19,res_no,atm_no,ang_no,i,idx,atm_typ_no,imax,ang_type
  integer :: gr_no, indx, prm_no
  character(len=4), dimension(4) :: atm, atom
!-----------------------------------------------------------------------
  open(unit, file = trim(infile_topo_eng))
  res_no = 0; atm_typ_no = 0 
  do  ! exit at "end" or "RES NTER"
     read(unit,*) string; backspace(unit)
!-----------------------------------------------------------------------
     if (string == 'END') then   ! complete structure of last residue
        eref%num_ang = ang_no; res_ref_eng(res_no) = eref
        eng_para%num_q_typ = atm_typ_no
        exit
     end if
!-----------------------------------------------------------------------
! main: actions to set up various residues                       RES CYS  
     if (string == 'RES') then
        if (res_no > 0) then  ! previous residue done, complete structure
           eref%num_ang = ang_no; res_ref_eng(res_no) = eref
           eng_para%num_q_typ = atm_typ_no
        end if
        read(unit,*) string, resname
        if (resname == 'NTER') exit
        call find_res(resname, res_no)
        ref = residue_ref(res_no); eref%num_atm = ref%num_atm
        atm_no = 0; ang_no = 0; eref%na_gr = 0
!         print*,resname
!-----------------------------------------------------------------------
! read atom: name   ATM N     1   NH1   -0.35000    1    ENDC    0
     else if (string == 'ATM') then
        atm_no=atm_no+1; atm_typ_no=atm_typ_no+1; num_prm=num_prm+1
        read(unit,*) string, atom_type; backspace(unit)
        call find_atm(ref,atom_type,idx)
        read(unit,*) string,atom_type,gr_no,atom_type, &
             eng_para%charge(atm_typ_no),indx, &
             new_numbering(res_no)%endc(indx), &
             new_numbering(res_no)%ring(indx)
        eng_para%q_type(1:2,atm_typ_no)=(/resname,ref%atom_type(atm_no)/)
        new_numbering(res_no)%atom(indx)=ref%atom_type(idx)
        call find_atm_cls(atom_type,eref%atm_cls(idx))
        eref%q_idx(idx) = atm_typ_no; eref%gr_no(idx) = gr_no
        eref%na_gr(gr_no) = eref%na_gr(gr_no) + 1 ! no of atoms in gr
        if (idx == 3) eref%num_gr = gr_no ! C is in the last group
!         write(*,'(i2,1x,a4,1x,a4)')&
!              atm_no,atom_type,eng_para%atom_cls(eref%atm_cls(idx))
!-----------------------------------------------------------------------
! angles (dih, imp, or bnd ang)
     else if (string=='DIH'.or.string=='IMP'.or.string=='ANG') then
        ang_no = ang_no + 1
        if (string=='DIH') then
           imax=4; ang_type=1
        else if (string=='IMP') then
           imax=4; ang_type=3
        else if (string=='ANG') then
           imax=3; ang_type=2
        end if
        eref%ang_type(ang_no)=ang_type; eref%na_ang(ang_no)=imax
        read(unit,*) string, (atm(i),i=1,imax); atom=atm
        prm_no = -1
        do i = 1, imax
           call find_atm(ref,atm(i),idx)
           eref%aa_idx(i,ang_no) = idx
! ang involving atoms in other res: dep on sequence, should be done separately
           if (idx<1 .or. idx>ref%num_atm) then
              prm_no = 0 
           else
              atm(i) = eng_para%atom_cls(eref%atm_cls(idx))
           end if
        end do
        if (prm_no /=0) call find_ang_prm(atm,ang_type,prm_no)
        eref%prm(ang_no) = prm_no
!         print*,resname,ang_type,atm(1:imax),atom(1:imax)
!-----------------------------------------------------------------------
     else
        read(unit,*)
     end if
  end do
  close(unit)
!-----------------------------------------------------------------------
end subroutine read_topology_eng
!-----------------------------------------------------------------------
subroutine setup_eng_para()
!-----------------------------------------------------------------------
  implicit none
  type(reference_residue_type) :: ref, next_ref, prev_ref
  type(ref_res_eng_type) :: eref, prev_eref, next_eref
  integer :: num_res,res_no,ang_no,bnd_no,idx,i,ang_type,imax, res_type
  character(len=4) :: string,resname,prev_atm,next_atm
  character(len=4), dimension(4) :: atm,atom
!-----------------------------------------------------------------------
  num_res=protein%num_res
  do res_no = 2, num_res-1
     res_type = protein%residue(res_no)%res_type
     ref=residue_ref(res_type)
     eref=res_ref_eng(res_type)
     resname = ref%residue_type
     prev_ref=residue_ref(protein%residue(res_no-1)%res_type)
     next_ref=residue_ref(protein%residue(res_no+1)%res_type)
     prev_eref=res_ref_eng(protein%residue(res_no-1)%res_type)
     next_eref=res_ref_eng(protein%residue(res_no+1)%res_type)
!-----------------------------------------------------------------------
!     print*,resname
     do ang_no = 1, eref%num_ang
        ang_type = eref%ang_type(ang_no)
        imax = eref%na_ang(ang_no)
        atm(1:imax) = ref%atom_type(eref%aa_idx(1:imax,ang_no));atom=atm
        if (eref%prm(ang_no) == 0) then
           do i = 1, imax
              call find_atm(ref,atm(i),idx)
              if (idx<1) then
                 prev_atm=atm(i); prev_atm=prev_atm(2:)
                 call find_atm(prev_ref,prev_atm,idx)
                 atm(i) = eng_para%atom_cls(prev_eref%atm_cls(idx))
              else if (idx>ref%num_atm) then
                 next_atm=atm(i); next_atm=next_atm(2:)
                 call find_atm(next_ref,next_atm,idx)
                 atm(i) = eng_para%atom_cls(next_eref%atm_cls(idx))
              else
                 atm(i) = eng_para%atom_cls(eref%atm_cls(idx))
              end if
           end do
           call find_ang_prm(atm,ang_type,ref_res_eng(res_no)%prm(ang_no))
        else
           ref_res_eng(res_no)%prm(ang_no) = eref%prm(ang_no)
        end if
!        i = res%prm(ang_no)
!        if (ang_type==1) then
!           write(*,'(i2,1x,4(a4,1x),i3,1x,4(a4,1x),3f8.2)')ang_type,atom(1:4),&
!                i,eng_para%atm_in_ang(1:4,i),eng_para%ang_para(1:3,i)
!        else if (ang_type==3) then
!           write(*,'(i2,1x,4(a4,1x),i3,1x,4(a4,1x),2f8.2)')ang_type,atom(1:4),&
!                i,eng_para%atm_in_ang(1:4,i),eng_para%ang_para(1:2,i)
!        else if (ang_type==2) then
!           write(*,'(i2,1x,3(a4,1x),i3,1x,3(a4,1x),2f8.2)')ang_type,atom(1:3),&
!                i,eng_para%atm_in_ang(1:3,i),eng_para%ang_para(1:2,i)
!        end if
     end do
  end do
  do res_no = 1, num_res
     ref=ref_res(res_no)
     eref=ref_res_eng(res_no)
     if (res_no < num_res) then
        next_ref=ref_res(res_no+1)
        next_eref=ref_res_eng(res_no+1)
     end if
     do bnd_no = 1, ref%num_bnd
        atm(1:2) = ref%atom_type(ref%b_idx(1:2,bnd_no));atom=atm
        do i = 1, 2
           call find_atm(ref,atm(i),idx)
           if (idx>ref%num_atm) then
              next_atm=atm(i); next_atm=next_atm(2:)
              call find_atm(next_ref,next_atm,idx)
              atm(i) = eng_para%atom_cls(next_eref%atm_cls(idx))
           else
              atm(i) = eng_para%atom_cls(eref%atm_cls(idx))
           end if
        end do
        call find_bnd_prm(atm,idx)
        ref_res_eng(res_no)%bnd_prm(bnd_no)=idx
!        write(*,'(4(a4,1x),2f8.2)')atom(1:2),&
!             eng_para%atm_in_bnd(1:2,idx),eng_para%bnd_para(2,idx),ref%b_len(bnd_no)
        if (abs(eng_para%bnd_para(2,idx)-ref%b_len(bnd_no))>1.0d-8) then
!           write(*,*) 'Error in assigning equilibrium bond lengths.'
!           write(*,*) res_no,ref%residue_type,bnd_no,atom(1:2), &
!                eng_para%bnd_para(2,idx),ref%b_len(bnd_no)
           ref%b_len(bnd_no) = eng_para%bnd_para(2,idx)
!           stop
        end if
     end do
  end do
!-----------------------------------------------------------------------
! terminal residues
  do res_no = 1, num_res, num_res-1
     if (res_no==1) then
        eref = Nter_ref_eng
        ref = Nter_ref
     else if (res_no==num_res) then
        eref = Cter_ref_eng
        ref = Cter_ref
     end if
!     print*,'ter',res_no
     do ang_no = 1, eref%num_ang
        ref_res_eng(res_no)%prm(ang_no)=eref%prm(ang_no)
!        i = eref%prm(ang_no)
!        ang_type = eref%ang_type(ang_no)
!        if (ang_type==1) then
!           write(*,'(i2,1x,4(a4,1x),i3,1x,4(a4,1x),3f8.2)')ang_type,&
!                ref%atom_type(eref%aa_idx(1:4,ang_no)), &
!                i,eng_para%atm_in_ang(1:4,i),eng_para%ang_para(1:3,i)
!        else if (ang_type==3) then
!           write(*,'(i2,1x,4(a4,1x),i3,1x,4(a4,1x),2f8.2)')ang_type,&
!                ref%atom_type(eref%aa_idx(1:4,ang_no)), &
!                i,eng_para%atm_in_ang(1:4,i),eng_para%ang_para(1:2,i)
!        else if (ang_type==2) then
!           write(*,'(i2,1x,3(a4,1x),i3,1x,3(a4,1x),2f8.2)')ang_type,&
!                ref%atom_type(eref%aa_idx(1:3,ang_no)), &
!                i,eng_para%atm_in_ang(1:3,i),eng_para%ang_para(1:2,i)
!        end if
     end do
  end do
end subroutine setup_eng_para
!-----------------------------------------------------------------------
subroutine setup_ter_ref_eng(ter_type)
!-----------------------------------------------------------------------
! set up Nter_ref_eng or Cter_ref_eng depending on the sequence
!-----------------------------------------------------------------------
  implicit none
  character(len=1), intent(in) :: ter_type
  type(reference_residue_type) :: ref, ter_ref, next_ref, prev_ref
  type(ref_res_eng_type) :: eref, prev_eref, next_eref
  character(len=4) :: resname, del_atm, ter_name, string, atom_type
  character(len=4) :: prev_atm, next_atm
  character(len=4), dimension(4) :: atm, atom
  character(len=4), dimension(4,max_ang_eng) :: atm_in_ang, atom_in_ang
  integer, dimension(max_atm) :: atm_cls, q_idx, group_no
  integer, dimension(max_ang_eng) :: angle_type,na_ang,delete_ang
  integer :: num_del_atm,num_add_atm,del_atm_idx,ang_no,ang_type,gr_no,na_gr
  integer :: num_del_ang,num_add_ang,imax,prev_res_type,next_res_type
  integer :: i,j,idx,unit=23, atm_typ_no, res_type, del

  if (ter_type == 'N') then
     res_type = protein%residue(1)%res_type 
     ter_ref = Nter_ref
     next_res_type = protein%residue(2)%res_type 
     next_ref=residue_ref(next_res_type)
     next_eref=res_ref_eng(next_res_type)
  end if
  if (ter_type == 'C') then
     res_type = protein%residue(protein%num_res)%res_type 
     ter_ref = Cter_ref
     prev_res_type=protein%residue(protein%num_res-1)%res_type
     prev_ref = residue_ref(prev_res_type)
     prev_eref = res_ref_eng(prev_res_type)
  end if
  ref = residue_ref(res_type); eref = res_ref_eng(res_type)
  resname = ref%residue_type
!-----------------------------------------------------------------------
! info needed for terminal residue
  if (ter_type == 'N') then
     if (resname == 'GLY') then
        del_atm = 'H000'; ter_name = 'GLYP'; num_del_atm=1; num_add_atm=3
     else if (resname == 'PRO') then
        del_atm = 'CD00'; ter_name = 'PROP' ! CD is not actually  deleted
        num_del_atm=0; num_add_atm=2
     else 
        del_atm = 'H000'; ter_name = 'NTER'; num_del_atm=1; num_add_atm=3
     end if
  else if (ter_type == 'C') then
     del_atm = 'O000'; ter_name = 'CTER'; num_del_atm=1; num_add_atm=2
  else
     write(*,*)'Unknown terminal type'; stop
  end if
!-----------------------------------------------------------------------
! save data before changing indices
  q_idx = eref%q_idx; atm_cls = eref%atm_cls; group_no = eref%gr_no
  angle_type=eref%ang_type; na_ang=eref%na_ang
  do i = 1, eref%num_ang
     imax = na_ang(i)
     atm_in_ang(1:imax,i) = ref%atom_type(eref%aa_idx(1:imax,i))
  end do
!-----------------------------------------------------------------------
  open(unit, file = trim(infile_topo_eng))
  do ! find the start of terminal residue type
     read(unit,*) string
     if (string == 'RES') then
        backspace(unit); read(unit,*) string, string
        if (string == ter_name) exit
     end if
  end do
!-----------------------------------------------------------------------
! change atom data
  if (ter_type == 'N') gr_no = eref%gr_no(1) ! group no for N
  if (ter_type == 'C') gr_no = eref%gr_no(3) ! group no for C
  call find_atm(ref,del_atm,del_atm_idx) ! find del_atm indx
  do i = 1, ref%num_atm  ! change idx 
     if (i/=del_atm_idx) then
        call find_atm(ter_ref,ref%atom_type(i),idx)
        eref%q_idx(idx)=q_idx(i)
        eref%atm_cls(idx)=atm_cls(i)
        eref%gr_no(idx)=group_no(i)
     end if
  end do
  eref%num_atm = ter_ref%num_atm
  atm_typ_no = eng_para%num_q_typ; na_gr = 0
  do ! change atom_cls and charge
     read(unit,*) string, atom_type; backspace(unit)
     if (string /= 'ATM') exit
     call find_atm(ter_ref,atom_type,idx)
     atm_typ_no = atm_typ_no + 1; na_gr = na_gr + 1; num_prm = num_prm + 1
     read(unit,*) string,atom_type,i,atom_type,eng_para%charge(atm_typ_no)
     eng_para%q_type(1:2,atm_typ_no) = &
          (/ter_type//resname(1:3),ter_ref%atom_type(idx)/)
     eref%q_idx(idx) = atm_typ_no; eref%gr_no(idx) = gr_no
     call find_atm_cls(atom_type,eref%atm_cls(idx))
  end do
  eng_para%num_q_typ = atm_typ_no; eref%na_gr(gr_no) = na_gr
!  print*,'done with atoms in ', ref%residue_type,' ', ter_type, 'ter'
!  do i = 1, eref%num_atm
!     write(*,'(i2,1x,a4,1x,a4,1x,f8.2)')i,ter_ref%atom_type(i),&
!          eng_para%atom_cls(eref%atm_cls(i)),eng_para%charge(eref%q_idx(i))
!  end do
!-----------------------------------------------------------------------
! change angle data
  delete_ang = 0; num_del_ang = 0
  do i = 1, eref%num_ang ! find ang to be deleted
     imax = na_ang(i); del = 0
     do j = 1, imax
        atom_type=atm_in_ang(j,i)
        if (num_del_atm/=0.and.atom_type==del_atm) then
           del = 1
        else
           call find_atm(ter_ref,atom_type,idx)
           if ((ter_type=='N'.and.idx<1) .or. &
                (ter_type=='C'.and.idx>ter_ref%num_atm)) then
              del = 1
           end if
        end if
     end do
     if (del == 1) then
        delete_ang(i) = 1; num_del_ang = num_del_ang+1
!        print*,'del ang',num_del_ang,(atm_in_ang(j,i),j=1,imax)
     end if
  end do
  idx = 0; num_add_ang = 0
  do ! add angles
     read(unit,*) string; backspace(unit)
     if (string/='DIH'.and.string/='IMP'.and.string/='ANG') exit
     idx = idx + 1; num_add_ang = num_add_ang + 1
     if (string=='DIH') then
        imax=4; eref%ang_type(idx)=1
     else if (string=='IMP') then
        imax=4; eref%ang_type(idx)=3
     else if (string=='ANG') then
        imax=3; eref%ang_type(idx)=2
     end if
     eref%na_ang(idx)=imax
     read(unit, *)string,(atom_in_ang(j,idx),j=1,imax)
  end do
  close(unit)
  do i = 1, eref%num_ang
     if (delete_ang(i)==0) then
        idx = idx + 1
        atom_in_ang(:,idx)=atm_in_ang(:,i)
        eref%ang_type(idx)=angle_type(i)
        eref%na_ang(idx)=na_ang(i)
     end if
  end do
  eref%num_ang = eref%num_ang+num_add_ang-num_del_ang ! update num_ang
!  print*, 'num ang in eng', eref%num_ang, idx
!-----------------------------------------------------------------------
  do ang_no = 1, eref%num_ang ! find indices: aa_idx, prm
     ang_type = eref%ang_type(ang_no)
     imax = eref%na_ang(ang_no)
     atm(1:imax)=atom_in_ang(1:imax,ang_no); atom = atm
     do i = 1, imax
        call find_atm(ter_ref,atm(i),idx)
        eref%aa_idx(i,ang_no) = idx
        if (ter_type=='N'.and.idx>ter_ref%num_atm) then
           next_atm=atm(i); next_atm=next_atm(2:) ! remove +
           call find_atm(next_ref,next_atm,idx)
           atm(i) = eng_para%atom_cls(next_eref%atm_cls(idx))
        else if (ter_type=='C'.and.idx<1) then
           prev_atm=atm(i); prev_atm=prev_atm(2:) ! remove -
           call find_atm(prev_ref,prev_atm,idx)
           atm(i) = eng_para%atom_cls(prev_eref%atm_cls(idx))
        else
           atm(i) = eng_para%atom_cls(eref%atm_cls(idx))
        end if
     end do
     call find_ang_prm(atm,ang_type,eref%prm(ang_no))
!     print*, ter_type, ang_no, atom(1:imax),' ', atm(1:imax)
!     idx = eref%prm(ang_no)
!     print*,'    ',eng_para%atm_in_ang(1:imax,idx),real(eng_para%ang_para(1:3,idx))
!     print*,''
  end do
!-----------------------------------------------------------------------
  if (ter_type == 'N') Nter_ref_eng = eref ! save eref 
  if (ter_type == 'C') Cter_ref_eng = eref ! save eref
end subroutine setup_ter_ref_eng
!-----------------------------------------------------------------------
subroutine initialize_ss_bnd()
!-----------------------------------------------------------------------
! set up Nter_ref or Cter_ref depending on the sequence
!-----------------------------------------------------------------------
  implicit none
  type(reference_residue_type) :: ref
  type(ref_res_eng_type) :: eref
  character(len=4) :: resname, string, atom_type
  integer :: atm_no, bnd_no, atm_typ_no, gr_no, ang_no, imax, ang_type
  integer :: idx, unit=19, i
  character(len=4), dimension(0:max_atm+1) :: atom_cls
  character(len=4), dimension(4) :: atm, atom
!-----------------------------------------------------------------------
  open(unit, file = infile_topo_geo)
  do ! find the start of the patch type
     read(unit,*) string
     if (string == 'RES') then
        backspace(unit); read(unit,*) string, string
        if (string == 'DISU') exit
     end if
  end do
!-----------------------------------------------------------------------
! atom
  atm_no = 0
  do 
     read(unit,*) string, atom_type; backspace(unit)
     if (string /= 'ATM') exit
     atm_no = atm_no + 1
     read(unit,*) string,ref%atom_type(atm_no),atom_cls(atm_no)
     call fill_atom_type4('DIS', ref%atom_type(atm_no)) !replace blanks by 0's
     call find_atm_cls(atom_cls(atm_no),eref%atm_cls(atm_no))
  end do
  ref%num_atm = atm_no
  ref%atom_type(0) = '1CA0'; ref%atom_type(atm_no+1) = '2CA0'
  atom_cls(0) = 'CH1E'; atom_cls(atm_no+1) = 'CH1E'
!-----------------------------------------------------------------------
! bond
  bnd_no = 0
  do 
     read(unit,*) string; backspace(unit)
     if (string /= 'BND') exit
     bnd_no = bnd_no + 1
     read(unit,*) string, (atm(i),i=1,2), ref%b_len(bnd_no)
     ! bnd parameter
     do i = 1, 2
        call find_atm(ref,atm(i),idx)
        atm(i) = atom_cls(idx)
     end do
     call find_bnd_prm(atm,eref%bnd_prm(bnd_no))
  end do
  ref%num_bnd = bnd_no
  ref%num_ang = 0
  close(unit)
!-----------------------------------------------------------------------
  open(unit, file = infile_topo_eng)
  do ! find the start of the patch type
     read(unit,*) string
     if (string == 'RES') then
        backspace(unit); read(unit,*) string, resname
        if (resname == 'DISU') exit
     end if
  end do
!-----------------------------------------------------------------------
! atom
  do 
     read(unit,*) string, atom_type; backspace(unit)
     if (string /= 'ATM') exit
     call find_atm(ref,atom_type,idx)
     if (idx <= ref%num_atm/2) then
        atm_typ_no=eng_para%num_q_typ+idx; num_prm=num_prm+1
     else 
        atm_typ_no=eng_para%num_q_typ-idx+ref%num_atm+1
     end if
     read(unit,*) string,atom_type,gr_no,atom_type, &
              eng_para%charge(atm_typ_no)
     eng_para%q_type(1:2,atm_typ_no)=(/'SCYS',atom_type/)
     eref%q_idx(idx) = atm_typ_no
  end do
  eref%num_gr = 2
  eng_para%num_q_typ = eng_para%num_q_typ + ref%num_atm/2
!-----------------------------------------------------------------------
  ang_no = 0
  do 
     read(unit,*) string; backspace(unit)
     if (string/='DIH'.and.string/='IMP'.and.string/='ANG') exit
     ang_no = ang_no + 1
     if (string=='DIH') then
        imax=4; ang_type=1
     else if (string=='IMP') then
        imax=4; ang_type=3
     else if (string=='ANG') then
        imax=3; ang_type=2
     end if
     eref%ang_type(ang_no)=ang_type; eref%na_ang(ang_no)=imax
     read(unit,*) string, (atm(i),i=1,imax); atom=atm
     do i = 1, imax
        call find_atm(ref,atm(i),idx)
        eref%aa_idx(i,ang_no) = idx
        atm(i) = atom_cls(idx)
     end do
     call find_ang_prm(atm,ang_type,eref%prm(ang_no))
     !         print*,resname,ang_type,atm(1:imax),atom(1:imax)
  end do
  eref%num_ang = ang_no
  close(unit)
!-----------------------------------------------------------------------
  ss_ref = ref; ss_ref_eng = eref
end subroutine initialize_ss_bnd
!-----------------------------------------------------------------------
subroutine setup_ss_bnd()
!-----------------------------------------------------------------------
  implicit none
  type(reference_residue_type) :: ref
  type(ref_res_eng_type) :: eref
  integer :: i, j, k, idx, ang_no
  character(len=4) :: atom
  
  do i = 1, protein%num_ss_bnd
     do j = 1, 2
        ref = ref_res(ref_ss(i)%res_no(j))
        eref = ref_res_eng(ref_ss(i)%res_no(j))
        call find_atm(ref,'SG00',ref_ss(i)%s_atm_no(j))
        do k = 1, ref%num_atm/2 ! change atom data
           atom = ss_ref%atom_type(k); atom = atom(2:)
           call find_atm(ref,atom,idx)
           eref%atm_cls(idx) = ss_ref_eng%atm_cls(k)
           eref%q_idx(idx) = ss_ref_eng%q_idx(k)
        end do
        if (j == 1) then ! add angles in the first residue
           ang_no = eref%num_ang
           do k = 1, ss_ref_eng%num_ang
              ang_no = ang_no + 1
              eref%ang_type(ang_no) = ss_ref_eng%ang_type(k)
              eref%na_ang(ang_no) = ss_ref_eng%na_ang(k)
              eref%aa_idx(:,ang_no) = ss_ref_eng%aa_idx(:,k) + 1000
              eref%prm(ang_no) = ss_ref_eng%prm(k)
           end do
           eref%num_ang = ang_no
        end if
        ref_res(ref_ss(i)%res_no(j)) = ref
        ref_res_eng(ref_ss(i)%res_no(j)) = eref
     end do
  end do
end subroutine setup_ss_bnd
!-----------------------------------------------------------------------
subroutine find_atm_cls(atm, index)
!-----------------------------------------------------------------------
! find index for atm als in eng_para
implicit none
   character(len=4), intent(in) :: atm
   integer, intent(out) :: index
   integer :: i
   do i = 1, eng_para%num_a_cls
      if (atm == eng_para%atom_cls(i)) then
         index = i; return
      end if
   end do
   write(*,*)'Error in finding atom cls no for ', atm
   stop
end subroutine find_atm_cls
!-----------------------------------------------------------------------
subroutine fill_atom_type3(atom_type)
!-----------------------------------------------------------------------
! fill atom_type with 0's if length of atom_type is less than 3
implicit none
   character(len=3), intent(inout) :: atom_type
   select case(len_trim(atom_type))
   case(1)
      atom_type = trim(atom_type)//'00'
   case(2)
      atom_type = trim(atom_type)//'0'
   case(3)
      atom_type = atom_type
   end select
end subroutine fill_atom_type3
!-----------------------------------------------------------------------
subroutine find_res(resname, res_no)
!-----------------------------------------------------------------------
! find res_no corresponding to resname
  implicit none
  character(len=4), intent(in) :: resname
  integer, intent(out) :: res_no
  integer :: i
  res_no = -1
  do i = 1, num_ref_res
     if (resname == residue_name(i)) then
        res_no = i; return
     end if
  end do
  write(*,*)'error: unknown residue name'
  stop
end subroutine find_res
!-----------------------------------------------------------------------
subroutine find_atm(ref, atm, index)
!-----------------------------------------------------------------------
! find index for atm in ref
implicit none
   type(reference_residue_type), intent(in) :: ref
   character(len=4), intent(in) :: atm
   integer, intent(out) :: index
   character(len=4) :: atom
   integer :: i
   atom = atm
   if (len_trim(atom) < 4) then
      call fill_atom_type4(ref%residue_type, atom)
   end if
   do i = -1, ref%num_atm+1
      if (atom == ref%atom_type(i)) then
         index = i; return
      end if
   end do
   write(*,*)'Error in finding atom no for ', atm
   stop
end subroutine find_atm
!-----------------------------------------------------------------------
subroutine find_bnd(ref,atm1,atm2,index,option)
!-----------------------------------------------------------------------
! find index for bnd (atm1,atm2) in ref
implicit none
   type(reference_residue_type), intent(in) :: ref
   character(len=4), intent(in) :: atm1, atm2
   integer, intent(out) :: index
   integer, intent(in) :: option
   integer :: i, head, tail
   call find_atm(ref, atm1, head)
   call find_atm(ref, atm2, tail)
   do i = -1, ref%num_bnd
      if (head == ref%b_idx(1,i) .and. tail == ref%b_idx(2,i)) then
         index = i; return
      end if
      if (head == ref%b_idx(2,i) .and. tail == ref%b_idx(1,i).and.i>=0) then
         index = ref%num_bnd+i+1; return
      end if
   end do
   if (option /= 1) then
      write(*,*)'Error in finding bond no for',atm1, atm2
      stop
   else
      index = -100; return
   end if
 end subroutine find_bnd
!-----------------------------------------------------------------------
subroutine find_ang(ref,ang_type,num_ang,atm,index, option)
!-----------------------------------------------------------------------
! find index for ang in ref, given ang type and atms in ang
! if option = 1, no warning if failed in finding ang. returns -100
implicit none
   type(reference_residue_type), intent(in) :: ref
   integer, intent(in) :: ang_type,num_ang, option
   character(len=4), dimension(:),intent(in) :: atm
   integer, intent(out) :: index
   integer :: i,j, imax
   integer, dimension(3) :: idx
   if (ang_type==1) imax=4; if (ang_type==2) imax=3
   if (ang_type==3 .and. option==1) then
      index = -100; return
   end if
   do i = 1, imax-1
      call find_bnd(ref,atm(i),atm(i+1),idx(i),1)
   end do
   do i = 0, num_ang !find ang corresponding to the bnd
      if (ref%ang_type(i)==ang_type) then
         do j = 1, imax-1
            if (idx(j) /= ref%a_idx(j,i)) goto 112
         end do
         index = i; return
      end if
112   continue
   end do
   if (option /= 1) then
      write(*,*) 'Error in finding ang no for', (atm(i),i=1,imax)
      stop
   else
      index = -100; return
   end if
 end subroutine find_ang
!-----------------------------------------------------------------------
subroutine fill_atom_type4(resname, atom_type)
!-----------------------------------------------------------------------
! fill atom_type with 0's if length of atom_type is less than 4
implicit none
   character(len=3), intent(in) :: resname
   character(len=4), intent(inout) :: atom_type
   if (resname == 'LYS' .and. atom_type(1:2) == 'HZ') then
      !insert '0' in the middle
      atom_type=atom_type(1:2)//'0'//atom_type(3:3) 
   else
      select case(len_trim(atom_type))
      case(1)
         atom_type = trim(atom_type)//'000'
      case(2)
         atom_type = trim(atom_type)//'00'
      case(3)
         atom_type = trim(atom_type)//'0'
      case(4)
         atom_type = atom_type
      end select
   end if
end subroutine fill_atom_type4
!-----------------------------------------------------------------------
subroutine find_ang_prm(atm,ang_type,idx)
!-----------------------------------------------------------------------
! find index for dih ang parameter in eng_para structure
  character(len=4), dimension(:), intent(in) :: atm
  integer, intent(in) :: ang_type
  integer, intent(out) :: idx
  character(len=4), dimension(4) :: atom
  integer :: i,imax
  if (ang_type == 1) then ! dih ang
     imax=4 
     do i = 1, eng_para%num_ang
        if (ang_type == eng_para%ang_type(i)) then
           atom(1:imax) = eng_para%atm_in_ang(1:imax,i)
           if (((atom(1)=='X'.and.atom(4)=='X') .and. &
                 ((atm(2)==atom(2).and.atm(3)==atom(3)) .or. &
                 (atm(2)==atom(3).and.atm(3)==atom(2)))) .or. &
                (atm(1)==atom(1).and.atm(2)==atom(2) .and. &
                 atm(3)==atom(3).and.atm(4)==atom(4)) .or. &
                (atm(1)==atom(4).and.atm(2)==atom(3) .and. &
                 atm(3)==atom(2).and.atm(4)==atom(1))) then
              idx = i; return
           end if
        end if
     end do
  else if (ang_type == 3) then ! imp ang
     imax=4
     do i = 1, eng_para%num_ang
        if (ang_type == eng_para%ang_type(i)) then
           atom(1:imax) = eng_para%atm_in_ang(1:imax,i)
           if (((atom(2)=='X'.and.atom(3)=='X') .and. &
                 ((atm(1)==atom(1).and.atm(4)==atom(4)) .or. &
                 (atm(1)==atom(4).and.atm(4)==atom(1)))) .or. &
                (atm(1)==atom(1).and.atm(2)==atom(2) .and. &
                 atm(3)==atom(3).and.atm(4)==atom(4)) .or. &
                (atm(1)==atom(4).and.atm(2)==atom(3) .and. &
                 atm(3)==atom(2).and.atm(4)==atom(1))) then
              idx = i; return
           end if
        end if
     end do
  else if (ang_type == 2) then
     imax=3
     do i = 1, eng_para%num_ang ! get bnd angles
        if (ang_type == eng_para%ang_type(i)) then
           atom(1:imax) = eng_para%atm_in_ang(1:imax,i)
           if (atm(2)==atom(2) .and. &
                ((atm(1)==atom(1).and.atm(3) == atom(3)) .or. &
                (atm(1)==atom(3).and.atm(3) == atom(1)))) then
              idx = i; return
           end if
        end if
     end do
  end if
   write(*,*)'Warning EEF1 can handle protein of residue size range 0 to',max_res 
  write(*,*)'failed in getting ang parm for',ang_type,atm(1:imax)
  stop
end subroutine find_ang_prm
!-----------------------------------------------------------------------
subroutine find_bnd_prm(atm,idx)
!-----------------------------------------------------------------------
! find index for bnd parameter in eng_para structure
  character(len=4), dimension(:), intent(in) :: atm
  integer, intent(out) :: idx
  character(len=4), dimension(2) :: atom
  integer :: i
  do i = 1, eng_para%num_bnd
     atom(1:2) = eng_para%atm_in_bnd(1:2,i)
     if ((atm(1)==atom(1).and.atm(2)==atom(2)) .or. &
          (atm(2)==atom(1).and.atm(1)==atom(2))) then
        idx = i; return
     end if
  end do
  write(*,*)'failed in getting ang parm for',ang_type,atm(1:imax)
  stop
end subroutine find_bnd_prm
!-----------------------------------------------------------------------
subroutine get_loca(resname,ter_type,loca,endc,ring)
!-----------------------------------------------------------------------
! returns data necessary for non-bonded list and 1,4 pair.
! numbering atoms in terms of position in the chain.   
!-----------------------------------------------------------------------
  implicit none
  character(len=4), intent(in) :: resname
  character(len=1), intent(in) :: ter_type
  character(len=4), dimension(:), intent(out) :: loca, endc
  integer, dimension(:), intent(out) :: ring
  integer :: na, res_no, i, idx
  call find_res(resname, res_no) 
  na = residue_ref(res_no)%num_atm
  loca(1:na) = new_numbering(res_no)%atom(1:na)
  endc(1:na) = new_numbering(res_no)%endc(1:na)
  ring(1:na) = new_numbering(res_no)%ring(1:na)
! correction for terminal residues
  if (ter_type=='N') then
     if (resname/='PRO') then
        do i = na,3,-1
           idx=i+2
           loca(idx)=loca(i); endc(idx)=endc(i); ring(idx)=ring(i)
        end do
        loca(1:4)=(/'HT10','N000','HT20','HT30'/)
        endc(1:4)=(/'ENDC','ENDC','HT20','HT30'/)
        ring(1:4)=(/     0,     0,     0,     0/)
     else
        do i = na,2,-1
           idx=i+2
           loca(idx)=loca(i); endc(idx)=endc(i); ring(idx)=ring(i)
        end do
        loca(1:3)=(/'HT10','N000','HT20'/)
        endc(1:3)=(/'ENDC','ENDC','HT20'/)
        ring(1:3)=(/     0,     1,     0/)
    end if
  else if (ter_type=='C') then
     loca(na)='OT10'; loca(na+1)='OT20'
     endc(na)='OT10'; endc(na+1)='OT20'; ring(na+1)=0
  end if
end subroutine get_loca
!-----------------------------------------------------------------------
subroutine find_num_string(long_string, num_string, strings)
!-----------------------------------------------------------------------
! scan long_string, get num of "words", and save in array strings
!-----------------------------------------------------------------------
  implicit none
  character(len=len_l), intent(in) :: long_string
  integer, intent(out) :: num_string
  character(len=len_s), dimension(:), intent(out) :: strings
  character(len=len_s) :: string
  integer :: i, n, k
  character(len=1) :: prev_char, cur_char
  n = len(long_string)
  if (n>len_l) then
     write(*,*) 'length of the long string is greater than maximum,',len_l
  end if
  num_string = 0; prev_char = ' '; string = ' '; k = 0
  do i = 1, n
     cur_char = long_string(i:i)
     if (cur_char == '!') exit
     if (prev_char==' ' .and. cur_char/=' ') then
        num_string = num_string + 1
        if (num_string>1) then
           strings(num_string-1)=string
           string=' '; k = 0
        end if
     end if
     if (cur_char /= ' ') then
        k = k + 1
        string(k:k) = cur_char
     end if
     prev_char = cur_char
  end do
  if (num_string>0) strings(num_string)=string
end subroutine find_num_string
!-----------------------------------------------------------------------
END MODULE in_out
!-----------------------------------------------------------------------
