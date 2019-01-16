!-----------------------------------------------------------------------
module geometry 
use in_out
contains
!-----------------------------------------------------------------------
subroutine fix_angles()
!-----------------------------------------------------------------------
! make list of fixed angles and angles for which angle change should be made
! 1. vary dih and bnd angles
! 2. vary dih angles only, and set bnd angles to ref ang
! 3. vary backbone dih angles only, and set other angs to ref 
! 4. vary phi and psi
!-----------------------------------------------------------------------
  implicit none
  integer :: i, j, k, ang_type, imax, bond, include_this_angle, ang_no
  integer :: num_ang, aux_ang, res_typ, num_res
  type(reference_residue_type) :: ref
  character(len=4) :: atom

  if (include_ang_type == 0) return

  num_res = protein%num_res

  ang_no = 0; aux_ang = 0
  do i = 1, num_res
     ref = ref_res(i)
     do j = 1, ref%num_ang
        ang_type = ref%ang_type(j); imax = ref%na_ang(j)-1
        include_this_angle = 0
        select case (include_ang_type)
        case(1) ! include all angles 
           include_this_angle = 1
        case(2) ! include dih angles only
           if (ang_type == 1) then
              include_this_angle = 1
           end if
        case(3) ! include backbone dihedrals only
           bond = ref%a_idx(imax,j)
           if (ang_type == 1 .and. &
                (bond==1.or.bond==2.or.bond==3) .and. &
                (i/=1.or.ref%a_idx(1,j)/=0) .and. &
                (i/=num_res.or.bond/=3)) then 
              include_this_angle = 1
           end if
        case(4) ! include phi and psi only
           bond = ref%a_idx(imax,j)
           if (ang_type == 1 .and. &
                (bond==2.or.bond==3).and. &
                (i/=1.or.ref%a_idx(1,j)/=0) .and. &
                (i/=num_res.or.bond/=3)) then 
              include_this_angle = 1
           end if
        case(5) ! side chain dihedral angles only
           bond = ref%a_idx(imax,j)
           atom = ref%atom_type(ref%b_idx(2,bond))
           if (ang_type == 1 .and. &
                atom /= 'N000' .and. atom /= 'CA00' .and. &
                atom /= 'C000' .and. atom /= 'O000' .and. &
                atom /= 'H000' .and. atom(1:2) /= 'OT' .and. &
                atom(1:2) /= 'HT' .and. atom /= '+N00') then
              include_this_angle = 1
           end if
        case default
           write(*,*)'Unknown angle type specifier.'
           stop
        end select
        if (i == 1) then ! exclude aux ang in first res
           do k = 1, imax
              if (ref%a_idx(k,j)<1) then
                 aux_ang = aux_ang + 1
                 if (include_ang_type == 2 .or. (aux_ang <= 3 &
                      .and. include_ang_type == 1)) then
                    include_this_angle = 0; exit
                 end if
              end if
           end do
        end if
        if (include_this_angle == 1) then
           ang_no = ang_no + 1
           ref_res(i)%ang_fixed(j) = 0
        else
           ref_res(i)%ang_fixed(j) = 1
        end if
     end do
  end do
  num_ang = ang_no

  if (include_ang_type == 3 .or. include_ang_type == 4) then
     ang_no = 0
     do i = 1, num_res
        ref = ref_res(i)
        ref_res(i)%ang_change = 0
        do j = 1, ref%num_ang
           if (ref%ang_fixed(j)==0) then
              ! delete first angle
              if (i == 1 .and. ang_no == 0) then
                 ref_res(i)%ang_fixed(j) = 1
                 goto 134
              end if
           end if
        end do
     end do
134  continue
  end if

  if (include_ang_type == 3 .or. include_ang_type == 4) then
     do i = 1, num_res
        ref = ref_res(i)
        ref_res(i)%ang_change = 0
        ref_res(i)%fixed_branch = 0
        do j = 1, ref%num_ang
           if (ref%ang_fixed(j)==0) then
              ang_type = ref%ang_type(j); bond = ref%a_idx(2,j)
              do k = 1, ref%num_ang
                 if (ref%ang_fixed(k) == 1 .and. &
                      ref%ang_type(k) == ang_type .and. &
                      ref%a_idx(2,k) == bond) then
                    ref_res(i)%ang_change(k) = j
                    ref_res(i)%fixed_branch(j) = k
                 end if
              end do
           end if
        end do
     end do
  end if

end subroutine fix_angles
!-----------------------------------------------------------------------
subroutine cartesian2internal(mode)
!-----------------------------------------------------------------------
! calc dih and bnd angles defined in "residue_ref" from bnd vectors
! set native values for fixed degrees of freedom
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: mode
  integer :: i, j, k, index, ang_type, imax, idx1, idx2
  real(dp), dimension(3,3) :: r
  real(dp), dimension(3) :: b
  real(dp) :: angle, small = 1.0d-10
  type(reference_residue_type) :: ref
  integer, dimension(2) :: ss_res, ss_atm

  if (mode > 0) then
     do i = 1, protein%num_res
        ref = ref_res(i)
        if (i < protein%num_res) then ! get coord for +N 
           protein%residue(i)%R(:,ref%num_atm+1) = protein%residue(i+1)%R(:,1)
        end if
        do j = 1, ref%num_bnd ! calc bond vectors
           idx1 = ref%b_idx(1,j); idx2 = ref%b_idx(2,j)
           b = protein%residue(i)%R(:,idx2) - protein%residue(i)%R(:,idx1)
           protein%residue(i)%b(:,j) = b
           protein%residue(i)%b_len(j) = sqrt(dot_product(b,b))
        end do
        if (i>1) then
           do j = -1, 0 ! (-C,N) & (-CA,-C) bnd
              protein%residue(i)%b(:,j) = protein%residue(i-1)%b(:,3+j)
              protein%residue(i)%b_len(j) = protein%residue(i-1)%b_len(3+j)
           end do
        else ! aux atms: set to arbitrary values
           protein%residue(i)%b(:,-1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
           protein%residue(i)%b_len(-1) = 1.0d0
           protein%residue(i)%b(:,0) = (/ 0.0d0, 1.0d0, 0.0d0 /)
           protein%residue(i)%b_len(0) = 1.0d0
        end if
     end do
     do i = 1, protein%num_ss_bnd
        ss_res(1:2) = ref_ss(i)%res_no(1:2)
        ss_atm(1:2) = ref_ss(i)%s_atm_no(1:2)
        b = protein%residue(ss_res(2))%R(:,ss_atm(2)) &
             - protein%residue(ss_res(1))%R(:,ss_atm(1))
        protein%ss_bnd(i)%b = b
        protein%ss_bnd(i)%b_len = sqrt(dot_product(b,b))
     end do
  end if

  do i = 1, protein%num_res
     ref = ref_res(i)
     do j = 1, ref%num_ang
        ang_type = ref%ang_type(j); imax = ref%na_ang(j)-1
        do k = 1, imax  ! get bond vectors consisting of the angle
           index = ref%a_idx(k,j)
!           if (index > 0 .or. i>1) then
              r(:,k) = protein%residue(i)%b(:,index)
!           else 
!              angle = ref%w0(j); goto 111 ! arbitrary aux ang
!           end if
           if (abs(sum(r(:,k)))<small) then !\ ang invloving 
              angle = ref%w0(j); goto 111   !/ undefined bnd
           end if
        end do
        call calc_angle(ang_type, r, angle)
        ref_res(i)%w0(j) = angle
111     protein%residue(i)%w(j) = angle
     end do

     do j = 1, ref%num_bnd
        b = protein%residue(i)%b(:,j)
        if (abs(sum(b))<small) then ! undefined bond
           protein%residue(i)%b_len(j) = ref%b_len(j)
        else
           protein%residue(i)%b_len(j) = sqrt(dot_product(b,b))
           ref_res(i)%b_len(j) = protein%residue(i)%b_len(j)
        end if
     end do

  end do
end subroutine cartesian2internal
!-----------------------------------------------------------------------
subroutine calc_angle(ang_type, r, angle)
!-----------------------------------------------------------------------
! dih or imp: r1=Rab, r2=Rbc, r3=Rcd : angle between planes abc and bcd
! bnd ang: r1=Rab, r2=Rbc: angle between -r1 and r2
!-----------------------------------------------------------------------
implicit none
   integer, intent(in) :: ang_type
   real(dp), dimension(:,:), intent(in) :: r
   real(dp), intent(out) :: angle
   real(dp), dimension(3) :: p, q, s
   real(dp) :: arg, one
   one = 1.0
   if (ang_type == 1 .or. ang_type ==3) then ! dihedral or improper
      call cross(r(:,1), r(:,2), p)
      call cross(r(:,2), r(:,3), q)
      call cross(r(:,3), r(:,1), s)
      arg = dot_product(p,q)/sqrt(dot_product(p,p)*dot_product(q,q))
      arg = sign(min(abs(arg),one),arg) ! to be sure abs(arg)<=1
      angle = sign(acos(arg), dot_product(s,r(:,2)))
   else if (ang_type == 2) then ! bond angle 
      p = -r(:,1); q = r(:,2)
      arg = dot_product(p,q)/sqrt(dot_product(p,p)*dot_product(q,q))
      arg = sign(min(abs(arg),one),arg) ! to be sure abs(arg)<=1
      angle = acos(arg)
   else
      write(*,*) 'Unknown angle type'
      stop
   end if
end subroutine calc_angle
!-----------------------------------------------------------------------
subroutine cross(p, q, s)
!-----------------------------------------------------------------------
implicit none
   real(dp), dimension(:), intent(in) :: p, q
   real(dp), dimension(:), intent(out) :: s
s(1)=p(2)*q(3)-p(3)*q(2);s(2)=p(3)*q(1)-p(1)*q(3);s(3)=p(1)*q(2)-p(2)*q(1)
end subroutine cross
!-----------------------------------------------------------------------
subroutine internal2cartesian_ww
!-----------------------------------------------------------------------
! calc cartesian coords from internal coords
! after change in backbone dihedral angles
!-----------------------------------------------------------------------
  implicit none
  type(residue_type) :: res
  type(reference_residue_type) :: ref
  integer :: i, j, k, bnd_no, ang_type, changed_ang
  integer, dimension(max_bnd) :: b_done
  real(dp), dimension(4) :: q
  real(dp), dimension(:,:), allocatable :: res_q
  real(dp), dimension(3) :: U, R_0, b
  real(dp) :: angle, angle_4,tan_w,tan_sqr,tan1, cosine,sine,p1,p2,p4
  real(dp) :: q1,q2,q3,q4
  integer, dimension(2) :: ss_res, ss_atm

  ! rotate some of fixed angles due to change in backbone angles
  do i = 1, protein%num_res
     ref = ref_res(i)
     do j = 1, ref%num_ang
        if (ref%ww_fixed(j)==1) then
           changed_ang = ref%ww_change(j)
           if (changed_ang /= 0) then
              protein%residue(i)%w(j) = protein_native%residue(i)%w(j) + &
                   protein%residue(i)%w(changed_ang) - protein_native%residue(i)%w(changed_ang)
           end if
        end if
     end do
  end do

  allocate(res_q(4,max_ang))
  q = unit_q; R_0 = origin ! initial q and and coord of first atm
  protein%residue(1)%R(:,:) = 0.0d0
  do i = 1, protein%num_res
     res = protein%residue(i) ! make a copy of residue(i)
     b_done = 0 ! bnd not calculated yet
     ref = ref_res(i)
     res%R(:,1) = R_0
     do k = 1, ref%num_dpn_ang(0) ! pass quaternion from prev res
        res_q(:,ref%i_dpn_ang(k,0)) = q
     end do
     do j = 1, ref%num_ang
        ang_type = ref%ang_type(j)
        angle = res%w(j)
        ! calc quaternion and q product
        angle_4  = 0.25*angle
        tan_w = tan(angle_4)
        tan_sqr = tan_w * tan_w; tan1 = 1.0d0 + tan_sqr
        cosine = (1.0d0 - tan_sqr)/tan1
        sine = 2.0d0*tan_w/tan1
        q1=res_q(1,j); q2=res_q(2,j); q3=res_q(3,j); q4=res_q(4,j)
        if (ang_type == 1) then
           p1 = cosine; p2 = sine
           q(1) = q1*p1 - q2*p2; q(2) = q1*p2 + p1*q2
           q(3) = p1*q3 + q4*p2; q(4) = p1*q4 - q3*p2
        else if (ang_type == 2) then
           p1 = sine; p4 = cosine
           q(1) = q1*p1 - q4*p4; q(2) = p1*q2 + q3*p4
           q(3) = p1*q3 - q2*p4; q(4) = p1*q4 + q1*p4
           res%sin_w(j) = 2.0d0*sine*cosine 
        end if
        res_q(:,j) = q
        do k = 1, ref%num_dpn_ang(j) ! pass q over to dpn ang
           res_q(:,ref%i_dpn_ang(k,j)) = q
        end do
        if (ref%num_dpn_bnd(j)>0) then ! 1/2*U
           q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4)
           U(1) = q1*q1-0.5d0+q2*q2 
           U(2) = q2*q3+q1*q4 
           U(3) = q2*q4-q1*q3
           bnd_no = ref%i_dpn_bnd(1,j)
           b = 2.0d0*res%b_len(bnd_no)*U
           res%b(:,bnd_no) = b; b_done(bnd_no) = 1
           res%R(:,ref%b_idx(2,bnd_no))=b+res%R(:,ref%b_idx(1,bnd_no))
        end if
     end do
     do j = 1, ref%num_bnd ! remaining bnds, e.g. bnd closing a ring
        if (b_done(j)==0) then
           b=res%R(:,ref%b_idx(2,j))-res%R(:,ref%b_idx(1,j))
           res%b(:,j)=b
           res%b_len(j)=sqrt(dot_product(b,b))
        end if
     end do
     if (i>1) then
        res%b(:,-1:0) = protein%residue(i-1)%b(:,2:3)!(-C,N)&(-CA,-C) bnd
        res%b_len(-1:0) = protein%residue(i-1)%b_len(2:3)
     end if
     protein%residue(i) = res ! save residue(i)
     q=res_q(:,ref%num_ang); R_0=res%R(:,ref%num_atm+1) !next q and R_0
  end do
  do i = 1, protein%num_ss_bnd      !ss bnds
     ss_res(1:2) = ref_ss(i)%res_no(1:2)
     ss_atm(1:2) = ref_ss(i)%s_atm_no(1:2)
     b = protein%residue(ss_res(2))%R(:,ss_atm(2)) &
          - protein%residue(ss_res(1))%R(:,ss_atm(1))
     protein%ss_bnd(i)%b = b
     protein%ss_bnd(i)%b_len = sqrt(dot_product(b,b))
  end do
  deallocate(res_q)
end subroutine internal2cartesian_ww
!-----------------------------------------------------------------------
subroutine internal2cartesian()
!-----------------------------------------------------------------------
! calc cartesian coords from internal coords 
!-----------------------------------------------------------------------
  implicit none
  type(residue_type) :: res
  type(reference_residue_type) :: ref
  integer :: i, j, k, bnd_no, ang_type, changed_ang
  integer, dimension(max_bnd) :: b_done
  real(dp), dimension(4) :: q
  real(dp), dimension(:,:), allocatable :: res_q
  real(dp), dimension(3) :: U, R_0, b
  real(dp) :: angle, angle_4,tan_w,tan_sqr,tan1, cosine,sine,p1,p2,p4
  real(dp) :: q1,q2,q3,q4
  integer, dimension(2) :: ss_res, ss_atm

  if (include_ang_type == 3 .or. include_ang_type == 4) then
     ! rotate some of fixed angles due to change in backbone angles
     do i = 1, protein%num_res
        ref = ref_res(i)
        do j = 1, ref%num_ang
           if (ref%ang_fixed(j)==1) then
              changed_ang = ref%ang_change(j)
              if (changed_ang /= 0) then
                 protein%residue(i)%w(j) = ref%w0(j) + &
             protein%residue(i)%w(changed_ang) - ref%w0(changed_ang)
              end if
           end if
        end do
     end do
  end if

  allocate(res_q(4,max_ang))
  q = unit_q; R_0 = origin ! initial q and and coord of first atm
  protein%residue(1)%R(:,:) = 0.0d0
  do i = 1, protein%num_res
     res = protein%residue(i) ! make a copy of residue(i)
     b_done = 0 ! bnd not calculated yet
     ref = ref_res(i)
     res%R(:,1) = R_0
     do k = 1, ref%num_dpn_ang(0) ! pass quaternion from prev res
        res_q(:,ref%i_dpn_ang(k,0)) = q
     end do
     do j = 1, ref%num_ang
        ang_type = ref%ang_type(j)
        angle = res%w(j)
        ! calc quaternion and q product
        angle_4  = 0.25*angle
        tan_w = tan(angle_4)
        tan_sqr = tan_w * tan_w; tan1 = 1.0d0 + tan_sqr
        cosine = (1.0d0 - tan_sqr)/tan1
        sine = 2.0d0*tan_w/tan1
        q1=res_q(1,j); q2=res_q(2,j); q3=res_q(3,j); q4=res_q(4,j)
        if (ang_type == 1) then
           p1 = cosine; p2 = sine
           q(1) = q1*p1 - q2*p2; q(2) = q1*p2 + p1*q2
           q(3) = p1*q3 + q4*p2; q(4) = p1*q4 - q3*p2
        else if (ang_type == 2) then
           p1 = sine; p4 = cosine
           q(1) = q1*p1 - q4*p4; q(2) = p1*q2 + q3*p4
           q(3) = p1*q3 - q2*p4; q(4) = p1*q4 + q1*p4
           res%sin_w(j) = 2.0d0*sine*cosine 
        end if
        res_q(:,j) = q
        do k = 1, ref%num_dpn_ang(j) ! pass q over to dpn ang
           res_q(:,ref%i_dpn_ang(k,j)) = q
        end do
        if (ref%num_dpn_bnd(j)>0) then ! 1/2*U
           q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4)
           U(1) = q1*q1-0.5d0+q2*q2 
           U(2) = q2*q3+q1*q4 
           U(3) = q2*q4-q1*q3
           bnd_no = ref%i_dpn_bnd(1,j)
           b = 2.0d0*res%b_len(bnd_no)*U
           res%b(:,bnd_no) = b; b_done(bnd_no) = 1
           res%R(:,ref%b_idx(2,bnd_no))=b+res%R(:,ref%b_idx(1,bnd_no))
        end if
     end do
     do j = 1, ref%num_bnd ! remaining bnds, e.g. bnd closing a ring
        if (b_done(j)==0) then
           b=res%R(:,ref%b_idx(2,j))-res%R(:,ref%b_idx(1,j))
           res%b(:,j)=b
           res%b_len(j)=sqrt(dot_product(b,b))
        end if
     end do
     if (i>1) then
        res%b(:,-1:0) = protein%residue(i-1)%b(:,2:3)!(-C,N)&(-CA,-C) bnd
        res%b_len(-1:0) = protein%residue(i-1)%b_len(2:3)
     end if
     protein%residue(i) = res ! save residue(i)
     q=res_q(:,ref%num_ang); R_0=res%R(:,ref%num_atm+1) !next q and R_0
  end do
  do i = 1, protein%num_ss_bnd      !ss bnds
     ss_res(1:2) = ref_ss(i)%res_no(1:2)
     ss_atm(1:2) = ref_ss(i)%s_atm_no(1:2)
     b = protein%residue(ss_res(2))%R(:,ss_atm(2)) &
          - protein%residue(ss_res(1))%R(:,ss_atm(1))
     protein%ss_bnd(i)%b = b
     protein%ss_bnd(i)%b_len = sqrt(dot_product(b,b))
  end do
  deallocate(res_q)
end subroutine internal2cartesian
!-----------------------------------------------------------------------
subroutine quaternion(axis_type, axis, quater_w, p)
!-----------------------------------------------------------------------
! calculate quaternion, given rotation axis and angle. 
! axis_type = 1 (axis on the x axis), 2 (z axis), 3 (3-d)
!-----------------------------------------------------------------------
implicit none
   integer, intent(in) :: axis_type
   real(dp), dimension(:), intent(in) :: axis
   real(dp), intent(in) :: quater_w
   real(dp), dimension(:), intent(out) :: p
   real(dp) :: tan_w, tan_sqr, tan1, cosine, sine
   tan_w = tan(quater_w)
   tan_sqr = tan_w * tan_w; tan1 = 1.0d0 + tan_sqr
   cosine = (1.0d0 - tan_sqr)/tan1
   sine = 2.0d0*tan_w/tan1
   p(1) = cosine
   if (axis_type == 1) then
      p(2) = sine; p(3:4) = 0.0d0
   else if (axis_type == 2) then
      p(2:3) = 0.0d0; p(4) = sine
   else if (axis_type == 3) then
      p(2:4) = axis(1:3) * sine
   else
      write(*,*) 'Unknown axis type.'
      stop
   end if
end subroutine quaternion
!-----------------------------------------------------------------------
subroutine q_product(axis_type, q, p, s)
!-----------------------------------------------------------------------
! calculates product of two quaternions p and q.
! axis_type = 1 (p(3:4)=0), 2 (p(2:3)=0), 3 
!-----------------------------------------------------------------------
implicit none
   integer, intent(in) :: axis_type
   real(dp), dimension(:), intent(in) :: p, q
   real(dp), dimension(:), intent(out) :: s
   real(dp), dimension(3) :: pi, qi, qcp
   real(dp) :: p1,p2,p4,q1,q2,q3,q4
   if (axis_type == 1) then ! p(3:4)=0
     p1=p(1); p2=p(2)
     q1=q(1); q2=q(2); q3=q(3); q4=q(4)
     s(1) = q1*p1 - q2*p2; s(2) = q1*p2 + p1*q2
     s(3) = p1*q3 + q4*p2; s(4) = p1*q4 - q3*p2
   else if (axis_type == 2) then ! p(2:3)=0
     p1=p(1); p4=p(4)
     q1=q(1); q2=q(2); q3=q(3); q4=q(4)
     s(1) = q1*p1 - q4*p4; s(2) = p1*q2 + q3*p4
     s(3) = p1*q3 - q2*p4; s(4) = p1*q4 + q1*p4
  else if (axis_type == 3) then
     qi = q(2:4);  pi = p(2:4)
     s(1) = q(1)*p(1) - dot_product(qi,pi)
     call cross(qi, pi, qcp)
     s(2:4) = q(1)*pi + p(1)*qi + qcp
  else
     write(*,*) 'Unknown axis type.'
     stop
  end if
end subroutine q_product
!-----------------------------------------------------------------------
subroutine rotation_matrix(vec_type, q, U)
!-----------------------------------------------------------------------
! constructs rotation matrix U from quaternion q.
! vec_type = 1 (vector rotated by U is on the x axis), 3 (3-d)
!-----------------------------------------------------------------------
implicit none
   integer, intent(in) :: vec_type
   real(dp), dimension(:), intent(in) :: q
   real(dp), dimension(:,:), intent(out) :: U
   real(dp) :: q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12,q13,q22,q23,q33  
   q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4)
   b0 = 2.0d0*q0; b1 = 2.0d0*q1
   q00 = b0*q0-1.0d0; q02 = b0*q2; q03 = b0*q3
   q11 = b1*q1;       q12 = b1*q2; q13 = b1*q3  
   if (vec_type == 1) then ! need only U(:,1)
      U(1,1) = q00+q11; U(2,1) = q12+q03; U(3,1) = q13-q02
   else if (vec_type == 3) then 
      b2 = 2.0d0*q2; b3 = 2.0d0*q3
      q01 = b0*q1; q22 = b2*q2; q23 = b2*q3; q33 = b3*q3 
      U(1,1) = q00+q11; U(1,2) = q12-q03; U(1,3) = q13+q02
      U(2,1) = q12+q03; U(2,2) = q00+q22; U(2,3) = q23-q01
      U(3,1) = q13-q02; U(3,2) = q23+q01; U(3,3) = q00+q33
   else
      write(*,*) 'Unknown vector type.'
      stop
   end if
end subroutine rotation_matrix
!-----------------------------------------------------------------------
end module geometry
!-----------------------------------------------------------------------
