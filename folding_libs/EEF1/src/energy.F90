!-----------------------------------------------------------------------
MODULE energy
!-----------------------------------------------------------------------
  use geometry
! lbfgs variables
  integer :: m_lbfgsb = 5 ! no of limited memory corrections stored
!-----------------------------------------------------------------------
! shifting function parameter
!  real(dp), parameter :: Rc_sqr = 64.0d0
! switching function parameter
  real(dp), parameter :: Ron=7.0d0, Roff=9.0d0, Rout = 10.0d0, &
       R_cufoff_slv_sqr=81.0d0
  real(dp) :: Ron_sqr, Roff_sqr, R3on_off, Ron_off_3, Rout_sqr, Rdiff_sqr
!-----------------------------------------------------------------------
  integer, dimension(:), allocatable :: num_pair, no_pair
  integer, dimension(:,:), allocatable :: i_14, pair_list, i_R, i_w, i_rr
  integer, dimension(:), allocatable :: LJ_type, ES_type, ipivot
  integer, dimension(:), allocatable :: start_rot, end_rot, invert_axis
  real(dp) :: Efac = 332.0716D0 ! prefactor multiplied in ES intxn
  real(dp) :: Eref_slv
  real(dp), dimension(:,:), allocatable :: R, R_gr, dfdr, R_gr_prev
  integer, dimension(:), allocatable :: Nagr
  integer, dimension(:,:), allocatable :: atm_gr, pair
  real(dp), parameter :: small = 1.0d-10
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
subroutine localmin_energy(f, iprint, status, mode)
!-----------------------------------------------------------------------
! This subroutine uses lbfgs-b for energy minimization
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: iprint, mode
  integer, intent(out) :: status
  real(dp), intent(out) :: f
  integer :: i, j, tnw, tnb, tnf, tnatm, min, tngr, iter, ii, unit = 31, npt, nfrm
  integer, parameter :: mmax=17 !max no of limited memory corrections.
  character(len=60) :: task, csave
  logical, dimension(4) :: lsave
  integer, dimension(44) :: isave
  real(dp), dimension(29) :: dsave
  integer, dimension(:), allocatable :: iwa, nbd
  real(dp), dimension(:), allocatable :: w_l, w_u, g, wa, w
  character(len=12) :: dcd_name
  type(protein_type) :: proteini
!-----------------------------------------------------------------------
  tnw = tot_num_ang
  tnf = tot_num_dof
  allocate(nbd(tnf), w_l(tnf), w_u(tnf), g(tnf), w(tnf))
  allocate(wa(2*mmax*tnf+4*tnf+12*mmax*mmax+12*mmax),iwa(3*tnf))
  if (mode == 0) then
     nbd(1:tnf) = 0 ! no bounds
  else
     nbd(1:tnw) = 0 ! no bounds
  endif

  N_update = 0; min = 1
  tngr = tot_num_gr
  allocate(R_gr_prev(3,tngr))

  if (mode == 0) then
     do i = 1, tnf, 3
        j = int((i-1)/3) + 1
        w(i:i+2) = protein%residue(i_rr(1,j))%R(1:3,i_rr(2,j))
     end do
  else
     ! get angles
     do i = 1, tnw
        w(i) = protein%residue(i_w(1,i))%w(i_w(2,i))
     end do
  endif

  iter = 0
  task = 'START'   ! start the iteration by initializing task.
!        ------- the beginning of the loop ----------
111 continue
!     This is the call to the L-BFGS-B code.
  call setulb(tnf,m_lbfgsb,w,w_l,w_u,nbd,f,g,factr_eng,pgtol_eng, &
       wa,iwa,task,iprint,csave,lsave,isave,dsave)
  if (isave(30) > max_iter_lbfgsb) then
     !write(*,*) 'Number of interation exceeded maximum', max_iter_lbfgsb
     !write(*,'(a6,2a12)') 'Tnf', 'Projg','F'
     !write(*,'(i6,e12.3,f12.3)') isave(34), dsave(13), dsave(2)
     status = 1
     goto 121
  end if

  if (task(1:2) == 'FG') then
!        the minimization routine has returned to request the
!        function f and gradient g values at the current x.
     if (mode == 0) then
        do i = 1, tnf, 3
           j = int((i-1)/3) + 1
           protein%residue(i_rr(1,j))%R(1:3,i_rr(2,j)) = w(i:i+2)
        end do
     else
        do i = 1, tnw ! change data in protein
           protein%residue(i_w(1,i))%w(i_w(2,i)) = w(i)
        end do
        call internal2cartesian() ! w -> R
     endif
     call energy_and_gradient(f, g, 0, min, 1, mode)
     !        go back to the minimization routine.
     goto 111
  endif
  if (task(1:5) == 'NEW_X')  goto 111
!        the minimization routine has returned with a new iterate,
!         and we have opted to continue the iteration.

!           ---------- the end of the loop -------------
 
!     If task is neither FG nor NEW_X we terminate execution.
  if (task(1:4)=='ABNO' .or. task(1:5)=='ERROR') then
     if (iprint >= 0) then
        write(*,*)'Abnormal termination or error in the lbfgs minimization.'
     end if
     status = 0
  else
     status = 1
  end if

121 continue
  close(unit)

! save the minimized configuration
  if (mode == 0) then
     do i = 1, tnf, 3
        j = int((i-1)/3) + 1
        protein%residue(i_rr(1,j))%R(1:3,i_rr(2,j)) = w(i:i+2)
     end do
  else
     do i = 1, tnw
        protein%residue(i_w(1,i))%w(i_w(2,i)) = w(i)
     end do
     call internal2cartesian()
  endif

  deallocate(R_gr_prev)
  deallocate(nbd, w_l, w_u, g, w, wa, iwa)

end subroutine localmin_energy
!-----------------------------------------------------------------------
subroutine initialize_energy(mode)
!-----------------------------------------------------------------------
! initialize for Cartesian min, if mode == 0, 
! and for internal min, if mode == 1
  implicit none
  integer, intent(in) :: mode
  integer :: num_res
  real(dp) :: Ron_off, Rdiff, E, x
  integer :: m1, m2, m3, n1, i, j, gr_no, na, tna, res_no, atm_no

  num_res = protein%num_res
  tot_num_atm = sum(ref_res(1:num_res)%num_atm)
     
  call make_ang_list(0) ! get num_tot_ang

  if (mode == 0) then ! CARTESIAN_MIN
     tot_num_dof = tot_num_atm*3
  else ! internal min
     tot_num_dof = tot_num_ang
  endif ! end CARTESIAN_MIN
  tot_num_ang_e = sum(ref_res_eng(1:num_res)%num_ang)

! allocate variables that depend only on the atom types and connectivity,
! and not on the conformation of protein
! **** variables ****
! i_R(1:3,i,p): res no,atm no,gr_no for the new atm i, protein p 
! LJ_type, ES_type: atm type in the purpose of LJ or ES intrxn 
! ipivot, start_rot, end_rot, invert_axis: needed for gradient calc
  m1 = tot_num_atm
  m2 = tot_num_ang
  m3 = tot_num_dof
  allocate(i_14(m1,m1), num_pair(m1), pair_list(m1,m1), R(3,m1), &
       i_R(3,m1), LJ_type(m1), ES_type(m1), no_pair(m1), pair(m1,m1), i_w(2,m2)) 
  if (mode == 0) then ! CARTESIAN_MIN
     allocate(i_rr(2,m3))
  else ! internal min
     allocate(invert_axis(m2), ipivot(m2), start_rot(m2),end_rot(m2))
  endif ! end CARTESIAN_MIN

! for switching function
  Ron_sqr = Ron*Ron; Roff_sqr = Roff*Roff 
  Ron_off = Roff_sqr - Ron_sqr
  Ron_off_3 = 1.0d0/(Ron_off**3); R3on_off = Roff_sqr - 3.0d0*Ron_sqr
! for group pair list update
  Rout_sqr = Rout*Rout; Rdiff = (Rout-Roff)/2; Rdiff_sqr = Rdiff*Rdiff

  call make_ang_list(1) 

  call make_lists(mode) 
! eef1 ref slv E
  E = 0.0d0
  do j = 1, tot_num_atm
     E = E + eng_para%slv_para(2,LJ_type(j))
  end do
  Eref_slv = E

! group no
  m3 = tot_num_gr
  allocate(Nagr(m3), atm_gr(5,m3), R_gr(3,m3))
  Nagr = 0
  do j = 1, tot_num_atm
     gr_no = i_R(3,j)
     na = Nagr(gr_no) + 1
     Nagr(gr_no) = na
     atm_gr(na,gr_no) = j
  end do

end subroutine initialize_energy
!-----------------------------------------------------------------------
subroutine finalize_energy(mode)
!-----------------------------------------------------------------------
! for Cartesian min, if mode == 0, 
! and for internal min, if mode == 1
  implicit none
  integer, intent(in) :: mode

  deallocate(i_14, num_pair, pair_list, R, i_R, LJ_type, ES_type)
  deallocate(no_pair, pair, i_w, Nagr, atm_gr, R_gr) 
  if (mode == 0) then ! CARTESIAN_MIN
     deallocate(i_rr)
  else
     deallocate(invert_axis, ipivot, start_rot, end_rot)
  endif

end subroutine finalize_energy
!-----------------------------------------------------------------------
subroutine energy_and_gradient(f, g, print_eng, min, calc_g, mode)
!-----------------------------------------------------------------------
! calculates energy and the gradient of energy wrt angles.
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: print_eng, min, calc_g, mode
  real(dp), intent(out) :: f
  real(dp), dimension(:), intent(out) :: g
  real(dp) :: Evdw, Ees, Eslv, Eb
  real(dp), dimension(3) :: Eint
  integer ::  tnatm, tngr
 
  tnatm = tot_num_atm 
  allocate(dfdr(3,tnatm))

  call pairlist_update(min)

  ! initialize energy derivatives
  dfdr = 0.0d0
  
  call bond_energy(Eb)

  call angle_energy(Eint)
 
  call non_bonded_energy(Evdw,Ees,Eslv)

  f = Evdw + Eb + sum(Eint) + Ees + Eslv

  if (print_eng == 1) then
     write(*,'(7a11)')'vdw','bnd','dih','ang','imp','es','slv'
     write(*,'(7f11.5)')Evdw,Eb,Eint(1:3),Ees,Eslv
     write(*,'(a11)')'Etot'
     write(*,'(f11.5)')f
  end if

  if (calc_g == 1) then
     call calc_gradient(g,mode)
  end if

  deallocate(dfdr)

end subroutine energy_and_gradient
!-----------------------------------------------------------------------
subroutine non_bonded_energy(Evdw,Ees,Eslv)
!-----------------------------------------------------------------------
  implicit none
  real(dp), intent(out) :: Evdw, Ees, Eslv
  real(dp), dimension(3) :: Ri, Rij
  integer :: i, j, k, idx, calc_r_sqr, calc_r, np
  integer :: LJi, LJj, ESi, ESj, gri, grj, Nai, Naj
  real(dp) :: fvdw, fes, fslv
  real(dp) :: r_sqr, epsilon_r6, r_abs, r_inv, r2, r6
  real(dp) :: xij,xji,fact_ij,fact_ji,L_i,L_j,tmpij,tmpji,tmp1,tmp2
  real(dp) :: df0, df0_slv
  real(dp) :: dr1, dr2, sw, sw_p, Rgr_sqr, sig2, fac, tmp
  real(dp), dimension(3) :: Rgr, df_sw_i, df_sw_j, df
  integer :: eps_idx, sig_idx, e14

  Eslv = Eref_slv; Evdw = 0.0d0; Ees = 0.0d0
  do i = 1, tot_num_atm-1
     Ri = R(:,i)
     gri = i_R(3,i)
     Nai = Nagr(gri)
     LJi = LJ_type(i)
     ESi = ES_type(i)
     do np = 1, no_pair(i)
        j = pair(np,i)
        grj = i_R(3,j)

        ! LJ and ES intrxn by groups 
        fvdw = 0.0d0; fes = 0.0d0; fslv = 0.0d0; df0 = 0.0d0
        calc_r_sqr = 0
        Rgr = R_gr(:,grj) - R_gr(:,gri) ! group distance
        Rgr_sqr = Rgr(1)*Rgr(1) + Rgr(2)*Rgr(2) + Rgr(3)*Rgr(3)
        if (Rgr_sqr <= Roff_sqr) then 

           ! index for 1,4 pair
           e14 = i_14(j,i)
           if (e14 == 1) then
              eps_idx = 3; sig_idx = 4 
           else
              eps_idx = 1; sig_idx = 2
           end if

           LJj = LJ_type(j)
           ESj = ES_type(j)

           ! LJ
           Rij = R(:,j) - Ri
           r_sqr = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)   
           calc_r_sqr = 1
           sig2 = eng_para%aux_para(sig_idx,LJi,LJj)
           r2 = sig2/r_sqr; r6 = r2*r2*r2
           epsilon_r6 = r6 * eng_para%aux_para(eps_idx,LJi,LJj)
           fvdw = epsilon_r6 * (2.0d0 - r6) ! epsilone is negative
           df0 = -12.0d0 * (fvdw - epsilon_r6)
           ! ES
           if (abs(eng_para%charge(ESi))>small .and. &
                abs(eng_para%charge(ESj))>small) then ! finite charge
              if (e14==1) then
                 fac = Efac*eng_para%E14fac
              else
                 fac = Efac
              end if
              fes = fac*eng_para%charge(ESi)*eng_para%charge(ESj) / r_sqr
              df0 = df0 - 2.0d0*fes
           end if

           ! EEF1 solvation
           calc_r = 0; df0_slv = 0
           tmp = eng_para%aux_para(5,LJi,LJj)
           if (abs(tmp)>small) then
              fact_ij = tmp
              r_abs = sqrt(r_sqr)
              r_inv = 1.0d0/r_abs
              calc_r = 1
              L_i = eng_para%slv_para(4,LJi) 
              xij = (r_abs-eng_para%slv_para(5,LJi))*L_i
              tmpij = fact_ij*exp(-xij*xij) 
              fslv = -tmpij
              df0_slv = tmpij*(xij*L_i+r_inv)
           end if
           tmp = eng_para%aux_para(5,LJj,LJi)
           if (abs(tmp)>small) then
              fact_ji = tmp
              if (calc_r == 0) then
                 r_abs = sqrt(r_sqr) 
                 r_inv = 1.0d0/r_abs
                 calc_r = 1
              end if
              L_j = eng_para%slv_para(4,LJj)
              xji = (r_abs-eng_para%slv_para(5,LJj))*L_j
              tmpji = fact_ji*exp(-xji*xji)
              fslv = fslv - tmpji
              df0_slv = df0_slv + tmpji*(xji*L_j+r_inv)
           end if
           if (calc_r == 1) then
              fslv = fslv/r_sqr; df0 = df0 + 2.0d0*df0_slv*r_inv
           end if

           ! Switching function 
           if (Rgr_sqr > Ron_sqr) then
              dr1 = Roff_sqr - Rgr_sqr
              dr2 = R3on_off + 2.0d0*Rgr_sqr
              tmp1 = dr1*Ron_off_3
              sw = dr1*tmp1*dr2
              sw_p = 4.0d0*tmp1*(dr1-dr2)
              df0 = df0*sw
              tmp1 = fvdw + fes + fslv
              fvdw = fvdw*sw; fes = fes*sw; fslv = fslv*sw
              tmp2 = tmp1*sw_p; Naj = Nagr(grj)
              df_sw_i = -tmp2/Nai*Rgr; df_sw_j = tmp2/Naj*Rgr
              do k = 1, Nai
                 idx = atm_gr(k,gri)
                 dfdr(:,idx) = dfdr(:,idx) + df_sw_i
              end do
              do k = 1, Naj
                 idx = atm_gr(k,grj)
                 dfdr(:,idx) = dfdr(:,idx) + df_sw_j
              end do
           end if

           df = df0/r_sqr*Rij
           Evdw = Evdw + fvdw; Ees = Ees + fes; Eslv = Eslv + fslv
           dfdr(:,i) = dfdr(:,i) - df; dfdr(:,j) = dfdr(:,j) + df
        end if

     end do ! j loop
  end do ! i loop

end subroutine non_bonded_energy
!-----------------------------------------------------------------------
subroutine non_bonded_energy_only(Evdw,Ees,Eslv)
!-----------------------------------------------------------------------
  implicit none
  real(dp), intent(out) :: Evdw, Ees, Eslv
  real(dp), dimension(3) :: Ri, Rij
  integer :: i, j, k, idx, calc_r_sqr, calc_r, np
  integer :: LJi, LJj, ESi, ESj, gri, grj, Nai, Naj
  real(dp) :: fvdw, fes, fslv
  real(dp) :: r_sqr, epsilon_r6, r_abs, r_inv, r2, r6
  real(dp) :: xij,xji,fact_ij,fact_ji,L_i,L_j,tmpij,tmpji,tmp1,tmp2
  real(dp) :: dr1, dr2, sw, sw_p, Rgr_sqr, sig2, fac, tmp
  real(dp), dimension(3) :: Rgr
  integer :: eps_idx, sig_idx, e14

  Eslv = Eref_slv; Evdw = 0.0d0; Ees = 0.0d0
  do i = 1, tot_num_atm-1
     Ri = R(:,i)
     gri = i_R(3,i)
     Nai = Nagr(gri)
     LJi = LJ_type(i)
     ESi = ES_type(i)
     do np = 1, no_pair(i)
        j = pair(np,i)
        grj = i_R(3,j)

        ! LJ and ES intrxn by groups 
        fvdw = 0.0d0; fes = 0.0d0; fslv = 0.0d0
        calc_r_sqr = 0
        Rgr = R_gr(:,grj) - R_gr(:,gri) ! group distance
        Rgr_sqr = Rgr(1)*Rgr(1) + Rgr(2)*Rgr(2) + Rgr(3)*Rgr(3)
        if (Rgr_sqr <= Roff_sqr) then 

           ! index for 1,4 pair
           e14 = i_14(j,i)
           if (e14 == 1) then
              eps_idx = 3; sig_idx = 4 
           else
              eps_idx = 1; sig_idx = 2
           end if

           LJj = LJ_type(j)
           ESj = ES_type(j)

           ! LJ
           Rij = R(:,j) - Ri
           r_sqr = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)   
           calc_r_sqr = 1
           sig2 = eng_para%aux_para(sig_idx,LJi,LJj)
           r2 = sig2/r_sqr; r6 = r2*r2*r2
           epsilon_r6 = r6 * eng_para%aux_para(eps_idx,LJi,LJj)
           fvdw = epsilon_r6 * (2.0d0 - r6) ! epsilone is negative
           ! ES
           if (abs(eng_para%charge(ESi))>small .and. &
                abs(eng_para%charge(ESj))>small) then ! finite charge
              if (e14==1) then
                 fac = Efac*eng_para%E14fac
              else
                 fac = Efac
              end if
              fes = fac*eng_para%charge(ESi)*eng_para%charge(ESj) / r_sqr
           end if

           ! EEF1 solvation
           calc_r = 0
           tmp = eng_para%aux_para(5,LJi,LJj)
           if (abs(tmp)>small) then
              fact_ij = tmp
              r_abs = sqrt(r_sqr)
              r_inv = 1.0d0/r_abs
              calc_r = 1
              L_i = eng_para%slv_para(4,LJi) 
              xij = (r_abs-eng_para%slv_para(5,LJi))*L_i
              tmpij = fact_ij*exp(-xij*xij) 
              fslv = -tmpij
           end if
           tmp = eng_para%aux_para(5,LJj,LJi)
           if (abs(tmp)>small) then
              fact_ji = tmp
              if (calc_r == 0) then
                 r_abs = sqrt(r_sqr) 
                 r_inv = 1.0d0/r_abs
                 calc_r = 1
              end if
              L_j = eng_para%slv_para(4,LJj)
              xji = (r_abs-eng_para%slv_para(5,LJj))*L_j
              tmpji = fact_ji*exp(-xji*xji)
              fslv = fslv - tmpji
           end if
           if (calc_r == 1) then
              fslv = fslv/r_sqr
           end if

           ! Switching function 
           if (Rgr_sqr > Ron_sqr) then
              dr1 = Roff_sqr - Rgr_sqr
              dr2 = R3on_off + 2.0d0*Rgr_sqr
              tmp1 = dr1*Ron_off_3
              sw = dr1*tmp1*dr2
              sw_p = 4.0d0*tmp1*(dr1-dr2)
              tmp1 = fvdw + fes + fslv
              fvdw = fvdw*sw; fes = fes*sw; fslv = fslv*sw
              tmp2 = tmp1*sw_p; Naj = Nagr(grj)
           end if
           Evdw = Evdw + fvdw; Ees = Ees + fes; Eslv = Eslv + fslv
        end if

     end do ! j loop
  end do ! i loop

end subroutine non_bonded_energy_only
!-----------------------------------------------------------------------
subroutine angle_energy(Eint)
!-----------------------------------------------------------------------
! calc f and df = dfi/d(theta)
!-----------------------------------------------------------------------
  implicit none
  real(dp), dimension(3), intent(out) :: Eint
  integer :: ang_type, prm_no, na, i, j, k, idx
  integer, dimension(2) :: indx, cor_ang
  real(dp) :: angle, fint, sg, arg, one, pq, sine
  real(dp), dimension(3,4) :: df_int
  real(dp), dimension(3,3) :: ra
  real(dp) :: n, dw, kw, tmp, w0, tan_w, tan_sqr, tan1, cosw, sinw
  real(dp) :: r23_sqr, r23_norm, p_norm_sqr, q_norm_sqr, p_dot_q, df0
  real(dp), dimension(3) :: p, q, s, ra1, ra2, ra3
  Eint = 0.0d0
  do i = 1, protein%num_res
     do j = 1, ref_res_eng(i)%num_ang
        ang_type = ref_res_eng(i)%ang_type(j)
        do k = 1, ref_res_eng(i)%na_ang(j)-1  ! get bond vectors
           indx(1:2) = ref_res_eng(i)%atm_idx(k:k+1,j)
           ra(:,k) = R(:,indx(2))-R(:,indx(1))
        end do
        ! calc angle
        if (ang_type == 1 .or. ang_type ==3) then ! torsion angle
           ra1=ra(:,1); ra2=ra(:,2); ra3=ra(:,3)
           p(1)=ra1(2)*ra2(3)-ra1(3)*ra2(2)
           p(2)=ra1(3)*ra2(1)-ra1(1)*ra2(3)
           p(3)=ra1(1)*ra2(2)-ra1(2)*ra2(1)
           q(1)=ra2(2)*ra3(3)-ra2(3)*ra3(2)
           q(2)=ra2(3)*ra3(1)-ra2(1)*ra3(3)
           q(3)=ra2(1)*ra3(2)-ra2(2)*ra3(1)
           s(1)=ra3(2)*ra1(3)-ra3(3)*ra1(2)
           s(2)=ra3(3)*ra1(1)-ra3(1)*ra1(3)
           s(3)=ra3(1)*ra1(2)-ra3(2)*ra1(1)
           sg = s(1)*ra2(1) + s(2)*ra2(2) + s(3)*ra2(3)
        else if (ang_type == 2) then ! bond angle 
           p = -ra(:,1); q = ra(:,2)
           sg = 1.0d0
        end if
        p_dot_q = p(1)*q(1) + p(2)*q(2) + p(3)*q(3)
        p_norm_sqr = p(1)*p(1) + p(2)*p(2) + p(3)*p(3)
        q_norm_sqr = q(1)*q(1) + q(2)*q(2) + q(3)*q(3)
        pq = sqrt(p_norm_sqr*q_norm_sqr)
        arg = p_dot_q/pq
        one = 1.0
        arg = sign(min(abs(arg),one),arg) ! to be sure abs(arg)<=1
        angle = sign(acos(arg), sg)

        ! energy
        prm_no = ref_res_eng(i)%prm(j)
        w0 = eng_para%ang_para(2,prm_no)
        kw = eng_para%ang_para(1,prm_no)
        if (ang_type == 1) then ! proper torsion
           n = eng_para%ang_para(3,prm_no)
           tan_w = tan(0.5d0*(n*angle - w0))
           tan_sqr = tan_w * tan_w; tan1 = 1.0d0 + tan_sqr
           cosw = (1.0d0 - tan_sqr)/tan1
           sinw = 2.0d0*tan_w/tan1
           fint = kw*(1.0d0+cosw)
           df0 = -kw*n*sinw
        else  ! bnd ang & improper torsion
           dw = angle - w0
           tmp = kw*dw
           fint = tmp*dw
           df0 = 2.0d0*tmp
        end if
        ! derivatives
        if (ang_type == 1 .or. ang_type == 3) then
           r23_sqr = ra2(1)*ra2(1) + ra2(2)*ra2(2) + ra2(3)*ra2(3) 
           r23_norm = sqrt(r23_sqr)
           df_int(:,1) = -p/p_norm_sqr*r23_norm
           df_int(:,4) =  q/q_norm_sqr*r23_norm
           df_int(:,2) = -df_int(:,1)*(1.0d0 + &
                (ra1(1)*ra2(1) + ra1(2)*ra2(2) + ra1(3)*ra2(3))/r23_sqr) &
                + df_int(:,4)* &
                (ra2(1)*ra3(1) + ra2(2)*ra3(2) + ra2(3)*ra3(3))/r23_sqr
           df_int(:,3) = -df_int(:,1)-df_int(:,4)-df_int(:,2)
           na = 4
        else
           sine = sin(angle)
           tmp = -pq*sine
           df_int(:,1) = (q - p/p_norm_sqr*p_dot_q)/tmp
           df_int(:,3) = (p - q/q_norm_sqr*p_dot_q)/tmp
           df_int(:,2) = -df_int(:,1)-df_int(:,3)
           na = 3
        end if
        df_int(:,1:na) = df0*df_int(:,1:na)
        Eint(ang_type) = Eint(ang_type) + fint
        do k = 1, na
           idx = ref_res_eng(i)%atm_idx(k,j)
           dfdr(:,idx) = dfdr(:,idx) + df_int(:,k)
        end do
     end do ! ang
  end do ! res
end subroutine Angle_energy
!-----------------------------------------------------------------------
subroutine angle_energy_only(Eint)
!-----------------------------------------------------------------------
  implicit none
  real(dp), dimension(3), intent(out) :: Eint
  integer :: ang_type, prm_no, na, i, j, k, idx
  integer, dimension(2) :: indx, cor_ang
  real(dp) :: angle, fint, sg, arg, one, pq, sine
  real(dp), dimension(3,4) :: df_int
  real(dp), dimension(3,3) :: ra
  real(dp) :: n, dw, kw, tmp, w0, tan_w, tan_sqr, tan1, cosw, sinw
  real(dp) :: r23_sqr, r23_norm, p_norm_sqr, q_norm_sqr, p_dot_q
  real(dp), dimension(3) :: p, q, s, ra1, ra2, ra3

  Eint(:) = 0.0d0
  do i = 1, protein%num_res
     do j = 1, ref_res_eng(i)%num_ang
        ang_type = ref_res_eng(i)%ang_type(j)
        do k = 1, ref_res_eng(i)%na_ang(j)-1  ! get bond vectors
           indx(1:2) = ref_res_eng(i)%atm_idx(k:k+1,j)
           ra(:,k) = R(:,indx(2))-R(:,indx(1))
        end do
        ! calc angle
        if (ang_type == 1 .or. ang_type ==3) then ! torsion angle
           ra1=ra(:,1); ra2=ra(:,2); ra3=ra(:,3)
           p(1)=ra1(2)*ra2(3)-ra1(3)*ra2(2)
           p(2)=ra1(3)*ra2(1)-ra1(1)*ra2(3)
           p(3)=ra1(1)*ra2(2)-ra1(2)*ra2(1)
           q(1)=ra2(2)*ra3(3)-ra2(3)*ra3(2)
           q(2)=ra2(3)*ra3(1)-ra2(1)*ra3(3)
           q(3)=ra2(1)*ra3(2)-ra2(2)*ra3(1)
           s(1)=ra3(2)*ra1(3)-ra3(3)*ra1(2)
           s(2)=ra3(3)*ra1(1)-ra3(1)*ra1(3)
           s(3)=ra3(1)*ra1(2)-ra3(2)*ra1(1)
           sg = s(1)*ra2(1) + s(2)*ra2(2) + s(3)*ra2(3)
        else if (ang_type == 2) then ! bond angle 
           p = -ra(:,1); q = ra(:,2)
           sg = 1.0d0
        end if
        p_dot_q = p(1)*q(1) + p(2)*q(2) + p(3)*q(3)
        p_norm_sqr = p(1)*p(1) + p(2)*p(2) + p(3)*p(3)
        q_norm_sqr = q(1)*q(1) + q(2)*q(2) + q(3)*q(3)
        pq = sqrt(p_norm_sqr*q_norm_sqr)
        arg = p_dot_q/pq
        one = 1.0
        arg = sign(min(abs(arg),one),arg) ! to be sure abs(arg)<=1
        angle = sign(acos(arg), sg)

        ! energy
        prm_no = ref_res_eng(i)%prm(j)
        w0 = eng_para%ang_para(2,prm_no)
        kw = eng_para%ang_para(1,prm_no)
        if (ang_type == 1) then ! proper torsion
           n = eng_para%ang_para(3,prm_no)
           tan_w = tan(0.5d0*(n*angle - w0))
           tan_sqr = tan_w * tan_w; tan1 = 1.0d0 + tan_sqr
           cosw = (1.0d0 - tan_sqr)/tan1
           sinw = 2.0d0*tan_w/tan1
           fint = kw*(1.0d0+cosw)
        else  ! bnd ang & improper torsion
           dw = angle - w0
           tmp = kw*dw
           fint = tmp*dw
        end if

        Eint(ang_type) = Eint(ang_type) + fint

     end do ! ang
  end do ! res
end subroutine Angle_energy_only
!-----------------------------------------------------------------------
subroutine bond_energy(Eb)
!-----------------------------------------------------------------------
! bond energy: not a constant because of bonds closing a ring
!-----------------------------------------------------------------------
  implicit none
  real(dp), intent(out) :: Eb
  integer :: i, j
  integer, dimension(2) :: indx
  real(dp), dimension(3) :: b, dEb
  real(dp) :: b_len, db, tmp
  integer :: prm_no, res_no, ang_no, ang_type, np, tnw, bnd_no
  
  Eb = 0.0d0
  do i = 1, protein%num_res
     do j = 1, ref_res(i)%num_bnd
        indx(1:2) = ref_res_eng(i)%b_idx(1:2,j)
        b = R(:,indx(2)) - R(:,indx(1))
        b_len = sqrt(b(1)*b(1)+ b(2)*b(2)+ b(3)*b(3))
        prm_no = ref_res_eng(i)%bnd_prm(j)
        db = b_len-eng_para%bnd_para(2,prm_no)
        tmp = eng_para%bnd_para(1,prm_no)*db
        Eb = Eb + tmp*db
        dEb = 2.0d0*tmp/b_len*b
        dfdr(:,indx(1)) = dfdr(:,indx(1)) - dEb
        dfdr(:,indx(2)) = dfdr(:,indx(2)) + dEb
     end do
  end do
! disulfide bond
  do i = 1, protein%num_ss_bnd
     indx(1:2) = ref_ss(i)%b_idx(1:2)
     b = R(:,indx(2)) - R(:,indx(1))
     b_len = sqrt(b(1)*b(1)+ b(2)*b(2)+ b(3)*b(3))
     prm_no = ss_ref_eng%bnd_prm(1)
     db = b_len-eng_para%bnd_para(2,prm_no)
     tmp = eng_para%bnd_para(1,prm_no)*db
     Eb = Eb + tmp*db
     dEb = 2.0d0*tmp/b_len*b
     dfdr(:,indx(1)) = dfdr(:,indx(1)) - dEb
     dfdr(:,indx(2)) = dfdr(:,indx(2)) + dEb
  end do
end subroutine bond_energy
!-----------------------------------------------------------------------
subroutine bond_energy_only(Eb)
!-----------------------------------------------------------------------
! bond energy: not a constant because of bonds closing a ring
!-----------------------------------------------------------------------
  implicit none
  real(dp), intent(out) :: Eb
  integer :: i, j
  integer, dimension(2) :: indx
  real(dp), dimension(3) :: b, dEb
  real(dp) :: b_len, db, tmp
  integer :: prm_no, res_no, ang_no, ang_type, np, tnw, bnd_no
  
  Eb = 0.0d0
  do i = 1, protein%num_res
     do j = 1, ref_res(i)%num_bnd
        indx(1:2) = ref_res_eng(i)%b_idx(1:2,j)
        b = R(:,indx(2)) - R(:,indx(1))
        b_len = sqrt(b(1)*b(1)+ b(2)*b(2)+ b(3)*b(3))
        prm_no = ref_res_eng(i)%bnd_prm(j)
        db = b_len-eng_para%bnd_para(2,prm_no)
        tmp = eng_para%bnd_para(1,prm_no)*db
        Eb = Eb + tmp*db
     end do
  end do
! disulfide bond
  do i = 1, protein%num_ss_bnd
     indx(1:2) = ref_ss(i)%b_idx(1:2)
     b = R(:,indx(2)) - R(:,indx(1))
     b_len = sqrt(b(1)*b(1)+ b(2)*b(2)+ b(3)*b(3))
     prm_no = ss_ref_eng%bnd_prm(1)
     db = b_len-eng_para%bnd_para(2,prm_no)
     tmp = eng_para%bnd_para(1,prm_no)*db
     Eb = Eb + tmp*db
  end do
end subroutine bond_energy_only
!-----------------------------------------------------------------------
subroutine calc_gradient(g, mode)
!-----------------------------------------------------------------------
! for Cartesian min, if mode == 0, 
! and for internal min, if mode == 1
  implicit none
  integer, intent(in) :: mode
  real(dp), dimension(:), intent(out) :: g
  real(dp) :: gi
  real(dp) :: Eb, norm, Rgr_sqr
  real(dp), dimension(3) :: r_pv, axis, dR, drdw, r1, r2, r_tr
  integer :: i, k, res_no, ang_no, ang_type, tnw, bnd_no, j
  integer, dimension(2) :: indx

  if (mode == 0) then ! CARTESIAN_MIN
     k = 0
     do i = 1, tot_num_atm
        g(k+1:k+3) = dfdr(1:3,i)
        k = k + 3
     end do
  else
     ! gradient wrt angles
     tnw = tot_num_ang
     do i = 1, tnw
        gi = 0.0d0
        r_pv = R(:,ipivot(i))
        ! axis of rotation
        res_no = i_w(1,i); ang_no = i_w(2,i)
        ang_type = ref_res(res_no)%ang_type(ang_no)
        if (ang_type == 1) then 
           indx(1) = ref_res(res_no)%a_idx(2,ang_no)        
           axis = protein%residue(res_no)%b(:,indx(1))
           norm = protein%residue(res_no)%b_len(indx(1))
        else 
           indx(1:2) = ref_res(res_no)%a_idx(1:2,ang_no)
           r1 = protein%residue(res_no)%b(:,indx(1))
           r2 = protein%residue(res_no)%b(:,indx(2))
           axis(1)=r1(3)*r2(2) - r1(2)*r2(3)
           axis(2)=r1(1)*r2(3) - r1(3)*r2(1)
           axis(3)=r1(2)*r2(1) - r1(1)*r2(2)
           norm = protein%residue(res_no)%b_len(indx(1))* &
                protein%residue(res_no)%b_len(indx(2))* &
                protein%residue(res_no)%sin_w(ang_no)
        end if
        if (invert_axis(i) == -1) axis = -axis
        axis=axis/norm
        do j = start_rot(i), end_rot(i)
           dR = R(:,j) - r_pv
           drdw(1)=axis(2)*dR(3)-axis(3)*dR(2)
           drdw(2)=axis(3)*dR(1)-axis(1)*dR(3)
           drdw(3)=axis(1)*dR(2)-axis(2)*dR(1)
           gi = gi + dfdr(1,j)*drdw(1) + dfdr(2,j)*drdw(2) + dfdr(3,j)*drdw(3)
        end do
        g(i) = gi
     end do
  endif
end subroutine calc_gradient
!-----------------------------------------------------------------------
subroutine pairlist_update(min)
  implicit none
  integer, intent(in) :: min
  integer :: i, j, gr_no, tnatm, tngr, m, gri, grj, update, np
  real(dp), dimension(3) :: Rgr, r1
  real(dp) :: Rgr_sqr
  integer, dimension(:,:), allocatable :: excl_gr

  tnatm = tot_num_atm 
  tngr = tot_num_gr

  allocate(excl_gr(tngr,tngr))

  R_gr(:,:) = 0.0d0
  do i = 1, tnatm
     R(:,i) = protein%residue(i_R(1,i))%R(:,i_R(2,i))
     gr_no = i_R(3,i)
     R_gr(:,gr_no) = R_gr(:,gr_no)+R(:,i)
  end do

  ! group coord
  do gr_no = 1, tngr
     R_gr(:,gr_no) = R_gr(:,gr_no)/Nagr(gr_no)
  end do

  ! make group cut-off list
  if (min==1) then
     update = 0
     if (N_update == 0) then
        update = 1
     else
        do i = 1, tngr
           r1 = R_gr(:,i) - R_gr_prev(:,i)
           if ((r1(1)*r1(1)+r1(2)*r1(2)+r1(3)*r1(3))> Rdiff_sqr) then
              update = 1
              exit
           end if
        end do
     end if
  else
     update = 1
  end if

  if (update == 1) then
     ! update list
     N_update = N_update + 1
     do i = 1, tngr
        r1 = R_gr(:,i)
        do j = i, tngr  
           Rgr = r1 - R_gr(:,j) ! group distance
           Rgr_sqr = Rgr(1)*Rgr(1) + Rgr(2)*Rgr(2) + Rgr(3)*Rgr(3)
           if (Rgr_sqr > Rout_sqr) then
              excl_gr(j,i) = 1; excl_gr(i,j) = 1
           else
              excl_gr(j,i) = 0; excl_gr(i,j) = 0
           end if
        end do
     end do

     do i = 1, tot_num_atm-1
        gri = i_R(3,i)
        m = 0
        do np = 1, num_pair(i)
           j = pair_list(np,i) 
           grj = i_R(3,j)
           if (excl_gr(gri,grj) == 0) then
              m = m + 1
              pair(m,i) = j 
           end if
        end do
        no_pair(i) = m
     end do

     if (min==1) R_gr_prev = R_gr
  end if
  deallocate(excl_gr)

end subroutine pairlist_update
!-----------------------------------------------------------------------
subroutine make_ang_list(option)
!-----------------------------------------------------------------------
! if option = 0, get total number of angles. If option = 1, get
! ang index i_w(1:2,i,p): res no & ang no for the new ang i, protein p
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: option
  type(reference_residue_type) :: ref
  integer :: res_no,i,j,ang_no

  ang_no = 0
  do res_no = 1, protein%num_res
     ref = ref_res(res_no)
     do i = 1, ref%num_ang
        if (ref%ang_fixed(i) == 0) then
           ang_no = ang_no + 1
           if (option == 1) i_w(1:2,ang_no)=(/ res_no, i /)
        end if
     end do
  end do
  if (option == 0) tot_num_ang = ang_no
end subroutine make_ang_list
!-----------------------------------------------------------------------
subroutine make_lists(mode)
!-----------------------------------------------------------------------
! get 1,4 neighbors, list of nonbonded atms, lists for gradient calc, 
! and LJ and ES atom types.
!-----------------------------------------------------------------------
! for Cartesian min, if mode == 0, 
! and for internal min, if mode == 1
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: mode
  type(reference_residue_type) :: ref
  type(ref_res_eng_type) :: eref
  integer :: num_res,loc_num,res_no,i,j,k,m,n,p,idx,n14,tnw
  integer :: coor_no,ii,jj,kk,mm,pp,tna,gr_no,atm_no, bb, del_branch
  integer :: ang_no,ang_type,imax,ip, sr, er, na_ang, bnd_no
  integer :: pivot_right,pivot_left, loc_CA, loc_C, add_branch
  real(dp), dimension(3,2) :: r
  integer, dimension(4) :: indx, aa_idx
  integer, dimension(4) :: atom, atm
  character(len=4) :: resname, atom_type, atm_type
  character(len=1) :: ter_type
  character(len=4), dimension(max_atm) :: loca, endc
  character(len=4), dimension(:), allocatable :: loca_exp
  integer, dimension(max_atm) :: ring
  integer, dimension(:,:), allocatable :: c_atm,i14,br1,br2,iR,loc_no,irr
  integer, dimension(:), allocatable :: end_no,ring_no,cn,n_br
  integer, dimension(2) :: ss_res, ss_atm

  num_res = protein%num_res
  tna=tot_num_atm
  allocate(c_atm(6,tna),end_no(tna),ring_no(tna),cn(tna),loca_exp(0:tna))
  allocate(i14(12,tna),br1(12,tna),br2(12,tna),iR(3,tna),n_br(tna))
  allocate(loc_no(-1:max_atm+1,num_res))
  if (mode == 0) then ! CARTESIAN_MIN
     allocate(irr(2,tot_num_dof))
  endif
! assign new number to all atoms in the protein
  loc_num=0; gr_no = 0 ! initialize location no. and group no.
  loca_exp(0) = '-C00'
  if (mode == 0) then ! CARTESIAN_MIN
     atm_no = 0
  endif
  do res_no = 1, num_res
     ref = ref_res(res_no); eref = ref_res_eng(res_no)
     resname = ref%residue_type
     if (res_no==1) then
        ter_type = 'N'
     else if (res_no==num_res) then
        ter_type = 'C'
     else
        ter_type = ' '
     end if
     call get_loca(resname,ter_type,loca,endc,ring)
     do i = 1, ref%num_atm
        ii=loc_num+i
        call find_atm(ref,loca(i),idx)
        loc_no(idx,res_no) = ii ! atom position
        iR(1:3,ii)=(/ res_no, idx, gr_no+eref%gr_no(idx) /)
        if (mode == 0) then ! CARTESIAN_MIN
           atm_no = atm_no + 1
           irr(1:2,atm_no)=(/ res_no, idx /)
        endif
        ring_no(ii) = ring(i) ! ring
        loca_exp(ii) = loca(i)
     end do
     loc_num = loc_num + ref%num_atm; gr_no = gr_no + eref%num_gr
     ! take care of -CA and -C
     if (res_no > 1) then
        loc_no(-1:0,res_no) = loc_no(2:3,res_no-1) 
     end if
  end do
  tot_num_gr = gr_no
! make bonded atm list, c_atm(j,i)
  loc_num=0 ! initialize location no.
  do res_no = 1, num_res
     ref = ref_res(res_no); eref = ref_res_eng(res_no)
     resname = ref%residue_type
     if (res_no==1) then
        ter_type = 'N'
     else if (res_no==num_res) then
        ter_type = 'C'
     else
        ter_type = ' '
     end if
     ! +N is N of the next residue
     if (res_no<num_res)loc_no(ref%num_atm+1,res_no)=loc_no(1,res_no+1)
     ! idx for new atm no in ang
!     print*,'ang in energy'
     do i = 1, eref%num_ang
        imax = eref%na_ang(i)
        if (eref%aa_idx(1,i)<1000) then
           ref_res_eng(res_no)%atm_idx(1:imax,i) = &
                loc_no(eref%aa_idx(1:imax,i),res_no)
        else
           atm(1:imax) = eref%aa_idx(1:imax,i) - 1000
           do j = 1, imax
              atom_type = ss_ref%atom_type(atm(j))
              atm_type = atom_type(2:)
              if (atom_type(1:1) == '1') then
                 call find_atm(ref,atm_type,idx)
                 ref_res_eng(res_no)%atm_idx(j,i) = loc_no(idx,res_no)
!                 print*,res_no,i,eref%ang_type(i),atom_type, &
!                      loca_exp(loc_no(idx,res_no))
              else
                 do k = 1, protein%num_ss_bnd
                    if (res_no == ref_ss(k)%res_no(1)) then
                       m = ref_ss(k)%res_no(2)
                       exit
                    end if
                 end do
                 call find_atm(ref_res(m),atm_type,idx)
                 ref_res_eng(res_no)%atm_idx(j,i) = loc_no(idx,m)
!                 print*,res_no,i,eref%ang_type(i),atom_type, &
!                      loca_exp(loc_no(idx,m))
              end if
           end do
        end if
!       print*,res_no,i,eref%ang_type(i),ref%atom_type(eref%aa_idx(1:imax,i)), &
!             ':', loca_exp(residue(res_no)%atm_idx(1:imax,i))
     end do
     if (mode == 1) then ! internal min
        ! angles in ref_res_eng -> corresponding ang in ref_res
        do i = 1, num_res
           do j = 1, ref_res_eng(i)%num_ang
              na_ang=ref_res_eng(i)%na_ang(j)
              indx(1:na_ang) = ref_res_eng(i)%aa_idx(1:na_ang,j)
              if (indx(1) >= 1000) then
                 ref_res_eng(i)%cor_ang(1,j) = 0
                 goto 102
              end if
              do k = 1, ref_res(i)%num_ang
                 if (ref_res(i)%na_ang(k)==na_ang) then
                    aa_idx(1) = ref_res(i)%b_idx &
                         (1,ref_res(i)%a_idx(1,k))
                    aa_idx(2:na_ang) = ref_res(i)%b_idx &
                         (2,ref_res(i)%a_idx(1:(na_ang-1),k))
                    do m = 1, na_ang
                       if (indx(m) /= aa_idx(m)) then
                          ref_res_eng(i)%cor_ang(1,j) = 0
                          goto 101
                       end if
                    end do
                    ref_res_eng(i)%cor_ang(1:2,j) = (/1,k/)
                    exit
101                 continue
                 end if
              end do
102           continue
           end do
        end do
     endif

     do i = 1, ref%num_bnd
        ref_res_eng(res_no)%b_idx(1:2,i)=loc_no(ref%b_idx(1:2,i),res_no)
     end do
     do i = 1, protein%num_ss_bnd
        ss_res(1:2) = ref_ss(i)%res_no(1:2)
        ss_atm(1:2) = ref_ss(i)%s_atm_no(1:2)
        ref_ss(i)%b_idx(1)=loc_no(ss_atm(1),ss_res(1))
        ref_ss(i)%b_idx(2)=loc_no(ss_atm(2),ss_res(2))
     end do
!     print*,'ang in energy done'
     call get_loca(resname,ter_type,loca,endc,ring)
     do i = 1, ref%num_atm
        ii=loc_num+i
        if (endc(i) == 'ENDC') then
           end_no(ii) = tna   ! end of branch
        else
           call find_atm(ref,endc(i),idx)
           end_no(ii)=loc_no(idx,res_no)
        end if
        coor_no=0 ! initialize no of coordinated atms
        do j = 1, ref%num_bnd
           indx(1:2) = loc_no(ref%b_idx(1:2,j),res_no)
           do k = 1, 2
              if (indx(k)==ii) then
                 coor_no = coor_no + 1
                 c_atm(coor_no,ii) = indx(3-k) ! bonded atm
              end if
           end do
        end do
        cn(ii) = coor_no ! no of bonded atm
     end do
     loc_num=loc_num+ref%num_atm
  end do
! find 1,4 neighbors and branches of nonbonded atms
  loc_num=0 ! initialize location no.
  do res_no = 1, num_res
     ref = ref_res(res_no)
     resname = ref%residue_type
     call find_atm(ref,'CA  ',idx)
     loc_CA = loc_no(idx,res_no) ! CA position
     call find_atm(ref,'C   ',idx)
     loc_C = loc_no(idx,res_no) ! C position
     do i = 1, ref%num_atm
        ii=loc_num+i !atm1
        n_br(ii)=0 ! no of 1,4 neighbors and branches
        do j = 1,cn(ii)
           jj=c_atm(j,ii) !atm2
           do k = 1,cn(jj)
              kk = c_atm(k,jj) !atm3
              if (kk/=ii) then 
                 do m = 1,cn(kk)
                    mm=c_atm(m,kk) !atm4
                    if (mm/=jj.and.mm>ii) then 
                       ! make sure that the list is not redundant
                       do n = 1, n_br(ii)
                          if (mm == br1(n,ii).and.i14(n,ii) == 1) goto 111
                       end do
                       if (ring_no(ii)==1.and.ring_no(mm)==1.and.iR(1,ii)==iR(1,mm)) then
                          ! if in the same ring, not a 1,4 neighbor
                          ! but the next atm could be start of the branch
                          do p = 1,cn(mm)
                             pp = c_atm(p,mm)
                             if (pp>ii.and.ring_no(pp)/=1) then
                                do n = 1, n_br(ii)
                                   if (pp == br1(n,ii)) goto 111
                                end do
                                n_br(ii)=n_br(ii)+1
                                i14(n_br(ii),ii)=0
                                br1(n_br(ii),ii)=pp
                                br2(n_br(ii),ii)=end_no(pp)
                                exit
                             end if
                          end do
                       else
                          n_br(ii)=n_br(ii)+1
                          i14(n_br(ii),ii)=1
                          br1(n_br(ii),ii)=mm
                          br2(n_br(ii),ii)=mm
                          if (end_no(mm)>mm) then
                             n_br(ii)=n_br(ii)+1
                             i14(n_br(ii),ii)=0
                             br1(n_br(ii),ii)=mm+1
                             br2(n_br(ii),ii)=end_no(mm)
                          end if
                       end if
                    end if
111                 continue
                 end do ! m (atm4)
              end if
           end do ! k (atm3)
        end do ! j (atm2)
        add_branch = 0 ! add branch for side chain atoms
        if (ii>loc_CA .and. ii<loc_C) then
           do j = 1, n_br(loc_CA)
              if (ii>=br1(j,loc_CA).and.ii<=br2(j,loc_CA)) add_branch = 1
           end do
        end if
        if (add_branch == 1) then
           n_br(ii) = n_br(ii)+1
           i14(n_br(ii),ii)=0
           br1(n_br(ii),ii)=loc_C
           br2(n_br(ii),ii)=end_no(loc_C)
        end if
     end do ! i (atm1)
     if (resname=='TRP') then ! special cases
        call find_atm(ref,'HE10',idx)
        ii = loc_no(idx,res_no)
        call find_atm(ref,'CE30',idx)
        do i = 1, n_br(ii)
           if (br1(i,ii) == loc_no(idx,res_no)) then
              call find_atm(ref,'CZ30',idx)
              br2(i,ii)=end_no(loc_no(idx,res_no))
           end if
        end do
     else if (resname=='HIS') then
        call find_atm(ref,'HD10',idx)
        ii = loc_no(idx,res_no)
        call find_atm(ref,'NE20',idx)
        del_branch = 0
        do i = 1, n_br(ii)
           if (del_branch == 1) then
              i14(i-1,ii) = i14(i,ii)
              br1(i-1,ii) = br1(i,ii); br2(i-1,ii) = br2(i,ii)
           end if
           if (br1(i,ii) == loc_no(idx,res_no) .and. i14(i,ii) == 0) then
              del_branch = 1
           end if
        end do
        n_br(ii) = n_br(ii) - 1
     end if
!     print*,'-------'
!     print*,res_no,resname
!     print*,'branches'
!     do i = 1, ref%num_atm
!        ii=loc_num+i
!        print*,iR(1,ii),loca_exp(ii), ii
!        do j=1,n_br(ii)
!            write(*,*)j,i14(j,ii),iR(1,br1(j,ii)),loca_exp(br1(j,ii)),br1(j,ii),':',&
!                iR(1,br2(j,ii)),loca_exp(br2(j,ii)),br2(j,ii)
!        end do
!     end do
     loc_num=loc_num+ref%num_atm
  end do ! res_no
!  n14 = 0
!  do i = 1, tna
!     do j = 1, n_br(i)
!        if (i14(j,i)==1) n14 = n14 + 1
!     end do
!  end do
!  write(*,*) 'Number of 1-4 pairs =',n14
  i_14(:,:) = -1
  do i = 1, tna-1
     m = 0
     do j = 1, n_br(i)
        do k = br1(j,i), br2(j,i)
           m = m + 1
           pair_list(m,i) = k
           i_14(k,i)=i14(j,i)
           i_14(i,k)=i14(j,i)
        end do
     end do
     num_pair(i) = m
  end do
!  num_br(1:tna)=n_br(1:tna)
!  i_14(:,1:tna)=i14(:,1:tna)
!  stt_br(:,1:tna)=br1(:,1:tna)
!  end_br(:,1:tna)=br2(:,1:tna)
  i_R(:,1:tna)=iR(:,1:tna)
  if (mode == 0) then ! CARTESIAN_MIN
     i_rr(:,1:tot_num_dof) = irr(:,1:tot_num_dof)
  endif
!-----------------------------------------------------------------------
  if (mode == 1) then ! internal min
     ! make lists for gradient calculation 
     tnw = tot_num_ang
     do ang_no = 1, tnw
        res_no = i_w(1,ang_no); i = i_w(2,ang_no)
        ref = ref_res(res_no)
        ang_type = ref%ang_type(i)
        imax = ref%na_ang(i)
        atom(1) = ref%b_idx(1,ref%a_idx(1,i))
        atom(2:imax) = ref%b_idx(2,ref%a_idx(1:(imax-1),i)) ! old numbering
        atm(1:imax) = loc_no(atom(1:imax),res_no)           ! new numbering
        bb = 0 ! not a backbone angle
        if (atom(imax)==1.or.atom(imax)==2.or.atom(imax)==3.or. &
             atom(imax)==ref%num_atm+1) bb=1 ! bb ang
        if (ang_type == 1) then
           pivot_right = atm(3); pivot_left = atm(2)
        else if (ang_type == 2) then
           pivot_right = atm(2); pivot_left = atm(2)
        end if
        if (include_ang_type == 3 .or.  include_ang_type == 4) then
           if ((end_no(atm(3))-atm(3))<=(end_no(atm(2)+1)-1)) then
              ipivot(ang_no) = pivot_right
              start_rot(ang_no) = atm(3)
              end_rot(ang_no) = end_no(atm(3))
              invert_axis(ang_no) = 1
           else
              ipivot(ang_no) = pivot_left
              start_rot(ang_no) = 1
              end_rot(ang_no) = end_no(atm(2)+1)
              invert_axis(ang_no) = -1
           end if
        else if (include_ang_type == 0 .and. ref%fixed_branch(i) /= 0) then
           if (bb==1) then
              if ((end_no(atm(3))-atm(3))<=(end_no(atm(2)+1)-1)) then
                 ipivot(ang_no) = pivot_right
                 start_rot(ang_no) = atm(3)
                 end_rot(ang_no) = end_no(atm(3))
                 invert_axis(ang_no) = 1
              else
                 ipivot(ang_no) = pivot_left
                 start_rot(ang_no) = 1
                 end_rot(ang_no) = end_no(atm(2)+1)
                 invert_axis(ang_no) = -1
              end if
           else 
              ipivot(ang_no) = pivot_right
              start_rot(ang_no) = atm(3)
              end_rot(ang_no) = end_no(atm(3))
              invert_axis(ang_no) = 1
           end if
        else
           if (bb==0.or.(end_no(atm(imax))-atm(imax))<=(atm(imax)-2)) then
              ipivot(ang_no) = pivot_right
              start_rot(ang_no) = atm(imax)
              end_rot(ang_no) = end_no(atm(imax))
              invert_axis(ang_no) = 1
           else
              ipivot(ang_no) = pivot_left
              start_rot(ang_no) = 1
              end_rot(ang_no) = atm(imax)-1
              invert_axis(ang_no) = -1
           end if
        end if
!     ip = ipivot(ang_no); sr = start_rot(ang_no); er = end_rot(ang_no)
!     print*,res_no,ang_no,ang_type,invert_axis(ang_no), &
!          loca_exp(atm(1:imax)),':',iR(1,ip),loca_exp(ip), &
!          ':',iR(1,sr),loca_exp(sr),':',iR(1,er),loca_exp(er)
     end do
  endif
!-----------------------------------------------------------------------
!get LJ type and ES type
  do i = 1, tna
     res_no = i_R(1,i); atm_no = i_R(2,i)
     eref = ref_res_eng(res_no)
     LJ_type(i) = eref%atm_cls(atm_no)
     ES_type(i) = eref%q_idx(atm_no)
  end do
!-----------------------------------------------------------------------
  deallocate(c_atm,i14,br1,br2,end_no,ring_no,cn,n_br,loca_exp,iR,loc_no)
  if (mode == 0) then ! CARTESIAN_MIN
     deallocate(irr)
  endif
end subroutine make_lists
!-----------------------------------------------------------------------
end MODULE energy
!-----------------------------------------------------------------------
