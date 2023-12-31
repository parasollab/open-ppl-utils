!-----------------------------------------------------------------------
MODULE rmsd
 use globals
CONTAINS
!-----------------------------------------------------------------------
subroutine calc_rmsd(protein1, protein2, rmsd_type, &
     option, U, error)
  implicit none
  integer, intent(in) :: option
  type(protein_type), intent(in) :: protein1, protein2
  character(len=2), intent(in) :: rmsd_type
  real(dp), dimension(:,:), intent(out) :: U
  real(dp), intent(out) :: error
  real(dp), dimension(:,:), allocatable :: coord1, coord2
  integer :: num_atm1, num_atm2, num_atm

  if (rmsd_type == 'bb') then
     num_atm = 3*protein1%num_res
  else if (rmsd_type == 'aa') then
     num_atm = sum(ref_res(1:protein1%num_res)%num_atm)
  else
     write(*,*)'Unknown rmsd type.'
     stop
  end if
  allocate(coord1(3,num_atm), coord2(3,num_atm))
  call get_coord(protein1,coord1,num_atm1,rmsd_type,1) 
  call get_coord(protein2,coord2,num_atm2,rmsd_type,1) 
  if (num_atm1 /= num_atm2) then
     write(*,*)'Error in rmsd. Number of atoms in two proteins are not the same'
     stop
  end if
  call ls_rmsd(num_atm1,coord1,coord2,option,U,error)
  deallocate(coord1, coord2)
  return
end subroutine calc_rmsd
!-----------------------------------------------------------------------
subroutine get_coord(protein0,coord,num_atm,coor_type,option)
!-----------------------------------------------------------------------
! option=0: get coord for atoms in the original pdb file
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: option
  integer, intent(out) :: num_atm
  character(len=2) :: coor_type
  character(len=4) :: atom
  type(protein_type), intent(in) :: protein0
  real(dp), dimension(:,:), intent(out) :: coord
  integer :: i, j, index
  index = 0
  if (coor_type == 'aa') then ! all atom coor
     do i = 1, protein0%num_res
        do j = 1, ref_res(i)%num_atm
           if (option/=0.or. &
                (option==0.and.protein0%residue(i)%atom_pdb(j)==1)) then
              index = index + 1
              coord(:,index)=protein0%residue(i)%R(:,j)
           end if
        end do
     end do
  else if (coor_type == 'bb') then ! backbone coor
     do i = 1, protein0%num_res
        do j = 1, ref_res(i)%num_atm
           atom = ref_res(i)%atom_type(j)
           if (atom=='N000'.or.atom=='CA00'.or.atom=='C000') then 
              if (option/=0.or. &
                   (option==0.and.protein0%residue(i)%atom_pdb(j)==1)) then
                 index = index + 1
                 coord(:,index)=protein0%residue(i)%R(:,j)
              end if
           end if
        end do
     end do
  else
     write(*,*) 'Unknown coordinate type:',coor_type
     stop
  end if
  num_atm = index

end subroutine get_coord
!-----------------------------------------------------------------------
subroutine ls_rmsd(n,coord1,coord2,option,U,error)
!-----------------------------------------------------------------------
!  This subroutine calculates the least square rmsd of the two coordinate
!  sets coord1(3,n) and coord2(3,n) using a method based on quaternion.
!  If option=1, then the geometric centers and best-fit rotation matrix
!  are returned.
!-----------------------------------------------------------------------
implicit none
   integer, intent(in) :: n, option
   real(dp), dimension(:,:), intent(in) :: coord1, coord2
   real(dp), dimension(:,:), intent(out) :: U
   real(dp), intent(out) :: error
   real(dp), dimension(3) :: x_center, y_center
   integer :: i, j
   real(dp), dimension(3,n) :: x, y
   real(dp), dimension(n) :: xi, yi
   real(dp) :: x_norm, y_norm, lambda
   real(dp), dimension(3,3) :: Rmatrix
   real(dp), dimension(4,4) :: S
   real(dp), dimension(4) :: q
! make copies of the original coordinates
      x(:,1:n)    = coord1(:,1:n);  y(:,1:n)      = coord2(:,1:n)
! calculate the barycenters, centroidal coordinates, and the norms
      x_norm      = 0.0d0;            y_norm      = 0.0d0
   do i = 1, 3
      xi          = x(i,:);           yi          = y(i,:)
      x_center(i) = sum(xi)/dble(n);  y_center(i) = sum(yi)/dble(n)
      xi          = xi - x_center(i); yi          = yi - y_center(i)
      x(i,:)      = xi;               y(i,:)      = yi
      x_norm = x_norm + dot_product(xi, xi)
      y_norm = y_norm + dot_product(yi, yi)
   end do
! calculate the R matrix
   do i = 1, 3
      do j = 1, 3
         Rmatrix(i,j) = dot_product(x(i,:),y(j,:))
      end do
   end do
! lower triangular part of the S matrix
   S(1, 1) = Rmatrix(1, 1) + Rmatrix(2, 2) + Rmatrix(3, 3)
   S(2, 1) = Rmatrix(2, 3) - Rmatrix(3, 2);  S(1, 2) = S(2, 1)
   S(2, 2) = Rmatrix(1, 1) - Rmatrix(2, 2) - Rmatrix(3, 3)
   S(3, 1) = Rmatrix(3, 1) - Rmatrix(1, 3);  S(1, 3) = S(3, 1)
   S(3, 2) = Rmatrix(1, 2) + Rmatrix(2, 1);  S(2, 3) = S(3, 2)
   S(3, 3) =-Rmatrix(1, 1) + Rmatrix(2, 2) - Rmatrix(3, 3)
   S(4, 1) = Rmatrix(1, 2) - Rmatrix(2, 1);  S(1, 4) = S(4, 1)
   S(4, 2) = Rmatrix(1, 3) + Rmatrix(3, 1);  S(2, 4) = S(4, 2)
   S(4, 3) = Rmatrix(2, 3) + Rmatrix(3, 2);  S(3, 4) = S(4, 3)
   S(4, 4) =-Rmatrix(1, 1) - Rmatrix(2, 2) + Rmatrix(3, 3) 
! Calculate eigenvalues and eigenvectors. For explanation about arguments,
! see the end of this subroutine.
   call dstmev(S, lambda, q)
   if (option .eq. 1) then
      call rotation_matrix(q, U)
   end if
! lambda is the largest eigenvalue
   error = sqrt(max(0.0d0,(0.5d0*(x_norm+y_norm)-lambda))/dble(n))
return
end subroutine ls_rmsd
!-----------------------------------------------------------------------
subroutine rotation_matrix(q, U)
!-----------------------------------------------------------------------
! This subroutine constructs rotation matrix U from quaternion q.
!-----------------------------------------------------------------------
implicit none
   real(dp), dimension(:), intent(in) :: q
   real(dp), dimension(:,:), intent(out) :: U
   real(dp) :: q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12,q13,q22,q23,q33
   q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4)
   b0 = 2.0d0*q0; b1 = 2.0d0*q1; b2 = 2.0d0*q2; b3 = 2.0d0*q3
   q00 = b0*q0-1.0d0; q01 = b0*q1; q02 = b0*q2; q03 = b0*q3
   q11 = b1*q1; q12 = b1*q2; q13 = b1*q3  
   q22 = b2*q2; q23 = b2*q3; q33 = b3*q3 
   U(1,1) = q00+q11; U(1,2) = q12-q03; U(1,3) = q13+q02
   U(2,1) = q12+q03; U(2,2) = q00+q22; U(2,3) = q23-q01
   U(3,1) = q13-q02; U(3,2) = q23+q01; U(3,3) = q00+q33
return
end subroutine rotation_matrix
!-----------------------------------------------------------------------
subroutine DSTMEV(A,lambda,evec)
implicit none
!-----------------------------------------------------------------------
! a simple subroutine to compute the leading eigenvalue and eigenvector
! of a symmetric, traceless 4x4 matrix A by an inverse power iteration:
! (1) the matrix is converted to tridiagonal form by 3 Givens
! rotations;  V*A*V' = T
! (2) Gershgorin's theorem is used to estimate a lower
! bound for the leading negative eigenvalue:
! lambda_1 > g=min(T11-t12,-t21+T22-t23,-t32+T33-t34,-t43+T44)
!          =
! where tij=abs(Tij)
! (3) Form the positive definite matrix 
!     B = T-gI
! (4) Use svd (algorithm svdcmp from "Numerical Recipes")
!     to compute eigenvalues and eigenvectors for SPD matrix B
! (5) Shift spectrum back and keep leading singular vector
!     and largest eigenvalue.
! (6) Convert eigenvector to original matrix A, through 
!     multiplication by V'.  
!-----------------------------------------------------------------------
real(dp), dimension(4,4):: A, T, V, SV
integer                         :: i
integer, dimension(1)           :: max_loc
real(dp), dimension(4)  :: evec, SW
real(dp), dimension(8)  :: rv1
real(dp)                :: lambda
!-----------------------------------------------------------------------
!(I).   Convert to tridiagonal form, keeping similarity transform
! (a product of 3 Givens rotations)
call givens4(A,T,V)
!-----------------------------------------------------------------------
!(II)Estimate lower bound of smallest eigenvalue by Gershgorin's theorem
lambda=min(T(1,1)-abs(T(1,2)),-abs(T(2,1))+T(2,2)-abs(T(2,3)),&
          -abs(T(3,2))+T(3,3)-abs(T(3,4)),-abs(T(4,3))+T(4,4))
!-----------------------------------------------------------------------
!(III). Form positive definite matrix     T <== lambda*I - T
do i = 1,4
   T(i,i) = T(i,i)-lambda
enddo
!-----------------------------------------------------------------------
!(IV). Compute singular values/vectors of SPD matrix B
call svdcmp(4,T,4,4,SW,SV,rv1)
!-----------------------------------------------------------------------
!(V). Shift spectrum back
max_loc = maxloc(SW) 
lambda = SW(max_loc(1)) + lambda
!lambda = SW(1) + lambda
!-----------------------------------------------------------------------
!(VI). Convert eigenvector to original coordinates: (V is transposed!)
evec = matmul(V,SV(:,max_loc(1)))
!write(*,*)'-----------------------------------------------------------'
!write(*,*) 'lambda = ', lambda,'  eigenvector:  '
!write(*,99) evec
!write(*,*)'-----------------------------------------------------------'
!-----------------------------------------------------------------------
99 format(1x,4(d19.13,1x))
end subroutine dstmev
!-----------------------------------------------------------------------
subroutine givens4(S,T,V)
!-----------------------------------------------------------------------
implicit none
real(dp), dimension(4,4), intent(in)  :: S
real(dp), dimension(4,4), intent(out) :: T,V
real(dp)  :: c1,c2,c3, s1,s2,s3, r1,r2,r3, c1c2, s1c2
!real(dp) :: pythag
! external        pythag
!performs givens rotations to reduce symmetric 4x4 matrix to tridiagonal
T=S; V = 0.d0
!-----------------------------------------------------------------------
!Zero out entries T(4,1) and T(1,4)
! compute cos and sin of rotation angle in the 3-4 plane
r1 = pythag(T(3,1),T(4,1))
if(r1 .ne. 0.d0) then
  c1 = T(3,1)/r1; s1 = T(4,1)/r1
  V(3,3) = c1  ; V(3,4) = s1
  V(4,3) =-s1  ; V(4,4) = c1
  T(3,1) = r1  ; T(4,1) = 0.d0
  T(3:4,2:4) = matmul(V(3:4,3:4),T(3:4,2:4))
  T(1:2,3:4) = transpose(T(3:4,1:2))
  T(3:4,3:4) = matmul(T(3:4,3:4),transpose(V(3:4,3:4)))
else
  c1 = 1.d0; s1 = 0.d0
endif
!-----------------------------------------------------------------------
!Zero out entries T(3,1) and T(1,3)
! compute cos and sin of rotation angle in the 2-3 plane
r2 = pythag(T(3,1), T(2,1))
if(r2 .ne. 0.d0) then
  c2 = T(2,1)/r2; s2 = T(3,1)/r2
  V(2,2) = c2  ; V(2,3) = s2
  V(3,2) =-s2  ; V(3,3) = c2
  T(2,1) = r2  ; T(3,1) = 0.d0
  T(2:3,2:4) = matmul(V(2:3,2:3),T(2:3,2:4))
  T(1,2:3)   = T(2:3,1);  T(4,2:3) = T(2:3,4)
  T(2:3,2:3) = matmul(T(2:3,2:3),transpose(V(2:3,2:3)))
else
  c2 = 1.d0; s2 = 0.d0
endif
!-----------------------------------------------------------------------
!Zero out entries T(4,2) and T(2,4)
! compute cos and sin of rotation angle in the 3-4 plane
r3 = pythag(T(4,2), T(3,2))
if(r3 .ne. 0.d0) then
  c3 = T(3,2)/r3; s3 = T(4,2)/r3
  V(3,3) = c3  ; V(3,4) = s3
  V(4,3) =-s3  ; V(4,4) = c3
  T(3,2) = r3  ; T(4,2) = 0.d0
  T(3:4,3:4) = matmul(V(3:4,3:4),T(3:4,3:4))
  T(1:2,3:4) = transpose(T(3:4,1:2))
  T(3:4,3:4) = matmul(T(3:4,3:4),transpose(V(3:4,3:4)))
else
  c3 = 1.d0; s3 = 0.d0
endif
!-----------------------------------------------------------------------
!Compute net rotation matrix (accumulate similarity for evec. computation)
! To save transposing later, This is the transpose!
V(1,1)=1.d0; V(1,2:4) = 0.d0; V(2:4,1) = 0.d0
V(2,2) = c2;  V(3,2) = c1*s2 ; V(4,2) = s1*s2; c1c2 = c1*c2; s1c2=s1*c2
V(2,3) = -s2*c3 ; V(3,3) = c1c2*c3-s1*s3 ; V(4,3) =  s1c2*c3+c1*s3
V(2,4) =  s2*s3 ; V(3,4) =-c1c2*s3-s1*c3 ; V(4,4) = -s1c2*s3+c1*c3
!-----------------------------------------------------------------------
!write(*,*) (V(1:4,i) - W(1:4,i),i=1,4)
end subroutine givens4
!-----------------------------------------------------------------------
SUBROUTINE svdcmp(mmax,a,m,n,w,v,rv1)
!-----------------------------------------------------------------------
  implicit none
  integer ::    mmax
  INTEGER ::  m,n
  real(dp) :: a(mmax,*),v(mmax,*),w(*),rv1(*)
!     USES pythag
  INTEGER :: i,its,j,jj,k,l,nm
  real(dp) :: anorm,c,f,g,h,s,scale,x,y,z, one!,pythag
!      external pythag
!11---------------------------------------------------------------------
      g = 0.0d0
      scale = 0.0d0
      anorm = 0.0d0
      do i = 1, n
        l=i+1
        rv1(i)=scale*g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if(i.le.m)then
          do k=i,m
            scale=scale+abs(a(k,i))
          end do
          if(scale.ne.0.0d0)then
            do k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
            end do
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do j=l,n 
              s=0.0d0
              do k=i,m
                s=s+a(k,i)*a(k,j)
              end do
              f=s/h
              do k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
              end do
            end do
            do k=i,m
              a(k,i)=scale*a(k,i)
            end do
          endif
        endif
        w(i)=scale *g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if((i.le.m).and.(i.ne.n))then
          do k=l,n
            scale=scale+abs(a(i,k))
          end do
          if(scale.ne.0.0d0)then
            do k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
            end do
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do k=l,n
              rv1(k)=a(i,k)/h
            end do
            do j=l,m
              s=0.0d0
              do k=l,n
                s=s+a(j,k)*a(i,k)
              end do
              do k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
              end do
            end do
            do k=l,n
              a(i,k)=scale*a(i,k)
            end do
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
      end do
!85---------------------------------------------------------------------
      do i = n, 1, -1
        if(i .lt. n) then
          if(g.ne.0.0d0)then
            do j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
            end do
            do j=l,n
              s=0.0d0
              do k=l,n
                s=s+a(i,k)*v(k,j)
              end do
              do k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
              end do
            end do
          endif
          do j=l,n
            v(i,j)=0.0d0
            v(j,i)=0.0d0
          end do
        endif
        v(i,i)=1.0d0
        g=rv1(i)
        l=i
      end do
!-----------------------------------------------------------------------
      do i = min(m,n), 1, -1
        l=i+1
        g=w(i)
        do j=l,n
          a(i,j)=0.0d0
        end do
        if(g.ne.0.0d0)then
          g=1.0d0/g
          do j=l,n
            s=0.0d0
            do k=l,m
              s=s+a(k,i)*a(k,j)
            end do
            f=(s/a(i,i))*g
            do k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
            end do        
          end do
          do j=i,m
            a(j,i)=a(j,i)*g
          end do
        else
          do j= i,m
            a(j,i)=0.0d0
          end do
        endif
        a(i,i)=a(i,i)+1.0d0
      end do
!-----------------------------------------------------------------------
      do k=n,1,-1
        do its=1,30
          do l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
          end do
1         c=0.0d0
          s=1.0d0
          do i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0d0/h
            c= (g*h)
            s=-(f*h)
            do j=1,m     
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
            end do
          end do
2         z=w(k)
          if(l .eq. k)then
            if(z.lt.0.0d0)then
              w(k)=-z
              do j=1,n
                v(j,k)=-v(j,k)
              end do
            endif
            goto 3
          endif
          if(its.eq.30) then
              write(*,*) 'no convergence in svdcmp'
              stop
          endif
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
          one = 1.0
          g=pythag(f,one)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0d0
          s=1.0d0
          do j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g    
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
            end do
            z=pythag(f,h)
            w(j)=z
            if(z.ne.0.0d0)then
              z=1.0d0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
            end do
          end do
          rv1(l)=0.0d0
          rv1(k)=f       
          w(k)=x
        end do
3       continue
      end do
!-----------------------------------------------------------------------
      return
    END subroutine svdcmp
!-----------------------------------------------------------------------
    FUNCTION pythag(a,b)
!-----------------------------------------------------------------------
      implicit none
      real(dp) :: a,b, pythag
      real(dp) :: absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*dsqrt(1.0d0+(absb/absa)**2)
      else
        if(absb.eq.0.0d0)then
          pythag=0.0d0
        else
          pythag=absb*dsqrt(1.0d0+(absa/absb)**2)
        endif
      endif
      return
    END function pythag
!-----------------------------------------------------------------------
end MODULE rmsd
!-----------------------------------------------------------------------
