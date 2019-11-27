subroutine blckmgs(m,n,ndim,darray)
!----------------------------------------------------------------------
! Purpose:
! Subroutine performs modified Gramm-Schmidt orthonormalization of a
! real matrix. Orthonormalization is done in-place, so the darray is
! overwritten on exit. Linearly dependent vectors are set to zero.
!
! Input:
! m      - Number of rows in the matrix darray
! n      - Number of columns in the matrix darray
! ndim   - First array dimension as declared in the calling routine
! darray - Array to be orthonormalized
!
! Output:
! darray - Orthonormalized array
!----------------------------------------------------------------------
  implicit none

! Input/Output:
  integer, intent(in) :: m,n,ndim
  real*8, dimension(ndim,n), intent(out) :: darray
!----------------------------------------------------------------------
! Local:
  integer :: ii,jj,kk,ll,ibsize,nblcks,istrt,jstrt,iend,ncol,ierr
  real*8 :: tmp
  real*8, dimension(:,:), allocatable :: smat
!----------------------------------------------------------------------
! Local parameters
!csg  real*8, parameter :: thr = epsilon(1.0d0) !Threshold for zero vectors
  real*8 thr

! External BLAS functions
  real*8 :: ddot
!----------------------------------------------------------------------

  thr = epsilon(1.0d0) !Threshold for zero vectors
!----------------------------------------------------------------------
! Block size optimized for Athlon 1200 MHz with 2.0GB memory for
! matrices up to 5000x5000
  ibsize = 60
!----------------------------------------------------------------------

!-----------------------------------------------------------------
! Allocate overlap matrix
!-----------------------------------------------------------------
  allocate(smat(ibsize,ibsize),stat=ierr)
  if(ierr /= 0)  stop 'Mamory allocation error in blckmgs'

!-----------------------------------------------------------------
! Calculate the number of blocks
!-----------------------------------------------------------------
  nblcks = (n+ibsize-1)/ibsize
  ibsize = min(n,ibsize)

!-----------------------------------------------------------------
! Orthogonalize the first block using modified schmidt
!-----------------------------------------------------------------
  do ii=1,ibsize

     tmp = ddot(m,darray(1,ii),1,darray(1,ii),1)

! Linear dependence
     if(tmp < thr) then
        darray(1:m,ii) = 0.0d0
        cycle
     end if

     tmp = 1.0d0/sqrt(tmp)
     call dscal(m,tmp,darray(1,ii),1)

     do jj=ii+1,ibsize
        tmp = ddot(m,darray(1,ii),1,darray(1,jj),1)
        call daxpy(m,-tmp,darray(1,ii),1,darray(1,jj),1)
     end do

  end do

!-----------------------------------------------------------------
! Loop over remaining blocks
!-----------------------------------------------------------------
  do ii=1,nblcks-1

! Initial and final column and number of columns in the block ii+1
     istrt = ii*ibsize+1
     iend  = (ii+1)*ibsize
     iend  = min(n,iend)
     ncol  = iend - istrt + 1

! Orthogonalize the block ii+1 against the previous ones
     do jj=1,ii

! Initial index of the block jj
        jstrt = (jj-1)*ibsize+1

        call dgemm('t','n',ibsize,ncol,m,1.0d0,darray(1,jstrt),ndim,darray(1,istrt),ndim,0.0d0,smat,ibsize)
        call dgemm('n','n',m,ncol,ibsize,-1.0d0,darray(1,jstrt),ndim,smat,ibsize,1.0d0,darray(1,istrt),ndim)

     end do

! Othogonalize vectors on the block ii+1 among themself using modified schmidt
     do kk=istrt,iend

        tmp = ddot(m,darray(1,kk),1,darray(1,kk),1)

! Linear dependence
        if(tmp < thr) then
           darray(1:m,kk) = 0.0d0
           cycle
        end if

        tmp = 1.0d0/sqrt(tmp)
        call dscal(m,tmp,darray(1,kk),1)

        do ll=kk+1,iend
           tmp = ddot(m,darray(1,kk),1,darray(1,ll),1)
           call daxpy(m,-tmp,darray(1,kk),1,darray(1,ll),1)
        end do

     end do

  end do

! Clean up
  deallocate(smat,stat=ierr)
  if(ierr /= 0) stop 'Mamory deallocation error in blckmgs'

end subroutine blckmgs
