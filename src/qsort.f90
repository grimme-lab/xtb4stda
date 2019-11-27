recursive subroutine qsort(a, first, last, ind)
  implicit none
  real*8 a(*), x, t
  integer ind(*)
  integer first, last
  integer i, j, ii

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     ii=ind(i); ind(i) = ind(j);  ind(j) = ii
     i=i+1
     j=j-1
  end do
  if (first < i-1) call qsort(a, first, i-1, ind)
  if (j+1 < last)  call qsort(a, j+1, last, ind)
end subroutine qsort

recursive subroutine qqsort(a, first, last)
  implicit none
  real*8 a(*), x, t
  integer first, last
  integer i, j, ii

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call qqsort(a, first, i-1)
  if (j+1 < last)  call qqsort(a, j+1, last)
end subroutine qqsort
