      subroutine fodenmak(uhf,nmo,eps,occ,efermi)
      implicit none
      integer nmo
      logical uhf
      real*8 efermi
      real*8 occ(*)
      real*8 eps(*)

      integer i
      real*8 inte

      inte=2.0d0
      if(uhf) inte=1.0d0

      do i=1,nmo
         if(eps(i)*27.2113957.le.efermi) then
            occ(i)=inte-occ(i)
         endif
      enddo

      end
