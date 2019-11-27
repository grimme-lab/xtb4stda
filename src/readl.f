c     *****************************************************************

      subroutine readl(a1,x,n)
      implicit real*8 (a-h,o-z)
      character*(*) a1
      dimension x(*)
      i=0
      is=1
  10  i=i+1
      x(i)=readaa(a1,is,ib,ie)
      if(ib.gt.0 .and. ie.gt.0) then
                                is=ie
                                goto 10
      endif
      n=i-1
      return
      end


      function readaa(a,istart,iend,iend2)
      implicit real*8 (a-h,o-z)
      real*8 readaa
      character*(*) a

      NINE=ICHAR('9')
      IZERO=ICHAR('0')
      MINUS=ICHAR('-')
      IDOT=ICHAR('.')
      ND=ICHAR('D')
      NE=ICHAR('E')
      IBL=ICHAR(' ')

      iend=0
      iend2=0
      idig=0
      c1=0
      c2=0
      one=1.d0
      x = 1.d0
      nl=len(a)
      do 10 j=istart,nl-1
         n=ichar(a(j:j))
         m=ichar(a(j+1:j+1))
         if(n.le.nine.and.n.ge.izero .or.n.eq.idot)goto 20
         if(n.eq.minus.and.(m.le.nine.and.m.ge.izero
     1 .or. m.eq.idot)) goto 20
   10 continue
      readaa=0.d0
      return
   20 continue
      iend=j
      do 30 i=j,nl
         n=ichar(a(i:i))
         if(n.le.nine.and.n.ge.izero) then
            idig=idig+1
            if (idig.gt.10) goto 60
            c1=c1*10+n-izero
         elseif(n.eq.minus.and.i.eq.j) then
            one=-1.d0
         elseif(n.eq.idot) then
            goto 40
         else
            goto 60
         endif
   30 continue
   40 continue
      idig=0
      do 50 ii=i+1,nl
         n=ichar(a(ii:ii))
         if(n.le.nine.and.n.ge.izero) then
            idig=idig+1
            if (idig.gt.10) goto 60
            c2=c2*10+n-izero
            x = x /10
         elseif(n.eq.minus.and.ii.eq.i) then
            x=-x
         else
            goto 60
         endif
   50 continue
c
c put the pieces together
c
   60 continue
      readaa= one * ( c1 + c2 * x)
      do 55 j=iend,nl
         n=ichar(a(j:j))
         iend2=j
         if(n.eq.ibl)return
   55 if(n.eq.nd .or. n.eq.ne)goto 57
      return

   57 c1=0.0d0
      one=1.0d0
      do 31 i=j+1,nl
         n=ichar(a(i:i))
         iend2=i
         if(n.eq.ibl)goto 70
         if(n.le.nine.and.n.ge.izero) c1=c1*10.0d0+n-izero
         if(n.eq.minus)one=-1.0d0
   31 continue
   61 continue
   70 readaa=readaa*10**(one*c1)
      return
      end
