c aoat,lao,fila in CAO basis
c aoat2,lao2,fila2 in AO basis

      subroutine dtrafo2(n,nbf,nao)
      implicit none
      integer n,nbf,nao
      include 'ehtcommon.fh'
      integer i,j,be,en

c     do i=1,n
c        write(*,*) '6d ',fila(1:2,i)
c     enddo

      j=0
      do i=1,nbf
         if(lao(i).ne.5)then
            j=j+1
            lao2(j)=lao(i)
            if(lao(i).gt.4) lao2(j)=lao2(j)-1
            aoat2 (j)=aoat (i)
            valao2(j)=valao(i)
            hdiag2(j)=hdiag(i)
         endif
      enddo

      do j=1,n
         be=100000
         en=-1
         do i=1,nao
            if(aoat2(i).eq.j.and.i.lt.be)be=i
            if(aoat2(i).eq.j.and.i.gt.en)en=i
         enddo
         fila2(1,j)=be
         fila2(2,j)=en
      enddo

c     do i=1,n
c        write(*,*) '5d ',fila(1:2,i)
c     enddo
c     write(*,'(''6d '',40i3)') lao (1:nbf)
c     write(*,'(''5d '',40i3)') lao2(1:nbf)


      end
