      subroutine atovlp(l,npri,nprj,alpa,alpb,conta,contb,ss)
      implicit none
      integer l,npri,nprj
      real*8 alpa(*),alpb(*)
      real*8 conta(*),contb(*)
      real*8 ss

      integer ii,jj
      real*8 ab,s00,sss,pi,ab05
      data pi/3.1415926535897932384626433832795029d0/

              SS=0.0d0
              do ii=1,npri
                 do jj=1,nprj
                    ab =1./(alpa(ii)+alpb(jj))
                    s00=(pi*ab)**1.50d00
                    if(l.eq.0)then
                       sss=s00
                    endif
                    if(l.eq.1)then
                       ab05=ab*0.5
                       sss=s00*ab05
                    endif
                    SS=SS+SSS*conta(ii)*contb(jj)
                 enddo
              enddo

      end
