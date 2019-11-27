CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ELEM(KEY1, NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) KEY1
      CHARACTER*2 ELEMNT(107),E

      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/

      nat=0
      e='  '
      do i=1,len(key1)
         if(key1(i:i).ne.' ') L=i
      enddo
      k=1
      DO J=1,L
         if (k.gt.2)exit
         N=ICHAR(key1(J:J))
         if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
            e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
            k=k+1
         endif
         if(n.ge.ichar('a') .and. n.le.ichar('z') )then
            e(k:k)=key1(j:j)
            k=k+1
         endif
      enddo

      DO I=1,107
         if(e.eq.elemnt(i))then
            NAT=I
            RETURN
         ENDIF
      ENDDO

      end

