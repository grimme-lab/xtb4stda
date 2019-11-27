C     *****************************************************************

      FUNCTION ASYM(I)
      CHARACTER*2 ASYM
      CHARACTER*2 ELEMNT(107), AS
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
      AS=ELEMNT(I)
      CALL UPPER(AS)
      ASYM=AS
      if(i.eq.103) asym='XX'
      RETURN
      END

      SUBROUTINE UPPER(AS)
      CHARACTER*2 AS
      NSP=ICHAR(' ')
      ND=ICHAR('A')-ICHAR('a')
      DO 10 I=1,2
         J=ICHAR(AS(I:I))
         IF(J.NE.NSP)J=J+ND
         AS(I:I)=CHAR(J)
  10  CONTINUE
      RETURN
      END

      SUBROUTINE UPPER10(AS)
      CHARACTER*(*) AS
      NSP=ICHAR(' ')
      ND=ICHAR('A')-ICHAR('a')
      DO 10 I=1,10
         J=ICHAR(AS(I:I))
         IF(J.NE.NSP)J=J+ND
         AS(I:I)=CHAR(J)
  10  CONTINUE
      RETURN
      END
