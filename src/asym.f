! This file is part of xtb4stda.
!
! Copyright (C) 2015-2019 Stefan Grimme
!
! xtb4stda is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb4stda is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb4stda.  If not, see <https://www.gnu.org/licenses/>.

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
