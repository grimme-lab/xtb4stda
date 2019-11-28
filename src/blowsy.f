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

      SUBROUTINE BLOWSY(ITY,A,B,N)
C
C     BLOW UP SYMMETRIC OR ANTISYMMETRIC MATRIX TO FULL SIZE
      REAL*8 A(*),B(N,N)
C
C     DETERMINE IF WE HAVE AN ANTISYMMETRIC INTEGRAL

      IF (ITY.EQ.-1) GOTO 99
      IJ=0
      DO 1 I=1,N
      DO 2 J=1,I-1
      IJ=IJ+1
      B(J,I)=A(IJ)
2     B(I,J)=A(IJ)
      IJ=IJ+1
      B(I,I)=A(IJ)
1     CONTINUE
      RETURN
99    IJ=0
      DO 11 I=1,N
      DO 12 J=1,I-1
      IJ=IJ+1
      B(J,I)=-A(IJ)
12    B(I,J)=A(IJ)
      IJ=IJ+1
      B(I,I)=0.D0
11    CONTINUE
      RETURN
      END

