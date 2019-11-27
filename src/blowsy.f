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

