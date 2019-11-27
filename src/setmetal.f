c which atoms are metals?
      subroutine setmetal
      implicit none
      include 'aoelementcommon.fh'

      metal=1
      metal(1:2)  =0
      metal(6:10) =0
      metal(14:18)=0
      metal(32:36)=0
      metal(50:54)=0
      metal(82:86)=0

      end

