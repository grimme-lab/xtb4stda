
c principal quantum number of valence shell

      integer function pqn(at)
      integer at

      if(at.le.2)then
         pqn=1
      elseif(at.le.10)then
         pqn=2
      elseif(at.le.18)then
         pqn=3
      elseif(at.le.36)then
         pqn=4
      elseif(at.le.54)then
         pqn=5
      else
         pqn=6
      endif

      end

      integer function ncore(at)
      integer at

      if(at.le.2)then
         ncore=0
      elseif(at.le.10)then
         ncore=2
      elseif(at.le.18)then
         ncore=10
      elseif(at.le.29)then   !zn
         ncore=18
      elseif(at.le.36)then
         ncore=28
      elseif(at.le.47)then
         ncore=36
      elseif(at.le.54)then
         ncore=46
      elseif(at.le.71)then
         ncore=54
      elseif(at.le.79)then
         ncore=68
      elseif(at.le.86)then
         ncore=78
      endif

      end
