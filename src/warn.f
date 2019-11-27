      subroutine warn(s)
      implicit none
      character*(*) s

      write(*,'(10x,65(''!''))')
      write(*,'(10x,65(''!''))')
      write(*,'(10x,''WARNING:'',a)')trim(s)
      write(*,'(10x,65(''!''))')
      write(*,'(10x,65(''!''))')

      end

      subroutine warn2(s1,s2)
      implicit none
      character*(*) s1,s2

      write(*,'(10x,45(''!''))')
      write(*,'(10x,45(''!''))')
      write(*,'(10x,a)')trim(s1)
      write(*,'(10x,a)')trim(s2)
      write(*,'(10x,45(''!''))')
      write(*,'(10x,45(''!''))')

      end
