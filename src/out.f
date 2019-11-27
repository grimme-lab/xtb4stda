      subroutine outvip(eneut,ekat,shift)
      implicit none
      real*8 eneut,ekat,shift
      real*8 ip,dum
      logical ex

      write(*,*)
      write(*,*)'vertical deltaSCC IP calculation'
      ip=ekat-eneut-shift
      write(*,'(''empirical IP/EA shift (eV):'',f10.4)')
     .            27.21139570d0*shift
      write(*,'(''delta SCC IP (eV)'',f10.4)') 27.21139570d0*ip

c compare with reference for fit
      inquire(file='.ip',exist=ex)
      if(ex)then
         open(unit=11,file='.ip')
         read(11,*) dum
         close(11)
         open(unit=12,file='.IP')
         write(12,'(2F12.6)') dum,ip*27.21139570d0
         close(12)
      endif

      end

      subroutine outvea(eneut,eani,shift)
      implicit none
      real*8 eneut,eani,shift
      real*8 ea,dum
      logical ex

      write(*,*)
      write(*,*)'vertical deltaSCC EA calculation'
      ea=eneut-eani-shift
      write(*,'(''empirical IP/EA shift (eV):'',f10.4)')
     .            27.21139570d0*shift
      write(*,'(''delta SCC EA (eV):'',f10.4)') 27.21139570d0*ea

c compare with reference for fit
      inquire(file='.ea',exist=ex)
      if(ex)then
         open(unit=11,file='.ea')
         read(11,*) dum
         close(11)
         open(unit=12,file='.EA')
         write(12,'(2F12.6)') dum,ea*27.21139570d0
         close(12)
      endif

      end
