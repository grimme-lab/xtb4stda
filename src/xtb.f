ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine xtb(n,at,xyz,q,z,cn,unsat,
     .               nel,nopen,ndim,eel,wb,kspd,gscal,tfac,kcn,
     .               lshift,mothr,
     .               pop,dipol,wbo,molden,wfnout,lmo,diff,iter,et)
      use gbobc, only: lgbsa,init_gbsa,compute_fgb,gsasa,ghb,lhb,
     .                 compute_gb_egrad,update_nnlist_gbsa

      implicit none
      integer n,at(n),nel,nopen,ndim,iter
      real*8 xyz(3,n),eel,wb(n,n),q(n),z(n),gscal,kspd(6),lshift,tfac,et
      real*8 cn(n),kcn,unsat(n)
      real*8 mothr,tt,nfoda,nfodb
      real*8 R27
      data R27/27.21139570d0/
      logical pop,dipol,wbo,lmo
      logical molden,wfnout,diff

      include 'ehtcommon.fh'
      include 'aoelementcommon.fh'

      real*8 ,allocatable ::X(:,:)
      real*8 ,allocatable ::P(:,:)
      real*8 ,allocatable ::S(:,:)
      real*8 ,allocatable ::H(:,:)
      real*8 ,allocatable ::T(:)
      real*8 ,allocatable ::Ptmp(:,:)
      real*8 ,allocatable ::Xcao(:,:)
      real*8 ,allocatable ::jab(:,:)
      real*8 ,allocatable ::aux(:)
      real*8 ,allocatable ::emo(:)
      real*8 ,allocatable ::focc(:),focca(:),foccb(:)
      real*8 ,allocatable ::xcen(:)
      real*8 ,allocatable ::pc(:,:)
      real*8 ,allocatable ::qpc(:)
      real*8 ,allocatable ::jpc(:,:)
      real*8 ,allocatable ::qmull(:)
      real*8 ,allocatable ::tmp(:)
      real*8 ,allocatable ::chir(:),cm5_true(:),cm5_appr(:)
      real*8 ,allocatable ::qlmom(:,:)
      real*8 ,allocatable ::fdum(:)
      integer,allocatable ::apc(:)
      integer,allocatable ::iwork(:)

      real*8 temp,xsum,eh1,xgab,rab,eold,dum,kmagic(4,4),gscalpc,xx(20)
      real*8 t0,t1,sum,rr,dip,scal,egap,hav,wbr(n,n),ga,gb
      real*8 efa,efb
      real*8 dipo(3)
      integer LWORK,LIWORK,INFO,IU,IFOUND
      integer npr,ii,jj,kk,i,j,m,k,iat,jat,lina,ibmax,ncao,jter
      integer mmm(20),ihomo,ishell,jshell,nao,np,ia,ndimv
      integer jmet,ihomoa,ihomob,ll,i1,i2,nn,lin
      data    mmm/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/
      character*1 flag
      character*2 asym
      character*128 atmp,ftmp
      logical ex
c OMP
      common /proc/ nproc
      integer nproc

c GBSA stuff
      real*8, allocatable :: fgb(:,:)
      real*8, allocatable :: fhb(:)
      real*8 :: gsolv,gborn,tgb,t8
c     for the CM5 charges
      real*8, allocatable :: cm5(:)
      real*8, allocatable :: fgba(:)

c the special J scaling for point charges
      gscalpc=1.0d0

      write(*,*)
c     write(*,*)'Natom : ',n
c     write(*,*)'Nel   : ',nel
      write(*,*)'Nbf   : ',ndim
      write(*,'('' T(el) : '',f7.1)')et


      np=0
      inquire(file='pcharge',exist=ex)
      if(ex)then
            open(unit=22,file='pcharge')
            write(*,*)'reading point charges ...'
            read(22,*) np
            allocate(pc(3,np),qpc(np),apc(np),jpc(np,n))
            do i=1,np
               read(22,*) qpc(i),pc(1:3,i),atmp
               call elem(atmp, j)
               apc(i)=j
            enddo
            close(22)
            write(*,*)'Npc   : ',np
      endif

c note: H is in eV!

      call cpu_time(t0)
      write(*,*)'calculating S/T integrals ...'
      allocate(T(ndim*(ndim+1)/2),X(ndim,ndim))
      call xstints(n,ndim,xyz,X,T)
      call cpu_time(t1)
      T = T * R27
      write(*,*) 'cpu time for ints ',t1-t0

c     call prmat(6,T,ndim,ndim,'T')
c     call prmat(6,X,ndim,ndim,'S')
c     stop

      call cpu_time(t0)
      ncao=ndim ! keep # of Cartesians for later
      call cao2sao(ndim,nao,X)

c d-functions exist
      if(nao.ne.ndim)then
        call cao2saop(ndim,nao,T)
        open(unit=22,file='xtm.temp',form='unformatted')
        write(22) X(1:nao,1:nao)
        write(22) T(1:nao*(nao+1)/2)
        deallocate(T,X)
        ndim=nao
      endif
      write(*,*)'Nao   : ',ndim

      liwork = 3 + 5*ndim
      lwork  = 1 + 6*ndim + 2*ndim**2
      liwork=8*ndim

      if(nao.ne.ncao)then
         allocate(H(ndim,ndim),aux(lwork),P(ndim,ndim),
     .            jab(n,n),iwork(liwork),emo(ndim),
     .            focca(ndim),foccb(ndim),
     .            focc(ndim),S(ndim,ndim),tmp(ndim),
     .            T(ndim*(ndim+1)/2),X(ndim,ndim))
         rewind(22)
         read(22) X(1:nao,1:nao)
         read(22) T(1:nao*(nao+1)/2)
         close(22,status='delete')
c        relocate fila,lao,hdiag
         call dtrafo2(n,ncao,nao)
      else
         lao2=lao
         aoat2=aoat
         valao2=valao
         fila2=fila
         hdiag2=hdiag
         allocate(H(ndim,ndim),aux(lwork),P(ndim,ndim),
     .            jab(n,n),iwork(liwork),emo(ndim),
     .            focca(ndim),foccb(ndim),
     .            focc(ndim),S(ndim,ndim),tmp(ndim))
      endif

      S = X

      allocate(qmull(n))

c     initialize the GBSA module (GBSA works with CM5 charges)
      if(lgbsa) then
        allocate(fgb(n,n),fhb(n),cm5(n))
        gsolv=0.d0
        gborn=0.d0
        gsasa=0.d0
        ghb=0.d0
        cm5=0.0d0
c       initialize the neighbor list
        call update_nnlist_gbsa(n,xyz)
c       initialize the fgb matrix (dielectric screening of the Coulomb potential)
        call compute_fgb(n,xyz,fgb,fhb)
c       initialize the CM5 charges computation
        call docm5(n,at,.false.,xyz,q,cm5)
      endif


c fill levels
      call occ(ndim,nel,nopen,ihomo,focc)

      do i=1,n
c atoms (no J scaling)
         do j=1,i
            jab(j,i)=xgab(n,xyz(1,j),xyz(1,i),at(j),at(i))*0.5d0
            jab(i,j)=jab(j,i)
          enddo
c point charges
          do j=1,np
            jpc(j,i)=xgab(n,pc(1,j),xyz(1,i),apc(j),at(i))*0.5d0*gscalpc
          enddo
      enddo


c iter counter for second SCC
      write(*,*)'making H0...'
      jter=0
      eold=0

      do i=1,4
         do j=1,4
            kmagic(i,j)=(kspd(i)+kspd(j))*0.5
         enddo
      enddo

      do ii=1,ndim
         iat=aoat2(ii)
         hdiag2(ii)=hdiag2(ii)*(1.0d0+0.01*kcn*cn(iat))
      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

4242  do ii=1,ndim
         ishell=mmm(lao2(ii))
         do jj=1,ii-1
           jshell=mmm(lao2(jj))
c          scaled <T>
           tt=-T(lin(jj,ii))*tfac
           hav=0.5*(hdiag2(ii)+hdiag2(jj))
c          Ryd-Ryd
           if(valao2(ii).eq.0.and.valao2(jj).eq.0) then
              H(jj,ii)=S(jj,ii)*kspd(5)*hav
           elseif (valao2(ii).eq.0.and.valao2(jj).ne.0) then
c          Ryd-val
              H(jj,ii)=S(jj,ii)*kspd(6)*hav
           elseif (valao2(ii).ne.0.and.valao2(jj).eq.0) then
c          Ryd-val
              H(jj,ii)=S(jj,ii)*kspd(6)*hav
           else
              H(jj,ii)=tt+S(jj,ii)*kmagic(jshell,ishell)*hav
           endif
           H(ii,jj)=0
           X(ii,jj)=0
         enddo
         H(ii,ii)=hdiag2(ii)
      enddo

      write(*,*)'making H1...'

      do i=1,ndim
         ii=aoat2(i)
         do j=1,i
            jj=aoat2(j)
            dum=S(j,i)
            if(abs(dum).gt.1.d-6)then
            eh1=0.0d0
            do kk=1,n
               eh1=eh1+q(kk)*(jab(kk,ii)+jab(kk,jj))
            enddo
            do kk=1,np
               eh1=eh1+qpc(kk)*(jpc(kk,ii)+jpc(kk,jj))
            enddo
            H(j,i)=H(j,i)-dum*eh1*gscal
            endif
         enddo
      enddo

c     add the gbsa SCC term
      if(lgbsa) then
       if(iter.eq.1.and.jter.lt.1) then ! only in first iteration, use CM5 charges from EN charges
         continue
       else
         cm5=q ! in all other cases, CM5 charges are already used
       endif
c      hbpow=2.d0*c3-1.d0
       do i=1,ndim
         ii=aoat2(i)
         do j=1,i
         jj=aoat2(j)
         dum=S(j,i)
c        GBSA SCC terms
         eh1=0.0d0
         do kk=1,n
            eh1=eh1+cm5(kk)*(fgb(kk,ii)+fgb(kk,jj))
         enddo
         t8=(fhb(ii)*cm5(ii)+fhb(jj)*cm5(jj))
         tgb=-dum*(0.5d0*eh1+t8)
         H(j,i)=H(j,i)+tgb
!         H(i,j)=H(j,i)
         enddo
       enddo
      endif


      call cpu_time(t1)
      write(*,*) 'cpu time for H    ',t1-t0

c     call prmat(6,x,ndim,ndim,'S')
c     call prmat(6,h,ndim,ndim,'H')

      call cpu_time(t0)
      write(*,*)'solving ...'

      call DSYGVD(1,'V','U',ndim,H,ndim,X,ndim,emo,aux,
     .            LWORK,IWORK,LIWORK,INFO)
      if(info.ne.0) stop 'diag error'
      call cpu_time(t1)
      write(*,*) 'cpu time for diag ',t1-t0
      write(*,*)

c     shift orbital energies in xTB runs (these are used in Fermi smearing)
      if(iter.eq.0)
     .  call shifteps(.false.,ndim,focc,focc,emo,emo,
     .                lshift,0.0d0,0.0d0)

      eel=0
      do i=1,ndim
         emo(i)=emo(i)/R27
         eel=eel+emo(i)*focc(i)
      enddo
      if(ihomo+1.le.ndim)then
      egap=(emo(ihomo+1)-emo(ihomo))*R27
      write(*,'('' gap (eV)           : '',f7.3)')
     .egap
      else
      egap=0
      endif

      write(*,'('' Koopmans IP (eV)   : '',F7.3 )')-emo(ihomo)*R27
      write(*,'('' Eel (Eh)           : '',F12.6)')eel

      if(egap.lt.0.2)then
         write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write(*,*)'!!!!!!!! WARNING: small HL-gap detected !!!!!!!!!'
         write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      endif

c     make density matrix and related stuff
      if(pop.or.dipol.or.wbo)then
c       Fermi smearing
        if(et.gt.0.1.and.ihomo+1.le.ndim)then
c       convert restricted occ first to alpha/beta
        call occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
        call FERMISMEAR
     .       (.true.,ndim,ihomoa,et,emo*R27,focca,nfoda,efa,ga)
        call FERMISMEAR
     .       (.true.,ndim,ihomob,et,emo*R27,foccb,nfodb,efb,gb)
        write(*,'('' Nfod               :'',f8.4)')nfoda+nfodb
        write(*,'('' Nfod/Nel           :'',f8.4)')
     .                                  (nfoda+nfodb)/float(nel)
        focc = focca + foccb
        endif
c restricted AO spin density
        if(nopen.gt.0)then
        tmp=0.0d0
        ndimv=0
        do k=1,ndim
         if(valao2(k).eq.1)ndimv=ndimv+1
         do i=1,ndim
            do j=1,ndim
               tmp(i)=tmp(i)+S(j,i)*H(j,k)*H(i,k)*(focca(k)-foccb(k))
            enddo
         enddo
        enddo
        open(unit=111,file='tmpxtb')
        write(111,*)ndimv
        do i=1,ndim
         if(valao2(i).eq.1)write(111,*)tmp(i)
        enddo
        close(111)
        endif
c     density matrix
        call dmat(ndim,focc,H,P)
      endif

      call preig(6,focc,27.213957d0,emo,
     .           max(ihomo-12,1),min(ihomo+11,ndim))

      write(*,*) '                ',jter+1,' SCC done.'
c     second SCC
      if(pop.and.iter.eq.1.and.jter.lt.1)then

      call mpop0(n,ndim,aoat2,S,P,qmull)

      qmull = z - qmull
      call docm5(n,at,.false.,xyz,qmull,q)
      X = S
      jter = jter +1
      goto 4242
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC END OF SCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      deallocate(T,X)

cccccccccccccccccccccccccccccccccccccccccccccc
c         localization
cccccccccccccccccccccccccccccccccccccccccccccc

      if(lmo) call local(n,at,ncao,ndim,ihomo,xyz,z,focc,S,H,emo)

cccccccccccccccccccccccccccccccccccccccccccccc
c         Population analysis
cccccccccccccccccccccccccccccccccccccccccccccc
      if(pop)then
        allocate(chir(n),cm5_true(n),cm5_appr(n),qlmom(3,n))

c       fitting stuff for charge TB
        inquire(file='charges',exist=ex)
        if(ex)then
         open(unit=22,file='charges')
         j=0
         do i=1,n
            j=j+1
            read(22,*,end=120) chir(i)
         enddo
120      close(22)
         if(j.ne.n) then
            ex=.false.
            goto 130
         endif
         call docm5(n,at,.false.,xyz,chir,cm5_true)
        endif
130     continue

c       Mulliken pop
        call mpop(n,ndim,aoat2,lao2,S,P,qmull,qlmom)

        q = z - qmull
        call docm5(n,at,.false.,xyz,q,cm5_appr)

        write(*,*)
        write(*,*)
        write(*,*) 'Mulliken/CM5 charges    n(s)   n(p)   n(d)   '
        do i=1,n
           write(*,'(i5,a3,2f8.4,1x,4f7.3)')
     .     i,asym(at(i)),q(i),cm5_appr(i),qlmom(1:3,i)
        enddo
        qmull=q      ! save

c       fitting stuff for charge TB
        if(ex)then
c        fitting stuff
         call execute_command_line('pwd > tmpxtb2')
         open(unit=43,file='tmpxtb2')
         read(43,'(a)') ftmp
         close(43,status='delete')
         open(unit=43,file='charges3')
         do i=1,n
           write(43,'(2F12.6,2x,a)')
     .     cm5_true(i)*10,cm5_appr(i)*10,trim(ftmp)
         enddo
         close(43)
        else
c       for QMDFF
         write(*,*) 'Hirshfeld charges for QMDFF written to <charges>'
         open(unit=43,file='charges')
         do i=1,n
           write(43,*) qmull(i)  ! QMDFF uses Hirshfeld to get CM5
         enddo
         close(43)
        endif

c       Loewdin pop
        allocate(Ptmp(ndim,ndim))
        call makel(ndim, S, H, Ptmp)
        q=0
        qlmom=0
        call lpop(n,ndim,aoat2,lao2,focc,Ptmp,1.0d0,q,qlmom)
        q = z - q
        write(*,*)
        write(*,*) 'Loewdin charges    n(s)   n(p)   n(d)'
        do i=1,n
           write(*,'(i5,a3,f8.4,1x,4f7.3)')
     .     i,asym(at(i)),q(i),qlmom(1:3,i)
        enddo

        if(diff)then
        write(*,*) 'because of diff. basis fct. Loewdin populations'
        write(*,*) 'are written to file <charges3>'
        call docm5(n,at,.false.,xyz,q,cm5_appr)
        open(unit=43,file='charges3')
        do i=1,n
           write(43,'(2F12.6,2x,a)')
     .     cm5_true(i)*10,cm5_appr(i)*10,trim(ftmp)
        enddo
        close(43)
        endif

        if(et.gt.1000)then
           allocate(fdum(ndim))
           fdum = focc
           call fodenmak(.false.,ndim,emo,fdum,efa)
           q=0
           qlmom=0
           call lpop(n,ndim,aoat2,lao2,fdum,Ptmp,1.0d0,q,qlmom)
           write(*,*)
           write(*,*) 'Loewdin FODpop     n(s)   n(p)   n(d)'
           open(unit=43,file='fod')
           do i=1,n
              write(*,'(i5,a3,f8.4,1x,4f7.3)')
     .        i,asym(at(i)),q(i),qlmom(1:3,i)
              write(43,'(F14.8)') q(i)
           enddo
           deallocate(fdum)
           write(*,*) 'FODpop written to file <fod>'
           close(43)
        endif
        deallocate(Ptmp)

c OUTPUT!
        q = cm5_appr

        inquire(file='orca.out',exist=ex)
        if(ex)then
         wbr=0
         open(unit=43,file='orca.out')
  10     read(43,'(a)',end=100) atmp
         if(index(atmp,'Mayer bond orders larger').ne.0) then
  20       read(43,'(a)',end=100) atmp
           do ll=1,len(atmp)
            if(atmp(ll:ll).eq.',') atmp(ll:ll)=' '
           enddo
           call readl(atmp,xx,nn)
           k=nn/3
           kk=0
           do i=1,k
            kk=kk+1
            i1=idint(xx(kk))+1
            if(i1.gt.n) then
               wbr=0
               goto 100
            endif
            kk=kk+1
            i2=idint(xx(kk))+1
            if(i2.gt.n) then
               wbr=0
               goto 100
            endif
            kk=kk+1
            wbr(i1,i2)=xx(kk)
            wbr(i2,i1)=xx(kk)
           enddo
           if(nn.gt.0) goto 20
         endif
         goto 10
 100     close(43)
        endif

      endif

cccccccccccccccccccccccccccccccccccccccccccccc
c           dipole moment
cccccccccccccccccccccccccccccccccccccccccccccc
      if(dipol)then
         call dipole(n,ncao,ndim,xyz,z,P,dip,dipo,.true.)
          inquire(file='dipole',exist=ex)
          if(ex)then
            open(unit=22,file='dipole')
            read(22,*)dum
            close(22)
            if(dip.gt.1.d-5.and.dum.gt.1.d-5)then
              open(unit=22,file='dipole2')
              write(22,*)2.5*dum,dip*2.5
              close(22)
            endif
          endif
      endif

cccccccccccccccccccccccccccccccccccccccccccccc
c            Wiberg BOs
cccccccccccccccccccccccccccccccccccccccccccccc
      if(wbo)
     .call wiberg(n,ndim,at,xyz,P,S,wbr,wb,ex)

cccccccccccccccccccccccccccccccccccccccccccccc
c write molden file
c only write MOs with eps < mothr Eh
cccccccccccccccccccccccccccccccccccccccccccccc

      if(wfnout)then
      write(*,*) 'writing mo output ...'
      write(*,*) 'molden style :',molden
! we need cartesian AOs, so transform X back to CAO basis
      if (ncao.ne.ndim) then
         allocate(Xcao(ncao,ndim))
         call sao2cao(ndim,H,ncao,Xcao)
         if(molden)then
         if(et.gt.1000)call fodenmak(.false.,ndim,emo,focc,efa)
         call printmold(n,ndim,ncao,xyz,at,Xcao,emo,focc,
     .                  mothr)
         else
         call printmos (n,ndim,ncao,xyz,at,Xcao,emo,focc,
     .                  mothr)
         endif
      else
         if(molden)then
         if(et.gt.1000)call fodenmak(.false.,ndim,emo,focc,efa)
         call printmold(n,ndim,ndim,xyz,at,H,emo,focc,
     .                  mothr)
         else
         call printmos (n,ndim,ndim,xyz,at,H,emo,focc,
     .                  mothr)
         endif
      endif
      endif

      if(lgbsa) deallocate(cm5,fgb,fhb)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real*8 function xgab(n,xyz1,xyz2,ati,atj)
      implicit none
      integer i,j,n,ati,atj
      real*8 xyz1(3),xyz2(3),rab
      include 'aoelementcommon.fh'

      rab=sqrt((xyz1(1)-xyz2(1))**2
     .        +(xyz1(2)-xyz2(2))**2
     .        +(xyz1(3)-xyz2(3))**2)

c makes no big difference, the one taken is sligthly better for charges

c     xgab=27.2113957/(rab+2./(gam(ati)+gam(atj)))

      xgab=27.2113957/(rab+1.00d0/sqrt(gam(ati)*gam(atj)) )

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine occ(ndim,nel,nopen,ihomo,focc)
      implicit none
      integer nel,nopen,ndim,ihomo
      real*8 focc(ndim)
      integer i,na,nb

      focc=0
c even nel
      if(mod(nel,2).eq.0)then
      ihomo=nel/2
      do i=1,ihomo
         focc(i)=2.0d0
      enddo
      if(2*ihomo.ne.nel) then
         ihomo=ihomo+1
         focc(ihomo)=1.0d0
         if(nopen.eq.0)nopen=1
      endif
      if(nopen.gt.1)then
         do i=1,nopen/2
            focc(ihomo-i+1)=focc(ihomo-i+1)-1.0
            focc(ihomo+i)=focc(ihomo+i)+1.0
         enddo
      endif
c odd nel
      else
      na=nel/2+(nopen-1)/2+1
      nb=nel/2-(nopen-1)/2
      do i=1,na
         focc(i)=focc(i)+1.
      enddo
      do i=1,nb
         focc(i)=focc(i)+1.
      enddo
      endif

      do i=1,ndim
         if(focc(i).gt.0.99) ihomo=i
      enddo

      end

