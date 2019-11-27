ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine uxtb(n,at,xyz,q,z,cn,sat,
     .               nel,nopen,ndim,eel,wb,kspd,gscal,tfac,kcn,
     .               lshift,lshifta,lshiftb,mothr,pop,dipol,wbo,
     .               molden,wfnout,lmo,diff,iter,et,epr)
      use gbobc, only: lgbsa,init_gbsa,compute_fgb,gsasa,ghb,lhb,
     .                 compute_gb_egrad,update_nnlist_gbsa
      implicit none
      integer n,at(n),nel,nopen,ndim,iter
      real*8 xyz(3,n),eel,wb(n,n),q(n),z(n),gscal,kspd(6),lshift,tfac
      real*8 cn(n),kcn,sat(n),et
      real*8 mothr,tt,lshifta,lshiftb
      logical pop,dipol,wbo
      logical molden,wfnout,diff
      logical epr,lmo

      include 'ehtcommon.fh'
      include 'aoelementcommon.fh'

      real*8 ,allocatable ::X(:,:)
      real*8 ,allocatable ::P(:,:),Pa(:,:),Pb(:,:),aospin(:)
      real*8 ,allocatable ::S(:,:)
      real*8 ,allocatable ::H(:,:)
      real*8 ,allocatable ::T(:)
      real*8 ,allocatable ::Ptmp(:,:)
      real*8 ,allocatable ::Xcaoa(:,:),Xcaob(:,:)
      real*8 ,allocatable ::jab(:,:)
      real*8 ,allocatable ::aux(:)
      real*8 ,allocatable ::Ha(:,:),Hb(:,:)
      real*8 ,allocatable ::emoa(:),emob(:)
      real*8 ,allocatable ::focca(:),foccb(:)
      real*8 ,allocatable ::pc(:,:)
      real*8 ,allocatable ::qpc(:)
      real*8 ,allocatable ::jpc(:,:)
      real*8 ,allocatable ::qmull(:)
      real*8 ,allocatable ::chir(:),cm5_true(:),cm5_appr(:)
      real*8 ,allocatable ::qlmom(:,:)
      real*8 ,allocatable ::fduma(:),fdumb(:)
      real*8 ,allocatable ::vector(:)
      integer,allocatable ::apc(:)
      integer,allocatable ::iwork(:)

      real*8 temp,xsum,eh1,xgab,rab,eold,dum,kmagic(4,4),gscalpc
      real*8 t0,t1,sum,rr,dip,scal,egapa,egapb,hav,f,ps,shift,ga,gb
      real*8 s2,nfoda,nfodb,wbr(n,n),xx(20),efa,efb,hfc(n),totspin
      real*8 dipo(3)
      integer LWORK,LIWORK,INFO
      integer npr,ii,jj,kk,i,j,m,k,iat,jat,lina,ibmax,ncao,jter,ib
      integer mmm(20),ihomoa,ihomob,ishell,jshell,nao,np,ia,mo
      integer jmet,lin,kshell,ndimv,ll,i1,i2,nn,maxiter
      data    mmm/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/
      character*1 flag
      character*2 asym
      character*128 atmp,ftmp
      logical ex,ex2

      real*8 R27
      data R27/27.21139570d0/
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

      maxiter=2

      totspin = float(nopen)*0.5

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
      if(nao.ne.ncao)then
         allocate(H(ndim,ndim),aux(lwork),Pa(ndim,ndim),Pb(ndim,ndim),
     .            jab(n,n),iwork(liwork),emoa(ndim),emob(ndim),
     .            focca(ndim),foccb(ndim),S(ndim,ndim),P(ndim,ndim),
     .            T(ndim*(ndim+1)/2),X(ndim,ndim),aospin(ndim),
     .            Ha(ndim,ndim),Hb(ndim,ndim))
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
         allocate(H(ndim,ndim),aux(lwork),Pa(ndim,ndim),Pb(ndim,ndim),
     .            jab(n,n),iwork(liwork),emoa(ndim),emob(ndim),
     .            focca(ndim),foccb(ndim),S(ndim,ndim),P(ndim,ndim),
     .            Ha(ndim,ndim),Hb(ndim,ndim),aospin(ndim))
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
      call occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)

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

c no spin polarization in first iter
      aospin=0
c read AO spin densities from file
c if this is the xTB run
      if(iter.eq.0) then
        write(*,*) 'reading <aospin>'
        open(unit=111,file='tmpxtb')
        read(111,*)ndimv
        allocate(vector(ndimv))
        do i=1,ndimv
         read(111,*)vector(i)
        enddo
        close(111,status='delete')
        j=0
        do i=1,ndim
         if(valao2(i).eq.1)then
            j=j+1
            aospin(i)=vector(j)
         endif
        enddo
        deallocate(vector)
c       write(*,*) 'aospin on file'
        call preig3(6,aospin,ndim)
      endif

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

4242  write(*,'(''================================================='')')
      write(*,'(20x,i2,''   SCC step'')')jter
      write(*,'(''================================================='')')
      do ii=1,ndim
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
c           Rybergs feel less atomic charge was tested: if(valao2(i).eq.0.or.valao2(j).eq.0) scal=gscal*1.00
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
       if(iter.eq.1.and.jter.le.1) then ! only in first two iterations (UHF case), use CM5 charges from EN charges
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
c     stop

      call cpu_time(t0)
      write(*,*)'solving ...'

      Ha=H
      Hb=H

c     spin-dependent part of Hamiltonian (wll is in eV)
      do i=1,ndim
         ii=aoat2(i)
         ishell=mmm(lao2(i))
         f=1.0d0+q(ii)*0.400
         do j=1,i
            dum=S(j,i)*0.5
            jj=aoat2(j)
            if(abs(dum).gt.1.d-6.and.ii.eq.jj)then
            jshell=mmm(lao2(j))
            eh1=0.0d0
            do k=fila2(1,ii),fila2(2,ii)
               kshell=mmm(lao2(k))
               eh1=eh1+aospin(k)*f*
     .         (wll(at(ii),lin(kshell,ishell))+
     .          wll(at(ii),lin(kshell,jshell)))
            enddo
            Ha(j,i)=Ha(j,i)+dum*eh1
            Hb(j,i)=Hb(j,i)-dum*eh1
            endif
         enddo
      enddo

c     call prmat(6,ha,ndim,ndim,'Ha')
c     call prmat(6,hb,ndim,ndim,'Ha')
c     alpha
      X=S
      call DSYGVD(1,'V','U',ndim,Ha,ndim,X,ndim,emoa,aux,
     .            LWORK,IWORK,LIWORK,INFO)
      if(info.ne.0) stop 'diag error'
c     beta
      X=S
      call DSYGVD(1,'V','U',ndim,Hb,ndim,X,ndim,emob,aux,
     .            LWORK,IWORK,LIWORK,INFO)
      if(info.ne.0) stop 'diag error'
      call cpu_time(t1)
      write(*,*) 'cpu time for diag ',t1-t0
      write(*,*)

c     shift orbital energies in xTB runs (these are used in Fermi smearing)
      if(iter.eq.0)
     .  call shifteps(.true.,ndim,focca,foccb,emoa,emob,
     .                lshift,lshifta,lshiftb)

      if(et.gt.0.1)then
      call FERMISMEAR(.true.,ndim,ihomoa,et,emoa,focca,nfoda,efa,ga)
      call FERMISMEAR(.true.,ndim,ihomob,et,emob,foccb,nfodb,efb,gb)
      endif

      eel=0
      do i=1,ndim
         emoa(i)=emoa(i)/R27
         emob(i)=emob(i)/R27
         eel=eel+emoa(i)*focca(i)+emob(i)*foccb(i)
      enddo
      if(ihomoa+1.le.ndim)then
      egapa=(emoa(ihomoa+1)-emoa(ihomoa))*R27
      egapb=(emob(ihomob+1)-emob(ihomob))*R27
      write(*,'('' a/b gap (eV)           : '',2f7.3)')
     .egapa,egapb
      else
      egapa=0
      egapb=0
      endif
      write(*,'('' Koopmans IP (eV)       : '',F7.3)')
     .                            -max(emoa(ihomoa),emob(ihomob))*R27
      write(*,'('' Eel (Eh)               : '',F12.6)')eel
      if(et.gt.0.1)then
      write(*,'('' Nfod                   :'',f8.4)')nfoda+nfodb
      write(*,'('' Nfod/Nel               :'',f8.4)')
     .                            (nfoda+nfodb)/float(nel)
      endif

      if(egapa.lt.0.2.or.egapb.lt.0.2)then
         write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write(*,*)'!!!!!!!! WARNING: small HL-gap detected !!!!!!!!!'
         write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      endif

c     make density matrix
      if(pop.or.dipol.or.wbo)then
      call dmat(ndim,focca,Ha,Pa)
      call dmat(ndim,foccb,Hb,Pb)
      ndimv=0
      do i=1,ndim
         if(valao2(i).eq.1)ndimv=ndimv+1
         aospin(i)=Pa(i,i)-Pb(i,i)
      enddo

      X = S
      P = Pa + Pb

      endif

      if(iter.eq.1.and.jter.eq.maxiter)then
        write(*,*) 'writing <aospin>'
        open(unit=111,file='tmpxtb')
        write(111,*)ndimv
        do i=1,ndim
         if(valao2(i).eq.1)write(111,*)aospin(i)
        enddo
        close(111)
      endif

      call preig(6,focca,R27,emoa,max(ihomoa-12,1),min(ihomoa+11,ndim))
      call preig(6,foccb,R27,emob,max(ihomob-12,1),min(ihomob+11,ndim))

      write(*,*) '                ',jter+1,' SCC done.'
c     next SCC
      if(pop.and.iter.eq.1.and.jter.lt.maxiter)then
         call mpop0(n,ndim,aoat2,S,P,qmull)
         qmull = z - qmull
c        update q for the first time in the second SCC i.e. do
c        in the first SCC just the aospin and start with this
c        and the EN the second SCC and one iteration i.e. three in total
         if(jter.gt.0) then
            write(*,*) 'updating charges'
            call docm5(n,at,.false.,xyz,qmull,q)
         endif
         jter = jter + 1
         goto 4242
      endif

      deallocate(T)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC END OF SCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c S^2
      s2=0
      do j=1,ihomob
         do ia=1,ndim
            dum=0.0d0
            do ib=1,ndim
               dum=dum+S(ib,ia)*Hb(ib,j)
            enddo
            X(ia,j)=dum
         enddo
      enddo
      do i=1,ihomoa
         do j=1,ihomob
         dum=0.0d0
            do ia=1,ndim
               dum=dum+Ha(ia,i)*X(ia,j)
            enddo
         s2=s2+dum**2
         enddo
      enddo

      deallocate(X)

      dum=ihomoa-ihomob
      write(*,'('' <S^2>           :'',f8.5)')
     .0.50d0*dum*(0.50d0*dum+1.0d0)+ihomob-s2

c     call edens(n,ncao,ndim,xyz,Pa-Pb)

cccccccccccccccccccccccccccccccccccccccccccccc
c         localization
cccccccccccccccccccccccccccccccccccccccccccccc

      if(lmo) then
         call local(n,at,ncao,ndim,ihomoa,xyz,z,focca,S,Ha,emoa)
         call local(n,at,ncao,ndim,ihomob,xyz,z,foccb,S,Hb,emob)
      endif

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
 120     close(22)
         if(j.ne.n) then
            ex=.false.
            goto 130
         endif
         call docm5(n,at,.false.,xyz,chir,cm5_true)
         call execute_command_line('pwd > tmpxtb2')
         open(unit=43,file='tmpxtb2')
         read(43,'(a)') ftmp
         close(43,status='delete')
        endif

 130    continue
        call mpop(n,ndim,aoat2,lao2,S,P,qmull,qlmom)
        q = z - qmull
        call docm5(n,at,.false.,xyz,q,cm5_appr)

        write(*,*)
        write(*,*)
        write(*,*) 'Mulliken/CM5 charges      n(s)   n(p)   n(d)   '
        do i=1,n
           write(*,'(i5,a3,2f8.4,1x,3f7.3)')
     .     i,asym(at(i)),q(i),cm5_appr(i),qlmom(1:3,i)
        enddo
        qmull=q      ! save

c       fitting stuff for charge TB
        if(ex)then
c        fitting stuff
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

        call mpop(n,ndim,aoat2,lao2,S,Pa-Pb,qmull,qlmom)

        write(*,*)
        write(*,*)'Pa-Pb analysis  n(s)   n(p)   n(d)       HFC (Mhz)'
        do i=1,n
           hfc(i) = qlmom(1,i) * mc(at(i)) * 1000. / totspin
           write(*,'(i5,a3,f7.3,3f7.3,2x,f8.2)')
     .     i,asym(at(i)),qmull(i),qlmom(1:3,i),hfc(i)
        enddo

c       Loewdin pop
        call umakel(ndim, S, Ha, Hb, Pa, Pb)
        q=0
        qlmom=0
        call lpop(n,ndim,aoat2,lao2,focca,Pa,1.0d0,q,qlmom)
        call lpop(n,ndim,aoat2,lao2,foccb,Pb,1.0d0,q,qlmom)
        q = z - q
        write(*,*)
        write(*,*) 'Loewdin charges    n(s)   n(p)   n(d)  '
        do i=1,n
           write(*,'(i5,a3,f8.4,1x,3f7.3)')
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

        q=0
        qlmom=0
        call lpop(n,ndim,aoat2,lao2,focca,Pa, 1.0d0,q,qlmom)
        call lpop(n,ndim,aoat2,lao2,foccb,Pb,-1.0d0,q,qlmom)
        write(*,*)
        write(*,*)'Pa-Pb analysis  n(s)   n(p)   n(d)       HFC (Mhz)'
        do i=1,n
           hfc(i) = qlmom(1,i) * mc(at(i)) * 1000. / totspin
           write(*,'(i5,a3,f7.3,3f7.3,2x,f8.2)')
     .     i,asym(at(i)),q(i),qlmom(1:3,i),
     .     hfc(i)
        enddo
        if(et.gt.1000)then
           allocate(fduma(ndim),fdumb(ndim))
           fduma= focca
           fdumb= foccb
           q=0
           qlmom=0
           call fodenmak(.true.,ndim,emoa,fduma,efa)
           call fodenmak(.true.,ndim,emob,fdumb,efb)
           call lpop(n,ndim,aoat2,lao2,fduma,Pa,1.0d0,q,qlmom)
           call lpop(n,ndim,aoat2,lao2,fdumb,Pb,1.0d0,q,qlmom)
           write(*,*)
           write(*,*) 'Loewdin FODpop     n(s)   n(p)   n(d)   '
           open(unit=43,file='fod')
           do i=1,n
              write(*,'(i5,a3,f8.4,1x,4f7.3)')
     .        i,asym(at(i)),q(i),qlmom(1:3,i)
              write(43,'(F14.8)') q(i)
           enddo
           deallocate(fduma,fdumb)
           write(*,*) 'FODpop written to file <fod>'
           close(43)
        endif
        deallocate(Pa,Pb)

        inquire(file='spin',exist=ex)
        if(ex)then
         open(unit=22,file='spin')
         do i=1,n
           read(22,*) chir(i)
         enddo
         close(22)
         open(unit=43,file='spin3')
         do i=1,n
           write(43,'(2F12.6,2x,a)')
     .     chir(i),q(i),trim(ftmp)
         enddo
         close(43)
        endif

c OUTPUT!
        q = cm5_appr

        inquire(file='orca.out',exist=ex)
        if(ex)then
         wbr=0
         chir=0
         open(unit=43,file='orca.out')
  10     read(43,'(a)',end=100) atmp
c        WBO
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
c        HFC
         if(index(atmp,'Nucleus ').ne.0) then
            call readl(atmp,xx,nn)
            i=idint(xx(1))+1
            if(i.gt.n) stop 'orca.out read error'
            do j=1,14
               read(43,'(a)') atmp
            enddo
            call readl(atmp,xx,nn)
            chir(i)=xx(nn)
         endif
c        next
         goto 10
 100     close(43)
         if(sum(abs(chir(1:n))).gt.1.d-6)then
           open(unit=22,file='hfc3')
           do i=1,n
           if(abs(chir(i)).lt.200) then
              write(22,'(2F14.6,2x,a)')
     .        chir(i), hfc(i), trim(ftmp)
           else
              write(22,'(2F14.6,2x,a)')
     .        chir(i)/5, hfc(i)/5, trim(ftmp)
           endif
           enddo
           close(22)
         endif
c       end orca.out read
        endif

c end pop
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
         allocate(Xcaoa(ncao,ndim),Xcaob(ncao,ndim))
         call sao2cao(ndim,Ha,ncao,Xcaoa)
         call sao2cao(ndim,Hb,ncao,Xcaob)
         if(molden)then
         if(et.gt.1000)call fodenmak(.true.,ndim,emoa,focca,efa)
         if(et.gt.1000)call fodenmak(.true.,ndim,emob,foccb,efb)
         call printumold(n,ndim,ncao,xyz,at,Xcaoa,Xcaob,emoa,emob,
     .                  focca,foccb,
     .                  mothr)
         else
         call printumos(n,ndim,ncao,xyz,at,Xcaoa,Xcaob,emoa,emob,
     .                  focca,foccb,
     .                  mothr)
         endif
      else
         if(molden)then
         if(et.gt.1000)call fodenmak(.true.,ndim,emoa,focca,efa)
         if(et.gt.1000)call fodenmak(.true.,ndim,emob,foccb,efb)
         call printumold(n,ndim,ndim,xyz,at,Ha,Hb,emoa,emob,
     .                  focca,foccb,
     .                  mothr)
         else
         call printumos(n,ndim,ndim,xyz,at,Ha,Hb,emoa,emob,
     .                  focca,foccb,
     .                  mothr)
         endif
      endif
      endif

      if(lgbsa) deallocate(cm5,fgb,fhb)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
      implicit none
      integer nel,nopen,ndim,ihomoa,ihomob
      real*8 focca(ndim), foccb(ndim)
      integer focc(ndim)
      integer i,na,nb,ihomo

      focc=0
      focca=0
      foccb=0
c even nel
      if(mod(nel,2).eq.0)then
      ihomo=nel/2
      do i=1,ihomo
         focc(i)=2
      enddo
      if(2*ihomo.ne.nel) then
         ihomo=ihomo+1
         focc(ihomo)=1
         if(nopen.eq.0)nopen=1
      endif
      if(nopen.gt.1)then
         do i=1,nopen/2
            focc(ihomo-i+1)=focc(ihomo-i+1)-1
            focc(ihomo+i)=focc(ihomo+i)+1
         enddo
      endif
c odd nel
      else
      na=nel/2+(nopen-1)/2+1
      nb=nel/2-(nopen-1)/2
      do i=1,na
         focc(i)=focc(i)+1
      enddo
      do i=1,nb
         focc(i)=focc(i)+1
      enddo
      endif

      do i=1,ndim
         if(focc(i).eq.2)then
            focca(i)=1.0d0
            foccb(i)=1.0d0
         endif
         if(focc(i).eq.1)focca(i)=1.0d0
      enddo

      ihomoa=0
      ihomob=0
      do i=1,ndim
         if(focca(i).gt.0.99) ihomoa=i
         if(foccb(i).gt.0.99) ihomob=i
      enddo

      if(ihomoa.lt.1) stop 'internal error in occu'
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
