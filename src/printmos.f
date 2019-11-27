ccccccccccccccccccccccccccccccccccccccccccc
!    write out unformatted sTDA input     c
ccccccccccccccccccccccccccccccccccccccccccc
! ncent  : # atoms
! nmo    : # MOs
! nbf    : # AOs
! nprims : # primitives (in total)
! xyz(4,ncent) : Cartesian coordinates & nuclear charge
! cont(nprims) : contraction coefficients of primitives
! alp(nprims) : exponents of primitives
! cmo(nbf,nmo) : LCAO-MO coefficients
! eval(nmo)    : orbital eigenvalues
! occ(nmo)     : occupation # of MO
! ipty(nprims) : angular momentum of primitive function
! ipao(nbf)    : # primitives in contracted AO
! ibf(ncent)   : # of contracted AOs on atom

      subroutine printmos(ncent,nmo,nbf,xyz,at,cmo,eval,
     .                    occ,mowrcut)

      include 'ehtcommon.fh'

      real*8, intent ( in ) :: xyz(3,ncent)
      real*8, intent ( in ) :: eval(nmo)
      real*8, intent ( in ) :: occ (nmo)
      real*8, intent ( in ) :: cmo(nbf,nmo)
      real*8, intent ( in ) :: mowrcut
      integer, intent( in ) :: at(ncent)
      integer, intent( in ) :: ncent,nmo,nbf
      ! temporary variables
      integer i,j,k,nprims,nmomax
      real*8 dum
      character*2 atyp

      nprims=0
      do i=1,nbf
        nprims=nprims+nprim(i)
      enddo

      iwfn=29
      open(unit=iwfn,file='wfn.xtb',form='unformatted',
     .     status='replace')


! only print out virtuals below cutoff
      nmomax=0
      do i=1,nmo
         if(eval(i).gt.mowrcut.and.nmomax.eq.0)nmomax=i-1
      enddo
      if(nmomax.eq.0) nmomax=nmo

                    !***********
                    ! RHF case *
                    !***********
! write dimensions
       write(iwfn)1
       write(iwfn)ncent,nbf,nmomax,nprims

! now write coordinates & atom symbol
       do i = 1,ncent
         call aasym(at(i),atyp)
         write(iwfn) atyp
       enddo

       do i = 1,ncent
         do j=1,3
            dum=xyz(j,i)
            write(iwfn) dum
         enddo
         write(iwfn) at(i)
       enddo

! Now print basis set data

! print ipty
       do i=1,nbf
          k = lao(i)
          do j=1,nprim(i)
             write(iwfn) k
          enddo
       enddo
! iaoat
       do i=1,nbf
          k=aoat(i)
          do j=1,nprim(i)
             write(iwfn) k
          enddo
       enddo
! ipao
       do i=1,nbf
          k=i
          do j=1,nprim(i)
             write(iwfn) k
          enddo
       enddo

! exponents and coefficients
         write(iwfn) alp(1:nprims)
         write(iwfn) cont(1:nprims)

! now the mo data
         write(iwfn) occ(1:nmomax)
         write(iwfn) eval(1:nmomax)

         write(iwfn) cmo(1:nbf,1:nmomax)
      close(iwfn)

      return
      end


C 'true' U version
      subroutine printumos(ncent,nmo,nbf,xyz,at,cmoa,cmob,evala,evalb,
     .                     occa,occb,mowrcut)

      include 'ehtcommon.fh'

      real*8, intent ( in ) :: xyz(3,ncent)
      real*8, intent ( in ) :: evala(nmo),evalb(nmo)
      real*8, intent ( in ) :: occa(nmo),occb(nmo)
      real*8, intent ( in ) :: cmoa(nbf,nmo),cmob(nbf,nmo)
      real*8, intent ( in ) :: mowrcut
      integer, intent( in ) :: at(ncent)
      integer, intent( in ) :: ncent,nmo,nbf
      ! temporary variables
      integer i,j,k,nprims,nmomax
      real*8 dum
      character*2 atyp

      nprims=0
      do i=1,nbf
        nprims=nprims+nprim(i)
      enddo

      iwfn=29
      open(unit=iwfn,file='wfn.xtb',form='unformatted',
     .      status='replace')

! only print out virtuals below cutoff
      nmomax=0
      do i=1,nmo
         if(evala(i).gt.mowrcut.and.nmomax.eq.0)nmomax=i-1
      enddo
      if(nmomax.eq.0) nmomax=nmo

! write dimensions
       write(iwfn)2
       write(iwfn)ncent,nbf,nmomax,nprims
! now write coordinates
       do i = 1,ncent
         call aasym(at(i),atyp)
         write(iwfn) atyp
       enddo
       do i = 1,ncent
         do j=1,3
            dum=xyz(j,i)
            write(iwfn) dum
         enddo
         write(iwfn) at(i)
       enddo

! Now print basis set data
! print ipty
       do i=1,nbf
          k = lao(i)
          do j=1,nprim(i)
             write(iwfn) k
          enddo
       enddo
! iaoat
       do i=1,nbf
          k=aoat(i)
          do j=1,nprim(i)
             write(iwfn) k
          enddo
       enddo
! ipao
       do i=1,nbf
          k=i
          do j=1,nprim(i)
             write(iwfn) k
          enddo
       enddo

! exponents and coefficients
       write(iwfn) alp(1:nprims)
       write(iwfn) cont(1:nprims)

! now the mo data
       write(iwfn) occa(1:nmomax)
       write(iwfn) evala(1:nmomax)
       write(iwfn) occb(1:nmomax)
       write(iwfn) evalb(1:nmomax)

       write(iwfn) cmoa(1:nbf,1:nmomax)
       write(iwfn) cmob(1:nbf,1:nmomax)

      close(iwfn)

      return
      end

