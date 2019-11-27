CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C compute coordination numbers by adding an inverse damping function
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine ncoord2(mzero,natoms,iz,xyz,cn,cn_thr)
      implicit none
      include 'aoelementcommon.fh'
      include 'd3common.fh'
      real*8 k1
      parameter (k1     =16)
      integer iz(*),natoms,i,max_elem
      real*8 xyz(3,*),cn(*),input
      real*8 cn_thr
      logical mzero

      integer iat
      real*8 dx,dy,dz,r,damp,xn,rr,rco,r2

      do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz
            if (r2.gt.cn_thr) cycle
            r=sqrt(r2)
c covalent distance in Bohr
            rco=rcov(iz(i))+rcov(iz(iat))
            rr=rco/r
c counting function exponential has a better long-range behavior than MHGs inverse damping
            damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
            xn=xn+damp
         endif
      enddo
      if(mzero.and.metal(iz(i)).eq.1) xn=0.0d0
      cn(i)=xn
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c coordination number tensor and p xyz exponent 'shifts'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine aniso(n,at,coord,pscale)
      implicit none
      integer at(n),n,i,j,k
      real*8 coord(3,*)
      real*8 t(6),eig(3),u(3,3),atmass,x,y,z,r,rr,dum,av
      real*8 pscale(3,n)
c     include 'aoelementcommon.fh'
      include 'd3common.fh'
      integer ncore

************************************************************************
*    matrix for moments of inertia is of form
*
*           |   y**2+z**2                         |
*           |    -y*x       z**2+x**2             | -i =0
*           |    -z*x        -z*y       x**2+y**2 |
*
************************************************************************
       pscale=0
       do i=1,n
          t=0
          do j=1,n
            if(i.eq.j)cycle
            x=coord(1,i)-coord(1,j)
            y=coord(2,i)-coord(2,j)
            z=coord(3,i)-coord(3,j)
c the contribution of atom j to the CN (='mass')
            r=sqrt(x*x+y*y+z*z)
            rr=(rcov(at(i))+rcov(at(j)))/r
            atmass=1.d0/(1.d0+exp(-12.0*(rr-1.0d0)))
            atmass=atmass*float(at(j)-ncore(at(j)))
c CN done
            t(1)=t(1)+atmass*(y**2+z**2)
            t(2)=t(2)-atmass*x*y
            t(3)=t(3)+atmass*(z**2+x**2)
            t(4)=t(4)-atmass*z*x
            t(5)=t(5)-atmass*y*z
            t(6)=t(6)+atmass*(x**2+y**2)
         enddo
         call rsp(t,3,3,eig,u)
         av=sum(eig)/3.+1.d-12
         eig=(av-eig)/av
c        write(*,'(2i3,3f8.3,2F12.5)')
c    .   i,at(i),eig,av
         do k=1,3
            do j=1,3
            pscale(j,i)=pscale(j,i)+eig(k)*u(j,k)**2
            enddo
         enddo
         write(*,'(5F14.6)') eig(1:3)
         write(*,'(i3,5F14.6)') i,av,pscale(1:3,i)
         call  prmat(6,u,3,3,'U')
      enddo

      end


! this computes a geometrical bond order sum for all atoms
! instead of CN
! use with parameters in file .param_geo_wbo
! radii in /stuff/ common
! read in main as:
!     open(unit=111,file='~/.param')
!     do i=1,86
!        read(111,*) radbo(i,1:3)
!     enddo
!     close(111)

c     subroutine geobo(natoms,iz,xyz,bo,dboij)
c     implicit none
c     integer, intent(in)  :: natoms,iz(natoms)
c     real*8,  intent(in)  :: xyz(3,natoms)
c     real*8,  intent(out) :: bo(natoms),dboij(3,natoms,natoms)
c     real*8 f1,bo_thr
c     parameter (f1     =16.0d0)
c     parameter (bo_thr=1000.d0)
c     integer i,j,k,l,max_elem
c     include 'stuff.fh'
c     real*8 input,expterm,dbo
c     integer iat,jat,nmax
c     real*8 dr(3),r,damp,xn,rr,rco,r2,ftmp

c     nmax=3 ! maximum number of possible bonds

c     bo=0.0d0
c     dboij=0.0d0
c     do iat=1,natoms
c       xn=0.0d0
c       do jat=1,natoms
c        if(jat.eq.iat)cycle
c           dr(1)=xyz(1,jat)-xyz(1,iat)
c           dr(2)=xyz(2,jat)-xyz(2,iat)
c           dr(3)=xyz(3,jat)-xyz(3,iat)
c           r2=sum(dr*dr)
c           if (r2.gt.bo_thr) cycle
c           r=sqrt(r2)
c covalent distance in Bohr
c           damp=0.0d0
c           dbo=0.0d0
c           ftmp=f1
c this loop goes over differently contracted rcovs (to get BO/CN > 1)
c           do l=1,nmax
c              rco=radbo(iz(iat),l)+radbo(iz(jat),l)
c              rr=rco/r
c              expterm=exp(-ftmp*(rr-1.d0))
c              damp=damp+1.d0/(1.d0+expterm)
c              dbo=dbo-ftmp*rco*expterm/
c    .           (r2*(expterm+1.d0)*(expterm+1.d0))
c              ftmp=ftmp*1.50d0 ! make damping more steep,
c           enddo
c           xn=xn+damp
c         ! now dBO/dxyz :  derivative of BO(j) for xyz-displacements of atom i
c           dboij(:,jat,iat)=-dbo*dr(:)/r
c           dboij(:,iat,iat)=dboij(:,iat,iat)-dbo*dr(:)/r
c----------------
c       enddo
c       bo(iat)=xn
c     enddo

c     end subroutine geobo


      subroutine d3modncoord(natoms,iz,xyz,bo,dboij)
      implicit none
      integer, intent(in)  :: natoms,iz(natoms)
      real*8,  intent(in)  :: xyz(3,natoms)
      real*8,  intent(out) :: bo(natoms),dboij(3,natoms,natoms)
c     include 'd3common.fh' ! rcov
c     include 'aoelementcommon.fh' ! EN
      real*8 f1,bo_thr,k4,k5,k6
      parameter (k4=4.10451d0)
      parameter (k5=19.08857d0)
      parameter (k6=2*11.28174d0**2)
      parameter (f1     =16.0d0)
      parameter (bo_thr=1600.d0)
      integer i,j,k,l,max_elem
      integer iat,jat
      real*8 dr(3),r,damp,xn,rr,rco,r2
      real*8 expterm,dbo,rcov(94)
c D3 radii
      data rcov/
     . 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865,
     . 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527,
     . 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820,
     . 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730,
     . 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923,
     . 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188,
     . 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349,
     . 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216,
     . 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717,
     . 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967,
     . 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625,
     . 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657,
     . 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833,
     . 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098,
     . 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878,
     . 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790,
     . 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584,
     . 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289,
     . 3.82984466, 3.85504098, 3.88023730, 3.90543362 /
      real*8 r0ab(94,94)

      call setr0ab(94,0.6d0,r0ab)

      dboij=0.0d0
      do iat=1,natoms
        xn=0.0d0
        do jat=1,natoms
         if(jat.eq.iat)cycle
            dr(1)=xyz(1,jat)-xyz(1,iat)
            dr(2)=xyz(2,jat)-xyz(2,iat)
            dr(3)=xyz(3,jat)-xyz(3,iat)
            r2=sum(dr*dr)
            if (r2.gt.bo_thr) cycle
            r=sqrt(r2)
cd4         den=k4*exp(-(abs((en(iz(iat))-en(iz(jat))))+ k5)**2/k6 )
! this loop goes over differently contracted rcovs (to get BO/CN > 1)
c              rco=rcov(iz(iat))+rcov(iz(jat))
               rco=r0ab(iz(iat),iz(jat))
               rr=rco/r
               expterm=exp(-f1*(rr-1.d0))
               damp=1.d0/(1.d0+expterm)
               dbo=-f1*rco*expterm/
     .           (r2*(expterm+1.d0)*(expterm+1.d0))
            xn=xn+damp
          ! now dBO/dxyz :  derivative of BO(j) for xyz-displacements of atom i
            dboij(:,jat,iat)=-dbo*dr(:)/r
            dboij(:,iat,iat)=dboij(:,iat,iat)-dbo*dr(:)/r
!----------------
        enddo
        bo(iat)=xn
      enddo

      end subroutine d3modncoord
