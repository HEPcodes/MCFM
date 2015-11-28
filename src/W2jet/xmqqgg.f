      subroutine xmqqgg(i,za,zb,xmatpc,xmatpv)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'nwz.f'
      double precision xmatpc,xmatpv,pp1,pp2,pm1,pm2,mp1,mp2,mm1,mm2
      double complex xla(-1:1,-1:1),yla(-1:1,-1:1),qla(-1:1,-1:1)
      double complex xal(-1:1,-1:1),yal(-1:1,-1:1),qal(-1:1,-1:1)
      integer j,k,i(6)
      logical first
      data first/.true./
      data pp1,pp2,pm1,pm2,mp1,mp2,mm1,mm2/8*0d0/

C----xla and yla are the two color-ordered amplitudes
      if (first) then
      first=.false.
         do j=-1,1
           do k=-1,1
           xla(j,k)=(0d0,0d0)
           yla(j,k)=(0d0,0d0)
           qla(j,k)=(0d0,0d0)
           xal(j,k)=(0d0,0d0)
           yal(j,k)=(0d0,0d0)
           qal(j,k)=(0d0,0d0)
           enddo
         enddo
      endif
      
      call subqcd(i(1),i(2),i(3),i(4),i(5),i(6),za,zb,xla)
      
C do the second ordering exchanging gluon momenta (5--6)

      call subqcd(i(1),i(2),i(3),i(4),i(6),i(5),za,zb,yla)

C---add up to give the symmetric combination of color amplitudes
C---ie QED-like piece
      qla(+1,+1)=xla(+1,+1)+yla(+1,+1)
      qla(+1,-1)=xla(+1,-1)+yla(-1,+1)
      qla(-1,+1)=xla(-1,+1)+yla(+1,-1)
      qla(-1,-1)=xla(-1,-1)+yla(-1,-1)

C---all the left-handed lepton couplings amplitudes
      pp1=abs(xla(+1,+1))**2+abs(yla(+1,+1))**2-abs(qla(+1,+1))**2*ninth
      mp1=abs(xla(+1,-1))**2+abs(yla(-1,+1))**2-abs(qla(+1,-1))**2*ninth
      pm1=abs(xla(-1,+1))**2+abs(yla(+1,-1))**2-abs(qla(-1,+1))**2*ninth
      mm1=abs(xla(-1,-1))**2+abs(yla(-1,-1))**2-abs(qla(-1,-1))**2*ninth

c      pp1=-abs(qla(+1,+1))**2*ninth
c      mp1=-abs(qla(+1,-1))**2*ninth
c      pm1=-abs(qla(-1,+1))**2*ninth
c      mm1=-abs(qla(-1,-1))**2*ninth

c      qla(+1,-1)=xla(+1,-1)+yla(+1,-1)
c      qla(-1,+1)=xla(-1,+1)+yla(-1,+1)
c      if ((dreal(za(5,6)*zb(6,5)) .lt. 30d0) .and.
c     . (dreal(za(1,5)*zb(5,1)) .gt. -30d0)) then
c      write(*,*) 'Perm ',i
c      write(*,*) 'x',abs(xla(+1,+1))**2+abs(xla(+1,-1))**2
c     .              +abs(xla(-1,+1))**2+abs(xla(-1,-1))**2
c      write(*,*) 'y',abs(yla(+1,+1))**2+abs(yla(+1,-1))**2
c     .              +abs(yla(-1,+1))**2+abs(yla(-1,-1))**2
c      write(*,*) 'q',abs(qla(+1,+1))**2+abs(qla(+1,-1))**2
c     .              +abs(qla(-1,+1))**2+abs(qla(-1,-1))**2
c      endif
  
C now exchange lepton momenta
      if (nwz .eq. 0) then
      call subqcd(i(1),i(2),i(4),i(3),i(6),i(5),za,zb,yal)
C do the second ordering exchanging gluon momenta (5-6)
      call subqcd(i(1),i(2),i(4),i(3),i(5),i(6),za,zb,xal)
      qal(+1,+1)=xal(+1,+1)+yal(+1,+1)
      qal(+1,-1)=xal(+1,-1)+yal(-1,+1)
      qal(-1,+1)=xal(-1,+1)+yal(+1,-1)
      qal(-1,-1)=xal(-1,-1)+yal(-1,-1)

C---all the right-handed lepton couplings amplitudes
      pp2=abs(xal(+1,+1))**2+abs(yal(+1,+1))**2-abs(qal(+1,+1))**2*ninth
      mp2=abs(xal(+1,-1))**2+abs(yal(-1,+1))**2-abs(qal(+1,-1))**2*ninth
      pm2=abs(xal(-1,+1))**2+abs(yal(+1,-1))**2-abs(qal(-1,+1))**2*ninth
      mm2=abs(xal(-1,-1))**2+abs(yal(-1,-1))**2-abs(qal(-1,-1))**2*ninth

      endif

C---xmatc left handed quark line, left-handed lepton
      xmatpc=V*XN/four*((pp1+pm1+mp1+mm1)+(pp2+pm2+mp2+mm2))/2d0
C---xmatv left handed quark line, right-handed lepton
      xmatpv=V*XN/four*((pp1+pm1+mp1+mm1)-(pp2+pm2+mp2+mm2))/2d0

      return
      end
