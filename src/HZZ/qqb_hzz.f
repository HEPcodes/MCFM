      subroutine qqb_hzz(p,msq)
      implicit none
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  W^- (e^-(p3)+nubar(p4)) + W^+ (nu(p5)+e^+(p6))
c---
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),s,s12
      double precision decay,gg,Asq
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

c---set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

c--- Old code
c      decay=4d0*gwsq*esq**2*(wmass/(one-xw))**2
c     . *(((l1*l2)**2+(r1*r2)**2)*s(3,5)*s(4,6)
c     .  +((r1*l2)**2+(r2*l1)**2)*s(3,6)*s(4,5))

c      aw=gwsq/(4d0*pi)

c      s12=s(1,2)
c      gg=aw*as**2*V/(18d0*pi)*(s12/wmass)**2
c---calculate propagators
c      fac=one/((s(5,6)-zmass**2)**2+(zmass*zwidth)**2)
c      fac=fac/((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
c      fac=fac/((s12-hmass**2)**2+(hmass*hwidth)**2)
c      msq(0,0)=avegg*fac*gg*decay

c--- New code (consistent with gg_h.f etc.)
      s12=s(1,2)

      decay=gwsq**3*zmass**2*4d0*xw**2/(one-xw)*
     . ( ((l1*l2)**2+(r1*r2)**2)*s(3,5)*s(4,6)
     .  +((r1*l2)**2+(r2*l1)**2)*s(3,6)*s(4,5))
      decay=decay/((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      decay=decay/((s(5,6)-zmass**2)**2+(zmass*zwidth)**2)
      decay=decay/((s12-hmass**2)**2+(hmass*hwidth)**2)

      Asq=(as/(3d0*pi))**2/vevsq
      gg=0.5d0*V*Asq*s12**2

c---calculate propagators
      msq(0,0)=avegg*gg*decay

      return
      end
