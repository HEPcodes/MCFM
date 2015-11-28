      subroutine qqb_hww(p,msq)
      implicit none
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  W^+ (nu(p3)+e^+(p4))+W^- (e^-(p5)+nubar(p6))
c---
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),s,s12
      double precision decay,aw,fac,gg
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

c---set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo


      decay=gwsq**3*wmass**2*s(3,5)*s(4,6)
      aw=gwsq/(4d0*pi)


      s12=s(1,2)
      gg=aw*as**2*V/(18d0*pi)*(s12/wmass)**2
c---calculate propagators
      fac=one/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
      fac=fac/((s(5,6)-wmass**2)**2+(wmass*wwidth)**2)
      fac=fac/((s12-hmass**2)**2+(hmass*hwidth)**2)
      msq(0,0)=avegg*fac*gg*decay

      return
      end
