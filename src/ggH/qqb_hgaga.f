      subroutine qqb_hgaga(p,msq)
      implicit none
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H --> gamma(p3)+gamma(p4)
c---
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),s,s12
      double precision decay,aw,gg
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

c---set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo



      decay=1d0
      s12=s(1,2)
 
c---calculate propagators and decay
c      decay=gwsq**3*wmass**2*s(3,5)*s(4,6)
c      decay=decay/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
c      decay=decay/((s(5,6)-wmass**2)**2+(wmass*wwidth)**2)
      decay=decay/((s12-hmass**2)**2+(hmass*hwidth)**2)

      aw=gwsq/(4d0*pi)
      gg=aw*as**2*V/(18d0*pi)*(s12/wmass)**2
      msq(0,0)=avegg*gg*decay

      return
      end
