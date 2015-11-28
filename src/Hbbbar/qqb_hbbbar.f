      subroutine qqb_hbbbar(p,msq)
      implicit none
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  b(p3)+bar(p4)
c---
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),s12
      double precision hdecay,aw,prop,gg

c---set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      s12=2d0*(p(1,4)*p(2,4)-p(1,1)*p(2,1)-p(1,2)*p(2,2)-p(1,3)*p(2,3))
      hdecay=gwsq*mbsq/(4d0*wmass**2)*2d0*(s12-4d0*mb**2)*xn 
      aw=gwsq/(4d0*pi)
C---gg is the matrix element squared before spin and color averaging
C   for gg-->H production
      gg=aw*as**2*V/(18d0*pi)*(s12/wmass)**2
c---calculate propagator
      prop=one/((s12-hmass**2)**2+(hmass*hwidth)**2)
      msq(0,0)=avegg*prop*gg*hdecay

      return
      end
