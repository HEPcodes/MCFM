      subroutine gg_h(p,msq)
      implicit none
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  b(p3)+b(p4))
c---
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),s,s12
      double precision hdecay,gg,Asq
      double precision aw
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

c---set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo


C   Deal with Higgs decay to b-bbar
      s12=s(1,2)
      hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s12-4d0*mb**2) 
      hdecay=hdecay/((s12-hmass**2)**2+(hmass*hwidth)**2)

      Asq=(as/(3d0*pi))**2/vevsq
      gg=0.5d0*Asq*V*s12**2
      msq(0,0)=avegg*gg*hdecay
      return
      end
