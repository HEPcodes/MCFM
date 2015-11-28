      subroutine qqb_w2jet_g(p,msq)
c----Matrix element squared for the process
c    parton1(-p1)+parton2(-p1)-->V(lepton(p3)+alepton(p4))
c    +parton(p5)+parton(p6)+parton(p7)
c-------------
      implicit none
      include 'constants.f'
      include 'masses.f'
      logical first
      integer n1,n2,nwz
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),gsq,as,
     . ason2pi,dot,s34,propw
      common/nwz/nwz
      common/qcdcouple/gsq,as,ason2pi
      data first/.true./
      

      if (first) then
      first=.false.
      call coupling
      endif
      call itwojet(p,nwz)      

      s34=2d0*dot(p,3,4)
      propw=s34**2/((s34-wmass**2)**2+(wmass*wwidth)**2)
      
      do n1=-nf,nf
      do n2=-nf,nf
      msq(n1,n2)=propw
      enddo
      enddo
      return
      
      end


