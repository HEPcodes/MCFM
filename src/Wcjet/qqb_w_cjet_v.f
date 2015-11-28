      subroutine qqb_w_cjet_v(p,msq)
      implicit none
************************************************************************
*     Author: J.M. Campbell                                            *
*     February, 2004.                                                  *
************************************************************************
c---- One-loop matrix element for W+c production, including c mass
c---- averaged over initial colours and spins
c for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4)) + cbar(p5)
c For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ c(p5) 
c----
      include 'constants.f'
      include 'ckm.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scheme.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      double precision virtsqwcg,virtqg,virtgq,virtqbg,virtgqb
      
      scheme='dred'

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      virtqg=virtsqwcg(1,2,3,4,5,p)
      virtgq=virtsqwcg(2,1,3,4,5,p)
      virtgqb=virtsqwcg(2,1,4,3,5,p)
      virtqbg=virtsqwcg(1,2,4,3,5,p)
      fac=ason2pi*gsq*gwsq**2*V*aveqg

c--- fill up lowest order matrix elements
c      call qqb_w_cjet(p,msq0)

c      taugs=-2d0*dot(p,1,2)
c      taucs=+2d0*dot(p,1,5)
c      taucg=+2d0*dot(p,2,5)
      
c      mQsq=mass2**2   
c--- overall factor in the virtual terms is (fourpi*mQsq)^(epsilon)      
c--- to match with our usual definition, we should thus multiply by
c--- (musq/mQsq)^(epsilon)
c      epin=epinv+lnrat(musq,mQsq)
c      epin2=epinv**2+epinv*lnrat(musq,mQsq)+half*lnrat(musq,mQsq)**2

c--- poles for (q,g)
c      poles_qg=-(
c     . +xn*(1.5d0*epin2-epin*lnrat(-taucg,mQsq)
c     .     -epin*lnrat(-taugs,mQsq)+0.5d0*epin
c     .     +lnrat(-taucg,mQsq)**2+0.5d0*lnrat(-taugs,mQsq)**2+1d0)
c     . -0.5d0/xn*(epin2-2d0*epin*lnrat(-taucs,mQsq)
c     .           +epin+2d0+2d0*lnrat(-taucs,mQsq)**2)
c     . -(2d0*tr/3d0*dfloat(nf)-11d0/6d0*xn-1.5d0*cf)*epin
c     . +(11d0/6d0*xn-2d0/3d0*tr*dfloat(nf)+1.5d0*cf)*lnrat(musq,mQsq))

c--- poles for (g,q)
c      poles_gq=-(
c     . +xn*(1.5d0*epin2-epin*lnrat(-taucs,mQsq)
c     .     -epin*lnrat(-taugs,mQsq)+0.5d0*epin
c     .     +lnrat(-taucs,mQsq)**2+0.5d0*lnrat(-taugs,mQsq)**2+1d0)
c     . -0.5d0/xn*(epin2-2d0*epin*lnrat(-taucg,mQsq)
c     .           +epin+2d0+2d0*lnrat(-taucg,mQsq)**2)
c     . -(2d0*tr/3d0*dfloat(nf)-11d0/6d0*xn-1.5d0*cf)*epin
c     . +(11d0/6d0*xn-2d0/3d0*tr*dfloat(nf)+1.5d0*cf)*lnrat(musq,mQsq))
 
      do j=-(nf-2),(nf-2)
      do k=-(nf-2),(nf-2)
      if ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=fac*Vsq(j,-4)*virtqg
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=fac*Vsq(j,+4)*virtqbg
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=fac*Vsq(-4,k)*virtgq
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=fac*Vsq(+4,k)*virtgqb
      endif

      enddo
      enddo
      
      return
      end
