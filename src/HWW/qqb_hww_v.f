      subroutine qqb_hww_v(p,msqv)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     . p(mxpart,4)
      integer j,k
      
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0d0
      enddo
      enddo

      call qqb_hww(p,msq)

c -sum of virtual diagram and result of integrating subtraction term
      msqv(0,0)=ason2pi*xn/3d0*(11d0+2d0*pi**2)*msq(0,0)

      return
      end
     
