      subroutine qqb_dirgam_gvec(p,n,in,msq)
C*********************************************************************** 
c     Author: R.K. Ellis                                               *
c     October, 2002.                                                   *
c     Matrix element for gamma production                              *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p4)              *
c     q(-p1)+qbar(-p2)--> gamma(p3)+ g(p4)                             *
C*********************************************************************** 
      implicit none
      include 'constants.f'
      integer j,k,in
C--in is the label of the parton dotted with n
      double precision msq(-nf:nf,-nf:nf),msqa(-nf:nf,-nf:nf),
     . p(mxpart,4),n(4),nDn,nDp4

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
 
      call qqb_dirgam(p,msqa)
      nDn=n(4)**2-n(1)**2-n(2)**2-n(3)**2
      nDp4=n(4)*p(in,4)-n(3)*p(in,3)-n(2)*p(in,2)-n(1)*p(in,1)

c--- appropriate scale is approx 1d-3*energy(incoming)
c--- so of order(1) for the Tevatron
      if (abs(nDp4).gt.1d-3*abs(p(1,4))) then 
         write(*,*) 'Error for in=',in
         write(*,*) 'cutoff',1d-3*abs(p(in,4))
         write(6,*) 'nDp4',nDp4
         call flush(6)
         stop
      endif


      do j=-nf,nf
      if (in .eq. 1) then
      msq(0,j)=-0.5d0*nDn*msqa(0,j)
      elseif (in .eq. 2) then
      msq(j,0)=-0.5d0*nDn*msqa(j,0)
      elseif (in .eq. 4) then      
      msq(j,-j)=-0.5d0*nDn*msqa(j,-j)
      elseif (in .eq. 5) then      
      msq(j,-j)=-0.5d0*nDn*msqa(j,-j)
      endif
      enddo
      return
      end
