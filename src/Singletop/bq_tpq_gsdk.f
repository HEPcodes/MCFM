      subroutine bq_tpq_gsdk(p,msqc)
      implicit none

c     Matrix element for t-bbar production with radiation in decay
c      b(-p1)+u(-p2)-->t(n(p3)+e^+(p4)+b(p5)+g(p7))+d(p6)
C     averaged(summed) over initial(final) colours and spins
c--- g(p7) represents a gluon 


      include 'constants.f'
      include 'masses.f'
      include 'ptilde.f'
      include 'qcdcouple.f'
      double precision msq(-nf:nf,-nf:nf),msqc(maxd,-nf:nf,-nf:nf),
     . p(mxpart,4),q(mxpart,4),omz,z,fac,ptDpg,pbDpg,ptDpb
      integer j,k
      logical incldip(0:maxd)
      common/incldip/incldip

      do j=-nf,nf
      do k=-nf,nf
      msqc(1,j,k)=0d0
      enddo
      enddo

      ndmax=1
      incldip(1)=.true.

      call wtransform(p,q,pbDpg,ptDpg,ptDpb)
      omz=ptDpg/(ptDpb+ptDpg-pbDpg)
      z=1d0-omz
      call bq_tpq(q,msq) 
      fac=gsq*cf
     . *(1d0/pbDpg*(2d0/omz-1d0-z)-(mt/ptDpg)**2)


      do j=-nf,nf
      do k=-nf,nf
      msqc(1,j,k)=fac*msq(j,k)
      enddo
      enddo
      return
      end

