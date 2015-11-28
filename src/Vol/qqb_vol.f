      subroutine qqb_vol(P,msq)
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'dprodx.f'
      
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf)
      integer j,k,N
      N=6
      call dotem(N,p,s)
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      msq(2,-1)=1d0
      return
      end

      subroutine qqb_vol_g(P,msq)
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'dprodx.f'
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf)
      integer j,k,N
      N=7
      call dotem(N,p,s)
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      msq(2,-1)=1d0
      return
      end

      subroutine qqb_vol_gs(P,msq)
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'ptilde.f'
      include 'dprodx.f'
      double precision P(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      integer j,k,nd,N
      N=7
      call dotem(N,p,s)
      do nd=1,maxd
      do j=-nf,nf
      do k=-nf,nf
      msq(nd,j,k)=0d0
      enddo
      enddo
      enddo
      ndmax=0
      
      return
      end
