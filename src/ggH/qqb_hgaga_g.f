      subroutine qqb_hgaga_g(p,msq)
      implicit none
c----NLO matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H(gamma(p3)+gamma(p4))+g(p5)
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j,k,iglue
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision sh,ss,tt,uu,s(mxpart,mxpart),decay
      double precision aw,qqb,qg,gq,gg
      parameter(iglue=5)
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      aw=gwsq/(4d0*pi)
      call dotem(iglue,p,s)

      decay=1d0
c--   calculate propagators
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)
      sh=s(1,2)+s(1,iglue)+s(2,iglue)

      gg=aw*as**3*4d0*V/9d0*xn*(sh**4+ss**4+tt**4+uu**4)
     . /(ss*tt*uu*wmass**2)
      qqb=aw*as**3*2d0*V/9d0*(tt**2+uu**2)/(ss*wmass**2)
      gq=-aw*as**3*2d0*V/9d0*(ss**2+tt**2)/(uu*wmass**2)
      qg=-aw*as**3*2d0*V/9d0*(ss**2+uu**2)/(tt*wmass**2)


      gg=avegg*decay*gg
      gq=aveqg*decay*gq
      qg=aveqg*decay*qg
      qqb=aveqq*decay*qqb

c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      if ((k .eq. -j) .and. (j .ne. 0)) then
      msq(j,k)=qqb
      elseif ((j .eq. 0) .and. (k .ne. 0)) then
      msq(j,k)=gq
      elseif ((j .ne. 0) .and. (k .eq. 0)) then
      msq(j,k)=qg
      elseif ((k .eq. 0) .and. (j .eq. 0)) then
      msq(j,k)=gg
      endif
      enddo
      enddo
      return
      end
