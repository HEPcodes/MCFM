      subroutine gq_qqqb(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2002.                                                     *
*     calculate the element squared                                    *
*     for the process                                                  *
c----My notation                                                       *
C     g(-p1) +q(-p2)=q(p3)+q(p4)+qb(p5)                                *
C     massless quarks                                                  *
C***********************************************************************
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'prods.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
     . ,wtgq,wtqg,wtga,wtag,qqqqg

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call dotem(5,p,s)

      if (s(1,3)*s(2,3)/s(1,2) .lt. 100d0) return
      if (s(1,4)*s(2,4)/s(1,2) .lt. 100d0) return
      if (s(1,5)*s(2,5)/s(1,2) .lt. 100d0) return
      if (s(3,4).lt. 100d0) return
      if (s(3,5).lt. 100d0) return
      if (s(4,5).lt. 100d0) return

      wtqg=aveqg*qqqqg(1,5,3,4,2)     
      wtgq=aveqg*qqqqg(2,5,3,4,1)     
      wtag=aveqg*qqqqg(3,4,2,5,1)     
      wtga=aveqg*qqqqg(3,4,1,5,2)     
      
      do j=-nf,nf
      do k=-nf,nf
      if ((j.eq.0) .and. (k.eq.+5)) msq(j,k)=wtgq
      if ((j.eq.0) .and. (k.eq.-5)) msq(j,k)=wtga
      if ((j.eq.+5) .and. (k.eq.0)) msq(j,k)=wtqg
      if ((j.eq.-5) .and. (k.eq.0)) msq(j,k)=wtag
      enddo
      enddo
      return
      end

      double precision function qqqqg(i1,i2,i3,i4,i5)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'prods.f'
      integer i1,i2,i3,i4,i5
      double precision ss,tt,uu,sp,tp,up,e1234,e13,e14,e23,e24,e12,e34

      ss=s(i1,i2)
      tt=s(i1,i3)
      uu=s(i1,i4)

      sp=s(i3,i4)
      tp=s(i2,i4)
      up=s(i2,i3)

      e12=2d0*s(i1,i2)/(s(i1,i5)*s(i2,i5))
      e34=2d0*s(i3,i4)/(s(i3,i5)*s(i4,i5))
      e13=2d0*s(i1,i3)/(s(i1,i5)*s(i3,i5))
      e14=2d0*s(i1,i4)/(s(i1,i5)*s(i4,i5))
      e23=2d0*s(i2,i3)/(s(i2,i5)*s(i3,i5))
      e24=2d0*s(i2,i4)/(s(i2,i5)*s(i4,i5))

      e1234=2d0*(e12+e34)-e13-e14-e23-e24
      qqqqg=-2d0*V*gsq**3*(
     .+(ss**2+sp**2+uu**2+up**2)/(2d0*tt*tp)*(2d0*cf*(e14+e23)+e1234/xn)      
     .+(ss**2+sp**2+tt**2+tp**2)/(2d0*uu*up)*(2d0*cf*(e13+e24)+e1234/xn)      
     .-2d0/xn*(ss**2+sp**2)*(ss*sp-tt*tp-uu*up)/(4d0*tt*tp*uu*up)
     .*(2d0*cf*(e12+e34)+e1234/xn))
      return
      end
