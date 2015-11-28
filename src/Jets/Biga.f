      double precision function Biga(i1,i2,i3,i4,i5)
      implicit none
C     Expressions taken from Ellis-Stirling-Webber Page 257
C     q_i(p1)+q_j(p2)->q_i(p3)+q_j(p4)+g(p5)
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
      double precision ss,sp,tt,tp,uu,up,
     .  eik12,eik34,eik13,eik14,eik23,eik24
c      double precision c1,c2
      eik12=s(i1,i2)/(s(i1,i5)*s(i2,i5))
      eik34=s(i3,i4)/(s(i3,i5)*s(i4,i5))
      eik14=s(i1,i4)/(s(i1,i5)*s(i4,i5))
      eik24=s(i2,i4)/(s(i2,i5)*s(i4,i5))
      eik13=s(i1,i3)/(s(i1,i5)*s(i3,i5))
      eik23=s(i2,i3)/(s(i2,i5)*s(i3,i5))
      ss=s(i1,i2)
      sp=s(i3,i4)
      tt=s(i1,i3)
      tp=s(i2,i4)
      uu=s(i1,i4)
      up=s(i2,i3)
c      c1=V**2/xn
c      c2=V/xn
c      Biga=c1*((uu+up)*(ss*sp+tt*tp-uu*up)
c     . +uu*(ss*tt+sp*tp)+up*(ss*tp+sp*tt))
c      Biga=Biga-c2*((ss+sp)*(ss*sp-tt*tp-uu*up)
c     . +2d0*tt*tp*(uu+up)+2d0*uu*up*(tt+tp))

c      Biga=2d0*Biga*(ss**2+sp**2+uu**2+up**2)
c     . /(tt*tp*s(i1,i5)*s(i2,i5)*s(i3,i5)*s(i4,i5))


      Biga=8d0*gsq**3*xn*CF*(ss**2+sp**2+uu**2+up**2)/(2*tt*tp)
     . *(2d0*CF*(eik14+eik23)
     . +(2d0*(eik12+eik34)-eik13-eik14-eik23-eik24)/xn)

      return
      end
