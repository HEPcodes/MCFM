      double precision function Bigb(i1,i2,i3,i4,i5)
      implicit none
C     Expressions taken from Berends,Kleiss,DeC,Gastmans,Wu
C     Physics Letters 103B 124, (1981)
C     q_i(p1)+q_i(p2)->q_i(p3)+q_i(p4)+g(p5)
      include 'constants.f'
      include 'sprods_com.f'
      include 'qcdcouple.f'
      integer i1,i2,i3,i4,i5
      double precision ss,sp,tt,tp,uu,up
      double precision e12,e34,e13,e14,e23,e24
      e12=2d0*s(i1,i2)/(s(i1,i5)*s(i2,i5))
      e34=2d0*s(i3,i4)/(s(i3,i5)*s(i4,i5))
      e14=2d0*s(i1,i4)/(s(i1,i5)*s(i4,i5))
      e24=2d0*s(i2,i4)/(s(i2,i5)*s(i4,i5))
      e13=2d0*s(i1,i3)/(s(i1,i5)*s(i3,i5))
      e23=2d0*s(i2,i3)/(s(i2,i5)*s(i3,i5))
      ss=s(i1,i2)
      sp=s(i3,i4)
      tt=s(i1,i3)
      tp=s(i2,i4)
      uu=s(i1,i4)
      up=s(i2,i3)
      Bigb=4d0*gsq**3*xn*CF
     . *(+(ss**2+sp**2+uu**2+up**2)/(2d0*tt*tp)
     . *(2d0*CF*(e14+e23)+(2d0*(e12+e34)-e13-e14-e23-e24)/xn)
     .   +(ss**2+sp**2+tt**2+tp**2)/(2d0*uu*up)
     . *(2d0*CF*(e13+e24)+(2d0*(e12+e34)-e13-e14-e23-e24)/xn)
     .   -2d0/xn*(ss**2+sp**2)*(ss*sp-tt*tp-uu*up)/(4d0*tt*tp*uu*up)
     . *(2d0*CF*(e12+e34)+(2d0*(e12+e34)-e13-e14-e23-e24)/xn))
      return
      end

