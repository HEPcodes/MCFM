      double precision function wbb(i1,i2,i4,i5,i6,i7) 
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'dprodx.f'
      integer i1,i2,i4,i5,i6,i7
      double precision p1Dp2,p1Dp4,p1Dp5,p1Dp6,p1Dp7
      double precision       p2Dp4,p2Dp5,p2Dp6,p2Dp7
      double precision             p4Dp5,p4Dp6,p4Dp7
      double precision                   p5Dp6,p5Dp7
      double precision                         p6Dp7
      double precision s167,s267,s45,xl,xr,xlr,prop
      p1Dp2=0.5d0*s(i1,i2)
      p1Dp4=0.5d0*s(i1,i4)
      p1Dp5=0.5d0*s(i1,i5)
      p1Dp6=0.5d0*s(i1,i6)
      p1Dp7=0.5d0*s(i1,i7)

      p2Dp4=0.5d0*s(i2,i4)
      p2Dp5=0.5d0*s(i2,i5)
      p2Dp6=0.5d0*s(i2,i6)
      p2Dp7=0.5d0*s(i2,i7)
      p4Dp5=0.5d0*s(i4,i5)
      p4Dp6=0.5d0*s(i4,i6)
      p4Dp7=0.5d0*s(i4,i7)
      p5Dp6=0.5d0*s(i5,i6)
      p5Dp7=0.5d0*s(i5,i7)
      p6Dp7=0.5d0*s(i6,i7)
      s167=s(i1,i6)+s(i1,i7)+s(i6,i7)
      s267=s(i2,i6)+s(i2,i7)+s(i6,i7)
      s45=s(i4,i5)+2d0*mb**2
      prop=s45**2*((s(i6,i7)-wmass**2)**2+(wmass*wwidth)**2)
      xl= 
     &+p2Dp6*p1Dp5/s267**2
     &*(p2Dp4*(p2Dp7+p6Dp7)+p4Dp6*(p2Dp7+p6Dp7)-p2Dp6*p4Dp7)
     &+p1Dp7*p2Dp4/s167**2
     &*(p1Dp5*(p1Dp6+p6Dp7)+p5Dp7*(p1Dp6+p6Dp7 )-p1Dp7*p5Dp6)

     &+(p1Dp5+p5Dp7)/s167/s267
     &*(p1Dp2*(p2Dp4*p6Dp7+p4Dp6*p6Dp7)-p1Dp6*(p2Dp4*p2Dp7+p2Dp7*p4Dp6))
 
     &+p1Dp7/s167/s267
     &*(p5Dp6*(p2Dp4*p2Dp7+p2Dp7*p4Dp6)-p2Dp5*(p2Dp4*p6Dp7+p4Dp6*p6Dp7))
 
     &+ p2Dp6/s167/s267
     &*(p1Dp6*p4Dp7*(p5Dp7+p1Dp5)-p1Dp4*(p5Dp7*p6Dp7+p1Dp5*p6Dp7)
     & +p1Dp7*(p4Dp5*p6Dp7-p4Dp7*p5Dp6-p1Dp5*p2Dp4))

      xr= 
     &+p2Dp6*p1Dp4/s267**2
     &*(p2Dp5*(p2Dp7+p6Dp7)+p5Dp6*(p2Dp7+p6Dp7)-p2Dp6*p5Dp7)
     &+p1Dp7*p2Dp5/s167**2
     &*(p1Dp4*(p1Dp6+p6Dp7)+p4Dp7*(p1Dp6+p6Dp7 )-p1Dp7*p4Dp6)

     &+(p1Dp4+p4Dp7)/s167/s267
     &*(p1Dp2*(p2Dp5*p6Dp7+p5Dp6*p6Dp7)-p1Dp6*(p2Dp5*p2Dp7+p2Dp7*p5Dp6))
 
     &+p1Dp7/s167/s267
     &*(p4Dp6*(p2Dp5*p2Dp7+p2Dp7*p5Dp6)-p2Dp4*(p2Dp5*p6Dp7+p5Dp6*p6Dp7))
 
     &+ p2Dp6/s167/s267
     &*(p1Dp6*p5Dp7*(p4Dp7+p1Dp4)-p1Dp5*(p4Dp7*p6Dp7+p1Dp4*p6Dp7)
     & +p1Dp7*(p4Dp5*p6Dp7-p5Dp7*p4Dp6-p1Dp4*p2Dp5))
 
 
      xlr=
     &+mb**2*p2Dp6/s267**2
     &*((p1Dp2+p1Dp6)*(p2Dp7+p6Dp7)-p2Dp6*p1Dp7)
     &+mb**2*p1Dp7/s167**2
     &*((p1Dp2+p2Dp7)*(p1Dp6+p6Dp7)-p2Dp6*p1Dp7)
     &+mb**2/s167/s267
     &* (+(p1Dp2*p6Dp7-p1Dp6*p2Dp7+p2Dp6*p1Dp7)*(p1Dp6+p2Dp7+p6Dp7)
     &   +(p1Dp2*p6Dp7-p1Dp6*p2Dp7-p2Dp6*p1Dp7)*p1Dp2)
     
      wbb=16d0*(xl+xr+xlr)/prop
      return
      end
