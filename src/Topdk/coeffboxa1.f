      subroutine coeffboxa1(mt,zp1,zp2,zp3,zp4,b0,bp,bm,bn,bc,bd,be,bg)
c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%
C-----This is an implementation of KM, Eq. (3.7) and Appendix A

C-----Coefficient for Box a1
      implicit none
      include 'constants.f'
      include 'eplog.f'
      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex b0,bp(5),bm(2),bn(2),bc,bd(2),be(2),bg(5)
      double complex b(8),p(0:8,5),m(0:8,5),n(0:8,2),
     & c(8),d(8,2),e(8,2),f(8),g(8,5)
      double precision t1,s,z1,z2,mt,mtsq,zt,zu,u1,y,x,Den,tt
      double precision uu,ddilog
      double complex cdot,wlog,wlogsq
      integer i,j

      s=2d0*dble(cdot(zp1,zp2))
      t1=2d0*dble(cdot(zp1,zp3))
      u1=-s-t1
      mtsq=mt**2
      tt=t1+mtsq
      uu=u1+mtsq
      z1=mtsq*s-t1**2
      z2=s+2d0*t1
      zt=2d0*mtsq+t1 
      zu=2d0*mtsq+u1
      y=sqrt(1d0-4d0*mtsq/s)
      x=(1d0-y)/(1d0+y)

      wlog=log(x)+im*pi
      wlogsq=log(x)**2+2d0*im*pi*log(x)

      Den=mtsq*s+s*t1+t1**2

      f(1)=wlog 
      f(2)=wlogsq 
      f(3)=-2d0*wlog*log(1d0+x)
     & +2d0*wlog*log(-t1/mtsq)-2d0*ddilog(-x)+3d0*pisqo6
      f(4)=pisqo6-ddilog(tt/mtsq) 
      f(5)=pisqo6 
      f(6)=log(-t1/mtsq) 
      f(7)=1d0 
      f(8)=wlog*log(1d0-x)+ddilog(x)

      b(1)=eplog*(1d0/y+y)/t1
      b(2)=(zu-y*u1)/2d0/Den
      b(3)=zt*zu/y/t1/Den
      b(4)=2d0*zu/Den
      b(5)=(4d0*y*u1-3d0*zu)/Den
      b(6)=czip
      b(7)=czip
      b(8)=-2d0*(1/y+y)/t1

      p(0,1)=(zt*Den+t1*z1*y**2)/Den**2
      p(1,1)=2d0*y*zt/Den
      p(3,1)=s*t1**2*y**3/Den**2
      p(4,1)=-2d0*zt*(2d0*mtsq/t1-s*t1*y**2/Den)/Den
      p(6,1)=2d0*zt**2/tt/Den
      p(7,1)=czip
      p(8,1)=czip
      p(2,1)=(p(0,1)+p(3,1))/2d0
      p(5,1)=-(3d0*p(0,1)+4d0*p(3,1))

      p(0,2)=4d0*mtsq*t1*(-z1*y**2+u1*z2*Den/s**2)/Den**3
      p(1,2)=-4d0*mtsq*y*(zt+z1/s)/Den**2
      p(3,2)=4d0*mtsq*t1*y*(Den*(2d0+t1/s)-s*t1*y**2)/Den**3
      p(4,2)=-8d0*mtsq*(s*y**2*(Den*tt+t1**2*zt)+Den*mtsq*z2*tt/t1)
     & /t1/Den**3
      p(6,2)=4d0*mtsq*(2/t1*Den+s*y**2-4d0*t1*u1/s-t1**2/tt)/Den**2
      p(7,2)=4d0*zt*u1/t1/s/Den
      p(8,2)=czip
      p(2,2)=(p(0,2)+p(3,2))/2d0
      p(5,2)=-(3d0*p(0,2)+4d0*p(3,2))

      p(0,3)=2d0*(2d0*mtsq*u1*z1*y**2-(2d0*mtsq*tt
     & -2d0*mtsq*u1*(u1**2+t1**2)/s**2
     & +t1**2*y**2)*Den)/Den**3
      p(1,3)=-4d0*y*(2d0*(2d0*mtsq-s)*zt*Den/s**2/y**2
     & -2d0*mtsq*s*y**2-3d0*(mtsq/s+1d0)*Den+2d0*mtsq*t1*z2/s)/Den**2
      p(3,3)=-2d0*zt*(t1*(y+2d0*mtsq*z2/s**2/y)*Den+2d0*mtsq*u1*z2*y)
     & /Den**3
      p(4,3)=-4d0*(2d0*mtsq**2*z2*s*y**2+(4d0*mtsq*tt*Den/t1**2
     & -zt*(2d0*mtsq**2/t1-t1)+2d0*mtsq*z2)*Den)/Den**3
      p(6,3)=4d0*zt*(2d0*(1/t1*Den-z1/s+2d0*z2)
     & +3d0*t1*u1/tt-t1**2*(mtsq+2d0*t1)/tt**2)/Den**2
      p(7,3)=4d0*zt*(u1/s-mtsq/tt)/t1/Den
      p(8,3)=czip
      p(2,3)=(p(0,3)+p(3,3))/2d0
      p(5,3)=-(3d0*p(0,3)+4d0*p(3,3))

      p(0,4)=2d0*(2d0*t1*(mtsq*t1**2-s*tt**2)*y**2-(6*mtsq/s*Den+t1*zt
     & +2d0*mtsq/s**2*t1**2*z2)*Den)/Den**3
      p(1,4)=-4d0*y*(zt*Den/s/y**2+3d0*tt**2-t1**2*(mtsq-t1)/s)/Den**2
      p(3,4)=2d0*t1*(4d0*mtsq/s**2/y*Den**2+2d0*mtsq*u1/s*y*Den
     & -t1*(3d0*s*tt+t1**2)*y**3)/Den**3
      p(4,4)=4d0*zt*((2d0*mtsq**2/t1-t1)*Den-2d0*s*t1*tt*y**2)/Den**3
      p(6,4)=-8*(zt/s+t1**2*y**2/Den)/Den
      p(7,4)=-4d0*zt/s/Den
      p(8,4)=czip
      p(2,4)=(p(0,4)+p(3,4))/2d0
      p(5,4)=-(3d0*p(0,4)+4d0*p(3,4))

      p(0,5)=4d0*mtsq*t1*(u1/s**2*z2*Den-z1*y**2)/Den**3
      p(1,5)=-4d0*mtsq*y*(zt+z1/s)/Den**2
      p(3,5)=4d0*mtsq*t1*y*((2d0+t1/s)*Den-s*t1*y**2)/Den**3
      p(4,5)=-8*mtsq*(tt/t1**2*(Den-t1*zt)*Den+s*t1*zt*y**2)/Den**3
      p(6,5)=4d0*mtsq*(2/t1*Den+s*y**2-4d0*t1*u1/s-t1**2/tt)/Den**2
      p(7,5)=4d0*zt*u1/t1/s/Den
      p(8,5)=czip
      p(2,5)=(p(0,5)+p(3,5))/2d0
      p(5,5)=-(3d0*p(0,5)+4d0*p(3,5))


      m(0,1)=(2d0*(tt-mtsq*u1/s)-t1**2*zt/Den)/Den
      m(1,1)=-2d0*(4d0*mtsq/s-3d0+t1*zt/Den)/s/y
      m(3,1)=t1*(4d0*mtsq*u1/s**2/y+(2d0*s*tt+t1**2)*y/Den)/Den
      m(4,1)=-2d0*(2d0*mtsq*tt/t1-zt*(2d0*s*tt+t1**2)/Den)/Den
      m(6,1)=2d0*(zt-2d0*t1*u1/s)/Den
      m(7,1)=czip
      m(8,1)=czip
      m(2,1)=(m(0,1)+m(3,1))/2d0
      m(5,1)=-(3d0*m(0,1)+4d0*m(3,1))

      m(0,2)=(-2d0*mtsq*(1d0-t1/s)+t1*u1*zt/Den)/Den
      m(1,2)=-2d0*(-y+t1*zu/y/Den)/s
      m(3,2)=-t1**2*(4d0*mtsq/s**2/y-u1*y/Den)/Den
      m(4,2)=-2d0*(2d0*mtsq*(mtsq/t1+2d0)-t1*u1*zt/Den)/Den
      m(6,2)=2d0*(2d0*Den/s+2d0*s+t1+t1**2/tt)/Den
      m(7,2)=czip
      m(8,2)=czip
      m(2,2)=(m(0,2)+m(3,2))/2d0
      m(5,2)=-(3d0*m(0,2)+4d0*m(3,2))


      n(0,1)=(2d0*mtsq*(1d0-t1/s-s*t1*y**2/Den)-t1**2*zu/Den)/Den
      n(1,1)=-2d0*(2d0-2d0*mtsq*t1/s*z2/Den+3d0*mtsq*s*y**2/Den)/s/y
      n(3,1)=(2d0*(tt+2d0*mtsq*t1**2/s**2)/y
     & -(2d0*mtsq*z1-t1**2*u1)*y/Den)/Den
      n(4,1)=2d0*(2d0*mtsq*(mtsq/t1+2d0)*Den
     & -2d0*mtsq*t1*s*y**2-t1**2*zu)/Den**2
      n(6,1)=-2d0*(6*mtsq+2d0*u1**2/s-5*mtsq*t1/tt)/Den
      n(7,1)=czip
      n(8,1)=czip
      n(2,1)=(n(0,1)+n(3,1))/2d0
      n(5,1)=-(3d0*n(0,1)+4d0*n(3,1)) 

      n(0,2)=t1*(2d0*mtsq/s-1d0-z1*y**2/Den+t1*zu/Den)/Den
      n(1,2)=-2d0*(2d0*t1/s-2d0*t1**2/s*zt/Den+3d0*s*tt*y**2/Den)/s/y
      n(3,2)=(2d0*(-2d0*mtsq*t1*u1/s**2-tt)/y
     & -(2d0*mtsq*s*tt-t1**2*zt)*y/Den)/Den
      n(4,2)=2d0*(2d0*mtsq*tt*Den/t1-2d0*s*t1*tt*y**2-t1**2*zt)/Den**2
      n(6,2)=-2d0*(3d0*zt+2d0*t1**2/s)/Den
      n(7,2)=czip
      n(8,2)=czip
      n(2,2)=(n(0,2)+n(3,2))/2d0
      n(5,2)=-(3d0*n(0,2)+4d0*n(3,2))

      c(1)=4d0/s/y
      c(3)=-(s*y+zt/y+2d0*t1*z2/s/y)/Den
      c(4)=2d0*(s+3d0*t1)/Den
      c(6)=czip
      c(7)=czip
      c(8)=czip
      c(2)=(c(3)+c(4)/2d0)/2d0
      c(5)=-(4d0*c(3)+3d0*c(4)/2d0)

      d(1,1)=4d0*zt/s/y/Den
      d(3,1)=-2d0*tt*(2d0/s/y+s*y/Den)/Den
      d(4,1)=4d0*tt*z2/Den**2
      d(6,1)=-4d0/Den
      d(7,1)=czip
      d(8,1)=czip
      d(2,1)=(d(3,1)+d(4,1)/2d0)/2d0
      d(5,1)=-(4d0*d(3,1)+3d0*d(4,1)/2d0)

      d(1,2)=4d0*zu/s/y/Den
      d(3,2)=2d0*(2d0*mtsq*z1-t1*u1*z2)/s/y/Den**2
      d(4,2)=4d0*t1*zu/Den**2
      d(6,2)=-4d0*mtsq/tt/Den
      d(7,2)=czip
      d(8,2)=czip
      d(2,2)=(d(3,2)+d(4,2)/2d0)/2d0
      d(5,2)=-(4d0*d(3,2)+3d0*d(4,2)/2d0)

      e(1,1)=d(1,2)
      e(2,1)=d(2,2)
      e(3,1)=d(3,2)
      e(4,1)=d(4,2)
      e(5,1)=d(5,2)
      e(6,1)=d(6,2)
      e(7,1)=d(7,2)
      e(8,1)=d(8,2)

      e(1,2)=d(1,1)
      e(2,2)=d(2,1)
      e(3,2)=d(3,1)
      e(4,2)=d(4,1)
      e(5,2)=d(5,1)
      e(6,2)=d(6,1)
      e(7,2)=d(7,1)
      e(8,2)=d(8,1)

      g(1,1)=2d0*(t1*y/Den-4/s/y)
      g(3,1)=t1*(4d0*zt*Den/s/y-t1*z2*y)/Den**2
      g(4,1)=-2d0*(6*tt*Den-s*t1**2*y**2)/Den**2
      g(6,1)=2d0*t1*zt/tt/Den
      g(7,1)=czip
      g(8,1)=czip
      g(2,1)=(g(3,1)+g(4,1)/2d0)/2d0
      g(5,1)=-(4d0*g(3,1)+3d0*g(4,1)/2d0)

      g(1,2)=4d0*y*(1d0+(mtsq*z2/s-u1)/s/y**2-2d0*t1**2*u1/s/Den)/Den
      g(3,2)=2d0*(-2d0*mtsq*t1*(6*Den-t1*z2)/s**2/y+(2d0*mtsq*s+t1**2)*
     & y+2d0*mtsq*t1**2*z2*y/Den)/Den**2
      g(4,2)=4d0*(2d0*mtsq*(mtsq-s)+zt**2-2d0*mtsq*s*t1**2*y**2/Den)
     &       /Den**2
      g(6,2)=4d0*(4d0*mtsq*s-t1*u1*(3d0+4d0*t1/s)
     & +3d0*mtsq*t1**2/tt)/Den**2
      g(7,2)=4d0*u1/s/Den
      g(8,2)=czip
      g(2,2)=(g(3,2)+g(4,2)/2d0)/2d0
      g(5,2)=-(4d0*g(3,2)+3d0*g(4,2)/2d0)

      g(1,3)=-4d0*(Den*(mtsq+6*mtsq*t1/s+s)/s/y+2d0*t1*uu*y)/Den**2
      g(3,3)=-2d0*((8d0*mtsq+t1)*Den/s/y+2d0*mtsq*t1**2*z2/s**2/y
     & -t1*(3d0*mtsq+t1**2/s)*
     & y+2d0*mtsq*t1*u1*z2*y/Den)/Den**2
      g(4,3)=4d0*(6*mtsq**2-t1*zt+2d0*mtsq*s*t1*u1*y**2/Den)/Den**2
      g(5,3)=-(4d0*g(3,3)+3d0*g(4,3)/2)
      g(6,3)=4d0*(4d0*Den+2d0*t1*(2d0+t1/s)*z2
     & -3d0*t1**2*z2/tt+t1**3*zt/tt**2)/Den**2
      g(7,3)=-4d0*(1d0-mtsq*u1/Den)/tt/s
      g(8,3)=czip
      g(2,3)=(g(3,3)+g(4,3)/2d0)/2d0

      g(1,4)=4d0*(mtsq*Den*(2d0*t1/s-1d0)/y+2d0*t1**3*y)/s/Den**2
      g(3,4)=2d0*(2d0*(2d0*mtsq*u1/s+z2)*Den/s/y
     & -2d0*mtsq*t1**2*z2/s**2/y
     & +t1*zt*y+2d0*t1**3*zt*y/Den)/Den**2
      g(4,4)=4d0*(2d0*tt*(3d0*mtsq-s)+t1**2
     & -2d0*s*t1**2*tt*y**2/Den)/Den**2
      g(5,4)=-(4d0*g(3,4)+3d0*g(4,4)/2)
      g(6,4)=8d0*(s*tt-2d0*t1**2*u1/s)/Den**2
      g(7,4)=-4d0*t1/s/Den
      g(8,4)=czip
      g(2,4)=(g(3,4)+g(4,4)/2d0)/2d0

      g(1,5)=4d0*(Den*(mtsq*z2/s-t1)/y-2d0*t1**2*u1*y)/s/Den**2
      g(3,5)=2d0*(2d0*u1*(2d0*mtsq-s)/s**2/y*Den+t1*(2d0*tt-t1*u1/s)*
     & y+t1**2*zt/s/y+2d0*mtsq*t1**2*z2*y/Den)/Den**2
      g(4,5)=4d0*(6d0*mtsq**2-2d0*mtsq*s-t1**2
     & +2d0*s*t1*tt*u1*y**2/Den)/Den**2
      g(5,5)=-(4d0*g(3,5)+3d0*g(4,5)/2)
      g(6,5)=4d0*(2d0*Den-t1*(1d0+4d0*t1/s)*u1+mtsq*t1**2/tt)/Den**2
      g(7,5)=4d0*u1/s/Den
      g(8,5)=czip
      g(2,5)=(g(3,5)+g(4,5)/2d0)/2d0

      b0=czip
      bc=czip
      do j=1,2
        bm(j)=czip
        bn(j)=czip
        bd(j)=czip
        be(j)=czip
      enddo
      do j=1,5
        bp(j)=czip
        bg(j)=czip
      enddo
      
      do i=1,8
        b0=b0+f(i)*b(i)
        bc=bc+f(i)*c(i)
        do j=1,2
          bm(j)=bm(j)+f(i)*m(i,j)
          bn(j)=bn(j)+f(i)*n(i,j)
          bd(j)=bd(j)+f(i)*d(i,j)
          be(j)=be(j)+f(i)*e(i,j)
        enddo
        do j=1,5
          bp(j)=bp(j)+f(i)*p(i,j)
          bg(j)=bg(j)+f(i)*g(i,j)
        enddo
      enddo

      return
      end
