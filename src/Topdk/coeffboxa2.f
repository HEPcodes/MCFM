      subroutine coeffboxa2(mt,zp1,zp2,zp3,zp4,b0,bp,bm,bn,bc,bd,be,bg)
c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%
C-----This is an implementation of KM, Eq. (3.9) and Appendix A

C-----Coefficient for Box a2
      implicit none
      include 'constants.f'
      include 'eplog.f'
      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex b0,bp(5),bm(2),bn(2),bc,bd(2),be(2),bg(5)
      double complex f(9),b(9),p(8,5),m(8,5),n(8,2),
     & c(8),d(8,2),e(8,2),g(8,5)
      double precision t1,s,z1,z2,mt,mtsq,zt,zu,u1,y,x,Den,tt,zeta2
      double precision uu,ddilog
      double complex cdot,slog,slogsq,wlogsq,fll,flld
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
      zeta2=pisqo6

      slog=log(s/mtsq)-im*pi
      slogsq=log(s/mtsq)**2-2d0*im*pi*log(s/mtsq)
      wlogsq=log(x)**2+2d0*im*pi*log(x)

      Den=mtsq*s+s*t1+t1**2

c--- Eq. (2.28)
      flld=4d0*ddilog(-x)+wlogsq+2d0*zeta2
c--- Eq. (2.37)
      fll=log(-t1/mtsq)**2+ddilog(tt/mtsq) 

      f(1)=slog
      f(2)=slogsq
      f(3)=slogsq+4d0*fll-4d0*slog*log(-t1/mtsq)+2d0*zeta2 
      f(4)=flld 
      f(5)=zeta2
      f(6)=log(-t1/mtsq) 
      f(7)=1d0
      f(8)=fll 
      f(9)=slog*log(-t1/mtsq)


      b(1)=-2d0*eplog/t1
      b(2)=czip
      b(3)=(s/2d0-t1)/Den
      b(4)=-(z2/2d0+zt)/y/Den
      b(5)=-10d0/t1
      b(6)=-2d0*zt/t1/tt
      b(7)=2d0*(eplog+epsqlog+2d0)/t1
      b(8)=-4d0/t1
      b(9)=4d0/t1

      p(1,1)=-eplog/s-2d0*(mtsq+t1**2/s)/Den
      p(2,1)=1d0/2d0/s
      p(3,1)=t1**2*(z2/2d0/Den-2d0/s)/Den
      p(4,1)=(mtsq/s*(s-4d0*t1)*Den+s*t1*zt/2d0-t1**3/s*z2)/y/Den**2
      p(5,1)=-4d0/s
      p(6,1)=-8d0*eplog/s-2d0*t1*(2d0-t1/tt)/Den
      p(7,1)=(5d0*epsqlog+2d0*eplog+4d0)/s
      p(8,1)=8d0/s

      p(1,2)=-4d0*mtsq*(2d0*t1/Den-1d0/s/y**2)/Den
      p(2,2)=czip
      p(3,2)=2d0*z1*mtsq*t1/Den**3
      p(4,2)=-2d0*mtsq*(2d0*mtsq/s**2/y**2+t1*(mtsq*(s+4d0*t1)+t1**2)
     & /Den**2)/y/Den
      p(5,2)=czip
      p(6,2)=4d0*(mtsq*(s+3d0*t1)+mtsq**2*t1/tt)/Den**2
      p(7,2)=-4d0*u1/s/Den
      p(8,2)=czip

      p(1,3)=4d0*(2d0*mtsq*u1/Den+(mtsq+(3d0+4d0*t1/s)*(2d0*mtsq-s))
     & /s/y**2)/Den
      p(2,3)=czip
      p(3,3)=-zt*(2d0*mtsq*s*u1/Den-t1)/Den**2
      p(4,3)=(4d0*mtsq/s*(mtsq+z2)/s/y**2-2*mtsq/s-2d0*mtsq*s/Den
     & +t1**2/Den-2d0*mtsq**2*s**2*y**2/Den**2)/y/Den
      p(5,3)=czip
      p(6,3)=4d0*(Den*t1**2/tt**2+2d0*t1**2/tt*(s-t1)+10d0*mtsq*s
     & +4d0*t1*zt-4d0*t1**2*u1/s)/Den**2
      p(7,3)=4d0*(2d0+t1/s-t1/tt)/Den
      p(8,3)=czip

      p(1,4)=4d0*((3d0*mtsq+2d0*t1)/y**2+2d0*t1**3/Den)/s/Den
      p(2,4)=czip
      p(3,4)=-t1*(2d0*mtsq+3d0*t1+2d0*t1**2*zt/Den)/Den**2
      p(4,4)=-(3d0*mtsq+t1+(3d0*mtsq+2d0*t1)/y**2
     & +t1*(mtsq*s-3d0*t1**2)/Den+2d0*t1**4*s*y**2/Den**2)/s/y/Den
      p(5,4)=czip
      p(6,4)=-8d0*t1**2*z2/s/Den**2
      p(7,4)=4d0*t1/s/Den
      p(8,4)=czip

      p(1,5)=-4d0*mtsq*(2d0*t1/Den-1d0/s/y**2)/Den
      p(2,5)=czip
      p(3,5)=2d0*mtsq*t1*z1/Den**3
      p(4,5)=-2d0*mtsq*(2d0*mtsq/s**2/y**2+t1*(mtsq*z2+t1*zt)/Den**2)
     & /y/Den
      p(5,5)=czip
      p(6,5)=4d0*(mtsq*s+4d0*mtsq*t1-t1**2+t1**3/tt)/Den**2
      p(7,5)=-4d0*u1/s/Den
      p(8,5)=czip


      m(1,1)=2d0*eplog/s+(tt-t1**2/s-4d0*mtsq*zt/s/y**2)/Den
      m(2,1)=-1d0/s
      m(3,1)=t1*(z2/s+t1*(2d0*mtsq+t1/2d0)/Den)/Den
      m(4,1)=(2d0*mtsq*Den/s/y**2-mtsq*t1*z2/s/y**2+2d0*mtsq*z2+t1
     & *(2d0*mtsq-u1)+3d0/2d0*s*t1**3*y**2/Den+t1**3*zt/Den)/s/y/Den
      m(5,1)=8d0/s
      m(6,1)=8*eplog/s+2d0*t1*(3d0+4d0*t1/s)/Den
      m(7,1)=-3d0*(2d0*epsqlog+eplog+2d0)/s
      m(8,1)=-8d0/s

      m(1,2)=2d0*eplog/s+(3d0*(mtsq+t1*u1/s)-4d0*mtsq*t1*z2/s**2/y**2)
     & /Den
      m(2,2)=-1d0/s
      m(3,2)=t1*(1d0+4d0*t1/s+t1*u1*(3d0/2d0+2d0*t1/s)/Den)/Den
      m(4,2)=-(tt-t1*zt/s-mtsq*t1*z2/s**2/y**2+t1**2*u1*z2/2d0/s/Den
     & +mtsq*s*t1*y**2/Den)/y/Den
      m(5,2)=8d0/s
      m(6,2)=2d0*(eplog*(1d0/t1+4d0/s)-2d0*(z1/t1-t1*z2/s)/Den
     & +mtsq*t1/tt/Den)
      m(7,2)=-(epsqlog/t1+6d0*epsqlog/s-2d0*eplog/t1
     & +3d0*eplog/s-2d0*(2d0/t1-3d0/s))
      m(8,2)=-2d0*(1d0/t1+4d0/s)


      n(1,1)=-2d0*eplog/s+(5d0*mtsq-2d0*t1+t1**2/s+4d0*zt/y**2
     & +t1*z2/s/y**2)/Den
      n(2,1)=1d0/s
      n(3,1)=-(2d0*mtsq-3d0*t1)/Den+t1**2*u1*(1d0/2d0-2d0*u1/s)/Den**2
      n(4,1)=-(4d0*zt*(s*zt-3d0*mtsq*t1)/s**2/y**2
     & -6d0*mtsq**2*z2/Den-mtsq*s*t1*y**2/Den-t1*(2d0*mtsq/s+1d0))
     & /2d0/y/Den
      n(5,1)=-m(5,2)
      n(6,1)=-m(6,2)-16d0*mtsq**2/tt/Den
      n(7,1)=-m(7,2)
      n(8,1)=-m(8,2)

      n(1,2)=-2d0*eplog/s-((8d0*mtsq+s)*zt/s/y**2-mtsq-t1**2/s)/Den
      n(2,2)=1d0/s
      n(3,2)=-(2d0/s-s*t1*tt/Den**2-t1**2*(2d0*mtsq+t1/2d0)/Den**2)
      n(4,2)=(3d0*mtsq*zt/s/y**2-t1*(2d0*tt/s+1d0-(2d0*mtsq-t1/2d0)
     & *zt/Den))/y/Den
      n(5,2)=-m(5,1)
      n(6,2)=-m(6,1)-16d0*mtsq/Den
      n(7,2)=-m(7,1)
      n(8,2)=-m(8,1)

      c(1)=-2d0/s/y**2
      c(2)=czip
      c(3)=-s/2/Den
      c(4)=1d0/2d0/s/y**3+z2/2d0/y/Den
      c(5)=czip
      c(6)=-2/tt
      c(7)=czip
      c(8)=czip

      d(1,1)=-8d0*zt/s/y**2/Den  
      d(2,1)=czip    
      d(3,1)=(-s*tt+t1**2)/Den**2
      d(4,1)=(z2/s/y**2 + 2d0*t1*zt/Den)/y/Den
      d(5,1)=czip
      d(6,1)=-8d0/Den
      d(7,1)=czip
      d(8,1)=czip

      d(1,2)=-8d0*zu/s/y**2/Den  
      d(2,2)=czip    
      d(3,2)=-(mtsq*s+t1*u1)/Den**2
      d(4,2)=z2*(-1d0/s/y**2 + 2d0*mtsq/Den)/y/Den
      d(5,2)=czip
      d(6,2)=-8d0*mtsq/tt/Den 
      d(7,2)=czip  
      d(8,2)=czip

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


      g(1,1)=2d0*zt/y**2/Den
      g(2,1)=czip
      g(3,1)=-s*t1**2/2d0/Den**2
      g(4,1)=-(2d0/s+mtsq*z2/s/y**2/Den+s*t1*zt/2d0/Den**2)/y
      g(5,1)=czip
      g(6,1)=-2d0*t1**2/tt/Den
      g(7,1)=czip
      g(8,1)=czip

      g(1,2)=-4d0*(2d0+t1/s-2d0*uu/s/y**2-3d0*mtsq*z2/s**2/y**4
     & -2d0*mtsq*t1*z2/s/y**2/Den)/Den
      g(2,2)=czip
      g(3,2)=s*(4d0*mtsq+t1*zt/s+2d0*t1**3*u1/s/Den)/Den**2
      g(4,2)=(2d0*(2d0*mtsq**2+2d0*mtsq*u1-t1*tt)/Den+
     & 2d0*(4d0*mtsq**2/s+mtsq+2d0*mtsq*t1/s-t1)/s/y**2
     & +t1**2*z2/s/y**2/Den-12d0*mtsq**2*z2/s**3/y**4
     & -2d0*t1**3*zu/Den**2)/y/Den
      g(5,2)=czip
      g(6,2)=4d0*(2d0*Den*(3d0+t1/s)+3d0*t1*u1+5d0*mtsq*t1**2/tt)/Den**2
      g(7,2)=4d0*zu/s/y**2/Den
      g(8,2)=czip

      g(1,3)=-4d0*(2d0+t1/s+(2d0*mtsq-u1)/s/y**2+3d0*mtsq*z2/s**2/y**4
     & -2d0*t1*uu*z2/s/y**2/Den)/Den
      g(2,3)=czip
      g(3,3)=s*(4d0*mtsq+t1*zt/s-2d0*mtsq*t1*u1/Den)/Den**2
      g(4,3)=(2d0*(2d0*mtsq**2-s*zt-2d0*t1**2+t1**3/s)/Den
     & +2d0*(4d0*mtsq**2/s+3d0*mtsq+16d0*mtsq*t1/s)/s/y**2
     & -t1**2*z2/s/y**2/Den+12d0*mtsq**2*z2/s**3/y**4
     & -2d0*t1**2*(t1*zu+z2*u1)/Den**2)/y/Den
      g(5,3)=czip
      g(6,3)=4d0*(2d0*(2d0+t1/s)+2d0*t1**2*uu/tt/Den-t1**2/tt**2)/Den
      g(7,3)=4d0*(t1/tt+zt/s/y**2)/Den
      g(8,3)=czip

      g(1,4)=4d0*(2d0*mtsq*(s-2d0*u1)/s**2/y**2-3d0*mtsq*z2/s**2/y**4
     & +2d0*t1**2*zt/s/y**2/Den)/Den
      g(2,4)=czip
      g(3,4)=-(2*tt*u1+t1**2+2*t1**4/Den)/Den**2
      g(4,4)=(3*t1**2/Den+(3*mtsq+18*mtsq*t1/s-2*s-8*t1)/s/y**2-
     & 8*mtsq**2*z1/s**2/y**2/Den+3*mtsq*z2/s**2/y**4+
     & 2*t1**4*z2/s/Den**2)/y/Den
      g(5,4)=czip
      g(6,4)=8d0*(2d0+t1/s+t1**2/Den)/Den
      g(7,4)=4d0*zt/s/y**2/Den
      g(8,4)=czip

      g(1,5)=4d0*(2d0*mtsq*z2/s+3d0*tt/y**2-6d0*mtsq*t1/s/y**2
     &      +2d0*t1**2*zu/Den)/s/y**2/Den
      g(2,5)=czip
      g(3,5)=2d0/Den+t1*zt/Den**2+2d0*t1**3*u1/Den**3
      g(4,5)=(2d0*(2d0*mtsq**2+t1**2*u1/s-mtsq*t1**2*z2/Den)/Den
     & +t1**2*z2/s/y**2/Den
     & -2d0*(2d0*mtsq*(mtsq-s)+(mtsq+s)*z2+12d0*mtsq**2*zt/s/y**2)
     & /s**2/y**2)/y/Den
      g(5,5)=czip
      g(6,5)=4d0*(-2d0*Den*u1/s+t1*z2-t1**3/tt)/Den**2
      g(7,5)=4d0*zu/s/y**2/Den
      g(8,5)=czip
      
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
c--- extra term for b0
      b0=b0+f(9)*b(9)
      
      return
      end
