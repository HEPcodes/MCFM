      subroutine coeffboxa4(mt,zp1,zp2,zp3,zp4,b0,bp,bm,bn,bc,bd,be,bg)
c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%
C-----This is an implementation of KM, Eq. (3.11) and Appendix A

C-----Coefficient for Box a4
      implicit none
      include 'constants.f'
      include 'eplog.f'
      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex b0,bp(5),bm(2),bn(2),bc,bd(2),be(2),bg(5)
      double complex f(9),b(9),p(9,5),m(9,5),n(9,2),
     & c(9),d(9,2),e(9,2),g(9,5)
      double precision t1,s,z1,mt,mtsq,zt,zu,u1,Den,tt,zeta2
      double precision uu,ddilog
      double complex cdot,fll,fllu,fdl,fdlu
      integer i,j
      
      s=2d0*dble(cdot(zp1,zp2))
      t1=2d0*dble(cdot(zp1,zp3))
      u1=-s-t1
      mtsq=mt**2
      tt=t1+mtsq
      uu=u1+mtsq
      z1=mtsq*s-t1**2
      zt=2d0*mtsq+t1 
      zu=2d0*mtsq+u1
      zeta2=pisqo6

c      slog=log(s/mtsq)-im*pi
c      slogsq=log(s/mtsq)**2-2d0*im*pi*log(s/mtsq)
c      wlogsq=log(x)**2+2d0*im*pi*log(x)

      Den=mtsq*s+s*t1+t1**2

c--- Eq. (2.31)
      fdl=zeta2-ddilog(tt/mtsq)
      fdlu=zeta2-ddilog(uu/mtsq)

c--- Eq. (2.37)
      fll=log(-t1/mtsq)**2+ddilog(tt/mtsq) 
      fllu=log(-u1/mtsq)**2+ddilog(uu/mtsq) 

      f(1)=fll 
      f(2)=fllu 
      f(3)=fdl+fdlu-fll-fllu+2d0*log(-t1/mtsq)*log(-u1/mtsq)-4d0*zeta2 
      f(4)=fdl 
      f(5)=zeta2 
      f(6)=log(-t1/mtsq) 
      f(7)=log(-u1/mtsq) 
      f(8)=log(-t1/mtsq)*log(-u1/mtsq) 
      f(9)=1d0


      b(1)=2d0*(2d0-u1/t1)/s
      b(2)=2d0*(2d0-t1/u1)/s
      b(3)=(2d0*mtsq-3d0*t1*u1/s+s)/Den
      b(4)=czip
      b(5)=8d0*(3d0/s - s/t1/u1)  
      b(6)=2d0*eplog/u1 + zt/t1/tt
      b(7)=2d0*eplog/t1 + zu/u1/uu
      b(8)=-b(5)/2d0
      b(9)=(epsqlog + eplog + 2d0)*s/t1/u1

      p(1,1)= - 2d0*(zu/u1 + 2d0)/s
      p(2,1)=-2d0*(zu - 2d0*t1)/s/u1
      p(3,1)= - 2d0*zu/s/u1 + 2d0*t1/Den - t1*u1**2*zt/s/Den**2
      p(4,1)=-4d0*mtsq/t1/u1
      p(5,1)=4d0*p(2,1)
      p(6,1)= - 4d0*eplog/u1 - 2d0*mtsq*zt/tt/Den
      p(7,1)=2d0*zu*(t1/u1+2d0*mtsq/uu)/Den
      p(8,1)=-2d0*p(2,1)
      p(9,1)=2d0/u1*(2d0 + epsqlog + eplog)

      p(1,2)=8d0*mtsq/s**2/u1
      p(2,2)=p(1,2)
      p(3,2)=4d0*mtsq*(2d0*mtsq*s*(Den/u1-t1)*Den + t1*u1**2*z1)
     & /s**2/Den**3
      p(4,2)=-8d0*mtsq/t1**2/u1
      p(5,2)=4d0*p(1,2)
      p(6,2)=-4d0*mtsq*( 2d0*(1d0/t1-1d0/s)*Den-u1
     & *(3d0*t1-u1)/s + mtsq*t1/tt )
     & /Den**2
      p(7,2)=-4d0*mtsq*(mtsq**2/uu+mtsq-t1**2/s+t1*u1/s)/Den**2
      p(8,2)=-2d0*p(1,2)
      p(9,2)=4d0*zt/t1/Den

      p(1,3)=4d0*(2d0*mtsq/u1 - 3d0)/s**2 
      p(2,3)=p(1,3)
      p(3,3)=p(1,3) - 2d0*( 2d0*t1*u1/Den+(4d0*mtsq*u1**3-t1*u1**2*zt)
     & /Den**2
     &  - 2d0*mtsq*t1*u1**3*(t1-u1)/Den**3)/s**2
      p(4,3)=-8d0*mtsq/t1**2/u1
      p(5,3)=4d0*p(1,3)
      p(6,3)=-4d0*(2d0*mtsq*((1/t1-1/s)*Den+u1*(t1+5d0*u1)/s)
     & +t1*u1*(5d0*mtsq**2+3d0*t1*u1)/s/tt 
     & + mtsq*t1**2*(mtsq**2+t1*u1)/s/tt**2)/Den**2
      p(7,3)=-4d0*(2d0*mtsq+5d0*u1+2d0*t1*u1*zu/Den)/s/Den
      p(8,3)=-2d0*p(1,3)
      p(9,3)=4d0*mtsq*zt/t1/tt/Den

      p(1,4)=4d0*(2d0*mtsq*t1/u1**2 + 4d0*mtsq/u1 + 3d0)/s**2  
      p(2,4)=p(1,4)
      p(3,4)=p(1,4) + 2d0*(t1*(4d0*mtsq*t1*(t1-u1)
     &     -u1**2*(2d0*mtsq+3d0*t1))*Den
     &     + 2d0*mtsq*t1**3*u1*(t1-u1))/s**2/Den**3
      p(4,4)=-8d0*mtsq/t1/u1**2
      p(5,4)=4d0*p(1,4)    
      p(6,4)=4d0*(2d0*mtsq+5d0*t1+2d0*t1*u1*zt/Den)/s/Den
      p(7,4)=4d0*(2d0*mtsq*((s/u1-1d0)*Den+t1*(u1+5d0*t1))
     & + u1*t1*(5d0*mtsq**2+3d0*t1*u1)/uu
     &  + mtsq*u1**2*(mtsq**2+t1*u1)/uu**2)/s/Den**2
      p(8,4)=-2d0*p(1,4)
      p(9,4)=-4d0*mtsq*zu/u1/uu/Den

      p(1,5)=8d0*mtsq*(t1 + 2d0*u1)/s**2/u1**2
      p(2,5)=p(1,5)
      p(3,5)=p(1,5) - 4d0*mtsq*(t1*u1*(t1-2d0*u1)*Den 
     & + t1**2*u1**2*(t1-u1))/s**2/Den**3
      p(4,5)=-8d0*mtsq/t1/u1**2        
      p(5,5)=4d0*p(1,5)
      p(6,5)=4d0*mtsq*(2d0*mtsq*s-u1**2 + t1**2*uu/tt 
     & + 2d0*mtsq*t1*u1/tt)/s/Den**2
      p(7,5)=4d0*mtsq*(-2d0*mtsq*s*(t1/u1+2d0)+3d0*t1**2 
     & - u1**2*(mtsq-t1)/uu)/s/Den**2
      p(8,5)=-2d0*p(1,5)
      p(9,5)=-4d0*zu/u1/Den


      m(1,1)=-2d0*(2d0*mtsq + 3d0*t1*u1/s)/s/u1
      m(2,1)=m(1,1) + 2d0/u1
      m(3,1)=(-2d0*(2d0*mtsq*s/u1+t1)+2d0*s*(2d0*mtsq*u1-s*t1)/Den
     & -t1**3*u1**2/Den**2)/s**2
      m(4,1)=-4d0*mtsq/t1/u1
      m(5,1)=4d0*m(1,1) + 16d0/u1
      m(6,1)=4d0*eplog/u1 - 2d0*(2d0*Den - t1*u1)/s/Den
      m(7,1)=2d0*eplog/u1 + 2d0*(2d0*Den + t1**2)/s/Den
      m(8,1)=-m(5,1)/2d0     
      m(9,1)=-3d0*epsqlog/u1

      m(1,2)=2*(-2d0*mtsq*s/u1-5d0*t1-u1+u1**2/t1)/s**2
      m(2,2)=2d0*(-2d0*mtsq/s - 1d0 + 3d0*t1**2/s**2)/u1
      m(3,2)=-(2d0*s*(2d0*mtsq/u1-1d0) 
     & + 2d0*t1*(3d0*mtsq*s+t1**2)/Den - t1**2*u1**3/Den**2)/s**2
      m(4,2)=m(4,1)
      m(5,2)=4d0*m(2,2)
      m(6,2)=2d0*eplog*(2d0/u1-1/t1) 
     & - 2d0*(2d0*mtsq*(2d0+u1/t1)+t1*(t1-2d0*u1)/s + t1**2/tt)/Den
      m(7,2)=2d0*(t1*(-2d0*mtsq*s/u1+t1-2d0*u1)/s + mtsq*u1/uu)/Den
      m(8,2)=-2d0*m(2,2)
      m(9,2)=-epsqlog*(2d0/u1-1/t1) + 2d0*s*eplog/t1/u1
     & + 4*s/t1/u1


      n(1,1)=6d0*(2d0*t1 + u1)/s**2
      n(2,1)=-2d0*(3d0*t1**2/s**2 - 1d0)/u1
      n(3,1)=(2d0*(2d0*t1+3*u1) + 2d0*t1*(t1**2+t1*u1+5d0*u1**2)/Den 
     & + 5d0*t1**2*u1**3/Den**2)/s**2
      n(4,1)=czip   
      n(5,1)=4d0*n(2,1)
      n(6,1)=-4d0*eplog/u1+2d0*(Den/s+2d0*u1**2/s+7d0*mtsq**2/tt)/Den
      n(7,1)=2d0*(t1*(2d0-t1/u1)*zu+6d0*mtsq*u1 - mtsq*s*u1/uu)/s/Den
      n(8,1)=-2d0*n(2,1)
      n(9,1)=2d0*(epsqlog + eplog + 2d0)/u1

      n(1,2)=6d0*t1/s**2
      n(2,2)=-2d0*(2d0/u1 - 3d0*t1/s**2)
      n(3,2)=-t1*(2d0 - 2d0*(t1**2+u1**2)/Den + t1*u1*(6d0*mtsq*s-t1*u1)
     & /Den**2)/s**2
      n(4,2)=czip    
      n(5,2)=4d0*n(2,2)
      n(6,2)=-4d0*eplog/u1 - 2d0*(-8d0*mtsq + 3d0*t1*u1/s)/Den
      n(7,2)=-2d0*(t1*(2d0*mtsq*t1/u1-4d0*mtsq+3*t1)
     & - 6d0*mtsq*u1*(mtsq-t1)/uu)/s/Den
      n(8,2)=-2d0*n(2,2) 
      n(9,2)=n(9,1)

      c(1)=czip
      c(2)=czip
      c(3)=(u1-t1)/Den 
      c(4)=czip
      c(5)=czip 
      c(6)=1d0/tt
      c(7)=-1d0/uu
      c(8)=czip
      c(9)=czip

      d(1,1)=czip 
      d(2,1)=czip
      d(3,1)=2d0*tt*u1/Den**2
      d(4,1)=czip 
      d(5,1)=czip 
      d(6,1)=4d0/Den
      d(7,1)=-4d0*mtsq/uu/Den 
      d(8,1)=czip 
      d(9,1)=czip

      d(1,2)=czip 
      d(2,2)=czip
      d(3,2)=-2d0*t1*uu/Den**2
      d(4,2)=czip
      d(5,2)=czip
      d(6,2)=4d0*mtsq/tt/Den  
      d(7,2)=-4d0/Den
      d(8,2)=czip
      d(9,2)=czip

      e(1,1)=d(1,2)
      e(2,1)=d(2,2)
      e(3,1)=d(3,2) + 2d0*mtsq*u1/Den**2
      e(4,1)=d(4,2)
      e(5,1)=d(5,2) 
      e(6,1)=2d0*d(6,2)
      e(7,1)=2d0*d(7,2)
      e(8,1)=d(8,2)
      e(9,1)=d(9,2)

      e(1,2)=d(1,1)
      e(2,2)=d(2,1)
      e(3,2)=d(3,1) - 2d0*mtsq*t1/Den**2
      e(4,2)=d(4,1)
      e(5,2)=d(5,1)
      e(6,2)=2d0*d(6,1)
      e(7,2)=2d0*d(7,1)
      e(8,2)=d(8,1)
      e(9,2)=d(9,1)


      g(1,1)=6d0/s
      g(2,1)=g(1,1)
      g(3,1)=(6d0 - 2d0*t1**2/Den - t1**2*u1**2/Den**2)/s
      g(4,1)=czip
      g(5,1)=4d0*g(1,1)
      g(6,1)=-2d0*mtsq*t1/tt/Den
      g(7,1)=-2d0*(t1 + 2d0*mtsq*u1/uu)/Den
      g(8,1)=-2d0*g(1,1)
      g(9,1)=czip

      g(1,2)=-12d0/s**2
      g(2,2)=g(1,2)
      g(3,2)=-2d0*(6d0-u1**2*(4d0*mtsq*s+5d0*t1**2)/Den**2 
     & - 2d0*t1**3*u1**3/Den**3)/s**2
      g(4,2)=czip  
      g(5,2)=4d0*g(1,2)
      g(6,2)=4d0*(5d0*mtsq**2*t1/tt + 3d0*u1*(2d0*mtsq-t1*u1/s))/Den**2
      g(7,2)=-4d0*(mtsq*(2d0*t1+3d0*u1)-t1*u1
     & *(t1+4d0*u1)/s + mtsq*u1**2/uu)/Den**2
      g(8,2)=-2d0*g(1,2)
      g(9,2)=4d0/Den

      g(1,3)=g(1,2)
      g(2,3)=g(1,2)
      g(3,3)=g(3,2) + 4d0*mtsq**2*s*u1/Den**3
      g(4,3)=czip      
      g(5,3)=g(5,2)
      g(6,3)=4d0*(u1*(4d0*Den+t1*u1)/s + 3d0*mtsq*t1*uu/tt 
     & + mtsq**2*t1**2/tt**2)/Den**2
      g(7,3)=4d0*u1*(2d0*t1*u1 - Den)/s/Den**2  
      g(8,3)=-2d0*g(1,2)  
      g(9,3)=4d0*mtsq/tt/Den

      g(1,4)=g(1,2)
      g(2,4)=g(1,2)
      g(3,4)=g(3,3) + 4d0*t1**2*u1*zu/Den**3
      g(4,4)=czip
      g(5,4)=g(5,2)
      g(6,4)=4d0*(3d0*t1+4d0*u1 + 2d0*t1**2*u1/Den)/s/Den
      g(7,4)=4d0*u1*(-3d0*t1**2/s+3d0*mtsq*t1/uu-mtsq**3/uu**2)/Den**2
      g(8,4)=-2d0*g(1,2)  
      g(9,4)=4d0*mtsq/uu/Den

      g(1,5)=g(1,2)
      g(2,5)=g(1,2)
      g(3,5)=g(3,4) - 4d0*mtsq**2*s*t1/Den**3
      g(4,5)=czip
      g(5,5)=g(5,2)
      g(6,5)=4d0*(mtsq**2*t1/tt + u1*(2d0*mtsq-3d0*t1*u1/s))/Den**2
      g(7,5)=4d0*(mtsq**2*u1/uu + t1*(2d0*mtsq-3d0*t1*u1/s))/Den**2
      g(8,5)=-2d0*g(1,2)
      g(9,5)=g(9,2)

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
      
      do i=1,9
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
