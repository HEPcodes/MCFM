      subroutine coeffboxb(m,zp1,zp2,zp3,zp4,sum)
c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%

C-----Coefficient for Box b
      implicit none
      include 'constants.f'
      include 'eplog.f'
      double precision hb0(6),hb(4,6),zeta2
      double complex zp1(4),zp2(4),zp3(4),zp4(4),cdot,wlogsq
      double complex fb(6),sum(0:6)
      double precision m,mtsq,s,t1,u1,uu,Denu,x,zt,zu,z1u,z2u,
     & beta,ddilog
      integer i,k

c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)
      mtsq=m**2
      s=2d0*dble(cdot(zp1,zp2))
      u1=2d0*dble(cdot(zp2,zp3))
      uu=u1+mtsq;
      t1=-s-u1 
      Denu=s*u1+u1**2+mtsq*s

c      z1=mtsq*s-t1**2 
c      z2=s+2*t1 
      z1u=mtsq*s-u1**2
      z2u=s+2d0*u1
      zt=2d0*mtsq+t1 
      zu=2d0*mtsq+u1

      beta=sqrt(1d0-4d0*mtsq/s)
      x=(1d0-beta)/(1d0+beta)
      wlogsq=log(x)**2+2d0*im*pi*log(x)
      zeta2=pi**2/6d0

C-----This is an implementation of KM, Eq. (4.16)
      fb(6)=1d0
      fb(5)=log(-u1/mtsq)**2+ddilog(uu/mtsq)
      fb(4)=4d0*ddilog(-x)+wlogsq+2d0*zeta2
      fb(3)=log(s/mtsq)-im*pi 
      fb(2)=log(-u1/mtsq) 
      fb(1)=2d0*zeta2+4d0*fb(5)
     &     -4d0*log(-u1/mtsq)*(log(s/mtsq)-im*pi)
     &     +(log(s/mtsq)**2-2d0*im*pi*log(s/mtsq)) 

      hb0(1)=(mtsq*s*t1/Denu-mtsq-4d0*u1-2d0*u1**2/s)/2d0/Denu
      hb0(2)=-2d0*(2d0*(eplog+1d0)*uu/s+1d0-mtsq*t1/Denu)/uu
      hb0(3)=2d0*(1d0+zt*t1/beta**2/Denu)/s
      hb0(4)=(s*beta**2-2d0*zu-2d0*Denu/s/beta**2
     & +4d0*mtsq*u1*zu/s**2/beta**2
     & +(mtsq*z1u+s*u1*t1*beta**2)/Denu)/2d0/beta/Denu
      hb0(5)=4d0/s
      hb0(6)=2d0*epsqlog/s

      hb(1,1)=(mtsq*s+u1*zu)/Denu**2
      hb(2,1)=4d0*(zu/uu+2d0*u1/s)/Denu
      hb(3,1)=-2d0*(1d0/beta**2+1d0)*z2u/s/Denu
      hb(4,1)=((-mtsq*s*beta**2+u1*z1u/s)/Denu
     & +(8d0*mtsq**2/s+u1)/s/beta**2)/beta/Denu

      hb(1,2)=z1u/4d0/Denu**2
      hb(2,2)=zu/uu/Denu
      hb(3,2)=-1d0/Denu
      hb(4,2)=(s*u1*beta**2/Denu-1d0)/4d0/beta/Denu

      hb(1,3)=u1/8d0/Denu
      hb(2,3)=0d0
      hb(3,3)=0d0
      hb(4,3)=zu/beta/8d0/Denu

      hb(1,4)=-u1**2/Denu**2
      hb(2,4)=4d0/denu
      hb(3,4)=4d0*zu/s/beta**2/Denu
      hb(4,4)=-(zu/beta**2+u1*z1u/Denu)/s/beta/Denu

      hb(1,5)=s*u1/2d0/Denu**2
      hb(2,5)=2d0*u1/uu/Denu
      hb(3,5)=2d0*z2u/s/beta**2/Denu
      hb(4,5)=-(z2u*u1/Denu+2d0*zu/s/beta**2)/2d0/beta/Denu

      hb(1,6)=s*u1/4d0/Denu**2
      hb(2,6)=u1/uu/Denu
      hb(3,6)=z2u/s/beta**2/Denu
      hb(4,6)=-(z2u*u1/Denu+2d0*zu/s/beta**2)/4d0/beta/Denu

      do k=0,6
      sum(k)=czip
      enddo

      do i=1,6
      sum(0)=sum(0)+fb(i)*hb0(i)
      enddo

      do k=1,6
      do i=1,4
      sum(k)=sum(k)+fb(i)*hb(i,k)
      enddo
      enddo

      return
      end
