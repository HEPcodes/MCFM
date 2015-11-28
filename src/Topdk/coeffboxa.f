      subroutine coeffboxa(m,zp1,zp2,zp3,zp4,sum)
c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%

C-----Coefficient for Box a
      implicit none
      include 'constants.f'
      include 'eplog.f'
      double precision ha0(6),ha(4,6),zeta2
      double complex zp1(4),zp2(4),zp3(4),zp4(4),cdot,fll,flld,wlogsq
      double complex fa(6),sum(0:6)
      double precision m,mtsq,s,t1,u1,tt,Den,x,z1,z2,zt,beta,ddilog
      integer i,k

c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)
      mtsq=m**2
      s=2d0*dble(cdot(zp1,zp2))
      t1=+2d0*dble(cdot(zp1,zp3))
      tt=t1+mtsq
      u1=-s-t1
      Den=mtsq*s-u1*t1
      
      z1=mtsq*s-t1**2
      z2=s+2d0*t1
      zt=2d0*mtsq+t1
c      zu=2*mtsq+u1

      beta=sqrt(1d0-4d0*mtsq/s)
      x=(1d0-beta)/(1d0+beta)      
      wlogsq=log(x)**2+2d0*im*pi*log(x)
      zeta2=pi**2/6d0

c      write(6,*) 't1,tt',t1,tt

c--- Eq. (2.28)
      flld=4d0*ddilog(-x)+wlogsq+2d0*zeta2
c--- Eq. (2.37)
      fll=log(-t1/mtsq)**2+ddilog(tt/mtsq) 

C-----This is an implementation of KM, Eq. (4.14)
      fa(1)=(log(s/mtsq)**2-2d0*im*pi*log(s/mtsq))
     &     +4d0*fll-4d0*log(-t1/mtsq)*(log(s/mtsq)-im*pi)
     &     +2d0*zeta2
      fa(2)=log(-t1/mtsq) 
      fa(3)=log(s/mtsq)-im*pi 
      fa(4)=flld
      fa(5)=fll
      fa(6)=1d0

C-----This is an implementation of KM, Eq. (B1)
      ha0(1)=(t1**3/Den+mtsq+2d0*t1**2/s)/2d0/Den
      ha0(2)=2d0*(2d0*eplog/s+2d0/s-t1/Den)
      ha0(3)=-2d0*(1d0+t1*zt/beta**2/Den)/s
      ha0(4)=-((mtsq*s*tt+t1**4/s)/Den-t1*(2d0*tt+zt/beta**2)/s)
     & /2d0/beta/Den
      ha0(5)=-4d0/s
      ha0(6)=-2d0*epsqlog/s

      ha(1,1)=2d0*t1*tt/Den**2
      ha(2,1)=8d0*t1/s/Den
      ha(3,1)=-8d0*(1d0-t1*zt/Den)/s**2/beta**2
      ha(4,1)=2d0*zt*(tt/Den+2d0*mtsq/s**2/beta**2)/beta/Den

      ha(1,2)=-z1/4d0/Den**2
      ha(2,2)=-zt/tt/Den
      ha(3,2)=1d0/Den
      ha(4,2)=(1d0-s*t1*beta**2/Den)/4d0/beta/Den

      ha(1,3)=t1/8d0/Den
      ha(2,3)=0d0
      ha(3,3)=0d0
      ha(4,3)=zt/beta/8d0/Den

      ha(1,4)=t1**2/Den**2
      ha(2,4)=-4d0/Den
      ha(3,4)=-4d0*zt/s/beta**2/Den
      ha(4,4)=(t1*z1/Den+zt/beta**2)/s/beta/Den

      ha(1,5)=s*t1/2d0/Den**2
      ha(2,5)=2d0*t1/tt/Den
      ha(3,5)=2d0*z2/s/beta**2/Den
      ha(4,5)=-(2d0*zt/s/beta**2+t1*z2/Den)/2d0/beta/Den

      ha(1,6)=s*t1/4d0/Den**2
      ha(2,6)=t1/tt/Den
      ha(3,6)=z2/s/beta**2/Den
      ha(4,6)=-(2d0*zt/s/beta**2+t1*z2/Den)/4d0/beta/Den

      do k=0,6
      sum(k)=czip
      enddo

      do i=1,6
      sum(0)=sum(0)+fa(i)*ha0(i)
      enddo

      do k=1,6
      do i=1,4
      sum(k)=sum(k)+fa(i)*ha(i,k)
      enddo
      enddo


      return
      end

