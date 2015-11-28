      double complex function KMqqb5f(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      implicit none
      include 'eplog.f'
      include 'scheme.f'
c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%
C-----This is an implementation of KM, Eq. (4.10)

      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex Ul(4),Vbl(4),Ubm(4),Vm(4),Bqq,wlogsq,flld
      double complex cdot,im,string1,string0
      double precision pi,Nc,mt,mtsq,beta,s,x,ddilog,zeta2
      parameter(Nc=3d0)
      mtsq=mt**2
      pi=2d0*dasin(1d0)
      im=dcmplx(0d0,1d0)
      zeta2=pi**2/6d0
      s=2d0*dble(cdot(zp1,zp2))
      beta=sqrt(1d0-4d0*mt**2/s)     
      x=(1d0-beta)/(1d0+beta)
      wlogsq=log(x)**2+2d0*im*pi*log(x)
C---Eq.2.28
      flld=4d0*ddilog(-x)+wlogsq+2d0*zeta2

c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)
      KMqqb5f=
     & (3d0*(eplog+2d0)-(log(s/mtsq)-im*pi)*(8d0*mtsq/s-1d0)/beta**2
     & +flld*mtsq/s/beta**3)
      if (scheme .eq. 'dred') then
        KMqqb5f=KMqqb5f+Nc/2d0/1.5d0
      endif
      KMqqb5f=KMqqb5f*Bqq(zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      KMqqb5f=KMqqb5f
     & -4d0*im*mt*string1(Vbl,zp3,Ul)*string0(Ubm,Vm)
     & *((log(s/mtsq)-im*pi)*(8d0*mtsq/s+1d0)/beta**2
     &                       -2d0-3d0*flld*mtsq/s/beta**3)/s**2/beta**2
      KMqqb5f=KMqqb5f*1.5d0
      return
      end
