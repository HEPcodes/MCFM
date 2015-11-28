      subroutine KMgg3h(mt,zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm,
     & KMgg3hab,KMgg3hba,KMgg3hdd)
      implicit none
      include 'constants.f'
      include 'eplog.f'
      include 'scheme.f'
c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%
C-----This is an implementation of KM Eq. (2.23)

      double complex zp1(4),zp2(4),zp3(4),zp4(4),lnrat,wlog
      double complex e1(4),e2(4),Ubm(4),Vm(4)
      double complex cdot
      double complex KMgg3hab,KMgg3hba,KMgg3hdd,Bs
      double precision mt,nlf,s,mtsq,beta,x
      parameter(nlf=5d0)
      s=2d0*dble(cdot(zp1,zp2))
      mtsq=mt**2
      beta=sqrt(1d0-4d0*mtsq/s)
      x=(1d0-beta)/(1d0+beta)
      wlog=log(x)+im*pi

c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)
      KMgg3hab=
     & ( xn*(5d0/3d0*(eplog-lnrat(-s,mtsq))+31d0/9d0)
     & -nlf*(2d0/3d0*(eplog-lnrat(-s,mtsq))+10d0/9d0))
     & -2d0/3d0*(eplog+5d0/3d0+4d0*mtsq/s
     & +(1d0+2d0*mtsq/s)*beta*wlog)
      if (scheme .eq. 'dred') then
        KMgg3hab=KMgg3hab-Nc/3d0
      endif
      
      KMgg3hab=KMgg3hab*Bs(zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm)
      KMgg3hba=-KMgg3hab
      KMgg3hdd=czip
      
      return
      end
