      double complex function KMqqb5g(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
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
C-----This is an implementation of KM, Eq. (4.2) (related to Eq. (2.23))

      double complex zp1(4),zp2(4),zp3(4),zp4(4),lnrat,wlog
      double complex Ul(4),Vbl(4),Ubm(4),Vm(4),Bqq,cdot
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
      KMqqb5g=
     & ( xn*(5d0/3d0*(eplog-lnrat(-s,mtsq))+31d0/9d0)
     & -nlf*(2d0/3d0*(eplog-lnrat(-s,mtsq))+10d0/9d0))
     & -2d0/3d0*(eplog+5d0/3d0+4d0*mtsq/s
     & +(1d0+2d0*mtsq/s)*beta*wlog)
      if (scheme .eq. 'dred') then
        KMqqb5g=KMqqb5g-Nc/3d0
      endif
      KMqqb5g=KMqqb5g*Bqq(zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)

      return
      end
