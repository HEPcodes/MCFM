      double complex function KMqqb5c(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
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
C-----This is an implementation of KM, Eq. (4.6)

      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex Ul(4),Vbl(4),Ubm(4),Vm(4),Bqq
      double complex cdot,im
      double precision cf,Nc,mt,s,mtsq,pi,zeta2
      parameter(cf=4d0/3d0,Nc=3d0,im=(0d0,1d0))
      s=2d0*dble(cdot(zp1,zp2))
      mtsq=mt**2
      pi=2d0*dasin(1d0)
      zeta2=pi**2/6d0

c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)
      KMqqb5c=
     &  (2d0*epsqlog+3d0*eplog-2d0*(log(s/mtsq)-im*pi)*eplog
     &  +(log(s/mtsq)**2-2d0*im*pi*log(s/mtsq))
     &  -3d0*(log(s/mtsq)-im*pi)-8d0*zeta2+8d0)/6d0
      if (scheme .eq. 'dred') then
        KMqqb5c=KMqqb5c+(cf-Nc/2d0)
      endif
      KMqqb5c=KMqqb5c*Bqq(zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)

      return
      end
