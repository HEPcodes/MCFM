      subroutine KMgg2d1(mt,zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm,
     & KMgg2d1ba,KMgg2d1dd)
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
C-----This is an implementation of KM Eq. (2.16)

      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex e1(4),e2(4),Ubm(4),Vm(4)
      double complex cdot,string2
      double complex KMgg2d1ba,KMgg2d1dd,Bt
      double precision mt,t1,tt,mtsq
      mtsq=mt**2
      t1=2d0*dble(cdot(zp1,zp3))
      tt=t1+mtsq

c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)
      KMgg2d1ba=
     & Cf*(-eplog-t1/tt+log(-t1/mtsq)*(4d0*t1/tt+t1**2/tt**2-4d0))
      if (scheme .eq. 'dred') then
        KMgg2d1ba=KMgg2d1ba-Cf
      endif      
      KMgg2d1ba=KMgg2d1ba*Bt(zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm)
      
      KMgg2d1ba=KMgg2d1ba-im*Cf*mt*string2(Ubm,e1,e2,Vm)*(
     & 1d0-2d0*log(-t1/mtsq)-log(-t1/mtsq)*t1/tt)/tt
      
      KMgg2d1dd=czip
      
      return
      end
