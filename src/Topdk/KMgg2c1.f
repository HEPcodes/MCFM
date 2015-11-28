      subroutine KMgg2c1(mt,zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm,
     & KMgg2c1ba,KMgg2c1dd)
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
C-----This is an implementation of KM Eq. (2.33)

      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex e1(4),e2(4),Ubm(4),Vm(4),fdl
      double complex cdot,string2,string1
      double complex KMgg2c1ba,KMgg2c1dd,Bt
      double precision mt,t1,tt,mtsq,zeta2,ddilog
      mtsq=mt**2
      t1=2d0*dble(cdot(zp1,zp3))
      tt=t1+mtsq
      zeta2=pisqo6
      
      fdl=zeta2-ddilog(tt/mtsq)

c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)
      KMgg2c1ba=
     & (-eplog+2d0*fdl*mtsq/t1-log(-t1/mtsq)*(6d0*mtsq+t1)/tt)/6d0
      if (scheme .eq. 'dred') then
        KMgg2c1ba=KMgg2c1ba+(Cf-Nc/2d0)
      endif      
      KMgg2c1ba=KMgg2c1ba*Bt(zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm)
      
      KMgg2c1ba=KMgg2c1ba+im*(
     & +cdot(e1,zp3)*string1(Ubm,e2,Vm)*(
     &  2d0*fdl*mtsq*tt/t1**2+log(-t1/mtsq)*(mtsq/tt-2d0*tt/t1)+1d0)
     & +mt*cdot(zp3,e1)*string2(Ubm,zp1,e2,Vm)*(
     &  log(-t1/mtsq)/tt-1d0/t1)
     & +mt*string2(Ubm,e1,e2,Vm)*log(-t1/mtsq))/3d0/tt
      
      KMgg2c1dd=czip
      
      return
      end
