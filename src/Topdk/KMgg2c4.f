      subroutine KMgg2c4(mt,zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm,
     & KMgg2c4ba,KMgg2c4dd)
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
C-----This is an implementation of KM Eq. (2.42)

      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex e1(4),e2(4),Ubm(4),Vm(4),fll
      double complex cdot,string2,string1,string0
      double complex KMgg2c4ba,KMgg2c4dd,Bt
      double precision mt,t1,tt,mtsq,ddilog
      mtsq=mt**2
      t1=2d0*dble(cdot(zp1,zp3))
      tt=t1+mtsq
      
c--- Eq. (2.37)
      fll=log(-t1/mtsq)**2+ddilog(tt/mtsq) 

c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)
      KMgg2c4ba=3d0*(
     &  3d0*eplog-epsqlog
     &  +(2d0*log(-t1/mtsq)-1d0)*eplog+4d0
     &  +6d0*log(-t1/mtsq)*mtsq/tt-2d0*fll)/2d0
      if (scheme .eq. 'dred') then
        KMgg2c4ba=KMgg2c4ba+Nc/2d0
      endif      
      KMgg2c4ba=KMgg2c4ba*Bt(zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm)
      
      KMgg2c4ba=KMgg2c4ba+3d0*im*(
     &  cdot(e2,zp4)*string1(Ubm,e1,Vm)*(
     &   +epsqlog-2d0*eplog-2d0*log(-t1/mtsq)*eplog
     &   +2d0*fll+2d0*log(-t1/mtsq)
     &    *mtsq*(mtsq/tt+1d0)/tt+2d0*mtsq/tt-6d0)/2d0/t1
     & -3d0*mt*string2(Ubm,e1,e2,Vm)*log(-t1/mtsq)/2d0/tt
     & -mt*cdot(e2,zp4)*string2(Ubm,e1,zp1,Vm)*(
     &   log(-t1/mtsq)*(mtsq+tt)+tt)/t1/tt**2
     & -2d0*mt*cdot(e1,zp3)*cdot(e2,zp4)*string0(Ubm,Vm)*(
     &   log(-t1/mtsq)*(mtsq+tt)+tt)/t1/tt**2)
      
      KMgg2c4dd=czip
      
      return
      end
