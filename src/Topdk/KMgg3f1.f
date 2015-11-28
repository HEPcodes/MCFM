      subroutine KMgg3f1(mt,zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm,
     & KMgg3f1ab,KMgg3f1ba,KMgg3f1dd)
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
C-----This is an implementation of KM Eq. (2.46)

      double complex zp1(4),zp2(4),zp3(4),zp4(4),wlog,wlogsq
      double complex e1(4),e2(4),Ubm(4),Vm(4)
      double complex cdot,string0
      double complex KMgg3f1ab,KMgg3f1ba,KMgg3f1dd,Bs
      double precision mt,s,t1,mtsq,beta,x,zeta2,ddilog
      mtsq=mt**2
      s=2d0*dble(cdot(zp1,zp2))
      t1=2d0*dble(cdot(zp1,zp3))
      zeta2=pisqo6
      beta=sqrt(1d0-4d0*mtsq/s)
      x=(1d0-beta)/(1d0+beta)
      wlog=log(x)+im*pi
      wlogsq=log(x)**2+2d0*im*pi*log(x)
            
c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)
      KMgg3f1ab=-s*(
     & +s*beta*eplog+(2d0*mtsq-s)*(2d0*wlog*eplog-4d0*ddilog(x)
     & -4d0*wlog*log(1d0-x)+6d0*wlog+wlogsq-8d0*zeta2)+3d0*wlog*s
     &             )/6d0/s**2/beta
      if (scheme .eq. 'dred') then
        KMgg3f1ab=KMgg3f1ab+(Cf-Nc/2d0)
      endif      
      KMgg3f1ab=KMgg3f1ab*Bs(zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm)
      
      KMgg3f1ab=KMgg3f1ab+2d0*im*(
     & mt*wlog*string0(Ubm,Vm)*(-cdot(e1,e2)*(s+2d0*t1)
     &                          -4d0*cdot(e1,zp3)*cdot(e2,zp4)
     &                          +4d0*cdot(e1,zp4)*cdot(e2,zp3))
     &                           )/6d0/s**2/beta
      
      KMgg3f1ba=-KMgg3f1ab
      KMgg3f1dd=czip
      
      return
      end
