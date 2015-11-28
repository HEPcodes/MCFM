      subroutine KMggtriHQ(mt,zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm,
     & KMggtriHQab,KMggtriHQba,KMggtriHQdd)
      implicit none
      include 'constants.f'
      include 'eplog.f'
c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%
C-----This is an implementation of KM Eq. (2.51)

      double complex zp1(4),zp2(4),zp3(4),zp4(4),wlog,wlogsq
      double complex e1(4),e2(4),Ubm(4),Vm(4)
      double complex cdot,string1
      double complex KMggtriHQab,KMggtriHQba,KMggtriHQdd,Bs
      double precision mt,s,mtsq,beta,x,zeta2
      mtsq=mt**2
      s=2d0*dble(cdot(zp1,zp2))
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
      KMggtriHQab=2d0*(
     & 3d0*eplog+3d0*wlog*(2d0*mtsq/s+1d0)*beta+5d0+12d0*mtsq/s)/9d0

      KMggtriHQab=KMggtriHQab*Bs(zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm)
      
      KMggtriHQab=KMggtriHQab+2d0*3d0*im*string1(Ubm,zp1,Vm)*(
     & +cdot(e1,e2)/s**2-2d0*cdot(e1,zp2)*cdot(e2,zp1)/s**3)*(
     &  3d0*(wlogsq+4d0*beta*wlog)*mtsq-18d0*zeta2*mtsq+24d0*mtsq+s)/9d0
      
      KMggtriHQba=-KMggtriHQab
      KMggtriHQdd=czip
      
      return
      end
