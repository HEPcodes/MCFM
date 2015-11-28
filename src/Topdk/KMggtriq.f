      subroutine KMggtriq(mt,zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm,
     & KMggtriqab,KMggtriqba,KMggtriqdd)
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
C-----This is an implementation of KM Eq. (2.50)

      double complex zp1(4),zp2(4),zp3(4),zp4(4),slog
      double complex e1(4),e2(4),Ubm(4),Vm(4)
      double complex cdot,string1
      double complex KMggtriqab,KMggtriqba,KMggtriqdd,Bs
      double precision mt,nlf,s,mtsq
      parameter(nlf=5d0)
      mtsq=mt**2
      s=2d0*dble(cdot(zp1,zp2))
      slog=log(s/mtsq)-im*pi
                  
c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)
      KMggtriqab=2d0*nlf*(
     & 3d0*eplog-3d0*slog+5d0)/9d0

      KMggtriqab=KMggtriqab*Bs(zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm)
      
      KMggtriqab=KMggtriqab+2d0*nlf*3d0*im*string1(Ubm,zp1,Vm)*(
     & +cdot(e1,e2)/s-2d0*cdot(e1,zp2)*cdot(e2,zp1)/s**2)/9d0
      
      KMggtriqba=-KMggtriqab
      KMggtriqdd=czip
      
      return
      end
