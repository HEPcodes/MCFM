      subroutine KMgg3j1(mt,zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm,
     & KMgg3j1ab,KMgg3j1ba,KMgg3j1dd)
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
C-----This is an implementation of KM, Eq. (2.22)
      double precision mt
      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex e1(4),e2(4),Ubm(4),Vm(4),Bs
      double complex KMgg3j1ab,KMgg3j1ba,KMgg3j1dd

      KMgg3j1ab=-2d0/3d0*eplog*Bs(zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm)
      
      KMgg3j1ba=-KMgg3j1ab
      KMgg3j1dd=czip

      return
      end
