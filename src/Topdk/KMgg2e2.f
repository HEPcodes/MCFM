      subroutine KMgg2e2(mt,zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm,
     & KMgg2e2ab,KMgg2e2ba,KMgg2e2dd)
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
C-----This is an implementation of KM, Eq. (2.20)
      double precision mt
      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex e1(4),e2(4),Ubm(4),Vm(4),Bt
      double complex KMgg2e2ab,KMgg2e2ba,KMgg2e2dd

      KMgg2e2ba=-2d0/3d0*eplog*Bt(zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm)
      
      KMgg2e2ab=czip
      KMgg2e2dd=czip

      return
      end
