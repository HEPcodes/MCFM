      double complex function KMqqb5h(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
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
C-----This is an implementation of KM, Eq. (4.4)

      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex Ul(4),Vbl(4),Ubm(4),Vm(4),Bqq
      double precision cf,mt
      parameter(cf=4d0/3d0)
      
      KMqqb5h=-cf*(3d0*eplog+4d0)
      if (scheme .eq. 'dred') then
        KMqqb5h=KMqqb5h-cf
      endif
      KMqqb5h=KMqqb5h*Bqq(zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)

      return
      end
