      subroutine KMggtrig(mt,zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm,
     & KMggtrigab,KMggtrigba,KMggtrigdd)
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
C-----This is an implementation of KM Eq. (2.49)

      double complex zp1(4),zp2(4),zp3(4),zp4(4),slog,slogsq
      double complex e1(4),e2(4),Ubm(4),Vm(4)
      double complex cdot,string1
      double complex KMggtrigab,KMggtrigba,KMggtrigdd,Bs
      double precision mt,s,mtsq,zeta2
      mtsq=mt**2
      s=2d0*dble(cdot(zp1,zp2))
      zeta2=pisqo6
      slog=log(s/mtsq)-im*pi
      slogsq=log(s/mtsq)**2-2d0*im*pi*log(s/mtsq)
                  
c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)
      KMggtrigab=Nc*(
     & +33d0*eplog-36d0*epsqlog-171d0*eplog+36d0*slog*eplog
     & +138d0*slog-18d0*slogsq+144d0*zeta2-284d0)/72d0
      if (scheme .eq. 'dred') then
        KMggtrigab=KMggtrigab+1d0
      endif      
      KMggtrigab=KMggtrigab*Bs(zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm)
      
      KMggtrigab=KMggtrigab-Nc*6d0*im*string1(Ubm,zp1,Vm)*(
     & +cdot(e1,e2)*(27d0*eplog-6d0*epsqlog-33d0*eplog
     &              +6d0*slog*eplog+6d0*slog-3d0*slogsq
     &              +24d0*zeta2-4d0)/s
     & -cdot(e1,zp2)*cdot(e2,zp1)*16d0/s**2)/72d0
      
      KMggtrigba=-KMggtrigab
      KMggtrigdd=czip
      
      return
      end
