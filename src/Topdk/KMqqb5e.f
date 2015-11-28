      double complex function KMqqb5e(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
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
C-----This is an implementation of KM, Eq. (4.9)

      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex Ul(4),Vbl(4),Ubm(4),Vm(4),Bqq,wlog,wlogsq
      double complex cdot,im,string1,string0
      double precision pi,cf,Nc,mt,beta,s,x,ddilog,zeta2
      parameter(cf=4d0/3d0,Nc=3d0)
      pi=2d0*dasin(1d0)
      im=dcmplx(0d0,1d0)
      zeta2=pi**2/6d0
      s=2d0*dble(cdot(zp1,zp2))
      beta=sqrt(1d0-4d0*mt**2/s)     
      x=(1d0-beta)/(1d0+beta)
      wlog=log(x)+im*pi
      wlogsq=log(x)**2+2d0*im*pi*log(x)
      
c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)
      KMqqb5e=
     &  (-eplog+3d0*wlog*beta+(1d0/beta+beta)*(eplog*wlog
     & -2d0*ddilog(x)-2d0*wlog*log(1d0-x)+0.5d0*wlogsq-4d0*zeta2))
      if (scheme .eq. 'dred') then
        KMqqb5e=KMqqb5e+(cf-Nc/2d0)*6d0
      endif
      KMqqb5e=KMqqb5e*Bqq(zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      KMqqb5e=KMqqb5e
     & +4d0*im*mt*string1(Vbl,zp3,Ul)*string0(Ubm,Vm)*wlog/s**2/beta
      KMqqb5e=KMqqb5e/6d0
      return
      end
