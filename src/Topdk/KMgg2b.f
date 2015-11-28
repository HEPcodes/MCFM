      subroutine KMgg2b(mt,zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm,
     & KMgg2bba,KMgg2bdd)
      implicit none
      include 'constants.f'
c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%
C-----This is an implementation of KM Eq. (2.29)

      double complex zp1(4),zp2(4),zp3(4),zp4(4),wlogsq,flld
      double complex e1(4),e2(4),Ubm(4),Vm(4)
      double complex cdot,string0,string1,string2
      double complex KMgg2bba,KMgg2bdd
      double precision mt,s,mtsq,beta,x,ddilog,zeta2
      mtsq=mt**2
      zeta2=pi**2/6d0
      s=2d0*dble(cdot(zp1,zp2))
      beta=sqrt(1d0-4d0*mtsq/s)
      x=(1d0-beta)/(1d0+beta)
      wlogsq=log(x)**2+2d0*im*pi*log(x)
C---Eq.2.28
      flld=4d0*ddilog(-x)+wlogsq+2d0*zeta2

c--- Note that analytic continuation of logarithms has been performed
c--- according to Eq. (5.1)
c--- Moreoever, in order to replace overall factor in Eq. (2.3) by
c--- the MCFM (and QCDLoop) standard, 1/Gamma(1-ep), we must also
c--- multiply by (1+ep^2*zeta2)

c--- NB: these diagrams are the same in both schemes
      KMgg2bba=im*xn*(
     & 4d0*(2d0*mt*cdot(e1,e2)*string0(Ubm,Vm)
     &     -(2d0*cdot(zp3,e2)+cdot(zp4,e2))*string1(Ubm,e1,Vm)
     &     +(cdot(zp3,e1)+2d0*cdot(zp4,e1))*string1(Ubm,e2,Vm))
     &    *(s*beta*(log(s/mtsq)-im*pi)-mtsq*flld)
     &-3d0*mt*s*string2(Ubm,e1,e2,Vm)
     &    *(4d0*beta*(log(s/mtsq)-im*pi)-flld)
     &               )/(4d0*s**2*beta**3)

      KMgg2bdd=-im*(
     & 2d0*(2d0*mt*cdot(e1,e2)*string0(Ubm,Vm)
     &     +(cdot(zp3,e2)-cdot(zp4,e2))*string1(Ubm,e1,Vm)
     &     +(cdot(zp3,e1)-cdot(zp4,e1))*string1(Ubm,e2,Vm))
     &    *(s*beta*(log(s/mtsq)-im*pi)-mtsq*flld)
     &-3d0*mt*s*beta**2*cdot(e1,e2)*flld*string0(Ubm,Vm)
     &               )/(4d0*s**2*beta**3)
                  
      return
      end
