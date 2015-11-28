      subroutine KMampsqggQQb(q1,q2,q3,q4,xmt,res0,res1)
      implicit none
c----notation
C     0--> g(q1)+g(q2)+Q(p3)+Qb(p4)
c---
c--- Computed using the helicity amplitudes from:
c---
c--- \bibitem{Korner:2002hy}
c--- J.~G.~Korner and Z.~Merebashvili,
c--- %``One-loop corrections to four-point functions with two external massive
c--- %fermions and two external massless partons,''
c--- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c--- [arXiv:hep-ph/0207054].
c--- %%CITATION = PHRVA,D66,054023;%%
      include 'constants.f'
      include 'scale.f'
      include 'epinv.f'
      include 'spinorsw.f'
      double precision q1(4),q2(4),q3(4),q4(4),xmt,res0,res1,nlf
      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex ze1(4),ze2(4)
      double complex Ubm(4),Vm(4)
      integer j,nu,j1,j2,j3,j4
      double complex renab,renba
      double complex Bs,Bt,Bu
      double complex d2bab,d2bba,D2bdd,D2bddu,D2c1ab,D2c1ba,D2c1dd,
     & D2c1ddu,D2c2ab,D2c2ba,D2c2dd,D2c2ddu,D2c3ab,D2c3ba,D2c3dd,
     & D2c3ddu,D2c4ab,D2c4ba,D2c4dd,D2c4ddu,D2d1ab,D2d1ba,D2d1dd,
     & D2d1ddu,D3f1ab,D3f1ba,D3f1dd,D3f2ab,D3f2ba,D3f2dd,D3hab,D3hba,
     & D3hdd,Da1ab,Da1ba,Da1dd,Da1ddu,Da2ab,Da2ba,Da2dd,Da2ddu,
     & Da34ab,Da34ba,Da34dd,Da34ddu,Dtrigab,Dtrigba,Dtrigdd,Dtrihqab,
     & Dtrihqba,Dtrihqdd,Dtriqab,Dtriqba,Dtriqdd
      double complex tab1,tab2,tab3,loab,loba,loqed,
     & KMggab,KMggba,KMggdd
      parameter(nlf=5d0)
      
c--- set up common variable to indicate KM routines should be used
c--- to calculate spinor strings
      spinorsw='KM'      

c--- these routines expect energy in the 1st component and, as usual,
c--- we modifiy this solution by exchanging (px <--> pz) and changing 
C--- the sign of py, (to remain in right handed system) [c.f. spinor_w.f]
      do j=1,4
      if (j .eq. 1) nu=4
      if (j .eq. 2) nu=3
      if (j .eq. 3) nu=2
      if (j .eq. 4) nu=1
      zp1(j)=dcmplx(q1(nu))
      zp2(j)=dcmplx(q2(nu))
      zp3(j)=dcmplx(q3(nu))
      zp4(j)=dcmplx(q4(nu))
      if (j .eq. 3) then
      zp1(j)=-zp1(j)
      zp2(j)=-zp2(j)
      zp3(j)=-zp3(j)
      zp4(j)=-zp4(j)
      endif
      enddo

c--- initialize to zero
      res0=0d0
      res1=0d0
      
      do j1=-1,1,2
      call vspinor(zp4,xmt,j1,Vm)

      do j2=-1,1,2
      call ubarspinor(zp3,xmt,j2,Ubm)

      do j3=-1,1,2
      call pol_mless(zp1,j3,ze1)

      do j4=-1,1,2
      call pol_mless(zp2,j4,ze2)

      tab1=Bs(zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm)
      tab3=Bt(zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm)
      tab2=Bu(zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm)

c--- t-channel diagrams
      call KMgg2b(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D2bba,D2bdd)
      call KMgg2c1(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D2c1ba,D2c1dd)
      call KMgg2c2(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D2c2ba,D2c2dd)
      call KMgg2c3(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D2c3ba,D2c3dd)
      call KMgg2c4(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D2c4ba,D2c4dd)
      call KMgg2d1(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D2d1ba,D2d1dd)

c-- u-channel diagrams obtained by symmetry, KM Eq. (2.5)
      call KMgg2b(xmt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & D2bab,D2bddu)
      call KMgg2c1(xmt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & D2c1ab,D2c1ddu)
      call KMgg2c2(xmt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & D2c2ab,D2c2ddu)
      call KMgg2c3(xmt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & D2c3ab,D2c3ddu)
      call KMgg2c4(xmt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & D2c4ab,D2c4ddu)
      call KMgg2d1(xmt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & D2d1ab,D2d1ddu)

c--- s-channel diagrams
      call KMgg3h(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D3hab,D3hba,D3hdd)
      call KMgg3f1(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D3f1ab,D3f1ba,D3f1dd)
      call KMgg3f2(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D3f2ab,D3f2ba,D3f2dd)
      call KMggtrig(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & Dtrigab,Dtrigba,Dtrigdd)
      call KMggtriq(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & Dtriqab,Dtriqba,Dtriqdd)
      call KMggtriHQ(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & DtriHQab,DtriHQba,DtriHQdd)

c--- boxes
      call KMgg2a(xmt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & Da1ba,Da1dd,Da2ba,Da2dd,Da34ba,Da34dd)

c-- Mt <-> Mu box symmetry, KM Eq. (2.5)
      call KMgg2a(xmt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & Da1ab,Da1ddu,Da2ab,Da2ddu,Da34ab,Da34ddu)

c--- put in overall wave function and charge renormalization by hand
      renab=tab2+tab1
      renba=tab3-tab1
      
      renab=renab*(
     & +2d0*((nlf/3d0-11d0/6d0*xn)*epinv+xn/6d0)
     & -(xn**2-1d0)/2d0/xn*(3d0*(epinv+log(musq/xmt**2))+5d0))
      renba=renba*(
     & +2d0*((nlf/3d0-11d0/6d0*xn)*epinv+xn/6d0)
     & -(xn**2-1d0)/2d0/xn*(3d0*(epinv+log(musq/xmt**2))+5d0))

c--- sum all diagrams (ab)
      KMggab=
     & +D3hab+D2bab+D2c1ab+D2c2ab+D2c3ab+D2c4ab+D3f1ab+D3f2ab
     & +Dtrigab+Dtriqab+DtriHQab+Da1ab+Da2ab+Da34ab
     & +D2d1ab+renab

c--- sum all diagrams (ba)
      KMggba=
     & +D3hba+D2bba+D2c1ba+D2c2ba+D2c3ba+D2c4ba+D3f1ba+D3f2ba
     & +Dtrigba+Dtriqba+DtriHQba+Da1ba+Da2ba+Da34ba
     & +D2d1ba+renba

c--- sum all diagrams (dd)
c--- note that D2bddu is excluded, since it is equal to D2bdd
      KMggdd=
     & +D3hdd+D2bdd+D2c1dd+D2c2dd+D2c3dd+D2c4dd+D3f1dd+D3f2dd
     & +Dtrigdd+Dtriqdd+DtriHQdd+Da1dd+Da2dd+Da34dd
     & +D2c1ddu+D2c2ddu+D2c3ddu+D2c4ddu+D2d1ddu
     & +Da1ddu+Da2ddu+Da34ddu
     & +D2d1dd
c--- apply additional factor of (2*Nc) that was present in the numerical
c--- routines, in order to use the same code for summing and squaring
      KMggdd=KMggdd*2d0*xn

      loab=tab2+tab1
      loba=tab3-tab1
      loqed=tab2+tab3

      res0=res0
     . +V*xn/4d0*(abs(loab)**2+abs(loba)**2-abs(loqed)**2/xn**2)

      res1=res1+V*xn/4d0*dreal(
     . +KMggab*dconjg(loab)+KMggba*dconjg(loba)
     . +(KMggdd-KMggab-KMggba)*dconjg(loqed)/xn**2)

      enddo
      enddo
      enddo
      enddo

      return
      end
