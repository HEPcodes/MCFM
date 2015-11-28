      subroutine KMwrap(p)
c--- Routine to take a normal MCFM phase space point and examine
c--- various amplitudes calculated from the paper:   
c--- \bibitem{Korner:2002hy}
c--- J.~G.~Korner and Z.~Merebashvili,
c--- %``One-loop corrections to four-point functions with two external massive
c--- %fermions and two external massless partons,''
c--- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c--- [arXiv:hep-ph/0207054].
c--- %%CITATION = PHRVA,D66,054023;%%
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'scale.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'spinorsw.f'
      integer h1,h2,h3,h4,j1,j2,j3,j4,nlf
      double precision p(mxpart,4),losq
      double precision beta,eta1(4),eta4(4)
      double complex  zp1(4),zp2(4),zp3(4),zp4(4),zeta1(4),zeta4(4),
     & ze1(4),ze2(4)
      double complex  Ubm(4),Vm(4),Vbl(4),Ul(4),string0,string1,temp,
     & temps,tempt,tempu
      double complex Da,Db,Dc,Dd,De,Df,Dg,Dh,Di
      double complex Bqq,KMqqb5c,KMqqb5d,KMqqb5e,KMqqb5f,KMqqb5g,
     & KMqqb5h,KMqqb5i
      double complex ren,renab,renba,dummy
      double complex Bs,Bt,Bu
      double complex D2e1ab,D2e1ba,D2e1dd,D2e2ab,D2e2ba,D2e2dd,
     & D3j1ab,D3j1ba,D3j1dd,D3j2ab,D3j2ba,D3j2dd,
     & D3hab,D3hba,D3hdd,D2bab,D2bba,D2bdd,
     & D2c1ab,D2c1ba,D2c1dd,D2c2ab,D2c2ba,D2c2dd,
     & D2c3ab,D2c3ba,D2c3dd,D2c4ab,D2c4ba,D2c4dd,
     & D2d1ab,D2d1ba,D2d1dd,
     & D3f1ab,D3f1ba,D3f1dd,D3f2ab,D3f2ba,D3f2dd,
     & Dtrigab,Dtrigba,Dtrigdd,Dtriqab,Dtriqba,Dtriqdd,
     & DtriHQab,DtriHQba,DtriHQdd,
     & Da1ba,Da1dd,Da2ba,Da2dd,Da34ba,Da34dd,
     & Da1ab,Da1ddu,Da2ab,Da2ddu,Da34ab,Da34ddu,
     & D2bddu,D2c1ddu,D2c2ddu,D2c3ddu,D2c4ddu,D2d1ddu
      double complex KMqqb(2,2,2),
     & KMggab(2,2,2,2),KMggba(2,2,2,2),KMggdd(2,2,2,2)
c--- definitions for numerical routines below here
      integer j,nu,ep
      double precision q1(4),q2(4),q3(4),q4(4)
      double precision k1(4),k2(4),k3(4),k4(4)
      double complex zk1(4),zk2(4),zk3(4),zk4(4)
      double complex Ubl(4),Vl(4)
      double complex Ublx(4),Vlx(4),Ubmx(4),Vmx(4)
      double complex lo1,
     . abtt1(-2:0),abtt2(-2:0),abtt3(-2:0),abtt4(-2:0),abtt5(-2:0),
     . abtt6(-2:0),abtt7(-2:0),abtt8(-2:0),abtt9(-2:0),abtt10(-2:0),
     . abtt11(-2:0)
      double complex ze3(4),ze4(4),
     . ab1(-2:0),ab2(-2:0),ab3(-2:0),ab4(-2:0),
     . ab5(-2:0),ab6(-2:0),ab7(-2:0),ab8(-2:0),ab9(-2:0),
     . ab10(-2:0),ab11(-2:0),ab12(-2:0),ab13(-2:0),ab14(-2:0),
     . ab15(-2:0),ab16(-2:0),ab17(-2:0),ab18(-2:0),ab19(-2:0),
     . ab20(-2:0),ab21(-2:0),ab22(-2:0),ab23(-2:0),ab24(-2:0),
     . ab25(-2:0),ab26(-2:0),ab27(-2:0),ab28(-2:0),ab29(-2:0),
     . ab30(-2:0),ab31(-2:0),ab32(-2:0),ab33(-2:0),ab34(-2:0),
     . ba1(-2:0),ba2(-2:0),ba3(-2:0),ba4(-2:0),
     . ba5(-2:0),ba6(-2:0),ba7(-2:0),ba8(-2:0),ba9(-2:0),
     . ba10(-2:0),ba11(-2:0),ba12(-2:0),ba13(-2:0),ba14(-2:0),
     . ba15(-2:0),ba16(-2:0),ba17(-2:0),ba18(-2:0),ba19(-2:0),
     . ba20(-2:0),ba21(-2:0),ba22(-2:0),ba23(-2:0),ba24(-2:0),
     . ba25(-2:0),ba26(-2:0),ba27(-2:0),ba28(-2:0),ba29(-2:0),
     . ba30(-2:0),ba31(-2:0),ba32(-2:0),ba33(-2:0),ba34(-2:0),
     . dd1(-2:0),dd2(-2:0),dd3(-2:0),dd4(-2:0),
     . dd5(-2:0),dd6(-2:0),dd7(-2:0),dd8(-2:0),dd9(-2:0),
     . dd10(-2:0),dd11(-2:0),dd12(-2:0),dd13(-2:0),dd14(-2:0),
     . dd15(-2:0),dd16(-2:0),dd17(-2:0),dd18(-2:0),dd19(-2:0),
     . dd20(-2:0),dd21(-2:0),dd22(-2:0),dd23(-2:0),dd24(-2:0),
     . dd25(-2:0),dd26(-2:0),dd27(-2:0),dd28(-2:0),dd29(-2:0),
     . dd30(-2:0),dd31(-2:0),dd32(-2:0),dd33(-2:0),dd34(-2:0),
     . tab1,tab2,tab3,
     . tba1,tba2,tba3,
     . tdd1,tdd2,tdd3
      double complex hoab(-2:0),hoba(-2:0),hodd(-2:0),ho(-2:0),
     . num,numab,numba,numdd

      logical first,showQ,showG
      data first/.true./
      parameter(nlf=5d0)
      
c--- flags to control diagram-by-diagram output
      showQ=.false.
      showG=.false.      
      
      scheme='dred'

C--- create complex momenta
c--- JMC note: swap 1 and 3 components in order to end up with the
c---           same spinors as are generated in spinor_w.f;
c---           change sign of 2nd component???
      zp1(1)=dcmplx(p(1,4))
      zp1(2)=dcmplx(p(1,3))
      zp1(3)=-dcmplx(p(1,2))
      zp1(4)=dcmplx(p(1,1))

      zp2(1)=dcmplx(p(2,4))
      zp2(2)=dcmplx(p(2,3))
      zp2(3)=-dcmplx(p(2,2))
      zp2(4)=dcmplx(p(2,1))

      zp3(1)=dcmplx(p(3,4))
      zp3(2)=dcmplx(p(3,3))
      zp3(3)=-dcmplx(p(3,2))
      zp3(4)=dcmplx(p(3,1))

      zp4(1)=dcmplx(p(4,4))
      zp4(2)=dcmplx(p(4,3))
      zp4(3)=-dcmplx(p(4,2))
      zp4(4)=dcmplx(p(4,1))

c--- set up common variable to indicate KM routines should be used
c--- to calculate spinor strings
      spinorsw='KM'      
      
      losq=0d0
      do h1=-1,1,2
      h2=h1
      do h3=-1,1,2
      do h4=-1,1,2
C---Thus for the moment 
      call vbar0spinor(zp2,h2,Vbl)
      call u0spinor(zp1,h1,Ul)
      call ubarspinor(zp3,mt,h3,Ubm)
      call vspinor(zp4,mt,h4,Vm)


      temp=Bqq(zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      call KMqqb5ab(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm,Da,Db)
      Dc=KMqqb5c(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      Dd=KMqqb5d(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      De=KMqqb5e(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      Df=KMqqb5f(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      Dg=KMqqb5g(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
c      Dh=KMqqb5h(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
c      Di=KMqqb5i(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      
c--- put in overall wave function and charge renormalization by hand
      ren=temp
      
      ren=ren*(
     & +2d0*((nlf/3d0-11d0/6d0*xn)*epinv+xn/6d0)
     & +2d0/3d0*(epinv+log(musq/mt**2))
     & -(xn**2-1d0)/2d0/xn*(3d0*(epinv+log(musq/mt**2))+5d0))
     
c--- get indices into array from helicity labels (-1 > 1, +1 -> 2)
      j1=(h1+3)/2
      j3=(h3+3)/2
      j4=(h4+3)/2
      
c--- sum all diagrams (ab)
      KMqqb(j1,j3,j4)=Da+Db+Dc+Dd+De+Df+Dg+ren
      
      if (showQ) then
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,Bqq',h1,h2,h3,h4,temp
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,Da ',h1,h2,h3,h4,Da/temp
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,Db ',h1,h2,h3,h4,Db/temp
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,Dc ',h1,h2,h3,h4,Dc/temp
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,Dd ',h1,h2,h3,h4,Dd/temp
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,De ',h1,h2,h3,h4,De/temp
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,Df ',h1,h2,h3,h4,Df/temp
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,Dg ',h1,h2,h3,h4,Dg/temp
c      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,Dh ',h1,h2,h3,h4,Dh/temp
c      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,Di ',h1,h2,h3,h4,Di/temp
      write(6,'(a16,4i3,2e20.12)')'h1,h2,h3,h4,ren',h1,h2,h3,h4,ren/temp
      endif

      losq=losq+abs(temp)**2
      
      enddo
      enddo
      enddo

      if (showQ) then
      write(6,*) 'Bqq squared',losq
      write(6,*)
      endif

      do h1=-1,1,2
      do h2=-1,1,2

      do h3=-1,1,2
      do h4=-1,1,2

      call pol_mless(zp1,h1,ze1)
      call pol_mless(zp2,h2,ze2)
      call ubarspinor(zp3,mt,h3,Ubm)
      call vspinor(zp4,mt,h4,Vm)

      temps=Bs(zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm)
      tempt=Bt(zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm)
      tempu=Bu(zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm)

c--- t-channel diagrams
      call KMgg2b(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D2bba,D2bdd)
      call KMgg2c1(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D2c1ba,D2c1dd)
      call KMgg2c2(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D2c2ba,D2c2dd)
      call KMgg2c3(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D2c3ba,D2c3dd)
      call KMgg2c4(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D2c4ba,D2c4dd)
      call KMgg2d1(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D2d1ba,D2d1dd)

c-- u-channel diagrams obtained by symmetry, KM Eq. (2.5)
      call KMgg2b(mt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & D2bab,D2bddu)
      call KMgg2c1(mt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & D2c1ab,D2c1ddu)
      call KMgg2c2(mt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & D2c2ab,D2c2ddu)
      call KMgg2c3(mt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & D2c3ab,D2c3ddu)
      call KMgg2c4(mt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & D2c4ab,D2c4ddu)
      call KMgg2d1(mt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & D2d1ab,D2d1ddu)

c--- s-channel diagrams
      call KMgg3h(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D3hab,D3hba,D3hdd)
      call KMgg3f1(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D3f1ab,D3f1ba,D3f1dd)
      call KMgg3f2(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & D3f2ab,D3f2ba,D3f2dd)
      call KMggtrig(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & Dtrigab,Dtrigba,Dtrigdd)
      call KMggtriq(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & Dtriqab,Dtriqba,Dtriqdd)
      call KMggtriHQ(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & DtriHQab,DtriHQba,DtriHQdd)

c--- boxes
      call KMgg2a(mt,zp1,zp2,zp3,zp4,ze2,ze1,Ubm,Vm,
     & Da1ba,Da1dd,Da2ba,Da2dd,Da34ba,Da34dd)

c-- Mt <-> Mu box symmetry, KM Eq. (2.5)
      call KMgg2a(mt,zp2,zp1,zp3,zp4,ze1,ze2,Ubm,Vm,
     & Da1ab,Da1ddu,Da2ab,Da2ddu,Da34ab,Da34ddu)

c--- put in overall wave function and charge renormalization by hand
      renab=tempu+temps
      renba=tempt-temps
      
      renab=renab*(
     & +2d0*((nlf/3d0-11d0/6d0*xn)*epinv+xn/6d0)
     & -(xn**2-1d0)/2d0/xn*(3d0*(epinv+log(musq/mt**2))+5d0))
      renba=renba*(
     & +2d0*((nlf/3d0-11d0/6d0*xn)*epinv+xn/6d0)
     & -(xn**2-1d0)/2d0/xn*(3d0*(epinv+log(musq/mt**2))+5d0))

c--- get indices into array from helicity labels (-1 > 1, +1 -> 2)
      j1=(h1+3)/2
      j2=(h2+3)/2
      j3=(h3+3)/2
      j4=(h4+3)/2
      
c--- sum all diagrams (ab)
      KMggab(j1,j2,j3,j4)=
     & +D3hab+D2bab+D2c1ab+D2c2ab+D2c3ab+D2c4ab+D3f1ab+D3f2ab
     & +Dtrigab+Dtriqab+DtriHQab+Da1ab+Da2ab+Da34ab
     & +D2d1ab+renab

c--- sum all diagrams (ba)
      KMggba(j1,j2,j3,j4)=
     & +D3hba+D2bba+D2c1ba+D2c2ba+D2c3ba+D2c4ba+D3f1ba+D3f2ba
     & +Dtrigba+Dtriqba+DtriHQba+Da1ba+Da2ba+Da34ba
     & +D2d1ba+renba

c--- sum all diagrams (dd)
c--- note that D2bddu is excluded, since it is equal to D2bdd
      KMggdd(j1,j2,j3,j4)=
     & +D3hdd+D2bdd+D2c1dd+D2c2dd+D2c3dd+D2c4dd+D3f1dd+D3f2dd
     & +Dtrigdd+Dtriqdd+DtriHQdd+Da1dd+Da2dd+Da34dd
     & +D2c1ddu+D2c2ddu+D2c3ddu+D2c4ddu+D2d1ddu
     & +Da1ddu+Da2ddu+Da34ddu
     & +D2d1dd
c--- apply additional factor of (2*Nc) that was present in the numerical
c--- routines, in order to use the same code for summing and squaring
      KMggdd(j1,j2,j3,j4)=KMggdd(j1,j2,j3,j4)*2d0*xn

      if (showG) then
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Bsab',h1,h2,h3,h4,+temps
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Btab',h1,h2,h3,h4,czip!tempt
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Buab',h1,h2,h3,h4,+tempu
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D2e1',h1,h2,h3,h4,
     & D2e1dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D2e2',h1,h2,h3,h4,
     & D2e2dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D3j1',h1,h2,h3,h4,
     & D3j1dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D3j2',h1,h2,h3,h4,
     & D3j2dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D3h ',h1,h2,h3,h4,
     & D3hdd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D2b ',h1,h2,h3,h4,
     & D2bdd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D2c1',h1,h2,h3,h4,
     & D2c1dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D2c2',h1,h2,h3,h4,
     & D2c2dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D2c3',h1,h2,h3,h4,
     & D2c3dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D2c4',h1,h2,h3,h4,
     & D2c4dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D3f1',h1,h2,h3,h4,
     & D3f1dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D3f2',h1,h2,h3,h4,
     & D3f2dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Dtrg',h1,h2,h3,h4,
     & Dtrigdd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Dtrq',h1,h2,h3,h4,
     & Dtriqdd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,DtrQ',h1,h2,h3,h4,
     & DtriHQdd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Da1 ',h1,h2,h3,h4,
     & Da1dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Da2 ',h1,h2,h3,h4,
     & Da2dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Da34',h1,h2,h3,h4,
     & Da34dd/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Ua1 ',h1,h2,h3,h4,
     & Da1ddu/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Ua2 ',h1,h2,h3,h4,
     & Da2ddu/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Ua34',h1,h2,h3,h4,
     & Da34ddu/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,U2b ',h1,h2,h3,h4,
     & D2bddu/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,U2c1',h1,h2,h3,h4,
     & D2c1ddu/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,U2c2',h1,h2,h3,h4,
     & D2c2ddu/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,U2c3',h1,h2,h3,h4,
     & D2c3ddu/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,U2c4',h1,h2,h3,h4,
     & D2c4ddu/temps*(2d0*xn)
      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,D2d1',h1,h2,h3,h4,
     & D2d1dd/temps*(2d0*xn)
c      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,ren ',h1,h2,h3,h4,
c     & renab/temps
c      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,sum2',h1,h2,h3,h4,
c     & (D2d1ab+renab)/temps
c      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Bsba',h1,h2,h3,h4,-temps
c      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Btba',h1,h2,h3,h4,+tempt
c      write(6,'(a18,4i3,2e20.12)')'h1,h2,h3,h4,Buba',h1,h2,h3,h4,czip!tempu
      write(6,*)
      endif
      
      enddo
      enddo
      enddo
      enddo
      
c--- now compute quantities using the numerical method
      spinorsw='xx'


      if (first) then
      first=.false.
      call pvsetmaxindex(3,3) ! Only need up to Cijk and Dijk
      call qlinit
      call pvarraysetup
      endif
      
      call pvsetmudim(scale)
      call pvclearcache


      do nu=1,4
      j=nu-1
      if (j .eq. 0) j=4
      j=nu
      q1(nu)=p(1,j)!/mt
      q2(nu)=p(2,j)!/mt
      q3(nu)=p(3,j)!/mt
      q4(nu)=p(4,j)!/mt
      enddo

C-----remember that diagrams are Qbar(k1)+Q(k2)+q(k3)+qbar(k4)
C     0--> qbar(q1)+q(q2)+Q(q3)+Qb(q4)
C     0--> qbar(k4)+q(k3)+Q(k2)+Qb(k1)
      do j=1,4
      k1(j)=q4(j)
      k2(j)=q3(j)
      k3(j)=q2(j)
      k4(j)=q1(j)
      zk1(j)=dcmplx(k1(j))
      zk2(j)=dcmplx(k2(j))
      zk3(j)=dcmplx(k3(j))
      zk4(j)=dcmplx(k4(j))
      enddo

      
      losq=0d0
      do h1=-1,1,2
      h2=-h1
      do h3=-1,1,2
      do h4=-1,1,2

      call Vspinor_w(zk1,mt,h4,Vm)
      call Ubarspinor_w(zk2,mt,h3,Ubm)
      call Ubarspinor_w(zk3,zip,h1,Ubl)
      call Vspinor_w(zk4,zip,h2,Vl)
      

      call QQqq1(Vm,Ubm,Ubl,Vl,k1,k2,k3,k4,lo1)
      call QQqql1(Vm,Ubm,Ubl,Vl,k1,k2,k3,k4,abtt1) ! S.E. (gluons)
      call QQqql2(Vm,Ubm,Ubl,Vl,k1,k2,k3,k4,abtt2) ! S.E. (light qrks)
      call QQqql3(Vm,Ubm,Ubl,Vl,k1,k2,k3,k4,abtt3) ! S.E. (ghosts)
      call QQqql4(Vm,Ubm,Ubl,Vl,k1,k2,k3,k4,abtt4) ! S.E. (heavy qrks)
      call QQqql5(Vm,Ubm,Ubl,Vl,k1,k2,k3,k4,abtt5) ! Triangle 5c
      call QQqql6(Vm,Ubm,Ubl,Vl,k1,k2,k3,k4,abtt6) ! Triangle 5d
      call QQqql7(Vm,Ubm,Ubl,Vl,k1,k2,k3,k4,abtt7) ! Triangle 5e
      call QQqql8(Vm,Ubm,Ubl,Vl,k1,k2,k3,k4,abtt8) ! Triangle 5f
      call QQqql9(Vm,Ubm,Ubl,Vl,k1,k2,k3,k4,abtt9) ! Box
      call QQqql10(Vm,Ubm,Ubl,Vl,k1,k2,k3,k4,abtt10) ! Box
      call QQqql11(Vm,Ubm,Ubl,Vl,k1,k2,k3,k4,abtt11) ! counter-terms

      if (showQ) then
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,lo1',h1,h2,h3,h4,lo1
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,d10',h1,h2,h3,h4,
     & (abtt10(-2)*epinv**2+abtt10(-1)*epinv+abtt10(0))/lo1
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,d9 ',h1,h2,h3,h4,
     & (abtt9(-2)*epinv**2+abtt9(-1)*epinv+abtt9(0))/lo1
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,d5 ',h1,h2,h3,h4,
     & (abtt5(-2)*epinv**2+abtt5(-1)*epinv+abtt5(0))/lo1
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,d6 ',h1,h2,h3,h4,
     & (abtt6(-2)*epinv**2+abtt6(-1)*epinv+abtt6(0))/lo1
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,d7 ',h1,h2,h3,h4,
     & (abtt7(-2)*epinv**2+abtt7(-1)*epinv+abtt7(0))/lo1
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,d8 ',h1,h2,h3,h4,
     & (abtt8(-2)*epinv**2+abtt8(-1)*epinv+abtt8(0))/lo1
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,d14',h1,h2,h3,h4,
     & (abtt1(-2)*epinv**2+abtt1(-1)*epinv+abtt1(0))/lo1
     &+(abtt2(-2)*epinv**2+abtt2(-1)*epinv+abtt2(0))/lo1
     &+(abtt3(-2)*epinv**2+abtt3(-1)*epinv+abtt3(0))/lo1
     &+(abtt4(-2)*epinv**2+abtt4(-1)*epinv+abtt4(0))/lo1
c      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,d1 ',h1,h2,h3,h4,
c     & (abtt1(-2)*epinv**2+abtt1(-1)*epinv+abtt1(0))/lo1
c      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,d2 ',h1,h2,h3,h4,
c     & (abtt2(-2)*epinv**2+abtt2(-1)*epinv+abtt2(0))/lo1
c      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,d3 ',h1,h2,h3,h4,
c     & (abtt3(-2)*epinv**2+abtt3(-1)*epinv+abtt3(0))/lo1
c      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,d4 ',h1,h2,h3,h4,
c     & (abtt4(-2)*epinv**2+abtt4(-1)*epinv+abtt4(0))/lo1
      write(6,'(a16,4i3,2e20.12)') 'h1,h2,h3,h4,d11',h1,h2,h3,h4,
     & (abtt11(-2)*epinv**2+abtt11(-1)*epinv+abtt11(0))/lo1
      pause
      endif
      losq=losq+abs(lo1)**2

c--- sum all diagrams
      do ep=-2,0
      ho(ep)=
     . +abtt1(ep)+abtt2(ep)+abtt3(ep)+abtt4(ep)+abtt5(ep)
     . +abtt6(ep)+abtt7(ep)+abtt8(ep)+abtt9(ep)+abtt10(ep)+abtt11(ep)
      enddo

c--- add up powers of epsilon
      num=ho(-2)*epinv**2+ho(-1)*epinv+ho(0)      

c--- get indices into array from helicity labels (-1 > 1, +1 -> 2)
      j1=(h1+3)/2
      j3=(h3+3)/2
      j4=(h4+3)/2
      
c--- print out comparison      
      write(6,*) 'KM, num (ab)',KMqqb(j1,j3,j4),num,
     & cdabs(KMqqb(j1,j3,j4)/num)

      enddo
      enddo
      enddo

      if (showQ) then
      write(6,*) 'lo1 squared',losq
      endif
      
      do h1=-1,1,2
      do h2=-1,1,2

      do h3=-1,1,2
      do h4=-1,1,2

      call Vspinor_w(zk1,mt,h4,Vm)
      call Ubarspinor_w(zk2,mt,h3,Ubm)
      call polarization(zk3,h1,ze3)
      call polarization(zk4,h2,ze4)
      
      call QQgg1(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,tab1,tba1,tdd1)
      call QQgg2(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,tab2,tba2,tdd2)
      call QQgg3(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,tab3,tba3,tdd3)

      call QQggl1(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab1,ba1,dd1) ! Bubble 3g2 -> trig 
      call QQggl2(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab2,ba2,dd2) ! Triangle 2b
      call QQggl3(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab3,ba3,dd3)
      call QQggl4(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab4,ba4,dd4)
      call QQggl5(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab5,ba5,dd5) ! S.E. (gluons)
      call QQggl6(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab6,ba6,dd6) ! S.E. (light qrks)
      call QQggl7(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab7,ba7,dd7) ! S.E. (ghosts)
      call QQggl8(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab8,ba8,dd8) ! S.E. (heavy qrks)
      call QQggl9(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab9,ba9,dd9) ! Bubble 2d1
      call QQggl10(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab10,ba10,dd10)

      call QQggl11(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab11,ba11,dd11) ! Triangle 3g1 (gluons) -> trig
      call QQggl12(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab12,ba12,dd12) ! Triangle 3g1 (light qrks) -> triq
      call QQggl13(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab13,ba13,dd13)
      call QQggl14(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab14,ba14,dd14) ! Triangle 3g1 (ghosts) -> trig
      call QQggl15(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab15,ba15,dd15) ! Triangle 3g1 (ghosts) -> trig
      call QQggl16(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab16,ba16,dd16) ! Triangle 3g1 (heavy qrks) -> triHQ
      call QQggl17(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab17,ba17,dd17)
      call QQggl18(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab18,ba18,dd18) ! Triangle 2c4 (u)
      call QQggl19(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab19,ba19,dd19) ! Triangle 2c2 (u)
      call QQggl20(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab20,ba20,dd20)

      call QQggl21(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab21,ba21,dd21)
      call QQggl22(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab22,ba22,dd22)
      call QQggl23(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab23,ba23,dd23)
      call QQggl24(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab24,ba24,dd24) ! Triangle 2c3 (u)
      call QQggl25(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab25,ba25,dd25) ! Triangle 2c1 (u)
      call QQggl26(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab26,ba26,dd26) ! Triangle 3f1
      call QQggl27(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab27,ba27,dd27) ! Triangle 3f2
      call QQggl28(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab28,ba28,dd28) ! Box 2a1
      call QQggl29(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab29,ba29,dd29) ! Box 2a2
      call QQggl30(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab30,ba30,dd30)

      call QQggl31(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab31,ba31,dd31)
      call QQggl32(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab32,ba32,dd32)
      call QQggl33(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab33,ba33,dd33)
      call QQggl34(Vm,Ubm,ze3,ze4,k1,k2,k3,k4,ab34,ba34,dd34)

      if (showG) then
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,tdd1',h1,h2,h3,h4,tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,tdd3',h1,h2,h3,h4,tab3
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,tdd2',h1,h2,h3,h4,tab2
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd5 ',h1,h2,h3,h4,
c     & (dd5(-2)*epinv**2+dd5(-1)*epinv+dd5(0))/tab1
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd6 ',h1,h2,h3,h4,
c     & (dd6(-2)*epinv**2+dd6(-1)*epinv+dd6(0))/tab1
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd7 ',h1,h2,h3,h4,
c     & (dd7(-2)*epinv**2+dd7(-1)*epinv+dd7(0))/tab1
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd8 ',h1,h2,h3,h4,
c     & (dd8(-2)*epinv**2+dd8(-1)*epinv+dd8(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd58',h1,h2,h3,h4,
     & (dd5(-2)*epinv**2+dd5(-1)*epinv+dd5(0))/tab1
     &+(dd6(-2)*epinv**2+dd6(-1)*epinv+dd6(0))/tab1
     &+(dd7(-2)*epinv**2+dd7(-1)*epinv+dd7(0))/tab1
     &+(dd8(-2)*epinv**2+dd8(-1)*epinv+dd8(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd2 ',h1,h2,h3,h4,
     & (dd2(-2)*epinv**2+dd2(-1)*epinv+dd2(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd25',h1,h2,h3,h4,
     & (dd25(-2)*epinv**2+dd25(-1)*epinv+dd25(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd19',h1,h2,h3,h4,
     & (dd19(-2)*epinv**2+dd19(-1)*epinv+dd19(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd24',h1,h2,h3,h4,
     & (dd24(-2)*epinv**2+dd24(-1)*epinv+dd24(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd18',h1,h2,h3,h4,
     & (dd18(-2)*epinv**2+dd18(-1)*epinv+dd18(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd26',h1,h2,h3,h4,
     & (dd26(-2)*epinv**2+dd26(-1)*epinv+dd26(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd27',h1,h2,h3,h4,
     & (dd27(-2)*epinv**2+dd27(-1)*epinv+dd27(0))/tab1
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd1 ',h1,h2,h3,h4,
c     & (dd1(-2)*epinv**2+dd1(-1)*epinv+dd1(0))/tab1
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd11',h1,h2,h3,h4,
c     & (dd11(-2)*epinv**2+dd11(-1)*epinv+dd11(0))/tab1
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd14',h1,h2,h3,h4,
c     & (dd14(-2)*epinv**2+dd14(-1)*epinv+dd14(0))/tab1
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd15',h1,h2,h3,h4,
c     & (dd15(-2)*epinv**2+dd15(-1)*epinv+dd15(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,ddtg',h1,h2,h3,h4,
     & (dd1(-2)*epinv**2+dd1(-1)*epinv+dd1(0))/tab1
     &+(dd11(-2)*epinv**2+dd11(-1)*epinv+dd11(0))/tab1
     &+(dd14(-2)*epinv**2+dd14(-1)*epinv+dd14(0))/tab1
     &+(dd15(-2)*epinv**2+dd15(-1)*epinv+dd15(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd12',h1,h2,h3,h4,
     & (dd12(-2)*epinv**2+dd12(-1)*epinv+dd12(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd16',h1,h2,h3,h4,
     & (dd16(-2)*epinv**2+dd16(-1)*epinv+dd16(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd28',h1,h2,h3,h4,
     & (dd28(-2)*epinv**2+dd28(-1)*epinv+dd28(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd29',h1,h2,h3,h4,
     & (dd29(-2)*epinv**2+dd29(-1)*epinv+dd29(0))/tab1
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd9 ',h1,h2,h3,h4,
c     & (dd9(-2)*epinv**2+dd9(-1)*epinv+dd9(0))/tab1
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd34',h1,h2,h3,h4,
c     & (dd34(-2)*epinv**2+dd34(-1)*epinv+dd34(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,sum2',h1,h2,h3,h4,
     & (dd9(-2)*epinv**2+dd9(-1)*epinv+dd9(0))/tab1
     &+(dd34(-2)*epinv**2+dd34(-1)*epinv+dd34(0))/tab1

      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd3 ',h1,h2,h3,h4,
     & (dd3(-2)*epinv**2+dd3(-1)*epinv+dd3(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd4 ',h1,h2,h3,h4,
     & (dd4(-2)*epinv**2+dd4(-1)*epinv+dd4(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd10',h1,h2,h3,h4,
     & (dd10(-2)*epinv**2+dd10(-1)*epinv+dd10(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd13',h1,h2,h3,h4,
     & (dd13(-2)*epinv**2+dd13(-1)*epinv+dd13(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd17',h1,h2,h3,h4,
     & (dd17(-2)*epinv**2+dd17(-1)*epinv+dd17(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd20',h1,h2,h3,h4,
     & (dd20(-2)*epinv**2+dd20(-1)*epinv+dd20(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd21',h1,h2,h3,h4,
     & (dd21(-2)*epinv**2+dd21(-1)*epinv+dd21(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd22',h1,h2,h3,h4,
     & (dd22(-2)*epinv**2+dd22(-1)*epinv+dd22(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd23',h1,h2,h3,h4,
     & (dd23(-2)*epinv**2+dd23(-1)*epinv+dd23(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd30',h1,h2,h3,h4,
     & (dd30(-2)*epinv**2+dd30(-1)*epinv+dd30(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd31',h1,h2,h3,h4,
     & (dd31(-2)*epinv**2+dd31(-1)*epinv+dd31(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd32',h1,h2,h3,h4,
     & (dd32(-2)*epinv**2+dd32(-1)*epinv+dd32(0))/tab1
      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,dd33',h1,h2,h3,h4,
     & (dd33(-2)*epinv**2+dd33(-1)*epinv+dd33(0))/tab1
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,tba1',h1,h2,h3,h4,tba1
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,tba3',h1,h2,h3,h4,tba3
c      write(6,'(a18,4i3,2e20.12)') 'h1,h2,h3,h4,tba2',h1,h2,h3,h4,tba2
      pause
      endif

c--- sum all diagrams
      do ep=-2,0
       hoab(ep)=
     . +ab1(ep)+ab2(ep)+ab3(ep)+ab4(ep)+ab5(ep)
     . +ab6(ep)+ab7(ep)+ab8(ep)+ab9(ep)+ab10(ep)
     . +ab11(ep)+ab12(ep)+ab13(ep)+ab14(ep)+ab15(ep)
     . +ab16(ep)+ab17(ep)+ab18(ep)+ab19(ep)+ab20(ep)
     . +ab21(ep)+ab22(ep)+ab23(ep)+ab24(ep)+ab25(ep)
     . +ab26(ep)+ab27(ep)+ab28(ep)+ab29(ep)+ab30(ep)
     . +ab31(ep)+ab32(ep)+ab33(ep)+ab34(ep)

       hoba(ep)=
     . +ba1(ep)+ba2(ep)+ba3(ep)+ba4(ep)+ba5(ep)
     . +ba6(ep)+ba7(ep)+ba8(ep)+ba9(ep)+ba10(ep)
     . +ba11(ep)+ba12(ep)+ba13(ep)+ba14(ep)+ba15(ep)
     . +ba16(ep)+ba17(ep)+ba18(ep)+ba19(ep)+ba20(ep)
     . +ba21(ep)+ba22(ep)+ba23(ep)+ba24(ep)+ba25(ep)
     . +ba26(ep)+ba27(ep)+ba28(ep)+ba29(ep)+ba30(ep)
     . +ba31(ep)+ba32(ep)+ba33(ep)+ba34(ep)

       hodd(ep)=
     . +dd1(ep)+dd2(ep)+dd3(ep)+dd4(ep)+dd5(ep)
     . +dd6(ep)+dd7(ep)+dd8(ep)+dd9(ep)+dd10(ep)
     . +dd11(ep)+dd12(ep)+dd13(ep)+dd14(ep)+dd15(ep)
     . +dd16(ep)+dd17(ep)+dd18(ep)+dd19(ep)+dd20(ep)
     . +dd21(ep)+dd22(ep)+dd23(ep)+dd24(ep)+dd25(ep)
     . +dd26(ep)+dd27(ep)+dd28(ep)+dd29(ep)+dd30(ep)
     . +dd31(ep)+dd32(ep)+dd33(ep)+dd34(ep)
      enddo

c--- add up powers of epsilon
      numab=hoab(-2)*epinv**2+hoab(-1)*epinv+hoab(0)      
      numba=hoba(-2)*epinv**2+hoba(-1)*epinv+hoba(0)      
      numdd=hodd(-2)*epinv**2+hodd(-1)*epinv+hodd(0)      

c--- get indices into array from helicity labels (-1 > 1, +1 -> 2)
      j1=(h1+3)/2
      j2=(h2+3)/2
      j3=(h3+3)/2
      j4=(h4+3)/2
      
c--- print out comparison      
c      write(6,*) 'KM, num (ab)',KMggab(j1,j2,j3,j4),numab,
c     & cdabs(KMggab(j1,j2,j3,j4)/numab)
c      write(6,*) 'KM, num (ba)',KMggba(j1,j2,j3,j4),
c     & numba,cdabs(KMggba(j1,j2,j3,j4)/numba)
c      write(6,*) 'KM, num (dd)',KMggdd(j1,j2,j3,j4),numdd,
c     & cdabs(KMggdd(j1,j2,j3,j4)/numdd)
      
      enddo
      enddo
      enddo
      enddo
      
      pause

      return
      end
      
