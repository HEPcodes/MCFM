      subroutine fundk(p,wtqqb,wtqbq,wtgg)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      double precision p(mxpart,4),ps(mxpart,4),fac,propw,
     . s34,s78,s35,s68,wtqbq,wtqqb,wtgg,
     . t1,t2,ro,tbar,ubar,b,xlp,xlm,
     . vlsm,vdt,vdw,vltm,vlwm,vdmp,vlpm,vlbl,vdmb,
     . f1,f2,f3,rmuom2,GRAMINV,bracks,brackt,brackw,ddilog,
     . AAAE,AAA0,AAA1,AAA2,AAA3,AAA4,AAA5,AAA6,XLF,
     . LOSQ,NLO0,NLO1,NLO2,NLO3,NLO4,NLO5,NLO6
      parameter(XLF=5d0)
      integer j,nu,kswap
      do kswap=1,2
      do nu=1,4
      if (kswap .eq. 1) then
      ps(1,nu)=p(1,nu)
      ps(2,nu)=p(2,nu)
      elseif (kswap .eq. 2) then
      ps(2,nu)=p(1,nu)
      ps(1,nu)=p(2,nu)
      endif
      do j=3,mxpart
      ps(j,nu)=zip
      enddo
      enddo
      
      s34=2d0*(p(3,4)*p(4,4)-p(3,3)*p(4,3)-p(3,2)*p(4,2)-p(3,1)*p(4,1))
      s78=2d0*(p(7,4)*p(8,4)-p(7,3)*p(8,3)-p(7,2)*p(8,2)-p(7,1)*p(8,1))
      s35=2d0*(p(3,4)*p(5,4)-p(3,3)*p(5,3)-p(3,2)*p(5,2)-p(3,1)*p(5,1))
      s68=2d0*(p(6,4)*p(8,4)-p(6,3)*p(8,3)-p(6,2)*p(8,2)-p(6,1)*p(8,1))
      propw=((s34-wmass**2)**2+(wmass*wwidth)**2)
     .     *((s78-wmass**2)**2+(wmass*wwidth)**2)
      fac=gwsq**4*s35*s68/((mt*twidth)**4*propw)
     . *gsq**2*ason2pi*V/4d0
c      we will have no further need for p3 and p5
c      we will have no further need for p6 and p8
      do nu=1,4
c t momentum
      ps(3,nu)=p(3,nu)+p(4,nu)+p(5,nu)
c t-bar momentum
      ps(4,nu)=p(6,nu)+p(7,nu)+p(8,nu)
c positron
      ps(5,nu)=p(4,nu)
c electron
      ps(6,nu)=p(7,nu)
      enddo
           
      call dotem(6,ps,s)
      ro=4d0*mt**2/s(1,2)
      rmuom2=2d0*log(scale/mt)
      t1=-s(2,4)/s(1,2)
      t2=-s(2,3)/s(1,2)
      tbar=t1-0.25*ro
      ubar=t2-0.25*ro
      b=sqrt(1.d0-ro)
      xlp=0.5*(1.d0+b)
      xlm=0.5*(1.d0-b)
      vlpm=log(xlp/xlm)
      vlsm=log(4.d0/ro)
      vltm=log(4.d0*t1/ro)
      vlwm=log(4.d0*t2/ro)
      vlbl=log(b/xlm)
      vdw=ddilog(1.d0-ro/(4.d0*t2))-0.5d0*vlwm**2
      vdt=ddilog(1.d0-ro/(4.d0*t1))-0.5d0*vltm**2
      vdmp=ddilog(-xlm/xlp)
      vdmb=ddilog(-xlm/b)+0.5d0*vlbl**2
      f1=(vlpm**2/2.d0-2*vdmb-pi**2/3.d0)/b
      f2=(-b*vlsm+vlpm**2/4.d0+vdmp+pi**2/12.d0)/b**3
      f3=(-b**3*vlsm-3*b*vlsm+0.75d0*vlpm**2
     & +3.d0*vdmp+pi**2/4.d0+2.d0*b**3)/b**5
      GRAMINV = 1d0/(t1*t2-0.25d0*ro)
      brackt=vlsm*vltm-vlsm**2/4d0-pisq/12d0+vdt
      brackw=vlsm*vlwm-vlsm**2/4d0-pisq/12d0+vdw
      bracks=vdmp+vlpm**2/4d0+pisq/12d0
      AAAE =  + b**(-2)*vlsm*xn**(-1) * (  - 1.D0/2.D0*t1**(-1)*GRAMINV
     &     - 2.D0*t1**(-1) - 2.D0*t1*GRAMINV + 2.D0*GRAMINV )
      AAAE = AAAE + b**(-1)*xn**(-1) * (  - 1.D0/2.D0*pisq )
      AAAE = AAAE + b*vlpm*TR * ( 2.D0 )
      AAAE = AAAE + b*vlpm*xn**(-1) * (  - 3.D0/2.D0 )
      AAAE = AAAE + b*xn**(-1)*bracks * ( 1.D0/2.D0*t1**(-1)*GRAMINV + 
     &    5.D0/32.D0*t1**(-1)*GRAMINV**2 - 1.D0/2.D0*t1**(-1) + 13.D0/4.
     &    D0*t1*GRAMINV**2 - 2.D0*t1**2*GRAMINV**2 + 1.D0/2.D0*t1**3*
     &    GRAMINV**2 - 3.D0/2.D0*GRAMINV**2 )
      AAAE = AAAE + b*xn**(-1) * (  - 1.D0/2.D0*pisq )
      AAAE = AAAE + b*xn*bracks * (  - 1.D0/4.D0*t1**(-1)*GRAMINV - 5.D0
     &    /64.D0*t1**(-1)*GRAMINV**2 + 1.D0/4.D0*t1**(-1) + 1.D0/4.D0*
     &    t1*GRAMINV - 9.D0/16.D0*t1*GRAMINV**2 + 1.D0/4.D0*t1**2*
     &    GRAMINV**2 - 1.D0/2.D0*GRAMINV + 3.D0/8.D0*GRAMINV**2 )
      AAAE = AAAE + b**3*vlpm*TR * (  - 2.D0/3.D0 )
      AAAE = AAAE + b**3*xn**(-1)*bracks * ( t1**(-1)*GRAMINV + 1.D0/16.
     &    D0*t1**(-1)*GRAMINV**2 - 1.D0/4.D0*t1*GRAMINV**2 - 1.D0/2.D0*
     &    GRAMINV**2 )
      AAAE = AAAE + b**3*xn*bracks * (  - 1.D0/4.D0*t1**(-1)*GRAMINV + 
     &    1.D0/32.D0*t1**(-1)*GRAMINV**2 + 1.D0/16.D0*t1*GRAMINV**2 + 1.
     &    D0/8.D0*GRAMINV**2 )
      AAAE = AAAE + b**5*xn**(-1)*bracks * (  - 7.D0/32.D0*t1**(-1)*
     &    GRAMINV**2 )
      AAAE = AAAE + b**5*xn*bracks * ( 3.D0/64.D0*t1**(-1)*GRAMINV**2 )
      AAAE = AAAE + vlsm*TR * ( 4.D0/3.D0*XLF )
      AAAE = AAAE + vlsm*xn**(-1) * (  - 3.D0/2.D0 + t1**(-1)*ro*
     &    GRAMINV + 5.D0/8.D0*t1**(-1)*ro*GRAMINV**2 + 2.D0*t1**(-1)*ro
     &     - 1.D0/2.D0*t1**(-1)*ro**2*GRAMINV + 11.D0/32.D0*t1**(-1)*
     &    ro**2*GRAMINV**2 - 1.D0/4.D0*t1**(-1)*ro**3*GRAMINV**2 + 3.D0
     &    *t1**(-1)*GRAMINV - 3.D0/2.D0*t1**(-1) - 7.D0/4.D0*t1*ro*
     &    GRAMINV**2 - 4.D0*t1*GRAMINV + 5.D0*t1*GRAMINV**2 - 1.D0/2.D0
     &    *t1**3*GRAMINV**2 - 2.D0*ro*GRAMINV - 5.D0/2.D0*ro*GRAMINV**2
     &     + 3.D0/2.D0*ro**2*GRAMINV**2 + 2.D0*GRAMINV - 2.D0*
     &    GRAMINV**2 )
      AAAE = AAAE + vlsm*xn * (  - 2.D0/3.D0 - 3.D0/4.D0*t1**(-1)*ro*
     &    GRAMINV + 1.D0/8.D0*t1**(-1)*ro*GRAMINV**2 - 1.D0/2.D0*
     &    t1**(-1)*ro + 1.D0/8.D0*t1**(-1)*ro**2*GRAMINV - 15.D0/64.D0*
     &    t1**(-1)*ro**2*GRAMINV**2 + 1.D0/16.D0*t1**(-1)*ro**3*
     &    GRAMINV**2 + 1.D0/2.D0*t1**(-1)*GRAMINV + 3.D0/4.D0*t1**(-1)
     &     + 5.D0/16.D0*t1*ro*GRAMINV**2 + 3.D0/4.D0*t1*GRAMINV - 1.D0/
     &    2.D0*t1*GRAMINV**2 - 1.D0/4.D0*t1**2*GRAMINV**2 + 9.D0/8.D0*
     &    ro*GRAMINV**2 - 3.D0/8.D0*ro**2*GRAMINV**2 + 1.D0/2.D0*
     &    GRAMINV - 1.D0/2.D0*GRAMINV**2 )
      AAAE = AAAE + vlsm**2*xn**(-1) * ( 1.D0/2.D0 )
      AAAE = AAAE + vltm*TBAR**(-1)*xn**(-1) * (  - 2.D0 - 1.D0/2.D0*
     &    t1**(-1)*ro*GRAMINV + t1**(-1)*ro + 1.D0/4.D0*t1**(-1)*ro**2*
     &    GRAMINV - 2.D0*t1**(-1) + 1.D0/2.D0*t1*ro*GRAMINV - 4.D0*t1
     &     - 3.D0/2.D0*ro*GRAMINV + ro + 2.D0*GRAMINV )
      AAAE = AAAE + vltm*TBAR**(-1)*xn * ( 1.D0 + 1.D0/4.D0*t1**(-1)*ro
     &    *GRAMINV - 1.D0/2.D0*t1**(-1)*ro - 1.D0/8.D0*t1**(-1)*ro**2*
     &    GRAMINV + t1**(-1) - 1.D0/4.D0*t1*ro*GRAMINV + 2.D0*t1 + 3.D0/
     &    4.D0*ro*GRAMINV - 1.D0/2.D0*ro - GRAMINV )
      AAAE = AAAE + vlwm*UBAR**(-1)*xn**(-1) * ( 2.D0 - t1**(-1)*ro - 1.
     &    D0/4.D0*t1**(-1)*ro**2*GRAMINV - t1*ro*GRAMINV - 4.D0*t1 + 1.D
     &    0/2.D0*t1**2*ro*GRAMINV + ro*GRAMINV - 1.D0/2.D0*ro + 1.D0/8.D
     &    0*ro**2*GRAMINV )
      AAAE = AAAE + vdw*xn**(-1) * ( 4.D0 )
      AAAE = AAAE + vdt*xn**(-1) * (  - 4.D0 )
      AAAE = AAAE + vdt*xn * ( 2.D0 )
      AAAE = AAAE + TR*rmuom2 * (  - 4.D0/3.D0*XLF )
      AAAE = AAAE + TR * (  - 20.D0/9.D0 - 20.D0/9.D0*XLF - 4.D0/3.D0*
     &    ro )
      AAAE = AAAE + rmuom2*xn * ( 11.D0/3.D0 )
      AAAE = AAAE + xn**(-1)*brackt * ( 4.D0 - 5.D0/2.D0*t1**(-1)*ro*
     &    GRAMINV + t1**(-1)*ro*GRAMINV**2 + t1**(-1)*ro + 3.D0/8.D0*
     &    t1**(-1)*ro**2*GRAMINV - 5.D0/8.D0*t1**(-1)*ro**2*GRAMINV**2
     &     + 1.D0/32.D0*t1**(-1)*ro**3*GRAMINV**2 + 4.D0*t1**(-1)*
     &    GRAMINV + t1*ro*GRAMINV - 5.D0/2.D0*t1*ro*GRAMINV**2 + 1.D0/8.
     &    D0*t1*ro**2*GRAMINV**2 - 2.D0*t1*GRAMINV + 2.D0*t1*GRAMINV**2
     &     + 1.D0/2.D0*t1**2*ro*GRAMINV**2 + 3.D0*ro*GRAMINV**2 + 2.D0*
     &    GRAMINV - 4.D0*GRAMINV**2 )
      AAAE = AAAE + xn**(-1)*brackw * (  - 4.D0 + 4.D0*t1**(-1)*ro*
     &    GRAMINV - t1**(-1)*ro*GRAMINV**2 - 1.D0/2.D0*t1**(-1)*ro - 1.D
     &    0/8.D0*t1**(-1)*ro**2*GRAMINV + 7.D0/8.D0*t1**(-1)*ro**2*
     &    GRAMINV**2 - 4.D0*t1**(-1)*GRAMINV + 2.D0*t1**(-1) + 5.D0*t1*
     &    ro*GRAMINV**2 + 1.D0/8.D0*t1*ro**2*GRAMINV**2 + 6.D0*t1*
     &    GRAMINV - 6.D0*t1*GRAMINV**2 - 5.D0/2.D0*t1**2*ro*GRAMINV**2
     &     + 2.D0*t1**2*GRAMINV**2 + 1.D0/2.D0*t1**3*ro*GRAMINV**2 - 2.D
     &    0*ro*GRAMINV - 3.D0*ro*GRAMINV**2 - 1.D0/2.D0*ro**2*
     &    GRAMINV**2 - 4.D0*GRAMINV + 4.D0*GRAMINV**2 )
      AAAE = AAAE + xn**(-1) * ( 6.D0 - 1.D0/3.D0*pisq )
      AAAE = AAAE + xn*brackt * (  - 2.D0 + 5.D0/4.D0*t1**(-1)*ro*
     &    GRAMINV - 1.D0/2.D0*t1**(-1)*ro*GRAMINV**2 - 1.D0/2.D0*
     &    t1**(-1)*ro - 3.D0/16.D0*t1**(-1)*ro**2*GRAMINV + 5.D0/16.D0*
     &    t1**(-1)*ro**2*GRAMINV**2 - 1.D0/64.D0*t1**(-1)*ro**3*
     &    GRAMINV**2 - 2.D0*t1**(-1)*GRAMINV - 1.D0/2.D0*t1*ro*GRAMINV
     &     + 5.D0/4.D0*t1*ro*GRAMINV**2 - 1.D0/16.D0*t1*ro**2*
     &    GRAMINV**2 + t1*GRAMINV - t1*GRAMINV**2 - 1.D0/4.D0*t1**2*ro*
     &    GRAMINV**2 - 3.D0/2.D0*ro*GRAMINV**2 - GRAMINV + 2.D0*
     &    GRAMINV**2 )
      AAAE = AAAE + xn * ( 31.D0/9.D0 - 1.D0/3.D0*pisq )
      AAAE = AAAE + f1*xn**(-1) * (  - 1.D0 + 1.D0/2.D0*ro )
      AAAE = AAAE + f2*xn**(-1) * (  - 5.D0/2.D0*t1**(-1)*ro*GRAMINV + 
     &    5.D0/8.D0*t1**(-1)*ro*GRAMINV**2 + 3.D0/2.D0*t1**(-1)*ro - 2.D
     &    0*t1**(-1)*ro**2*GRAMINV - 9.D0/32.D0*t1**(-1)*ro**2*
     &    GRAMINV**2 - 2.D0*t1**(-1)*ro**2 + 1.D0/2.D0*t1**(-1)*ro**3*
     &    GRAMINV - 19.D0/32.D0*t1**(-1)*ro**3*GRAMINV**2 + 1.D0/4.D0*
     &    t1**(-1)*ro**4*GRAMINV**2 + 5.D0/2.D0*t1**(-1)*GRAMINV - 11.D0
     &    /2.D0*t1**(-1) + 2.D0*t1*ro*GRAMINV - 27.D0/4.D0*t1*ro*
     &    GRAMINV**2 + 7.D0/4.D0*t1*ro**2*GRAMINV**2 - 8.D0*t1*GRAMINV
     &     + 5.D0*t1*GRAMINV**2 + 1.D0/2.D0*t1**3*ro*GRAMINV**2 - 1.D0/
     &    2.D0*t1**3*GRAMINV**2 - 2.D0*ro*GRAMINV - 1.D0/2.D0*ro*
     &    GRAMINV**2 + 2.D0*ro**2*GRAMINV + 4.D0*ro**2*GRAMINV**2 - 3.D0
     &    /2.D0*ro**3*GRAMINV**2 + 6.D0*GRAMINV - 2.D0*GRAMINV**2 )
      AAAE = AAAE + f2*xn * ( 2.D0 - t1**(-1)*ro*GRAMINV + 1.D0/8.D0*
     &    t1**(-1)*ro*GRAMINV**2 - 3.D0/2.D0*t1**(-1)*ro + 13.D0/16.D0*
     &    t1**(-1)*ro**2*GRAMINV - 23.D0/64.D0*t1**(-1)*ro**2*
     &    GRAMINV**2 + 1.D0/2.D0*t1**(-1)*ro**2 - 1.D0/8.D0*t1**(-1)*
     &    ro**3*GRAMINV + 19.D0/64.D0*t1**(-1)*ro**3*GRAMINV**2 - 1.D0/
     &    16.D0*t1**(-1)*ro**4*GRAMINV**2 + 1.D0/2.D0*t1**(-1)*GRAMINV
     &     + 7.D0/4.D0*t1**(-1) - 1.D0/2.D0*t1*ro*GRAMINV + 13.D0/16.D0
     &    *t1*ro*GRAMINV**2 - 5.D0/16.D0*t1*ro**2*GRAMINV**2 + 3.D0/4.D0
     &    *t1*GRAMINV - 1.D0/2.D0*t1*GRAMINV**2 + 1.D0/4.D0*t1**2*ro*
     &    GRAMINV**2 - 1.D0/4.D0*t1**2*GRAMINV**2 + 13.D0/8.D0*ro*
     &    GRAMINV**2 - 1.D0/2.D0*ro - 3.D0/2.D0*ro**2*GRAMINV**2 + 3.D0/
     &    8.D0*ro**3*GRAMINV**2 - 1.D0/2.D0*GRAMINV - 1.D0/2.D0*
     &    GRAMINV**2 )

      AAA0 =  + b**(-2)*vlsm*xn**(-1) * (  - 1.D0/2.D0*t1**(-1)*GRAMINV
     &     - 2.D0*t1**(-1) - 2.D0*t1*GRAMINV + 2.D0*GRAMINV )
      AAA0 = AAA0 + b**(-1)*xn**(-1) * (  - 1.D0/2.D0*pisq )
      AAA0 = AAA0 + b*vlpm*TR * ( 2.D0 )
      AAA0 = AAA0 + b*vlpm*xn**(-1) * (  - 3.D0/2.D0 )
      AAA0 = AAA0 + b*xn**(-1)*bracks * ( 1.D0/2.D0*t1**(-1)*GRAMINV + 
     &    5.D0/32.D0*t1**(-1)*GRAMINV**2 - 1.D0/2.D0*t1**(-1) + 13.D0/4.
     &    D0*t1*GRAMINV**2 - 2.D0*t1**2*GRAMINV**2 + 1.D0/2.D0*t1**3*
     &    GRAMINV**2 - 3.D0/2.D0*GRAMINV**2 )
      AAA0 = AAA0 + b*xn**(-1) * (  - 1.D0/2.D0*pisq )
      AAA0 = AAA0 + b*xn*bracks * (  - 1.D0/4.D0*t1**(-1)*GRAMINV - 5.D0
     &    /64.D0*t1**(-1)*GRAMINV**2 + 1.D0/4.D0*t1**(-1) + 1.D0/4.D0*
     &    t1*GRAMINV - 9.D0/16.D0*t1*GRAMINV**2 + 1.D0/4.D0*t1**2*
     &    GRAMINV**2 - 1.D0/2.D0*GRAMINV + 3.D0/8.D0*GRAMINV**2 )
      AAA0 = AAA0 + b**3*vlpm*TR * (  - 2.D0/3.D0 )
      AAA0 = AAA0 + b**3*xn**(-1)*bracks * ( t1**(-1)*GRAMINV + 1.D0/16.
     &    D0*t1**(-1)*GRAMINV**2 - 1.D0/4.D0*t1*GRAMINV**2 - 1.D0/2.D0*
     &    GRAMINV**2 )
      AAA0 = AAA0 + b**3*xn*bracks * (  - 1.D0/4.D0*t1**(-1)*GRAMINV + 
     &    1.D0/32.D0*t1**(-1)*GRAMINV**2 + 1.D0/16.D0*t1*GRAMINV**2 + 1.
     &    D0/8.D0*GRAMINV**2 )
      AAA0 = AAA0 + b**5*xn**(-1)*bracks * (  - 7.D0/32.D0*t1**(-1)*
     &    GRAMINV**2 )
      AAA0 = AAA0 + b**5*xn*bracks * ( 3.D0/64.D0*t1**(-1)*GRAMINV**2 )
      AAA0 = AAA0 + vlsm*TR * ( 4.D0/3.D0*XLF )
      AAA0 = AAA0 + vlsm*xn**(-1) * (  - 3.D0/2.D0 + t1**(-1)*ro*
     &    GRAMINV + 5.D0/8.D0*t1**(-1)*ro*GRAMINV**2 + 2.D0*t1**(-1)*ro
     &     - 1.D0/2.D0*t1**(-1)*ro**2*GRAMINV + 11.D0/32.D0*t1**(-1)*
     &    ro**2*GRAMINV**2 - 1.D0/4.D0*t1**(-1)*ro**3*GRAMINV**2 + 3.D0
     &    *t1**(-1)*GRAMINV - 3.D0/2.D0*t1**(-1) - 7.D0/4.D0*t1*ro*
     &    GRAMINV**2 - 4.D0*t1*GRAMINV + 5.D0*t1*GRAMINV**2 - 1.D0/2.D0
     &    *t1**3*GRAMINV**2 - 2.D0*ro*GRAMINV - 5.D0/2.D0*ro*GRAMINV**2
     &     + 3.D0/2.D0*ro**2*GRAMINV**2 + 2.D0*GRAMINV - 2.D0*
     &    GRAMINV**2 )
      AAA0 = AAA0 + vlsm*xn * (  - 2.D0/3.D0 - 3.D0/4.D0*t1**(-1)*ro*
     &    GRAMINV + 1.D0/8.D0*t1**(-1)*ro*GRAMINV**2 - 1.D0/2.D0*
     &    t1**(-1)*ro + 1.D0/8.D0*t1**(-1)*ro**2*GRAMINV - 15.D0/64.D0*
     &    t1**(-1)*ro**2*GRAMINV**2 + 1.D0/16.D0*t1**(-1)*ro**3*
     &    GRAMINV**2 + 1.D0/2.D0*t1**(-1)*GRAMINV + 3.D0/4.D0*t1**(-1)
     &     + 5.D0/16.D0*t1*ro*GRAMINV**2 + 3.D0/4.D0*t1*GRAMINV - 1.D0/
     &    2.D0*t1*GRAMINV**2 - 1.D0/4.D0*t1**2*GRAMINV**2 + 9.D0/8.D0*
     &    ro*GRAMINV**2 - 3.D0/8.D0*ro**2*GRAMINV**2 + 1.D0/2.D0*
     &    GRAMINV - 1.D0/2.D0*GRAMINV**2 )
      AAA0 = AAA0 + vlsm**2*xn**(-1) * ( 1.D0/2.D0 )
      AAA0 = AAA0 + vltm*TBAR**(-1)*xn**(-1) * (  - 2.D0 - 1.D0/2.D0*
     &    t1**(-1)*ro*GRAMINV + t1**(-1)*ro + 1.D0/4.D0*t1**(-1)*ro**2*
     &    GRAMINV - 2.D0*t1**(-1) + 1.D0/2.D0*t1*ro*GRAMINV - 4.D0*t1
     &     - 3.D0/2.D0*ro*GRAMINV + ro + 2.D0*GRAMINV )
      AAA0 = AAA0 + vltm*TBAR**(-1)*xn * ( 1.D0 + 1.D0/4.D0*t1**(-1)*ro
     &    *GRAMINV - 1.D0/2.D0*t1**(-1)*ro - 1.D0/8.D0*t1**(-1)*ro**2*
     &    GRAMINV + t1**(-1) - 1.D0/4.D0*t1*ro*GRAMINV + 2.D0*t1 + 3.D0/
     &    4.D0*ro*GRAMINV - 1.D0/2.D0*ro - GRAMINV )
      AAA0 = AAA0 + vlwm*UBAR**(-1)*xn**(-1) * ( 2.D0 - t1**(-1)*ro - 1.
     &    D0/4.D0*t1**(-1)*ro**2*GRAMINV - t1*ro*GRAMINV - 4.D0*t1 + 1.D
     &    0/2.D0*t1**2*ro*GRAMINV + ro*GRAMINV - 1.D0/2.D0*ro + 1.D0/8.D
     &    0*ro**2*GRAMINV )
      AAA0 = AAA0 + vdw*xn**(-1) * ( 4.D0 )
      AAA0 = AAA0 + vdt*xn**(-1) * (  - 4.D0 )
      AAA0 = AAA0 + vdt*xn * ( 2.D0 )
      AAA0 = AAA0 + TR*rmuom2 * (  - 4.D0/3.D0*XLF )
      AAA0 = AAA0 + TR * (  - 20.D0/9.D0 - 20.D0/9.D0*XLF - 4.D0/3.D0*
     &    ro )
      AAA0 = AAA0 + rmuom2*xn * ( 11.D0/3.D0 )
      AAA0 = AAA0 + xn**(-1)*brackt * ( 4.D0 - 5.D0/2.D0*t1**(-1)*ro*
     &    GRAMINV + t1**(-1)*ro*GRAMINV**2 + t1**(-1)*ro + 3.D0/8.D0*
     &    t1**(-1)*ro**2*GRAMINV - 5.D0/8.D0*t1**(-1)*ro**2*GRAMINV**2
     &     + 1.D0/32.D0*t1**(-1)*ro**3*GRAMINV**2 + 4.D0*t1**(-1)*
     &    GRAMINV + t1*ro*GRAMINV - 5.D0/2.D0*t1*ro*GRAMINV**2 + 1.D0/8.
     &    D0*t1*ro**2*GRAMINV**2 - 2.D0*t1*GRAMINV + 2.D0*t1*GRAMINV**2
     &     + 1.D0/2.D0*t1**2*ro*GRAMINV**2 + 3.D0*ro*GRAMINV**2 + 2.D0*
     &    GRAMINV - 4.D0*GRAMINV**2 )
      AAA0 = AAA0 + xn**(-1)*brackw * (  - 4.D0 + 4.D0*t1**(-1)*ro*
     &    GRAMINV - t1**(-1)*ro*GRAMINV**2 - 1.D0/2.D0*t1**(-1)*ro - 1.D
     &    0/8.D0*t1**(-1)*ro**2*GRAMINV + 7.D0/8.D0*t1**(-1)*ro**2*
     &    GRAMINV**2 - 4.D0*t1**(-1)*GRAMINV + 2.D0*t1**(-1) + 5.D0*t1*
     &    ro*GRAMINV**2 + 1.D0/8.D0*t1*ro**2*GRAMINV**2 + 6.D0*t1*
     &    GRAMINV - 6.D0*t1*GRAMINV**2 - 5.D0/2.D0*t1**2*ro*GRAMINV**2
     &     + 2.D0*t1**2*GRAMINV**2 + 1.D0/2.D0*t1**3*ro*GRAMINV**2 - 2.D
     &    0*ro*GRAMINV - 3.D0*ro*GRAMINV**2 - 1.D0/2.D0*ro**2*
     &    GRAMINV**2 - 4.D0*GRAMINV + 4.D0*GRAMINV**2 )
      AAA0 = AAA0 + xn**(-1) * ( 6.D0 - 1.D0/3.D0*pisq )
      AAA0 = AAA0 + xn*brackt * (  - 2.D0 + 5.D0/4.D0*t1**(-1)*ro*
     &    GRAMINV - 1.D0/2.D0*t1**(-1)*ro*GRAMINV**2 - 1.D0/2.D0*
     &    t1**(-1)*ro - 3.D0/16.D0*t1**(-1)*ro**2*GRAMINV + 5.D0/16.D0*
     &    t1**(-1)*ro**2*GRAMINV**2 - 1.D0/64.D0*t1**(-1)*ro**3*
     &    GRAMINV**2 - 2.D0*t1**(-1)*GRAMINV - 1.D0/2.D0*t1*ro*GRAMINV
     &     + 5.D0/4.D0*t1*ro*GRAMINV**2 - 1.D0/16.D0*t1*ro**2*
     &    GRAMINV**2 + t1*GRAMINV - t1*GRAMINV**2 - 1.D0/4.D0*t1**2*ro*
     &    GRAMINV**2 - 3.D0/2.D0*ro*GRAMINV**2 - GRAMINV + 2.D0*
     &    GRAMINV**2 )
      AAA0 = AAA0 + xn * ( 31.D0/9.D0 - 1.D0/3.D0*pisq )
      AAA0 = AAA0 + f1*xn**(-1) * (  - 1.D0 + 1.D0/2.D0*ro )
      AAA0 = AAA0 + f2*xn**(-1) * (  - 5.D0/2.D0*t1**(-1)*ro*GRAMINV + 
     &    5.D0/8.D0*t1**(-1)*ro*GRAMINV**2 + 3.D0/2.D0*t1**(-1)*ro - 2.D
     &    0*t1**(-1)*ro**2*GRAMINV - 9.D0/32.D0*t1**(-1)*ro**2*
     &    GRAMINV**2 - 2.D0*t1**(-1)*ro**2 + 1.D0/2.D0*t1**(-1)*ro**3*
     &    GRAMINV - 19.D0/32.D0*t1**(-1)*ro**3*GRAMINV**2 + 1.D0/4.D0*
     &    t1**(-1)*ro**4*GRAMINV**2 + 5.D0/2.D0*t1**(-1)*GRAMINV - 11.D0
     &    /2.D0*t1**(-1) + 2.D0*t1*ro*GRAMINV - 27.D0/4.D0*t1*ro*
     &    GRAMINV**2 + 7.D0/4.D0*t1*ro**2*GRAMINV**2 - 8.D0*t1*GRAMINV
     &     + 5.D0*t1*GRAMINV**2 + 1.D0/2.D0*t1**3*ro*GRAMINV**2 - 1.D0/
     &    2.D0*t1**3*GRAMINV**2 - 2.D0*ro*GRAMINV - 1.D0/2.D0*ro*
     &    GRAMINV**2 + 2.D0*ro**2*GRAMINV + 4.D0*ro**2*GRAMINV**2 - 3.D0
     &    /2.D0*ro**3*GRAMINV**2 + 6.D0*GRAMINV - 2.D0*GRAMINV**2 )
      AAA0 = AAA0 + f2*xn * ( 2.D0 - t1**(-1)*ro*GRAMINV + 1.D0/8.D0*
     &    t1**(-1)*ro*GRAMINV**2 - 3.D0/2.D0*t1**(-1)*ro + 13.D0/16.D0*
     &    t1**(-1)*ro**2*GRAMINV - 23.D0/64.D0*t1**(-1)*ro**2*
     &    GRAMINV**2 + 1.D0/2.D0*t1**(-1)*ro**2 - 1.D0/8.D0*t1**(-1)*
     &    ro**3*GRAMINV + 19.D0/64.D0*t1**(-1)*ro**3*GRAMINV**2 - 1.D0/
     &    16.D0*t1**(-1)*ro**4*GRAMINV**2 + 1.D0/2.D0*t1**(-1)*GRAMINV
     &     + 7.D0/4.D0*t1**(-1) - 1.D0/2.D0*t1*ro*GRAMINV + 13.D0/16.D0
     &    *t1*ro*GRAMINV**2 - 5.D0/16.D0*t1*ro**2*GRAMINV**2 + 3.D0/4.D0
     &    *t1*GRAMINV - 1.D0/2.D0*t1*GRAMINV**2 + 1.D0/4.D0*t1**2*ro*
     &    GRAMINV**2 - 1.D0/4.D0*t1**2*GRAMINV**2 + 13.D0/8.D0*ro*
     &    GRAMINV**2 - 1.D0/2.D0*ro - 3.D0/2.D0*ro**2*GRAMINV**2 + 3.D0/
     &    8.D0*ro**3*GRAMINV**2 - 1.D0/2.D0*GRAMINV - 1.D0/2.D0*
     &    GRAMINV**2 )

      AAA1 =  + vlsm*xn**(-1) * ( 4.D0*t1**(-1)*ro*GRAMINV + t1**(-1)*
     &    ro*GRAMINV**2 + t1**(-1)*ro**2*GRAMINV**2 + 4.D0*t1**(-1)*
     &    GRAMINV + 4.D0*t1*GRAMINV**2 + 4.D0*t1**2*GRAMINV**2 + 4.D0*
     &    ro*GRAMINV - 9.D0*ro*GRAMINV**2 + 2.D0*ro**2*GRAMINV**2 - 12.D
     &    0*GRAMINV )
      AAA1 = AAA1 + vlsm*xn * (  - 4.D0*t1**(-1)*ro*GRAMINV + t1**(-1)*
     &    ro*GRAMINV**2 - 3.D0/4.D0*t1**(-1)*ro**2*GRAMINV**2 + 4.D0*
     &    t1**(-1)*GRAMINV - 4.D0*t1**(-1) - ro*GRAMINV + 4.D0*ro*
     &    GRAMINV**2 - 1.D0/2.D0*ro**2*GRAMINV**2 + 8.D0*GRAMINV - 4.D0
     &    *GRAMINV**2 )
      AAA1 = AAA1 + vltm*TBAR**(-1)*xn**(-1) * (  - 8.D0 - 2.D0*
     &    t1**(-1)*ro*GRAMINV + 2.D0*t1**(-1)*ro + 1.D0/2.D0*t1**(-1)*
     &    ro**2*GRAMINV - 8.D0*t1**(-1) - 4.D0*ro*GRAMINV + 8.D0*
     &    GRAMINV )
      AAA1 = AAA1 + vltm*TBAR**(-1)*xn * ( 4.D0 + t1**(-1)*ro*GRAMINV
     &     - t1**(-1)*ro - 1.D0/4.D0*t1**(-1)*ro**2*GRAMINV + 4.D0*
     &    t1**(-1) + 2.D0*ro*GRAMINV - 4.D0*GRAMINV )
      AAA1 = AAA1 + vlwm*UBAR**(-1)*xn**(-1) * (  - 8.D0 - 2.D0*
     &    t1**(-1)*ro - 1.D0/2.D0*t1**(-1)*ro**2*GRAMINV - 4.D0*t1*
     &    GRAMINV + 4.D0*GRAMINV )
      AAA1 = AAA1 + xn**(-1)*brackt * (  - 2.D0*t1**(-1)*ro*GRAMINV + 2.
     &    D0*t1**(-1)*ro*GRAMINV**2 - 1.D0/2.D0*t1**(-1)*ro**2*
     &    GRAMINV**2 + 8.D0*t1**(-1)*GRAMINV + 4.D0*ro*GRAMINV**2 + 8.D0
     &    *GRAMINV - 8.D0*GRAMINV**2 )
      AAA1 = AAA1 + xn**(-1)*brackw * ( 2.D0*t1**(-1)*ro*GRAMINV - 
     &    t1**(-1)*ro*GRAMINV**2 + 1.D0/2.D0*t1**(-1)*ro**2*GRAMINV**2
     &     - 4.D0*t1**(-1)*GRAMINV + 4.D0*t1*GRAMINV**2 - 4.D0*t1**2*
     &    GRAMINV**2 - ro*GRAMINV**2 )
      AAA1 = AAA1 + xn*brackt * ( t1**(-1)*ro*GRAMINV - t1**(-1)*ro*
     &    GRAMINV**2 + 1.D0/4.D0*t1**(-1)*ro**2*GRAMINV**2 - 4.D0*
     &    t1**(-1)*GRAMINV - 2.D0*ro*GRAMINV**2 - 4.D0*GRAMINV + 4.D0*
     &    GRAMINV**2 )
      AAA1 = AAA1 + f2*xn**(-1) * ( t1**(-1)*ro*GRAMINV**2 - 4.D0*
     &    t1**(-1)*ro**2*GRAMINV - t1**(-1)*ro**3*GRAMINV**2 + 4.D0*
     &    t1**(-1)*GRAMINV - 4.D0*t1*ro*GRAMINV**2 + 4.D0*t1*GRAMINV**2
     &     - 4.D0*t1**2*ro*GRAMINV**2 + 4.D0*t1**2*GRAMINV**2 + 12.D0*
     &    ro*GRAMINV - 9.D0*ro*GRAMINV**2 - 4.D0*ro**2*GRAMINV + 11.D0*
     &    ro**2*GRAMINV**2 - 2.D0*ro**3*GRAMINV**2 - 8.D0*GRAMINV )
      AAA1 = AAA1 + f2*xn * (  - 7.D0*t1**(-1)*ro*GRAMINV + t1**(-1)*ro
     &    *GRAMINV**2 + 2.D0*t1**(-1)*ro + 7.D0/2.D0*t1**(-1)*ro**2*
     &    GRAMINV - 7.D0/4.D0*t1**(-1)*ro**2*GRAMINV**2 + 3.D0/4.D0*
     &    t1**(-1)*ro**3*GRAMINV**2 + 4.D0*t1**(-1)*GRAMINV - 6.D0*ro*
     &    GRAMINV + 8.D0*ro*GRAMINV**2 + ro**2*GRAMINV - 9.D0/2.D0*
     &    ro**2*GRAMINV**2 + 1.D0/2.D0*ro**3*GRAMINV**2 + 4.D0*GRAMINV
     &     - 4.D0*GRAMINV**2 )

      AAA2 =  + vlsm*xn**(-1) * ( 2.D0*t1**(-1)*ro*GRAMINV - 5.D0/4.D0*
     &    t1**(-1)*ro*GRAMINV**2 + 1.D0/2.D0*t1**(-1)*ro**2*GRAMINV**2
     &     - 5.D0*t1**(-1)*GRAMINV - 2.D0*t1*GRAMINV**2 - t1**2*
     &    GRAMINV**2 - 5.D0/4.D0*ro*GRAMINV**2 - GRAMINV + 4.D0*
     &    GRAMINV**2 )
      AAA2 = AAA2 + vlsm*xn * (  - 1.D0/2.D0*t1**(-1)*ro*GRAMINV + 1.D0/
     &    4.D0*t1**(-1)*ro*GRAMINV**2 - 1.D0/8.D0*t1**(-1)*ro**2*
     &    GRAMINV**2 + t1**(-1)*GRAMINV + 1.D0/2.D0*t1*GRAMINV**2 + 1.D0
     &    /2.D0*ro*GRAMINV**2 + GRAMINV - GRAMINV**2 )
      AAA2 = AAA2 + vltm*TBAR**(-1)*xn**(-1) * (  - 1.D0/4.D0*t1**(-1)*
     &    ro*GRAMINV - t1**(-1) - 1.D0/2.D0*ro*GRAMINV + GRAMINV )
      AAA2 = AAA2 + vltm*TBAR**(-1)*xn * ( 1.D0/8.D0*t1**(-1)*ro*
     &    GRAMINV + 1.D0/2.D0*t1**(-1) + 1.D0/4.D0*ro*GRAMINV - 1.D0/2.D
     &    0*GRAMINV )
      AAA2 = AAA2 + vlwm*UBAR**(-1)*xn**(-1) * ( t1*GRAMINV + 1.D0/2.D0
     &    *ro*GRAMINV - GRAMINV )
      AAA2 = AAA2 + xn**(-1)*brackt * (  - t1*GRAMINV**2 + 1.D0/2.D0*ro
     &    *GRAMINV**2 + GRAMINV )
      AAA2 = AAA2 + xn**(-1)*brackw * ( 1.D0/4.D0*t1**(-1)*ro*
     &    GRAMINV**2 + t1**(-1)*GRAMINV - t1*GRAMINV**2 + t1**2*
     &    GRAMINV**2 - 1.D0/4.D0*ro*GRAMINV**2 )
      AAA2 = AAA2 + xn*brackt * ( 1.D0/2.D0*t1*GRAMINV**2 - 1.D0/4.D0*
     &    ro*GRAMINV**2 - 1.D0/2.D0*GRAMINV )
      AAA2 = AAA2 + f2*xn**(-1) * ( 7.D0*t1**(-1)*ro*GRAMINV - 5.D0/4.D0
     &    *t1**(-1)*ro*GRAMINV**2 - 2.D0*t1**(-1)*ro**2*GRAMINV + 7.D0/
     &    4.D0*t1**(-1)*ro**2*GRAMINV**2 - 1.D0/2.D0*t1**(-1)*ro**3*
     &    GRAMINV**2 - 5.D0*t1**(-1)*GRAMINV + 2.D0*t1*ro*GRAMINV**2 - 
     &    2.D0*t1*GRAMINV**2 + t1**2*ro*GRAMINV**2 - t1**2*GRAMINV**2
     &     + ro*GRAMINV - 21.D0/4.D0*ro*GRAMINV**2 + 5.D0/4.D0*ro**2*
     &    GRAMINV**2 - GRAMINV + 4.D0*GRAMINV**2 )
      AAA2 = AAA2 + f2*xn * (  - 3.D0/2.D0*t1**(-1)*ro*GRAMINV + 1.D0/4.
     &    D0*t1**(-1)*ro*GRAMINV**2 + 1.D0/2.D0*t1**(-1)*ro**2*GRAMINV
     &     - 3.D0/8.D0*t1**(-1)*ro**2*GRAMINV**2 + 1.D0/8.D0*t1**(-1)*
     &    ro**3*GRAMINV**2 + t1**(-1)*GRAMINV - 1.D0/2.D0*t1*ro*
     &    GRAMINV**2 + 1.D0/2.D0*t1*GRAMINV**2 - 1.D0/2.D0*ro*GRAMINV
     &     + 3.D0/2.D0*ro*GRAMINV**2 - 1.D0/2.D0*ro**2*GRAMINV**2 + 1.D0
     &    /2.D0*GRAMINV - GRAMINV**2 )

      AAA3 =  + vlsm*xn**(-1) * ( 1.D0/2.D0*ro*GRAMINV - 1.D0/2.D0*
     &    GRAMINV )
      AAA3 = AAA3 + vlsm*xn * (  - 1.D0/16.D0*t1**(-1)*ro*GRAMINV - 1.D0
     &    /4.D0*t1**(-1) - 1.D0/8.D0*ro*GRAMINV + 1.D0/4.D0*GRAMINV )
      AAA3 = AAA3 + xn**(-1)*brackt * (  - 1.D0/8.D0*t1**(-1)*ro*
     &    GRAMINV - 1.D0/2.D0*t1**(-1) + 1.D0/2.D0*GRAMINV )
      AAA3 = AAA3 + xn**(-1)*brackw * ( 1.D0/8.D0*t1**(-1)*ro*GRAMINV
     &     + 1.D0/2.D0*t1**(-1) )
      AAA3 = AAA3 + xn*brackt * ( 1.D0/16.D0*t1**(-1)*ro*GRAMINV + 1.D0/
     &    4.D0*t1**(-1) - 1.D0/4.D0*GRAMINV )
      AAA3 = AAA3 + f2*xn**(-1) * ( ro*GRAMINV - 1.D0/2.D0*ro**2*
     &    GRAMINV - 1.D0/2.D0*GRAMINV )
      AAA3 = AAA3 + f2*xn * (  - 1.D0/16.D0*t1**(-1)*ro*GRAMINV + 1.D0/
     &    4.D0*t1**(-1)*ro + 1.D0/16.D0*t1**(-1)*ro**2*GRAMINV - 1.D0/4.
     &    D0*t1**(-1) - 3.D0/8.D0*ro*GRAMINV + 1.D0/8.D0*ro**2*GRAMINV
     &     + 1.D0/4.D0*GRAMINV )

      AAA4 =  + b**(-1)*vlpm*xn**(-1) * (  - 2.D0 )
      AAA4 = AAA4 + vlsm*xn**(-1) * (  - 6.D0*t1**(-1)*ro*GRAMINV + 3.D0
     &    *t1**(-1)*ro*GRAMINV**2 - t1**(-1)*ro**2*GRAMINV**2 + 12.D0*
     &    t1**(-1)*GRAMINV - 8.D0*t1**(-1) + 4.D0*t1**2*GRAMINV**2 + 3.D
     &    0*ro*GRAMINV**2 + 8.D0*GRAMINV - 8.D0*GRAMINV**2 )
      AAA4 = AAA4 + vlsm*xn * ( 3.D0/2.D0*t1**(-1)*ro*GRAMINV - 1.D0/2.D
     &    0*t1**(-1)*ro*GRAMINV**2 + 1.D0/4.D0*t1**(-1)*ro**2*
     &    GRAMINV**2 - 2.D0*t1**(-1)*GRAMINV + 2.D0*t1**(-1) - 3.D0/2.D0
     &    *ro*GRAMINV**2 - 4.D0*GRAMINV + 2.D0*GRAMINV**2 )
      AAA4 = AAA4 + vltm*TBAR**(-1)*xn**(-1) * ( t1**(-1)*ro*GRAMINV + 
     &    4.D0*t1**(-1) + ro*GRAMINV - 4.D0*GRAMINV )
      AAA4 = AAA4 + vltm*TBAR**(-1)*xn * (  - 1.D0/2.D0*t1**(-1)*ro*
     &    GRAMINV - 2.D0*t1**(-1) - 1.D0/2.D0*ro*GRAMINV + 2.D0*GRAMINV
     &     )
      AAA4 = AAA4 + vlwm*UBAR**(-1)*xn**(-1) * (  - 4.D0*t1*GRAMINV - 
     &    ro*GRAMINV + 4.D0*GRAMINV )
      AAA4 = AAA4 + xn**(-1)*brackt * (  - t1**(-1)*ro*GRAMINV**2 - 4.D0
     &    *t1**(-1)*GRAMINV - ro*GRAMINV**2 - 4.D0*GRAMINV + 4.D0*
     &    GRAMINV**2 )
      AAA4 = AAA4 + xn**(-1)*brackw * ( 8.D0*t1*GRAMINV**2 - 4.D0*t1**2
     &    *GRAMINV**2 - 4.D0*GRAMINV**2 )
      AAA4 = AAA4 + xn*brackt * ( 1.D0/2.D0*t1**(-1)*ro*GRAMINV**2 + 2.D
     &    0*t1**(-1)*GRAMINV + 1.D0/2.D0*ro*GRAMINV**2 + 2.D0*GRAMINV
     &     - 2.D0*GRAMINV**2 )
      AAA4 = AAA4 + f2*xn**(-1) * (  - 16.D0*t1**(-1)*ro*GRAMINV + 3.D0
     &    *t1**(-1)*ro*GRAMINV**2 + 8.D0*t1**(-1)*ro + 6.D0*t1**(-1)*
     &    ro**2*GRAMINV - 4.D0*t1**(-1)*ro**2*GRAMINV**2 + t1**(-1)*
     &    ro**3*GRAMINV**2 + 12.D0*t1**(-1)*GRAMINV - 4.D0*t1**2*ro*
     &    GRAMINV**2 + 4.D0*t1**2*GRAMINV**2 - 8.D0*ro*GRAMINV + 11.D0*
     &    ro*GRAMINV**2 - 3.D0*ro**2*GRAMINV**2 + 4.D0*GRAMINV - 8.D0*
     &    GRAMINV**2 )
      AAA4 = AAA4 + f2*xn * (  - 6.D0 + 3.D0*t1**(-1)*ro*GRAMINV - 1.D0/
     &    2.D0*t1**(-1)*ro*GRAMINV**2 - 2.D0*t1**(-1)*ro - 3.D0/2.D0*
     &    t1**(-1)*ro**2*GRAMINV + 3.D0/4.D0*t1**(-1)*ro**2*GRAMINV**2
     &     - 1.D0/4.D0*t1**(-1)*ro**3*GRAMINV**2 - 2.D0*t1**(-1)*
     &    GRAMINV + 3.D0*ro*GRAMINV - 7.D0/2.D0*ro*GRAMINV**2 + 3.D0/2.D
     &    0*ro**2*GRAMINV**2 - 2.D0*GRAMINV + 2.D0*GRAMINV**2 )
      AAA4 = AAA4 + f3*xn * ( 2.D0 )

      AAA5 =  + vlsm*xn**(-1) * (  - 1.D0/2.D0*t1**(-1)*ro*GRAMINV**2
     &     - 2.D0*t1**(-1)*GRAMINV - 2.D0*t1**2*GRAMINV**2 + 3.D0/2.D0*
     &    ro*GRAMINV**2 + 2.D0*GRAMINV )
      AAA5 = AAA5 + vlsm*xn * ( t1*GRAMINV**2 - 1.D0/2.D0*ro*GRAMINV**2
     &     - GRAMINV )
      AAA5 = AAA5 + vltm*TBAR**(-1)*xn**(-1) * ( 1.D0/2.D0*t1**(-1)*ro*
     &    GRAMINV + 2.D0*t1**(-1) - 2.D0*GRAMINV )
      AAA5 = AAA5 + vltm*TBAR**(-1)*xn * (  - 1.D0/4.D0*t1**(-1)*ro*
     &    GRAMINV - t1**(-1) + GRAMINV )
      AAA5 = AAA5 + vlwm*UBAR**(-1)*xn**(-1) * ( 2.D0*t1*GRAMINV - 2.D0
     &    *GRAMINV )
      AAA5 = AAA5 + xn**(-1)*brackt * (  - t1**(-1)*ro*GRAMINV**2 - 4.D0
     &    *t1**(-1)*GRAMINV - 2.D0*t1*GRAMINV**2 + 4.D0*GRAMINV**2 )
      AAA5 = AAA5 + xn**(-1)*brackw * (  - 1.D0/2.D0*t1**(-1)*ro*
     &    GRAMINV**2 - 2.D0*t1**(-1)*GRAMINV - 6.D0*t1*GRAMINV**2 + 2.D0
     &    *t1**2*GRAMINV**2 + 1.D0/2.D0*ro*GRAMINV**2 + 2.D0*GRAMINV + 
     &    4.D0*GRAMINV**2 )
      AAA5 = AAA5 + xn*brackt * ( 1.D0/2.D0*t1**(-1)*ro*GRAMINV**2 + 2.D
     &    0*t1**(-1)*GRAMINV + t1*GRAMINV**2 - 2.D0*GRAMINV**2 )
      AAA5 = AAA5 + f2*xn**(-1) * ( 2.D0*t1**(-1)*ro*GRAMINV - 1.D0/2.D0
     &    *t1**(-1)*ro*GRAMINV**2 + 1.D0/2.D0*t1**(-1)*ro**2*GRAMINV**2
     &     - 2.D0*t1**(-1)*GRAMINV + 2.D0*t1**2*ro*GRAMINV**2 - 2.D0*
     &    t1**2*GRAMINV**2 - 2.D0*ro*GRAMINV + 3.D0/2.D0*ro*GRAMINV**2
     &     - 3.D0/2.D0*ro**2*GRAMINV**2 + 2.D0*GRAMINV )
      AAA5 = AAA5 + f2*xn * (  - 1.D0/2.D0*t1**(-1)*ro*GRAMINV - 2.D0*
     &    t1**(-1) - t1*ro*GRAMINV**2 + t1*GRAMINV**2 + ro*GRAMINV - 1.D
     &    0/2.D0*ro*GRAMINV**2 + 1.D0/2.D0*ro**2*GRAMINV**2 )

      AAA6 =  + vlsm*xn**(-1) * ( 1.D0/4.D0*t1**(-1)*ro*GRAMINV**2 + 
     &    t1**(-1)*GRAMINV + t1**2*GRAMINV**2 - 3.D0/4.D0*ro*GRAMINV**2
     &     - GRAMINV )
      AAA6 = AAA6 + vlsm*xn * (  - 1.D0/2.D0*t1*GRAMINV**2 + 1.D0/4.D0*
     &    ro*GRAMINV**2 + 1.D0/2.D0*GRAMINV )
      AAA6 = AAA6 + vltm*TBAR**(-1)*xn**(-1) * (  - 1.D0/4.D0*t1**(-1)*
     &    ro*GRAMINV - t1**(-1) + GRAMINV )
      AAA6 = AAA6 + vltm*TBAR**(-1)*xn * ( 1.D0/8.D0*t1**(-1)*ro*
     &    GRAMINV + 1.D0/2.D0*t1**(-1) - 1.D0/2.D0*GRAMINV )
      AAA6 = AAA6 + vlwm*UBAR**(-1)*xn**(-1) * (  - t1*GRAMINV + 
     &    GRAMINV )
      AAA6 = AAA6 + xn**(-1)*brackt * ( 1.D0/2.D0*t1**(-1)*ro*
     &    GRAMINV**2 + 2.D0*t1**(-1)*GRAMINV + t1*GRAMINV**2 - 2.D0*
     &    GRAMINV**2 )
      AAA6 = AAA6 + xn**(-1)*brackw * ( 1.D0/4.D0*t1**(-1)*ro*
     &    GRAMINV**2 + t1**(-1)*GRAMINV + 3.D0*t1*GRAMINV**2 - t1**2*
     &    GRAMINV**2 - 1.D0/4.D0*ro*GRAMINV**2 - GRAMINV - 2.D0*
     &    GRAMINV**2 )
      AAA6 = AAA6 + xn*brackt * (  - 1.D0/4.D0*t1**(-1)*ro*GRAMINV**2
     &     - t1**(-1)*GRAMINV - 1.D0/2.D0*t1*GRAMINV**2 + GRAMINV**2 )
      AAA6 = AAA6 + f2*xn**(-1) * (  - t1**(-1)*ro*GRAMINV + 1.D0/4.D0*
     &    t1**(-1)*ro*GRAMINV**2 - 1.D0/4.D0*t1**(-1)*ro**2*GRAMINV**2
     &     + t1**(-1)*GRAMINV - t1**2*ro*GRAMINV**2 + t1**2*GRAMINV**2
     &     + ro*GRAMINV - 3.D0/4.D0*ro*GRAMINV**2 + 3.D0/4.D0*ro**2*
     &    GRAMINV**2 - GRAMINV )
      AAA6 = AAA6 + f2*xn * ( 1.D0/4.D0*t1**(-1)*ro*GRAMINV + t1**(-1)
     &     + 1.D0/2.D0*t1*ro*GRAMINV**2 - 1.D0/2.D0*t1*GRAMINV**2 - 1.D0
     &    /2.D0*ro*GRAMINV + 1.D0/4.D0*ro*GRAMINV**2 - 1.D0/4.D0*ro**2*
     &    GRAMINV**2 )

      call qbasis(t1,ro,NLO0,NLO1,NLO2,NLO3,NLO4,NLO5,NLO6)
C This is a place-holder for LOSQ in n-dimensions
      EPINV=0d0
      if (kswap .eq. 1) then
      LOSQ=1d0
      wtqqb=fac*(AAAE*LOSQ+AAA0*NLO0
     . +AAA1*NLO1+AAA2*NLO2+AAA3*NLO3
     . +AAA4*NLO4+AAA5*NLO5+AAA6*NLO6)
      elseif (kswap .eq. 2) then
      LOSQ=1d0
      wtqbq=fac*(AAAE*LOSQ+AAA0*NLO0
     . +AAA1*NLO1+AAA2*NLO2+AAA3*NLO3
     . +AAA4*NLO4+AAA5*NLO5+AAA6*NLO6)
      endif
      enddo
      return
      end
