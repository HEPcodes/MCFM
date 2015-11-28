      subroutine qqb_QQb_v(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2002.                                                     *
*     Virtual corrections                                              *
*     Calculates the element squared for the process                   *
*                                                                      *
*     q(-p1)+qbar(-p2) -> t(p3)+t~(p4)                                 *
*                                                                      * 
************************************************************************
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'scale.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),t1,ro,mass
      double precision qqsym,qqasy,ggsym,qqsdel,ggdel,qqanas,ggnas
      double precision ss,xm2,xlf,xmu,rmuom2,qqadel,qqsnas
      integer naem,nbeam1,nbeam2
      logical nascheck
      common/nascheck/nascheck
      common/para1/ss,xm2,xlf,xmu,rmuom2,naem,nbeam1,nbeam2

      ss=s(1,2)
      mass=mt
      xm2=mass**2

      xlf=dfloat(nf)
      xmu=scale
      rmuom2=dlog(xmu**2/xm2)
      naem=0
      ro=4d0*xm2/s(1,2)
      t1=-s(1,3)/s(1,2)                                                         

      nascheck=.false.

c--- Performs comparison with Nason if desired
      if (nascheck) then            
        qqsnas=qqsdel(t1,ro)  
        qqanas=qqadel(t1,ro)  
        ggnas=V*ggdel(t1,ro)  
        write(6,*) 'Results from NDE'
        write(6,*) 'qqsnas',qqsnas
        write(6,*) 'qqanas',qqanas
        write(6,*) 'qqsumnas',qqsnas+qqanas
        write(6,*) 'qqdifnas',qqsnas-qqanas
        write(6,*) 'ggnas',ggnas
      endif

      call dotem(4,p,s)
      call virteval(t1,ro,qqsym,qqasy,ggsym)  
      
      if (nascheck) then
        write(6,*) 'New results'
        write(6,*) 'qqsym',qqsym/gsq**3/8d0
        write(6,*) 'qqasy',qqasy/gsq**3/8d0
        write(6,*) 'qqsum',(qqsym+qqasy)/gsq**3/8d0
        write(6,*) 'qqdif',(qqsym-qqasy)/gsq**3/8d0
        write(6,*) 'ggsym',ggsym
        pause
      endif
      
      do j=-nf,nf
      k=-j
      if     (j .eq. 0) then
      msq(0,0)=avegg*gsq**3*ggsym/8d0/pi**2
      elseif (j .gt. 0) then
      msq(j,k)=aveqq*(qqsym+qqasy)/8d0/pi**2
      elseif (j .lt. 0) then
      msq(j,k)=aveqq*(qqsym-qqasy)/8d0/pi**2
      endif
      enddo

      return
      end


      subroutine virteval(t1,ro,qqsym,qqasy,ggsym)  
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'scale.f'
      external ddilog
      double precision t1,t2,ro,tbar,ubar,b,xlp,xlm,vlpm,vlsm,vltm,vlwm
      double precision vlbl,vdmp,vdmb,f1,f2,f3,f4,f5t1,f5t2,
     . qqQQv_0,qqQQv_1,ddilog,mass
      double precision vdt,vdw,xlf,rmuom2
      double precision qqss,qqsv,qqas,qqav
      double precision ggs
      double precision ggv1,ggv2,ggv3,ggv4,ggv5,ggv6,
     . ggv7,ggv8,ggv9,ggv10,ggv11
      double precision ggQQv,ggQQov,ggQQsv,qqsym,qqasy,ggsym,
     . ggQQv_0,ggQQov_0,ggQQsv_0,ggQQv_1,ggQQov_1,ggQQsv_1,
     . ggQQv_2,ggQQov_2,ggQQsv_2
      logical nascheck
      common/nascheck/nascheck
      
C---we arrive at the dred scheme by taking a HV result and applying
C---a finite renormalization on the 
C---namely a factor of as/4/pi*cf/2 for every massles quark leg (amplitude)
C---namely a factor of as/4/pi*xn/6 for every massles gluon leg (amplitude)
C---namely a factor of as/2/pi*cf/2 for every massles quark leg (square)
C---namely a factor of as/2/pi*xn/6 for every massles gluon leg (square)
      scheme='dred'
     
C Notation for logarithms and dilogarithms.
C Integrals were in unphysical region, but contination has now 
C been performed; 
C In unphysical region,        In physical region after continuation.
C vltm=log(at/m^2),
C vlpm=log(-lp/lm),            log((1+b)/(1-b))
C vlsm=log(as/m^2),            log(s/m^2)
C vlwm=log(aw/m^2),
C vlbl=log(-b/lm),             log(2*b/(1-b))

C vdw=li[2]((aw-m^2)/aw)-1/2*vlwm^2,
C vdt=li[2]((at-m^2)/at)-1/2*vltm^2,
C vdmp=li[2](-lm/lp),
C vdmb=li[2](-lm/b)+1/2*vlbl^2,

C as=-s
C lp=(1+b)/2,
C lm=(1-b)/2,
C at=t1*s=m^2-t,
C aw=t2*s=m^2-u,
C b=sqrt(1-ro),
C m=sqrt(ro*s)/2
C TBAR=-T/S=t1-1/4*ro
C UBAR=-U/S=t2-1/4*ro

      mass=mt
      xlf=dfloat(nf)
      rmuom2=2d0*dlog(scale/mass)
      t2=1d0-t1
      tbar=t1-0.25d0*ro
      ubar=t2-0.25d0*ro
      b=dsqrt(1d0-ro)   
      xlp=0.5d0*(1d0+b)
      xlm=0.5d0*(1d0-b)
      vlpm=dlog(xlp/xlm)     
      vlsm=dlog(4d0/ro)      
      vltm=dlog(4d0*t1/ro)      
      vlwm=dlog(4d0*t2/ro)      
      vlbl=dlog(b/xlm)       
      vdw=ddilog(1d0-ro/(4d0*t2))-0.5d0*vlwm**2      
      vdt=ddilog(1d0-ro/(4d0*t1))-0.5d0*vltm**2
      vdmp=ddilog(-xlm/xlp)    
      vdmb=ddilog(-xlm/b)+0.5d0*vlbl**2  

c--- Q-Qbar and Qbar-Q contributions

      f1=(vlpm**2/2d0-2d0*vdmb-pisq/3d0)/b
      f2=(-b*vlsm+vlpm**2/4d0+vdmp+pisq/12d0)/b**3
      f3=(-b**3*vlsm-3*b*vlsm+0.75d0*vlpm**2
     . +3d0*vdmp+pisq/4d0+2d0*b**3)/b**5
      f4 = (vlpm**2/4d0+vdmp+pi**2/12d0)/b
      f5t1 = (vltm**2+vdt+pi**2/6d0)/t1**3
      f5t2 = (vlwm**2+vdw+pi**2/6d0)/t2**3
      
c---  Singular parts
c---  qqQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps)
      qqQQv_0=gsq**2*V*(2d0*t1**2+2d0*t2**2+ro)
      if (scheme .eq. 'tH-V') then
      qqQQv_1=gsq**2*V*(-2d0)
      else
      qqQQv_1=0d0
      endif      
      
c--- These are the singular pieces, written in such a way that
c--- the limit EPINV -> 0 is smooth
      qqss=
     .  gsq/V*qqQQv_0*( 
     . -2d0*V/2/xn*EPINV**2-V/xn*EPINV
     . -3d0*V/2/xn*EPINV+xn*(vltm+vlwm)*EPINV
     . -(0.5d0*vlpm/b*(1d0+b**2)+vlsm)/xn*EPINV)
     . +gsq/V*qqQQv_1*(
     . -2d0*V/2/xn*EPINV-V/xn
     . -3d0*V/2/xn+xn*(vltm+vlwm)
     . -(0.5d0*vlpm/b*(1d0+b**2)+vlsm)/xn)
     
      qqas=gsq/V*qqQQv_0*(XN4/xn*(vltm-vlwm))*EPINV
     .    +gsq/V*qqQQv_1*(XN4/xn*(vltm-vlwm))

      if (nascheck) then
        qqss=0d0
        qqas=0d0
      endif

C---  Finite symmetric part
       qqsv=
     . +gsq/V*qqQQv_0*(xn*(22d0/9d0-0.5d0*vlsm*(vltm+vlwm)
     . +0.25d0*vlsm**2+11d0/3d0*rmuom2+11d0/6d0*vlsm-pisq/6d0)
     . +1d0/xn*(+6d0-pisq/3d0-1.5d0*vlsm+0.5d0*vlsm**2-b*vlpm
     . -0.5d0*(1d0+b**2)*(pisq/b+f1))
     . +TR*XLF*(-20d0/9d0+4d0/3d0*vlsm-4d0/3d0*rmuom2)
     . +TR*(-20d0/9d0+4d0/3d0*vlpm*b*(1d0+ro/2d0)-4d0/3d0*ro))
     . +gsq**3*0.5d0*xn*(t1-t2)**2*f3
     . +gsq**3*f2*(xn*(-5d0*(t1**2+t2**2)+2d0
     . +b**2*(0.5d0+6d0*t1**2+6d0*t2**2)-b**4))
     . -gsq**3*0.5d0/xn*vlpm/b*((t1-t2)**2+b**2)
     . +gsq**3*xn*((t2-t1)*(vdt-vdw)+0.5d0*ro*(vdt+vdw)
     . +0.5d0*vlsm**2*(t1**2 +t2**2)+1d0-pisq/12d0*ro
     . -vlsm*(vltm+vlwm)*(t1**2+t2**2)-vlsm*(vltm-vlwm)*(t1-t2)
     . -(vltm/TBAR+vlwm/UBAR)*(0.5d0*ro-t1*t2)
     . -vlsm*(1.5d0+2d0*ro))

c--- extra finite terms in DR scheme
      if (scheme .eq. 'dred') then
        qqsv=qqsv+0.5d0*(xn-1d0/xn)*gsq*qqQQv_0/V
      endif
        
C--- Finite antisymmetric part
      qqav=gsq/V*qqQQv_0*(-XN4/2d0/xn*vlsm*(vltm-vlwm))
     . +gsq**3*XN4/xn*(-(t1-t2)*(vdt+vdw)+0.5d0*ro*(vdt-vdw)
     . +(pisq/6d0+vlsm*(3d0-ro)+0.5d0*vlsm**2-vlsm*(vltm+vlwm))*(t1-t2)
     . -vlsm*(vltm-vlwm)*(t1**2+t2**2)
     . -(vltm/tbar-vlwm/ubar)*(0.5d0*ro-t1*t2)
     . +f2*(t1-t2)*(-1d0+2d0*b**2+b**4))

      qqsym=V*(qqss+qqsv)
      qqasy=V*(qqas+qqav)

c--- g-g contribution

c---  Singular parts
c---  ggQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps) and _2 for O(eps**2)
      ggQQv_0=2d0/xn*(V/t1/t2-2d0*xn**2)
     . *((t1**2+t2**2)+ro-ro**2/4d0/t1/t2)
      ggQQov_0=2d0/xn*((xn**2-2d0)/t1/t2-2d0*xn**2)
     . *((t1**2+t2**2)+ro-ro**2/4d0/t1/t2)
      ggQQsv_0=2d0/xn*(1d0/t1/t2)
     . *((t1**2+t2**2)+ro-ro**2/4d0/t1/t2)
      if (scheme .eq. 'tH-V') then
      ggQQv_1=2d0/xn*(V/t1/t2-2d0*xn**2)
     . *(-(t1**2+t2**2)-1d0)
      ggQQv_2=2d0/xn*(V/t1/t2-2d0*xn**2)
      ggQQov_1=2d0/xn*((xn**2-2d0)/t1/t2-2d0*xn**2)
     . *(-(t1**2+t2**2)-1d0)
      ggQQov_2=2d0/xn*((xn**2-2d0)/t1/t2-2d0*xn**2)
      ggQQsv_1=2d0/xn*(1d0/t1/t2)
     . *(-(t1**2+t2**2)-1d0)
      ggQQsv_2=2d0/xn*(1d0/t1/t2)
      else
      ggQQv_1=0d0
      ggQQv_2=0d0
      ggQQov_1=0d0
      ggQQov_2=0d0
      ggQQsv_1=0d0
      ggQQsv_2=0d0
      endif      

c--- These are the singular pieces, written in such a way that
c--- the limit EPINV -> 0 is smooth
      ggs=
     .  (-2d0*xn)*(EPINV**2*ggQQv_0+EPINV*ggQQv_1+ggQQv_2)
     . +(2d0*(2d0*TR/3d0*XLF-11d0/6d0*xn)-V/xn)*(EPINV*ggQQv_0+ggQQv_1)
     . +(EPINV*ggQQsv_0+ggQQsv_1)*(1d0+b**2)/b*vlpm*V/2d0/xn
     . +(EPINV*ggQQov_0+ggQQov_1)*(vlsm*xn-(1d0+b**2)/b*vlpm/xn/2d0)
     . +2d0*(EPINV*ggQQsv_0+ggQQsv_1)*(vlsm-vltm-vlwm)*xn
     . +4d0*xn**2*EPINV*vltm
     . *(-ro**2/4d0/t1**2+t2/t1+ro*t2/t1-2d0*t2**2)
     . +4d0*xn**2*EPINV*vlwm
     . *(-2d0*t1**2-ro**2/4d0/t2**2+t1/t2+ro*t1/t2)
      if (scheme .eq. 'tH-V') then
c--- These are the extra finite pieces (akin to _1) for the last 4 lines
c--- and they appear to have already been added in the finite pieces
c--- below, so we subtract it later there
      ggs=ggs+(
     . +4d0*xn**2*vltm
     . *(-2d0*t2/t1+2d0*t2**2)
     . +4d0*xn**2*vlwm
     . *(-2d0*t1/t2+2d0*t1**2))
      endif
      
      if (nascheck) then
        ggs=0d0
      endif
      
c--- extra finite terms in DR scheme
      if (scheme .eq. 'dred') then
        ggs=ggs+xn/3d0*ggQQv_0
      endif
        
C--- overall factor of V removed
      ggQQv=2d0/xn*(V/t1/t2-2d0*xn**2) 
     & *(t1**2+t2**2+ro-ro**2/(4d0*t1*t2))
      ggQQov=2d0/xn*((xn**2-2d0)/(t1*t2)-2d0*xn**2)
     & *(t1**2+t2**2+ro-ro**2/(4d0*t1*t2))
      ggQQsv=2d0/(xn*t1*t2)
     & *(t1**2+t2**2+ro-ro**2/(4d0*t1*t2))

      ggv1 = f3*(-2*t2**2+1/(t1*t2)/4.0-2*t1**2)*xn**2+f4*(24*t1*t2-ro**
     1   2/(t1*t2)/2.0+ro/(t1*t2)+11.0/(4.0*t1*t2)-4*ro-17)*xn**2+f2*(-2
     2   0*t1*t2+(-3.0)/(2.0*t1*t2)+11)*xn**2-2*ggqqv*rmuom2*((-11.0)*xn
     3   /6.0+2.0*tr*xlf/3.0)+(b**2+1)*ggqqsv*pi**2*v/(b*xn)/2.0-(b**2+1
     4   )*ggqqov*pi**2/(b*xn)/2.0-f5t2*ro**2*v/xn**2/4.0-f5t1*ro**2*v/x
     5   n**2/4.0+f4*((1-ro/2.0)*ro**2*(t2**2+t1**2)/(t1**2*t2**2)+(-ro*
     6   *3+5*ro**2-ro-6)/(t1*t2)+4)/xn**2+f1*(-(1-ro/2.0)*ro**2*(t2**2+
     7   t1**2)/(t1**2*t2**2)/2.0+(ro**3/2.0-2*ro**2+ro+2)/(t1*t2)+2*ro-
     8   4)/xn**2+f1*(4*t2**2+4*ro*t1*t2-(1-ro/2.0)*ro**2/(t1*t2)+4*t1**
     9   2-2*ro**2+2*ro)+f4*(-ro**3/(t1*t2)+4*ro**2/(t1*t2)-4/(t1*t2)-2*
     :   ro**2-2*ro+6)
      ggv2 = (-ro**3/t2/2d0+2*ro**2/t2-ro/t2-2/t2-ro**3/t2**2/2d0+ro**2/
     1   t2**2-ro**3/t1/2d0+3*ro**2/t1-4/t1+2)*vlpm*vlwm/(b*xn**2)+(-ro*
     2   *3/t2/2d0+3*ro**2/t2-4/t2-ro**3/t1/2d0+2*ro**2/t1-ro/t1-2/t1-ro
     3   **3/t1**2/2d0+ro**2/t1**2+2)*vlpm*vltm/(b*xn**2)+(ro**3/t2/2d0+
     4   (-5d0)*ro**2/(2d0*t2)+ro/t2/2d0+3/t2+ro**3/t2**2/4d0-ro**2/t2**
     5   2/2d0+ro**3/t1/2d0+(-5d0)*ro**2/(2d0*t1)+ro/t1/2d0+3/t1+ro**3/t
     6   1**2/4d0-ro**2/t1**2/2d0-2)*vlpm*vlsm/(b*xn**2)+(ro/(t1*t2)-4)*
     7   vlpm/(b*xn**2)
      ggv3 = (ro**3/t2-ro**2/t2+ro**3/t1-ro**2/t1-4*ro**3+4*ro**2)*tr*vl
     1   pm*xn/b+(-ro**3/t2/2d0+ro**2/t2-2*t1-ro**3/t1/2d0+3*ro**2/t1-4/
     2   t1-ro**2-ro+4)*vlpm*vlwm/b+(-ro**3/t2/2d0+3*ro**2/t2-4/t2+2*t1-
     3   ro**3/t1/2d0+ro**2/t1-ro**2-ro+2)*vlpm*vltm/b+(ro**3/t2/2d0-2*r
     4   o**2/t2+2/t2+ro**3/t1/2d0-2*ro**2/t1+2/t1+ro**2+ro-3)*vlpm*vlsm
     5   /b+(ro**2/t2/2d0-ro/t2/4d0-8*ro*t1**2+12*t1**2+8*ro*t1-12*t1+ro
     6   **2/t1/2d0-ro/t1/4d0-2*ro**2+ro+1)*vlpm/b
      ggv4 = (5*t1**2-2*ro*t1-8*t1+7d0*ro**2/8d0+3*ro-1)*vltm*xn**2/tbar
     1   +(-ro*t1**2/4d0+(-3d0)*ro**2*t1/8d0+ro*t1/4d0+3d0*ro**3/32d0-ro
     2   **2/8d0)*vltm*xn**2/tbar**2+(8*t1**2-16*t1+2*ro/t1-4/t1-ro**2/t
     3   1**2+16)*vltm*xn**2+((-14d0)*ro/(3d0*t2)+(-7d0)/(2d0*t2)+ro**2/
     4   t2**2+16*t1**2-16*t1+(-14d0)*ro/(3d0*t1)+(-7d0)/(2d0*t1)+ro**2/
     5   t1**2+26d0*ro/3d0+15)*xn**2-2*ro**2/t2+17d0*ro/(2d0*t2)+8/t2-2*
     6   ro**2/t2**2+ro/t2**2-16*t1**2+16*t1-2*ro**2/t1+17d0*ro/(2d0*t1)
     7   +8/t1-2*ro**2/t1**2+ro/t1**2-10*ro-19
      ggv5 = (-2*ro/t2-2/t2+ro**2/t2**2+8*t1**2-8*t1+2*ro+6)*vlsm*vlwm*x
     1   n**2+(2*ro/t2-4/t2-ro**2/t2**2+8*t1**2+8)*vlwm*xn**2+(8*t1**2-8
     2   *t1-2*ro/t1-2/t1+ro**2/t1**2+2*ro+6)*vlsm*vltm*xn**2+(-ro/t2/2.
     3   0-1/t2/2d0-ro/t1/2d0-1/t1/2d0+ro+1)*vlsm**2*xn**2+(ro/t2/2d0-1/
     4   t2/2d0-8*t1**2+8*t1+ro/t1/2d0-1/t1/2d0-2*ro+2)*vlsm*xn**2
      ggv6 = (5*t1**2+2*ro*t1-2*t1+7d0*ro**2/8d0+ro-4)*vlwm*xn**2/ubar+(
     1   -ro*t1**2/4d0+3d0*ro**2*t1/8d0+ro*t1/4d0+3d0*ro**3/32d0-ro**2/2
     2   d0)*vlwm*xn**2/ubar**2+(2*ro/t2+2/t2-8*t1-2*ro+2)*vdw*xn**2+(8*
     3   t1+2*ro/t1+2/t1-2*ro-6)*vdt*xn**2+(ro/2d0-ro*t1/4d0)*xn**2/ubar
     4   +(ro*t1/4d0+ro/4d0)*xn**2/tbar+pi**2*(ro/t2/6d0+1/t2/6d0-ro**2/
     5   t2**2/12d0+(-4d0)*t1**2/3d0+4d0*t1/3d0+ro/t1/6d0+1/t1/6d0-ro**2
     6   /t1**2/12d0-ro/3d0-1)*xn**2+(-ro**2/t2+ro/t2-2/t2-ro**2/t1+3*ro
     7   /t1+4/t1-ro**2/t1**2+ro/t1**2)*vltm/xn**2+(2*ro**2/t2-5*ro/t2-5
     8   /t2+ro**2/t2**2-ro/t2**2+2*ro**2/t1-5*ro/t1-5/t1+ro**2/t1**2-ro
     9   /t1**2+6)/xn**2
      ggv7 = (t1+4/t1+ro-4)*vlwm/(ubar*xn**2)+(3d0*ro*t1/4d0-2/t1+ro**2/
     1   4d0-ro/4d0+2)*vlwm/(ubar**2*xn**2)+(-ro**2/t2+3*ro/t2+4/t2-ro**
     2   2/t2**2+ro/t2**2-ro**2/t1+ro/t1-2/t1)*vlwm/xn**2+((-3d0)*ro**2/
     3   (4d0*t2)+ro/t2+4/t2+(-3d0)*ro**2/(4d0*t1)-ro/t1+2/t1+ro**2/t1**
     4   2/4d0-2)*vltm**2/xn**2+(4/t2-t1+ro-3)*vltm/(tbar*xn**2)+(-2/t2+
     5   (-3d0)*ro*t1/4d0+ro**2/4d0+ro/2d0+2)*vltm/(tbar**2*xn**2)+(-ro*
     6   *2/t2/4d0+3d0/(2d0*t2)-ro**2/t1/4d0+3d0/(2d0*t1)-1)*vlpm**2/xn*
     7   *2
      ggv8 = ((-3d0)*ro**2/(4d0*t2)-ro/t2+2/t2+ro**2/t2**2/4d0+(-3d0)*ro
     1   **2/(4d0*t1)+ro/t1+4/t1-2)*vlwm**2/xn**2+((-3d0)*ro**2/(4d0*t2)
     2   -ro/t2+2/t2+ro**2/t2**2/4d0+(-3d0)*ro**2/(4d0*t1)+ro/t1+4/t1-2)
     3   *vdw/xn**2+((-3d0)*ro**2/(4d0*t2)+ro/t2+4/t2+(-3d0)*ro**2/(4d0*
     4   t1)-ro/t1+2/t1+ro**2/t1**2/4d0-2)*vdt/xn**2+(2/t1-ro/4d0-2)/(ub
     5   ar*xn**2)+(2/t2-ro/4d0-2)/(tbar*xn**2)+pi**2*(-1/t2/2d0+ro**2/t
     6   2**2/24d0-1/t1/2d0+ro**2/t1**2/24d0+1d0/3d0)/xn**2
      ggv9 = (-2*ro**2/t2+4*ro/t2+4/t2-2*ro**2/t2**2-2*ro**2/t1+4*ro/t1+
     1   4/t1-2*ro**2/t1**2-8)*vltm*vlwm+(-ro**2/t2+2*ro/t2+2/t2-2*t1+ro
     2   /t1-2/t1-2*ro+2)*vltm**2+(-12/t2-t1**2+3*ro*t1+9*t1+(-7d0)*ro**
     3   2/8d0-4*ro+12)*vltm/tbar+(2/t2+(-3d0)*ro*t1**2/4d0+5d0*ro**2*t1
     4   /8d0+5d0*ro*t1/2d0+(-3d0)*ro**3/32d0+(-5d0)*ro**2/8d0-ro/2d0-2)
     5   *vltm/tbar**2+(ro**2/t2+10/t2+ro**2/t1-4*ro/t1-8/t1+2*ro**2/t1*
     6   *2-ro/t1**2)*vltm+(-ro**2/t2/8d0+ro/t2/4d0+1/t2-ro**2/t1/8d0+ro
     7   /t1/4d0+1/t1-ro+(-3d0)/2d0)*vlpm**2
      ggv10 = (ro/t2-2/t2+2*t1-ro**2/t1+2*ro/t1+2/t1-2*ro)*vlwm**2+(-t1*
     1   *2-3*ro*t1-7*t1-12/t1+(-7d0)*ro**2/8d0-ro+20)*vlwm/ubar+((-3d0)
     2   *ro*t1**2/4d0+(-5d0)*ro**2*t1/8d0-ro*t1+2/t1+(-3d0)*ro**3/32d0+
     3   5d0*ro/4d0-2)*vlwm/ubar**2+(ro**2/t2-4*ro/t2-8/t2+2*ro**2/t2**2
     4   -ro/t2**2+ro**2/t1+10/t1)*vlwm+(ro**2/t2-ro/t2-4/t2+2*t1-2*ro+4
     5   )*vdw+(-2*t1+ro**2/t1-ro/t1-4/t1-2*ro+6)*vdt
      ggv11 = (ro/t2/3d0+ro/t1/3d0+(-4d0)*ro/3d0)*tr*xlf*xn+(ro**2/t2/4.
     1   0+ro**2/t1/4d0-ro**2)*tr*vlpm**2*xn+(2*ro**2/t2+ro/t2/3d0+2*ro*
     2   *2/t1+ro/t1/3d0-8*ro**2+(-4d0)*ro/3d0)*tr*xn+pi**2*(-ro**2/t2/4
     3   d0-ro**2/t1/4d0+ro**2)*tr*xn+((-3d0)*ro*t1/4d0-2/t1-ro**2/4d0+3
     4   d0*ro/4d0+2)/ubar+(-2/t2+3d0*ro*t1/4d0-ro**2/4d0+2)/tbar+pi**2*
     5   (-ro**2/t2/24d0+ro/t2/4d0-1/t2+ro**2/t2**2/3d0-ro**2/t1/24d0+ro
     6   /t1/4d0-1/t1+ro**2/t1**2/3d0+ro/3d0+11d0/6d0)

c--- This is subtracting the extra finite piece mentioned above
      ggv11=ggv11-(
     . +4d0*xn**2*vltm
     . *(-2d0*t2/t1+2d0*t2**2)
     . +4d0*xn**2*vlwm
     . *(-2d0*t1/t2+2d0*t1**2))

C---replace the overall factor of V which was removed
      ggsym=V*(ggs+ggv1+ggv2+ggv3+ggv4+ggv5
     &     +ggv6+ggv7+ggv8+ggv9+ggv10+ggv11)
      return        
      end 

