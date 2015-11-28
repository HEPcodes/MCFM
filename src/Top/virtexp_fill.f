      subroutine virtexp_fill(t1_passed,t2_passed,mass,s12)
************************************************************************
*   Subroutine to fill the virtexp common block                        *
************************************************************************
      implicit none
      double precision t1_passed,t2_passed,mass,s12,ddilog
      include 'constants.f'
      include 'scale.f'
      include 'virtexp.f'
      
      t1=t1_passed
      t2=t2_passed
      xlf=dfloat(nf)
      ro=4d0*mass**2/s12
      omro=1d0-ro
      rmuom2=2d0*dlog(scale/mass)
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
      f1=(vlpm**2/2d0-2*vdmb-pi**2/3d0)/b
      f2=(-b*vlsm+vlpm**2/4d0+vdmp+pi**2/12d0)/b**3
      f3=(-b**3*vlsm-3d0*b*vlsm+0.75d0*vlpm**2
     & +3d0*vdmp+pi**2/4d0+2d0*b**3)/b**5
      brackt=vlsm*vltm-vlsm**2/4d0-pisq/12d0+vdt
      brackw=vlsm*vlwm-vlsm**2/4d0-pisq/12d0+vdw
      bracks=vdmp+vlpm**2/4d0+pisq/12d0

      INVG = 1d0/(t1*t2-0.25d0*ro)
      
      return
      end
      
