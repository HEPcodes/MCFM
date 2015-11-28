      function qqsdel(t1,ro)  
      implicit double precision (a-h,o-z)
      common/para1/s,xm2,xlf,xmu,rmuom2,naem,nbeam1,nbeam2
      external ddilog
      parameter (v=8.,xn=3.,xn4=5.,tr=0.5)
      data pi/3.141592653589793d0/
      
      t2=1.d0-t1
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
      
      
      srlgpr=log((1-b)/(1+b)*t2/t1)*log((1+b)/(1-b)*t2/t1)
      srl21p=ddilog(1.d0-(1.d0+b)/(1.d0-b)*t1/t2)
      srl21m=ddilog(1.d0-(1.d0-b)/(1.d0+b)*t1/t2)
      srl22p=ddilog(1.d0-(1.d0+b)/(1.d0-b)*t2/t1)
      srl22m=ddilog(1.d0-(1.d0-b)/(1.d0+b)*t2/t1)
C      
      srl12=log(t1/t2)
      srlgro=log(4.d0*t1*t2/ro)
      srlg1=log((1.d0+b)/(1.d0-b))/b
      srl2l=(ddilog(-4.d0*b/(1.d0-b)**2)-ddilog(4.d0*b/(1.+b)**2))/b
      srl212=0.5*log(4.d0*t1*t2/ro)**2+ddilog(1.d0-ro/4.d0/t1/t2)
      sl3515=(srlgpr+b**2*srlg1**2+srl21p+srl21m)
      sl3525=(srlgpr+b**2*srlg1**2+srl22p+srl22m)
C  -------------------------------------
      qqqqv=2.*t1**2+2.*t2**2+ro

      qqsv = qqqqv*(tr*(4.0*vlsm/3.0+(-4.0)*rmuom2/3.0+(-20.0)/9.0)*xlf+
     1   xn*(-vlsm*(vlwm+vltm)/2.0+vlsm**2/4.0+11.0*vlsm/6.0+11.0*rmuom2
     2   /3.0-pi**2/6.0+22.0/9.0)+(vlsm**2/2.0+(-3.0)*vlsm/2.0-b*vlpm-(b
     3   **2+1)*(pi**2/b+f1)/2.0-pi**2/3.0+6)/xn+tr*(4.0*b*(ro/2.0+1)*vl
     4   pm/3.0+(-4.0)*ro/3.0+(-20.0)/9.0))+xn*(-(ro/2.0-t1*t2)*(vlwm/ub
     5   ar+vltm/tbar)-(t2**2+t1**2)*vlsm*(vlwm+vltm)-(t1-t2)*vlsm*(vltm
     6   -vlwm)+(t2**2+t1**2)*vlsm**2/2.0-(2*ro+3.0/2.0)*vlsm+ro*(vdw+vd
     7   t)/2.0+(t2-t1)*(vdt-vdw)-pi**2*ro/12.0+1)-((t1-t2)**2+b**2)*vlp
     8   m/(b*xn)/2.0+f2*xn*(b**2*(6*t2**2+6*t1**2+1.0/2.0)-5*(t2**2+t1*
     9   *2)-b**4+2)+f3*xn*(t1-t2)**2/2.0
      qqsrv = qqqqv*(-2*(-vlsm*vlwm-vlsm*vltm+vlsm**2+srl212/2.0)/xn+rmu
     1   om2*v*(vlwm+vltm-2*vlsm+(-3.0)/2.0)/xn-2*v*(-vlsm**2+vlsm+(-1.0
     2   )/2.0)/xn-2*(1-ro/2.0)*(vlpm*vlsm+b*srl2l/4.0)/(b*xn)+(1-ro/2.0
     3   )*v*vlpm/(b*xn)+xn*(sl3525+sl3515)/2.0)

c--debug
      qqsrv=0d0
c--debug
      qqsdel=qqsv+qqsrv
      if (naem.eq.0) return
      vt1=vltm-vlsm
      vt2=vlwm-vlsm
      cdel=v/2./xn*(-9-2.*pi**2/3.+vt1**2+vt2**2+1.5*vt1+1.5*vt2)
      qqsdel=qqsdel-cdel*qqqqv
 31   continue      
      return        
      end 

      function qqadel(t1,ro)  
      implicit double precision (a-h,o-z)
      common/para1/s,xm2,xlf,xmu,rmuom2,naem,nbeam1,nbeam2
      external ddilog
      parameter (v=8.,xn=3.,xn4=5.,tr=0.5)
      data pi/3.141592653589793d0/
      
      t2=1.d0-t1
      tbar=t1-0.25*ro
      ubar=t2-0.25*ro
      b=sqrt(1.d0-ro)   
      xlp=0.5*(1.d0+b)
      xlm=0.5*(1.d0-b)
      vlpm=log(xlp/xlm)     
      vlsm=log(4.d0/ro)      
      vltm=log(4.d0*t1/ro)      
      vlwm=log(4.d0*t2/ro)
      vdw=ddilog(1.d0-ro/(4.d0*t2))-0.5d0*vlwm**2      
      vdt=ddilog(1.d0-ro/(4.d0*t1))-0.5d0*vltm**2
      vdmp=ddilog(-xlm/xlp)    
      f2=(-b*vlsm+vlpm**2/4.d0+vdmp+pi**2/12.d0)/b**3
            
      srlgpr=log((1-b)/(1+b)*t2/t1)*log((1+b)/(1-b)*t2/t1)
      srl21p=ddilog(1.d0-(1.d0+b)/(1.d0-b)*t1/t2)
      srl21m=ddilog(1.d0-(1.d0-b)/(1.d0+b)*t1/t2)
      srl22p=ddilog(1.d0-(1.d0+b)/(1.d0-b)*t2/t1)
      srl22m=ddilog(1.d0-(1.d0-b)/(1.d0+b)*t2/t1)
C      
      srlg1=log((1.d0+b)/(1.d0-b))/b
      sl3515=(srlgpr+b**2*srlg1**2+srl21p+srl21m)
      sl3525=(srlgpr+b**2*srlg1**2+srl22p+srl22m)
C  -------------------------------------
      qqqqv=2.*t1**2+2.*t2**2+ro

      qqav = ((t1-t2)*(-vlsm*(vlwm+vltm)+vlsm**2/2.0+(3-ro)*vlsm+pi**2/6
     1   .0)-(ro/2.0-t1*t2)*(vltm/tbar-vlwm/ubar)-(t2**2+t1**2)*vlsm*(vl
     2   tm-vlwm)-(t1-t2)*(vdw+vdt)+ro*(vdt-vdw)/2.0+(b**4+2*b**2-1)*f2*
     3   (t1-t2))*xn4/xn-qqqqv*vlsm*(vltm-vlwm)*xn4/xn/2.0
      qqarv = qqqqv*(2*vlsm*(vltm-vlwm)-sl3525/2.0+sl3515/2.0)*xn4/xn
c--debug
      qqarv=0d0
c--debug
      qqadel=qqav+qqarv
      return        
      end 

      subroutine qqplss(t1,t2,ro,qqsrp,qqarp,qqsrl,qqarl)        
      implicit double precision (a-h,o-z)
      common/para1/s,xm2,xlf,xmu,rmuom2,naem,nbeam1,nbeam2
      parameter (v=8.,xn=3.,xn4=5.)

C     Soft limit

      b=sqrt(1.d0-ro)
      vlsm=log(4.d0/ro)
      srl12=log(t1/t2)
      srlgro=log(4.d0*t1*t2/ro)
      srlg1=log((1.d0+b)/(1.d0-b))/b
      vomt1=log(t2)
      vomt2=log(t1)
      qqqqv=2.*t1**2+2.*t2**2+ro
      qqsrp = qqqqv*v*(4.*vlsm-2.*rmuom2-2.)/xn
     & -qqqqv*2.*srlg1*(1.-ro/2.0)/xn
     & +2.*qqqqv*srlgro/xn
     & -float(naem)*qqqqv*v/2./xn*(-2.*vomt1-3./2.)
     & -float(naem)*qqqqv*v/2./xn*(-2.*vomt2-3./2.)
      qqarp = 2.*qqqqv*srl12*xn4/xn
      qqsrl=v/2./xn*(8.*qqQQv-4.*float(naem)*qqQQv)
      qqarl=0.

      return
      end

      subroutine qqplus(t1,t2,ro,qqsrp,qqarp,qqsrl,qqarl)        
      implicit double precision (a-h,o-z)
      common/para1/s,xm2,xlf,xmu,rmuom2,naem,nbeam1,nbeam2
      parameter (v=8.,xn=3.,xn4=5.)
      
      omt1=1.d0-t1
      omt2=1.d0-t2
      omtx=t1+t2
      tx=1.d0-omtx

      dlam=sqrt(omtx**2-ro)
      ro45=4.d0*tx+ro
      vlsm=log(4.d0/ro)      
C
      rlgxro=log((4.d0*tx+ro)/ro)/tx    
      rlg11=log(t1/(omt1))                           
      rlg22=log(t2/(omt2))                           
      rlgro=log(4.*omt1*omt2/ro)

C     rl35=log((2.*tx+ro-2.*dlam-2.)/(2.*tx+ro+2.*dlam-2.))/dlam
      rl35=log(((1.+dlam)**2-tx**2)/((1.-dlam)**2-tx**2))/dlam

      qqqqv=2.*t1**2+2.*t2**2+ro
      qqqqv1=2.*t1**2+2.*omt1**2+ro*omt1/t2
      qqqqv2=2.*omt2**2+2.*t2**2+ro*omt2/t1

      vomt1=log(omt1)
      vomt2=log(omt2)
      
      qqsrp = qqqqv*v*(4.*vlsm-2.*rlgxro*tx-2.*rmuom2-2.)/xn
     & -qqqqv*rl35*(1.-ro/2.0)/xn
     & +2.*qqqqv*rlgro/xn+xn*qqqqv*(rlg22+rlg11)
     & -float(naem)*qqqqv1*v/2./xn*(-2.*vomt1-3./2.)
     & -float(naem)*qqqqv2*v/2./xn*(-2.*vomt2-3./2.)
      qqarp = qqqqv*(rlg11-rlg22)*xn4/xn

      qqsrl=v/2./xn*(8.*qqQQv-2.*float(naem)*(qqQQv1+qqQQv2))
      qqarl=0.

      return
      end

      function qqsrs(t1,t2,ro)   
      implicit double precision (a-h,o-z)
      common/para1/s,xm2,xlf,xmu,rmuom2,naem,nbeam1,nbeam2
      parameter (v=8.,xn=3.,xn4=5.)

      omtx=t1+t2
      tx=1.d0-omtx
      omt1=1.d0-t1
      omt2=1.d0-t2

      b=sqrt(1.d0-ro)
      dlam=sqrt(omtx**2-ro)
      ro45=4.d0*tx+ro
      vlsm=log(4.d0/ro)
      vtx=log(tx)
C     
      rlgxro=log((4.d0*tx+ro)/ro)/tx    
      rlg12=log(t1/(omt2))/tx                        
      rlg21=log(t2/(omt1))/tx                        
      rlg11=log(t1/(omt1))                           
      rlg22=log(t2/(omt2))                           
      rlgro=log(4.*omt1*omt2/ro)
c
c
      rl34=log((1.-(tx+dlam)**2)/(1.-(tx-dlam)**2))/dlam/tx
      rl35=log(((1.+dlam)**2-tx**2)/((1.-dlam)**2-tx**2))/dlam
C      

      qqsrs1 = rlg12*(4*t2**3/t1-8*t2**2/t1-ro*t2**2/t1**2+4*t2**2+4*t1*
     1   t2+6*t2/t1+2*ro*t2/t1**2-4*t2+4*t1**2-4*t1-2/t1-ro/t1**2+ro+2)/
     2   xn+xn*rlg12*(-4*t2**3/t1+8*t2**2/t1+ro*t2**2/t1**2-4*t1*t2-6*t2
     3   /t1-2*ro*t2/t1**2-2*t2-4*t1**2/omtx-2*t1/omtx+6*t1+2/t1+ro/t1**
     4   2-ro/omtx/2.0-1/omtx+2)+rlg21*(4*t2**2+4*t1*t2-4*t2+4*t1**3/t2-
     5   8*t1**2/t2+6*t1/t2-2/t2-ro*t1**2/t2**2+2*ro*t1/t2**2-ro/t2**2+4
     6   *t1**2-4*t1+ro+2)/xn+xn*rlg21*(-4*t1*t2+2*t2-4*t1**3/t2+8*t1**2
     7   /t2-6*t1/t2+2/t2+ro*t1**2/t2**2-2*ro*t1/t2**2+ro/t2**2-4*t1**2/
     8   omtx+2*t1/omtx+2*t1-ro/omtx/2.0-1/omtx)
      qqsrs2 = xn*rlgxro*(4*t2**2-t2/omt1+2*t2+4*t1**2-t1/omt2+2*t1+3.0*
     1   ro/2.0)+rlgxro*(-4*t2**2+t2/omt1-4*t1**2+t1/omt2-2*ro-2)/xn+xn*
     2   rlg22*(2*t2+4*t1**2/omtx+2*t1/omtx-2*t1+ro/omtx/2.0+1/omtx-1)+x
     3   n*rlg11*(2*t2+4*t1**2/omtx-2*t1/omtx-2*t1+ro/omtx/2.0+1/omtx+1)
     4   +4*rlgro/xn
      qqsrs3 = xn*(4*t2**2/t1-4*t2/t1-ro*t2/t1**2+4*t2+4*t1**2/t2-4*t1/t
     1   2+ro/t2+2/t2-ro*t1/t2**2+ro/t2**2+4*t1+ro/t1+2/t1+ro/t1**2-2/om
     2   t2-2/omt1+8)*vtx+(-4*t2**2/t1+4*t2/t1+ro*t2/t1**2-4*t2-4*t1**2/
     3   t2+4*t1/t2-ro/t2-2/t2+ro*t1/t2**2-ro/t2**2-4*t1-ro/t1-2/t1-ro/t
     4   1**2+2/omt2+2/omt1-8)*vtx/xn
      qqsrs4 = xn*rl34*((-3.0)*t2**7/2.0+(-9.0)*t1*t2**6/2.0+15.0*t2**6/
     1   4.0+(-3.0)*t1**2*t2**5/2.0+15.0*t1*t2**5/2.0+(-9.0)*t2**5/4.0+1
     2   5.0*t1**3*t2**4/2.0+(-15.0)*t1**2*t2**4/4.0+(-9.0)*t1*t2**4/4.0
     3   +3.0*t2**4/4.0+15.0*t1**4*t2**3/2.0-15*t1**3*t2**3+9.0*t1**2*t2
     4   **3/2.0+(-3.0)*t1**5*t2**2/2.0+(-15.0)*t1**4*t2**2/4.0+9.0*t1**
     5   3*t2**2/2.0+(-3.0)*t1**2*t2**2/2.0+(-9.0)*t1**6*t2/2.0+15.0*t1*
     6   *5*t2/2.0+(-9.0)*t1**4*t2/4.0+(-3.0)*t1**7/2.0+15.0*t1**6/4.0+(
     7   -9.0)*t1**5/4.0+3.0*t1**4/4.0)/dlam**4+xn*rl34*(5.0*t2**5/2.0+9
     8   .0*t1*t2**4/2.0+(-19.0)*t2**4/4.0+t1**2*t2**3-5*t1*t2**3+3.0*t2
     9   **3/4.0+t1**3*t2**2-t1**2*t2**2/2.0+9.0*t1*t2**2/4.0+9.0*t1**4*
     :   t2/2.0-5*t1**3*t2+9.0*t1**2*t2/4.0-t1*t2+5.0*t1**5/2.0+(-19.0)*
     ;   t1**4/4.0+3.0*t1**3/4.0)/dlam**2+xn*rl34*(-t2**3-t1*t2**2+5*t2*
     <   *2-t1**2*t2-t1*t2-ro*t2/2.0+(-7.0)*t2/2.0-t1**3-3*t1**2/omtx+5*
     =   t1**2-ro*t1/2.0-t1/2.0+ro/omtx/4.0+3.0*ro/4.0+9.0/4.0)+rl34*(2*
     >   t2+2*t1-2)/xn
      qqsrs5 = xn*(4*t2**2/t1-4*t2/t1-ro*t2/t1**2+4*t2+4*t1**2/t2-4*t1/t
     1   2+ro/t2+2/t2-ro*t1/t2**2+ro/t2**2+4*t1+ro/t1+2/t1+ro/t1**2-2/om
     2   t2-2/omt1+8)*vlsm+(-4*t2**2/t1+4*t2/t1+ro*t2/t1**2-4*t2-4*t1**2
     3   /t2+4*t1/t2-ro/t2-2/t2+ro*t1/t2**2-ro/t2**2-4*t1-ro/t1-2/t1-ro/
     4   t1**2+2/omt2+2/omt1-8)*vlsm/xn+rmuom2*(2*t2**2/t1-2*t2/t1-ro*t2
     5   /t1**2/2.0+2*t2+2*t1**2/t2-2*t1/t2+ro/t2/2.0+1/t2-ro*t1/t2**2/2
     6   .0+ro/t2**2/2.0+2*t1+ro/t1/2.0+1/t1+ro/t1**2/2.0-1/omt2-1/omt1+
     7   4)/xn+xn*rmuom2*(-2*t2**2/t1+2*t2/t1+ro*t2/t1**2/2.0-2*t2-2*t1*
     8   *2/t2+2*t1/t2-ro/t2/2.0-1/t2+ro*t1/t2**2/2.0-ro/t2**2/2.0-2*t1-
     9   ro/t1/2.0-1/t1-ro/t1**2/2.0+1/omt2+1/omt1-4)
      qqsrs6 = rl35*((-3.0)*t2**7/4.0+(-9.0)*t1*t2**6/4.0+9.0*t2**6/4.0+
     1   (-3.0)*t1**2*t2**5/4.0+9.0*t1*t2**5/2.0+(-3.0)*t2**5/4.0+15.0*t
     2   1**3*t2**4/4.0+(-9.0)*t1**2*t2**4/4.0+(-3.0)*t1*t2**4/4.0+(-3.0
     3   )*t2**4/4.0+15.0*t1**4*t2**3/4.0-9*t1**3*t2**3+3.0*t1**2*t2**3/
     4   2.0+(-3.0)*t1**5*t2**2/4.0+(-9.0)*t1**4*t2**2/4.0+3.0*t1**3*t2*
     5   *2/2.0+3.0*t1**2*t2**2/2.0+(-9.0)*t1**6*t2/4.0+9.0*t1**5*t2/2.0
     6   +(-3.0)*t1**4*t2/4.0+(-3.0)*t1**7/4.0+9.0*t1**6/4.0+(-3.0)*t1**
     7   5/4.0+(-3.0)*t1**4/4.0)/(dlam**4*xn)+xn*rl35*((-3.0)*t2**6/4.0+
     8   (-3.0)*t1*t2**5/2.0+3.0*t2**5/4.0+3.0*t1**2*t2**4/4.0+3.0*t1*t2
     9   **4/4.0+3*t1**3*t2**3+(-3.0)*t1**2*t2**3/2.0+3.0*t1**4*t2**2/4.
     :   0+(-3.0)*t1**3*t2**2/2.0+(-3.0)*t1**5*t2/2.0+3.0*t1**4*t2/4.0+(
     ;   -3.0)*t1**6/4.0+3.0*t1**5/4.0)/dlam**4
      qqsrsa = +rl35*(9.0*t2**5/4.0+13.0
     <   *t1*t2**4/4.0+(-17.0)*t2**4/4.0+(-3.0)*t1**2*t2**3/2.0-3*t1*t2*
     =   *3+t2**3/2.0+(-3.0)*t1**3*t2**2/2.0+5.0*t1**2*t2**2/2.0+t1*t2**
     >   2/2.0+t2**2/2.0+13.0*t1**4*t2/4.0-3*t1**3*t2+t1**2*t2/2.0+9.0*t
     ?   1**5/4.0+(-17.0)*t1**4/4.0+t1**3/2.0+t1**2/2.0)/(dlam**2*xn)+xn
     @   *rl35*(5.0*t2**4/4.0+t1*t2**3+(-5.0)*t2**3/4.0-t1**2*t2**2/2.0+
     1   t1*t2**2/4.0+t2**2/2.0+t1**3*t2+t1**2*t2/4.0-t1*t2+5.0*t1**4/4.
     2   0+(-5.0)*t1**3/4.0+t1**2/2.0)/dlam**2+rl35*((-3.0)*t2**3/2.0+t1
     3   *t2**2/2.0+2*t2**2+t1**2*t2/2.0-t1*t2+ro*t2/4.0+(-3.0)*t2/4.0+(
     4   -3.0)*t1**3/2.0+2*t1**2+ro*t1/4.0+(-3.0)*t1/4.0+7.0*ro/4.0+(-3.
     5   0)/4.0)/xn+xn*rl35*(-t2**2/2.0+3.0*t2/2.0+3*t1**2/omtx-t1**2/2.
     6   0+(-3.0)*t1/2.0-ro/omtx/4.0-ro/4.0+(-1.0)/2.0)
      qqsrs7 = xn*(ro*t2/t1**2/2.0+2*t2+ro/t2/2.0+1/t2+ro*t1/t2**2/2.0-r
     1   o/t2**2/2.0+2*t1+ro/t1/2.0+1/t1-ro/t1**2/2.0-6)+(-ro*t2/t1**2/2
     2   .0+t2-ro/t2/2.0-1/t2-ro*t1/t2**2/2.0+ro/t2**2/2.0+t1-ro/t1/2.0-
     3   1/t1+ro/t1**2/2.0+1/omt2+1/omt1+3)/xn
      qqsrs8 = xn*(6*t2**5+6*t1*t2**4-6*t2**4-12*t1**2*t2**3-12*t1**3*t2
     1   **2+12*t1**2*t2**2-3*t2**2+6*t1**4*t2+6*t1*t2-6*t2+6*t1**5-6*t1
     2   **4-3*t1**2+18*t1-12)/dlam**4+(-3*t2**5-3*t1*t2**4+3*t2**4+6*t1
     3   **2*t2**3+3*t2**3+6*t1**3*t2**2-6*t1**2*t2**2-3*t1*t2**2+3*t2**
     4   2-3*t1**4*t2-3*t1**2*t2-6*t1*t2+6*t2-3*t1**5+3*t1**4+3*t1**3+3*
     5   t1**2-18*t1+12)/(dlam**4*xn)+xn*(8*t2**2-8*t2+8*t1**2-24*t1+16)
     6   /ro45**2+(24*t1**2*t2-48*t1*t2+24*t2+24*t1**3-96*t1**2+120*t1-4
     7   8)/(dlam**4*xn*ro45)+xn*(-24*t1**2*t2+48*t1*t2-24*t2-24*t1**3+9
     8   6*t1**2-120*t1+48)/(dlam**4*ro45)+xn*(12*t2-12*t1+12)/ro45
      qqsrs9 = (6*t2**3-2*t1*t2**2-2*t2**2-2*t1**2*t2-3*t2+6*t1**3-2*t1*
     1   *2+5*t1-5)/(dlam**2*xn)+xn*(-8*t2**3+4*t2**2+6*t2-8*t1**3+4*t1*
     2   *2-18*t1+17)/dlam**2+(-8*t2**2+8*t2-8*t1**2+24*t1-16)/(xn*ro45*
     3   *2)+xn*(24*t1**2*t2-64*t1*t2+40*t2+24*t1**3-96*t1**2+144*t1-72)
     4   /(dlam**2*ro45)+xn*(16*t1**2*t2-32*t1*t2+16*t2+16*t1**3-64*t1**
     5   2+80*t1-32)/(dlam**2*ro45**2)+(-8*t1**2*t2+16*t1*t2-8*t2-8*t1**
     6   3+32*t1**2-48*t1+24)/(dlam**2*xn*ro45)+(-16*t1**2*t2+32*t1*t2-1
     7   6*t2-16*t1**3+64*t1**2-80*t1+32)/(dlam**2*xn*ro45**2)+(-4*t2+4*
     8   t1-4)/(xn*ro45)
      
      qqsrs=qqsrs1+qqsrs2+qqsrs3+qqsrs4+qqsrs5
     & +qqsrs6+qqsrsa+qqsrs7+qqsrs8+qqsrs9
      if (naem.eq.0) goto 33
      z1=t2/omt1
      omz1=1.-z1
      qqqqv1=2.*t1**2+2.*omt1**2+ro/z1
      cqoz1=v/2./xn*(omz1*log(omz1)
     & -(1.+z1**2)/omz1*log(z1)+3./2.+2.*z1)/z1

      z2=t1/omt2
      omz2=1.-z2
      qqqqv2=2.*omt2**2+2.*t2**2+ro/z2
      cqoz2=v/2./xn*(omz2*log(omz2)
     & -(1.+z2**2)/omz2*log(z2)+3./2.+2.*z2)/z2
      qqsrs=qqsrs-cqoz1/omt1*qqQQv1-cqoz2/omt2*qqQQv2
 33   continue
      return       
      end

      function qqars(t1,t2,ro)   
      implicit double precision (a-h,o-z)
      common/para1/s,xm2,xlf,xmu,rmuom2,naem,nbeam1,nbeam2
      parameter (v=8.,xn=3.,xn4=5.)

      omtx=t1+t2
      tx=1.-omtx
      omt1=1.-t1
      omt2=1.-t2

      b=sqrt(1.-ro)
      dlam=sqrt(omtx**2-ro)
C
      rlgxro=log((4.*tx+ro)/ro)/tx    
      rlg12=log(t1/(omt2))/tx                        
      rlg21=log(t2/(omt1))/tx                        
      rlg11=log(t1/(omt1))                           
      rlg22=log(t2/(omt2))                           
      rlgro=log(4*(omt1)*(omt2)/ro)
C
C     rl34=log((2.*tx**2+2.*dlam*tx-2.*tx-ro)
C    & /(2.*tx**2-2.*dlam*tx-2.*tx-ro))/dlam/tx
C     rl35=log((2.*tx+ro-2.*dlam-2.)/(2.*tx+ro+2.*dlam-2.))/dlam

      rl34=log((1.-(tx+dlam)**2)/(1.-(tx-dlam)**2))/dlam/tx
      rl35=log(((1.+dlam)**2-tx**2)/((1.-dlam)**2-tx**2))/dlam
C      
      qqars = (-rlg22-rlg21+rlg12+rlg11)*(2*(t2**2+t1**2)+ro/2.0+1)*xn4/
     1   (xn*omtx)+rl34*(t1**2-t2**2)*(-(t2+t1)**2+t2+t1-1)*xn4/(dlam**2
     2   *xn)/2.0+rl35*(t1**2-t2**2)*(1-(t2+t1)**2)*xn4/(dlam**2*xn)/2.0
     3   +2*(rlg21*t1-rlg12*t2)*xn4/(xn*omtx)+rl34*(t1-t2)*(t2+t1-ro/omt
     4   x+1/omtx+1)*xn4/xn/2.0+rl35*(t1-t2)*(t2+t1+(ro-1)/omtx+2)*xn4/x
     5   n/2.0+(rlg22+rlg11)*(t2-t1)*xn4/(xn*omtx)+rlgxro*(t1-t2)*xn4/xn
     6   -2*(t1-t2)*xn4/(dlam**2*xn)-(1/omt1-1/omt2)*xn4/xn
      return       
      end

