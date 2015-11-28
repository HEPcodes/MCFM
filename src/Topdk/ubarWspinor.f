      subroutine ubarWspinor(pn,pe,pb,xmt,xmw,ft)
      implicit none
C-----Spinor for outgoing on-shell top decaying to nu(pn)+e^+(pe)+b(pb)
c----- Note that (potentially scaled) top and W masses are now passed in
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      double precision m,pn(4),pe(4),pb(4),pw(4),sw,snb,xmt,xmw
      double complex zpe(4),zpt(4),f(4),ft(4),ifac
      integer i,j,nu

C----setup complex momenta
      do nu=1,4
      zpe(nu)=dcmplx(pe(nu))
      pw(nu)=pn(nu)+pe(nu)
      zpt(nu)=dcmplx(pn(nu)+pe(nu)+pb(nu))
      enddo

      snb=2d0*(pn(1)*pb(1)-pn(4)*pb(4)-pn(2)*pb(2)-pn(3)*pb(3))
      sw=pw(1)**2-pw(4)**2-pw(2)**2-pw(3)**2
C-----Top quark is on shell
      ifac=dcmplx(sw-xmw**2,xmw*wwidth)
     . *dcmplx(zip,xmt*twidth)/dsqrt(snb)
C---setup spinor for positron
      call ubar0spinor(zpe,1,f)
C---multiply with ptslash
      call spb(f,zpt,ft) 

      do j=1,4
      ft(j)=gwsq*(ft(j)+xmt*f(j))/ifac
      enddo

        
      return
      end
      

