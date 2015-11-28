      double precision function deltarj(i,j,p,pjet)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pjet(mxpart,4),phi1,phi2,etarap,dphi
      integer i,j
      
      phi1=atan2(p(i,1),p(i,2))
      phi2=atan2(pjet(j,1),pjet(j,2))
      dphi=phi1-phi2
      if (dphi .gt. pi) dphi=twopi-dphi
      if (dphi .lt. -pi) dphi=twopi+dphi
      deltarj=(etarap(i,p)-etarap(j,pjet))**2+dphi**2
      deltarj=dsqrt(deltarj)
      return
      end
