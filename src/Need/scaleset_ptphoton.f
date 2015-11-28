      subroutine scaleset_ptphoton(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c---  pt(photon)
      implicit none
      include 'constants.f'
      include 'process.f'
      integer n2,n3
      double precision p(mxpart,4),mu0,pt
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3

      if    ((case .eq. 'Wgamma') .or.
     &       (case .eq. 'Zgamma')) then
        mu0=pt(5,p)
      elseif((case .eq. 'dirgam')) then
        mu0=pt(3,p)
      else
        write(6,*) 'dynamicscale pt(photon)'//
     &             ' not supported for this process.'
	stop
      endif
      
      return
      end
      
