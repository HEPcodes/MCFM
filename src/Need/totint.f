      double precision function totint(r,wgt)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      double precision r(mxdim),virtint,realint,wgt
      totint=virtint(r,wgt)+realint(r,wgt)
      return
      end
