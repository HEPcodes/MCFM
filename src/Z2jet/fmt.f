      double complex function fmt(s12,s34,s56)
      implicit none
      include 'constants.f'
      include 'masses.f'
      double precision s12,s34,s56
      fmt=(1d0+(2d0*s34+s12+s56)/15d0/mt**2)/(24d0*mt**2)
      return
      end

