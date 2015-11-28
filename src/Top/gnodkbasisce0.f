      double precision function gnodkbasisce0()
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      gnodkbasisce0=  + ro*xnsq * ( 4.D0*t1**(-1) - 4.D0*t2**(-1) )
      gnodkbasisce0 = gnodkbasisce0 + ro**2*xnsq * (  - t1**(-2) + 
     &    t2**(-2) )
      gnodkbasisce0 = gnodkbasisce0 + xnsq * (  - 8.D0 + 4.D0*t1**(-1)
     &     + 16.D0*t1 - 4.D0*t2**(-1) )

      return
      end
