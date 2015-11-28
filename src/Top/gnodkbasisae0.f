      double precision function gnodkbasisae0()
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      gnodkbasisae0=  + ro*xnsq * (  - 8.D0 + 4.D0*t1**(-1) + 4.D0*
     &    t2**(-1) )
      gnodkbasisae0 = gnodkbasisae0 + ro * (  - 4.D0*t1**(-1) - 4.D0*
     &    t2**(-1) )
      gnodkbasisae0 = gnodkbasisae0 + ro**2*xnsq * (  - t1**(-2) - 
     &    t2**(-2) )
      gnodkbasisae0 = gnodkbasisae0 + ro**2 * ( t1**(-2) + 2.D0*
     &    t1**(-1) + t2**(-2) + 2.D0*t2**(-1) )
      gnodkbasisae0 = gnodkbasisae0 + xnsq * (  - 16.D0 + 4.D0*t1**(-1)
     &     + 16.D0*t1 - 16.D0*t1**2 + 4.D0*t2**(-1) )
      gnodkbasisae0 = gnodkbasisae0 + 8.D0 - 4.D0*t1**(-1) - 4.D0*
     &    t2**(-1)

      return
      end
