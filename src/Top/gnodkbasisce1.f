      double precision function gnodkbasisce1()
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      gnodkbasisce1=  + xnsq * ( 8.D0 - 8.D0*t1**(-1) - 16.D0*t1 + 8.D0
     &    *t2**(-1) )

      return
      end
