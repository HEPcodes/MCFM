      double precision function gdkbasisce2()
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      gdkbasisce2=  + ro*xnsq * ( 4.D0*s(1,2)*s(5,6)*t1**(-1) - 4.D0*s(
     &    1,2)*s(5,6)*t2**(-1) - 4.D0*s(1,5)*s(2,6)*t1**(-2) + 4.D0*s(1
     &    ,5)*s(2,6)*t1**(-1) + 4.D0*s(1,5)*s(2,6)*t2**(-1) - 4.D0*s(1,
     &    6)*s(2,5)*t1**(-1) + 4.D0*s(1,6)*s(2,5)*t2**(-2) - 4.D0*s(1,6
     &    )*s(2,5)*t2**(-1) )
      gdkbasisce2 = gdkbasisce2 + xnsq * (  - 16.D0*s(1,6)*s(3,5)*
     &    t1**(-1) + 16.D0*s(1,6)*s(3,5)*t2**(-1) - 16.D0*s(2,6)*s(3,5)
     &    *t1**(-1) + 16.D0*s(2,6)*s(3,5)*t2**(-1) - 16.D0*s(3,5)*s(3,6
     &    )*t1**(-1) + 16.D0*s(3,5)*s(3,6)*t2**(-1) )

      gdkbasisce2= gdkbasisce2/8d0
      return
      end
