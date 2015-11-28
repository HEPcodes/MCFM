      double precision function gdkbasis2()
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'virtexp.f'

      gdkbasis2 =  + ro * (  - 3.D0*s(1,2)*s(5,6) + 6.D0*s(1,2)*s(5,6)*
     &    t1 - 4.D0*s(1,2)*s(5,6)*t1**2 + s(1,5)*s(1,6) - 2.D0*s(1,5)*
     &    s(1,6)*t1**(-1) - 3.D0*s(1,5)*s(2,6) + 3.D0*s(1,5)*s(2,6)*t1
     &     + s(1,5)*s(3,6) - 2.D0*s(1,5)*s(3,6)*t1**(-1) + 4.D0*s(1,6)*
     &    s(2,5) - 3.D0*s(1,6)*s(2,5)*t1 + 7.D0*s(1,6)*s(3,5) - 5.D0*s(
     &    1,6)*s(3,5)*t1**(-1) + s(1,6)*s(3,5)*t2**(-1) + 2.D0*s(2,5)*
     &    s(2,6) + s(2,5)*s(3,6) + 7.D0*s(2,6)*s(3,5) - 3.D0*s(2,6)*s(3
     &    ,5)*t1**(-1) + s(2,6)*s(3,5)*t2**(-1) + 8.D0*s(3,5)*s(3,6) - 
     &    5.D0*s(3,5)*s(3,6)*t1**(-1) + s(3,5)*s(3,6)*t2**(-1) )
      gdkbasis2 = gdkbasis2 + ro**2 * (  - 2.D0*s(1,2)*s(5,6) + 3.D0/2.D
     &    0*s(1,2)*s(5,6)*t1**(-1) + 1.D0/2.D0*s(1,5)*s(1,6)*t1**(-2)
     &     + 3.D0/2.D0*s(1,5)*s(2,6)*t1**(-1) + 1.D0/2.D0*s(1,5)*s(3,6)
     &    *t1**(-2) - 2.D0*s(1,6)*s(2,5)*t1**(-1) - 1.D0/2.D0*s(1,6)*s(
     &    2,5)*t2**(-1) + s(1,6)*s(3,5)*t1**(-2) - 1.D0/2.D0*s(1,6)*s(3
     &    ,5)*t1**(-1) - 1.D0/2.D0*s(1,6)*s(3,5)*t2**(-1) - 1.D0/2.D0*
     &    s(2,5)*s(2,6)*t1**(-1) - 1.D0/2.D0*s(2,5)*s(2,6)*t2**(-1) - 1.
     &    D0/2.D0*s(2,5)*s(3,6)*t1**(-1) - 1.D0/2.D0*s(2,5)*s(3,6)*
     &    t2**(-1) + 1.D0/2.D0*s(2,6)*s(3,5)*t1**(-2) - s(2,6)*s(3,5)*
     &    t1**(-1) - s(2,6)*s(3,5)*t2**(-1) + s(3,5)*s(3,6)*t1**(-2) - 
     &    s(3,5)*s(3,6)*t1**(-1) - s(3,5)*s(3,6)*t2**(-1) )
      gdkbasis2 = gdkbasis2 + ro**3 * (  - 1.D0/4.D0*s(1,2)*s(5,6)*
     &    t1**(-2) + 1.D0/4.D0*s(1,2)*s(5,6)*t1**(-1) + 1.D0/4.D0*s(1,2
     &    )*s(5,6)*t2**(-1) - 1.D0/4.D0*s(1,5)*s(2,6)*t1**(-2) - 1.D0/4.
     &    D0*s(1,5)*s(2,6)*t1**(-1) - 1.D0/4.D0*s(1,5)*s(2,6)*t2**(-1)
     &     + 1.D0/4.D0*s(1,6)*s(2,5)*t1**(-2) + 1.D0/4.D0*s(1,6)*s(2,5)
     &    *t1**(-1) + 1.D0/4.D0*s(1,6)*s(2,5)*t2**(-1) )
      gdkbasis2 = gdkbasis2 + 2.D0*s(1,5)*s(1,6) - 4.D0*s(1,5)*s(1,6)*
     &    t1 + 2.D0*s(1,5)*s(1,6)*t1**2 + 2.D0*s(1,5)*s(2,6) - 4.D0*s(1
     &    ,5)*s(2,6)*t1 + 2.D0*s(1,5)*s(2,6)*t1**2 + 2.D0*s(1,5)*s(3,6)
     &     - 4.D0*s(1,5)*s(3,6)*t1 + 2.D0*s(1,5)*s(3,6)*t1**2 - 2.D0*s(
     &    1,6)*s(2,5)*t1 + 2.D0*s(1,6)*s(2,5)*t1**2 + 6.D0*s(1,6)*s(3,5
     &    ) - 24.D0*s(1,6)*s(3,5)*t1 + 14.D0*s(1,6)*s(3,5)*t1**2 - 2.D0
     &    *s(2,5)*s(2,6)*t1 + 2.D0*s(2,5)*s(2,6)*t1**2 - 2.D0*s(2,5)*s(
     &    3,6)*t1 + 2.D0*s(2,5)*s(3,6)*t1**2 + 6.D0*s(2,6)*s(3,5) - 22.D
     &    0*s(2,6)*s(3,5)*t1 + 14.D0*s(2,6)*s(3,5)*t1**2 + 6.D0*s(3,5)*
     &    s(3,6) - 24.D0*s(3,5)*s(3,6)*t1 + 16.D0*s(3,5)*s(3,6)*t1**2

      return
      end
