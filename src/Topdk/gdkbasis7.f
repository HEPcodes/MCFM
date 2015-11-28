      double precision function gdkbasis7()
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'virtexp.f'

      gdkbasis7 =  + ro * ( 2.D0*s(1,2)*s(5,6)*t1**(-1) - s(1,5)*s(1,6)
     &    *t1**(-2) - 2.D0*s(1,5)*s(1,6)*t1**(-1) - 2.D0*s(1,5)*s(2,6)
     &     - s(1,5)*s(2,6)*t1**(-2) + 2.D0*s(1,5)*s(2,6)*t1**(-1) - s(1
     &    ,5)*s(3,6)*t1**(-2) + 2.D0*s(1,5)*s(3,6)*t1**(-1) + 2.D0*s(1,
     &    6)*s(2,5) + s(1,6)*s(2,5)*t1**(-1) + s(1,6)*s(2,5)*t2**(-1)
     &     + 2.D0*s(1,6)*s(3,5)*t1**(-2) + s(1,6)*s(3,5)*t1**(-1) + 3.D0
     &    *s(1,6)*s(3,5)*t2**(-1) + s(2,5)*s(2,6)*t1**(-1) + s(2,5)*s(2
     &    ,6)*t2**(-1) + 3.D0*s(2,5)*s(3,6)*t1**(-1) + s(2,5)*s(3,6)*
     &    t2**(-1) + s(2,6)*s(3,5)*t1**(-2) + 2.D0*s(2,6)*s(3,5)*
     &    t2**(-1) + 2.D0*s(3,5)*s(3,6)*t1**(-2) + 2.D0*s(3,5)*s(3,6)*
     &    t1**(-1) + 2.D0*s(3,5)*s(3,6)*t2**(-1) )
      gdkbasis7 = gdkbasis7 + ro**2 * (  - 1.D0/2.D0*s(1,2)*s(5,6)*
     &    t1**(-2) - 1.D0/2.D0*s(1,2)*s(5,6)*t1**(-1) - 1.D0/2.D0*s(1,2
     &    )*s(5,6)*t2**(-1) + s(1,5)*s(1,6)*t1**(-2) + s(1,5)*s(1,6)*
     &    t1**(-1) + s(1,5)*s(1,6)*t2**(-1) - 1.D0/2.D0*s(1,5)*s(2,6)*
     &    t1**(-2) + 1.D0/2.D0*s(1,5)*s(2,6)*t1**(-1) + 1.D0/2.D0*s(1,5
     &    )*s(2,6)*t2**(-1) + 1.D0/2.D0*s(1,6)*s(2,5)*t1**(-2) - 1.D0/2.
     &    D0*s(1,6)*s(2,5)*t1**(-1) - 1.D0/2.D0*s(1,6)*s(2,5)*t2**(-1)
     &     )
      gdkbasis7 = gdkbasis7 - 8.D0*s(1,5)*s(1,6) + 4.D0*s(1,5)*s(1,6)*
     &    t1**(-1) + 4.D0*s(1,5)*s(1,6)*t1 - 8.D0*s(1,5)*s(2,6) + 4.D0*
     &    s(1,5)*s(2,6)*t1**(-1) + 4.D0*s(1,5)*s(2,6)*t1 - 8.D0*s(1,5)*
     &    s(3,6) + 4.D0*s(1,5)*s(3,6)*t1**(-1) + 4.D0*s(1,5)*s(3,6)*t1
     &     - 4.D0*s(1,6)*s(2,5) + 4.D0*s(1,6)*s(2,5)*t1 - 4.D0*s(1,6)*
     &    s(3,5)*t1**(-1) - 4.D0*s(1,6)*s(3,5)*t1 - 4.D0*s(2,5)*s(2,6)
     &     + 4.D0*s(2,5)*s(2,6)*t1 - 4.D0*s(2,5)*s(3,6) + 4.D0*s(2,5)*
     &    s(3,6)*t1 + 4.D0*s(2,6)*s(3,5) - 4.D0*s(2,6)*s(3,5)*t1**(-1)
     &     - 4.D0*s(2,6)*s(3,5)*t1 - 4.D0*s(3,5)*s(3,6)*t1**(-1)

      return
      end
