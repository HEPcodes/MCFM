      double precision function qdkbasis6(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qdkbasis6=  + ro * ( 2.D0*S(i1,5)*S(i1,6) - 4.D0*S(i1,5)*S(i1,6)*
     &    t1 + 2.D0*S(i1,5)*S(i1,6)*t1**2 + 2.D0*S(i1,5)*S(i2,6) - 4.D0
     &    *S(i1,5)*S(i2,6)*t1 + 2.D0*S(i1,5)*S(i2,6)*t1**2 - 2.D0*S(i1,
     &    6)*S(i2,5)*t1 + 2.D0*S(i1,6)*S(i2,5)*t1**2 - 2.D0*S(i2,5)*S(
     &    i2,6)*t1 + 2.D0*S(i2,5)*S(i2,6)*t1**2 - 2.D0*S(3,5)*S(i1,6)
     &     + 6.D0*S(3,5)*S(i1,6)*t1**2 - 2.D0*S(3,5)*S(i2,6) + 2.D0*S(3
     &    ,5)*S(i2,6)*t1 + 6.D0*S(3,5)*S(i2,6)*t1**2 - 2.D0*S(3,5)*S(3,
     &    6) + 8.D0*S(3,5)*S(3,6)*t1**2 + 2.D0*S(3,6)*S(i1,5) - 4.D0*S(
     &    3,6)*S(i1,5)*t1 + 2.D0*S(3,6)*S(i1,5)*t1**2 - 2.D0*S(3,6)*S(
     &    i2,5)*t1 + 2.D0*S(3,6)*S(i2,5)*t1**2 )
      qdkbasis6 = qdkbasis6 + ro**2 * ( S(i1,5)*S(i1,6)*t1 + 2.D0*S(i1,
     &    5)*S(i2,6)*t1 + S(i2,5)*S(i2,6)*t1 + 2.D0*S(3,5)*S(i1,6) - S(
     &    3,5)*S(i1,6)*t1 + S(3,5)*S(i2,6) - S(3,5)*S(i2,6)*t1 + 2.D0*
     &    S(3,5)*S(3,6) + S(3,6)*S(i1,5) + S(3,6)*S(i1,5)*t1 + S(3,6)*
     &    S(i2,5)*t1 + 2.D0*S(5,6)*S(i1,i2)*t1 - 2.D0*S(5,6)*S(i1,i2)*
     &    t1**2 )
      qdkbasis6 = qdkbasis6 + ro**3 * (  - 1.D0/2.D0*S(i1,5)*S(i2,6) + 
     &    1.D0/2.D0*S(i1,6)*S(i2,5) - 1.D0/2.D0*S(5,6)*S(i1,i2) )

      return
      end
