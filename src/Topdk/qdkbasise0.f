      double precision function qdkbasise0(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qdkbasise0=  + ro * (  - 2.D0*S(i1,5)*S(i1,6) - 2.D0*S(i1,5)*S(i2
     &    ,6) - 2.D0*S(i1,6)*S(i2,5) - 2.D0*S(i2,5)*S(i2,6) - 4.D0*S(1,
     &    2)*S(5,6)*t1 + 4.D0*S(1,2)*S(5,6)*t1**2 + 4.D0*S(1,5)*S(2,6)
     &     - 4.D0*S(1,5)*S(2,6)*t1 + 4.D0*S(1,6)*S(2,5)*t1 - 2.D0*S(3,5
     &    )*S(i1,6) - 2.D0*S(3,5)*S(i2,6) - 4.D0*S(3,5)*S(3,6) - 2.D0*
     &    S(3,6)*S(i1,5) - 2.D0*S(3,6)*S(i2,5) )
      qdkbasise0 = qdkbasise0 + ro**2 * ( S(1,2)*S(5,6) )
      qdkbasise0 = qdkbasise0 - 8.D0*S(3,5)*S(i1,6) + 16.D0*S(3,5)*S(i1
     &    ,6)*t1 - 16.D0*S(3,5)*S(i1,6)*t1**2 - 8.D0*S(3,5)*S(i2,6) + 
     &    16.D0*S(3,5)*S(i2,6)*t1 - 16.D0*S(3,5)*S(i2,6)*t1**2 - 8.D0*
     &    S(3,5)*S(3,6) + 16.D0*S(3,5)*S(3,6)*t1 - 16.D0*S(3,5)*S(3,6)*
     &    t1**2

      return
      end
