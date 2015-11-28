      double precision function qdkbasis3(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qdkbasis3=  + ro * (  - 20.D0*S(i1,5)*S(i1,6) - 12.D0*S(i1,5)*S(
     &    i1,6)*t1 + 20.D0*S(i1,5)*S(i2,6) - 52.D0*S(i1,5)*S(i2,6)*t1
     &     - 8.D0*S(i1,6)*S(i2,5) + 28.D0*S(i1,6)*S(i2,5)*t1 - 8.D0*S(
     &    i2,5)*S(i2,6) - 12.D0*S(i2,5)*S(i2,6)*t1 - 32.D0*S(3,5)*S(i1,
     &    6) + 12.D0*S(3,5)*S(i1,6)*t1 - 20.D0*S(3,5)*S(i2,6) + 12.D0*
     &    S(3,5)*S(i2,6)*t1 - 40.D0*S(3,5)*S(3,6) - 20.D0*S(3,6)*S(i1,5
     &    ) - 12.D0*S(3,6)*S(i1,5)*t1 - 8.D0*S(3,6)*S(i2,5) - 12.D0*S(3
     &    ,6)*S(i2,5)*t1 - 40.D0*S(5,6)*S(i1,i2)*t1 + 40.D0*S(5,6)*S(i1
     &    ,i2)*t1**2 )
      qdkbasis3 = qdkbasis3 + ro**2 * ( 6.D0*S(i1,5)*S(i2,6) - 6.D0*S(
     &    i1,6)*S(i2,5) + 10.D0*S(5,6)*S(i1,i2) )
      qdkbasis3 = qdkbasis3 - 32.D0*S(3,5)*S(i1,6) + 64.D0*S(3,5)*S(i1,
     &    6)*t1 - 160.D0*S(3,5)*S(i1,6)*t1**2 - 32.D0*S(3,5)*S(i2,6) + 
     &    64.D0*S(3,5)*S(i2,6)*t1 - 160.D0*S(3,5)*S(i2,6)*t1**2 - 32.D0
     &    *S(3,5)*S(3,6) + 64.D0*S(3,5)*S(3,6)*t1 - 160.D0*S(3,5)*S(3,6
     &    )*t1**2

      return
      end
