      double precision function qdkbasis1(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qdkbasis1=  + ro * (  - 1.D0/2.D0*S(i1,5)*S(i1,6) + 1.D0/2.D0*S(
     &    i1,5)*S(i1,6)*t1 + 1.D0/2.D0*S(i1,5)*S(i2,6) - 5.D0/2.D0*S(i1
     &    ,5)*S(i2,6)*t1 + 2.D0*S(i1,5)*S(i2,6)*t1**2 + 3.D0/2.D0*S(i1,
     &    6)*S(i2,5)*t1 - 2.D0*S(i1,6)*S(i2,5)*t1**2 + 1.D0/2.D0*S(i2,5
     &    )*S(i2,6)*t1 - S(3,5)*S(i1,6) + 3.D0/2.D0*S(3,5)*S(i1,6)*t1
     &     - 1.D0/2.D0*S(3,5)*S(i2,6) + 3.D0/2.D0*S(3,5)*S(i2,6)*t1 - 
     &    S(3,5)*S(3,6) + 2.D0*S(3,5)*S(3,6)*t1 - 1.D0/2.D0*S(3,6)*S(i1
     &    ,5) + 1.D0/2.D0*S(3,6)*S(i1,5)*t1 + 1.D0/2.D0*S(3,6)*S(i2,5)*
     &    t1 - S(5,6)*S(i1,i2)*t1 + 3.D0*S(5,6)*S(i1,i2)*t1**2 - 2.D0*
     &    S(5,6)*S(i1,i2)*t1**3 )
      qdkbasis1 = qdkbasis1 + ro**2 * ( 1.D0/4.D0*S(i1,5)*S(i2,6) - 1.D0
     &    /4.D0*S(i1,6)*S(i2,5) + 1.D0/4.D0*S(5,6)*S(i1,i2) - 1.D0/2.D0
     &    *S(5,6)*S(i1,i2)*t1 )
      qdkbasis1 = qdkbasis1 + 4.D0*S(3,5)*S(i1,6)*t1 - 12.D0*S(3,5)*S(
     &    i1,6)*t1**2 + 8.D0*S(3,5)*S(i1,6)*t1**3 + 4.D0*S(3,5)*S(i2,6)
     &    *t1 - 12.D0*S(3,5)*S(i2,6)*t1**2 + 8.D0*S(3,5)*S(i2,6)*t1**3
     &     + 4.D0*S(3,5)*S(3,6)*t1 - 12.D0*S(3,5)*S(3,6)*t1**2 + 8.D0*
     &    S(3,5)*S(3,6)*t1**3

      return
      end
