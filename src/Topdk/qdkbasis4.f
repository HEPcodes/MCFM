      double precision function qdkbasis4(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qdkbasis4=  + ro * (  - 1.D0/2.D0*S(i1,5)*S(i1,6) + 3.D0/2.D0*S(
     &    i1,5)*S(i1,6)*t1 - S(i1,5)*S(i1,6)*t1**2 - 1.D0/2.D0*S(i1,5)*
     &    S(i2,6) + 3.D0/2.D0*S(i1,5)*S(i2,6)*t1 - S(i1,5)*S(i2,6)*
     &    t1**2 + 1.D0/2.D0*S(i1,6)*S(i2,5)*t1 - S(i1,6)*S(i2,5)*t1**2
     &     + 1.D0/2.D0*S(i2,5)*S(i2,6)*t1 - S(i2,5)*S(i2,6)*t1**2 - S(3
     &    ,5)*S(i1,6) + 7.D0/2.D0*S(3,5)*S(i1,6)*t1 - 3.D0*S(3,5)*S(i1,
     &    6)*t1**2 - 1.D0/2.D0*S(3,5)*S(i2,6) + 5.D0/2.D0*S(3,5)*S(i2,6
     &    )*t1 - 3.D0*S(3,5)*S(i2,6)*t1**2 - S(3,5)*S(3,6) + 4.D0*S(3,5
     &    )*S(3,6)*t1 - 4.D0*S(3,5)*S(3,6)*t1**2 - 1.D0/2.D0*S(3,6)*S(
     &    i1,5) + 3.D0/2.D0*S(3,6)*S(i1,5)*t1 - S(3,6)*S(i1,5)*t1**2 + 
     &    1.D0/2.D0*S(3,6)*S(i2,5)*t1 - S(3,6)*S(i2,5)*t1**2 )
      qdkbasis4 = qdkbasis4 + ro**2 * (  - 1.D0/4.D0*S(i1,5)*S(i1,6) - 
     &    1.D0/2.D0*S(i1,5)*S(i2,6)*t1 - 1.D0/2.D0*S(i1,6)*S(i2,5) + 1.D
     &    0/2.D0*S(i1,6)*S(i2,5)*t1 - 1.D0/4.D0*S(i2,5)*S(i2,6) - 1.D0/
     &    2.D0*S(3,5)*S(i1,6) - 1.D0/2.D0*S(3,5)*S(i2,6) - S(3,5)*S(3,6
     &    ) - 1.D0/2.D0*S(3,6)*S(i1,5) - 1.D0/2.D0*S(3,6)*S(i2,5) - S(5
     &    ,6)*S(i1,i2)*t1 + S(5,6)*S(i1,i2)*t1**2 )
      qdkbasis4 = qdkbasis4 + ro**3 * ( 1.D0/4.D0*S(5,6)*S(i1,i2) )

      return
      end
