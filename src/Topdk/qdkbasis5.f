      double precision function qdkbasis5(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qdkbasis5=  + ro * ( S(i1,5)*S(i1,6) - 2.D0*S(i1,5)*S(i1,6)*t1 + 
     &    S(i1,5)*S(i1,6)*t1**2 + S(i1,5)*S(i2,6) - 2.D0*S(i1,5)*S(i2,6
     &    )*t1 + S(i1,5)*S(i2,6)*t1**2 - S(i1,6)*S(i2,5)*t1 + S(i1,6)*
     &    S(i2,5)*t1**2 - S(i2,5)*S(i2,6)*t1 + S(i2,5)*S(i2,6)*t1**2 - 
     &    S(3,5)*S(i1,6) - S(3,5)*S(i1,6)*t1**2 - S(3,5)*S(i2,6) + S(3,
     &    5)*S(i2,6)*t1 - S(3,5)*S(i2,6)*t1**2 - S(3,5)*S(3,6) + S(3,6)
     &    *S(i1,5) - 2.D0*S(3,6)*S(i1,5)*t1 + S(3,6)*S(i1,5)*t1**2 - S(
     &    3,6)*S(i2,5)*t1 + S(3,6)*S(i2,5)*t1**2 )
      qdkbasis5 = qdkbasis5 + ro**2 * (  - 1.D0/2.D0*S(i1,5)*S(i1,6) + 
     &    1.D0/2.D0*S(i1,5)*S(i2,6) - 1.D0/2.D0*S(i1,5)*S(i2,6)*t1 + 1.D
     &    0/2.D0*S(i1,6)*S(i2,5)*t1 )

      return
      end
