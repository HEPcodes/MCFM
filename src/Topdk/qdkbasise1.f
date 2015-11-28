      double precision function qdkbasise1(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qdkbasise1=  + ro * (  - 2.D0*S(1,2)*S(5,6) + 2.D0*S(1,5)*S(1,6)
     &     + 2.D0*S(1,5)*S(2,6) + 2.D0*S(1,6)*S(2,5) + 2.D0*S(2,5)*S(2,
     &    6) )
      qdkbasise1 = qdkbasise1 + 8.D0*S(3,5)*S(i1,6) + 8.D0*S(3,5)*S(i2,
     &    6) + 8.D0*S(3,5)*S(3,6)

      return
      end
