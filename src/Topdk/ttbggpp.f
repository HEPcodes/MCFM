      double complex function ttbggpp(i1,i2,i3,i4,i5,i6,i7,i8,i0)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer i1,i2,i3,i4,i5,i6,i7,i8,i0
      double precision  s12,s1345,s1678
      s12=s(i1,i2)
      s1345=s(i1,i3)+s(i1,i5)
      s1678=s(i6,i1)+s(i8,i1)
      ttbggpp =  + s1345**(-1) * ( 1/(za(i1,i0))/(za(i2,i0))*za(i3,i0)
     &    **2*za(i6,i7)*zb(i1,i3)*zb(i2,i6)*zb(i3,i4) + 1/(za(i1,i0))/(
     &    za(i2,i0))*za(i3,i0)*za(i5,i0)*za(i6,i7)*zb(i1,i5)*zb(i2,i6)*
     &    zb(i3,i4) + 1/(za(i1,i0))/(za(i2,i0))*za(i3,i0)*za(i7,i0)*zb(
     &    i1,i2)*zb(i3,i4)*mt**2 + 1/(za(i1,i0))/(za(i2,i0))*za(i3,i0)*
     &    za(i7,i0)*zb(i1,i4)*zb(i2,i3)*mt**2 + 1/(za(i1,i0))/(za(i2,i0
     &    ))*za(i5,i0)*za(i7,i0)*zb(i1,i4)*zb(i2,i5)*mt**2 - 1/(za(i2,
     &    i0))*za(i7,i0)*zb(i1,i2)*zb(i1,i4)*mt**2 )
      ttbggpp = ttbggpp + 1/(za(i1,i0))/(za(i1,i2))*za(i3,i0)*za(i6,i7)
     & *zb(i2,i6)*zb(i3,i4) + 1/(za(i1,i0))/(za(i1,i2))*za(i7,i0)*zb(i2
     &    ,i4)*mt**2 + 1/(za(i1,i2))/(za(i2,i0))*za(i3,i0)*za(i6,i7)*
     &    zb(i1,i6)*zb(i3,i4) + 1/(za(i1,i2))/(za(i2,i0))*za(i7,i0)*zb(
     &    i1,i4)*mt**2

      return
      end
