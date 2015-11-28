      double complex function ttbggmm(i1,i2,i3,i4,i5,i6,i7,i8,i0)
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
      ttbggmm =  + s1345**(-1) * ( 1/(zb(i1,i0))/(zb(i2,i0))*za(i1,i2)*
     &    za(i6,i7)*zb(i4,i0)*zb(i6,i0)*mt**2 + 1/(zb(i1,i0))/(zb(i2,i0
     &    ))*za(i1,i3)*za(i2,i3)*za(i6,i7)*zb(i3,i0)*zb(i3,i4)*zb(i6,i0
     &    ) + 1/(zb(i1,i0))/(zb(i2,i0))*za(i1,i3)*za(i2,i5)*za(i6,i7)*
     &    zb(i3,i4)*zb(i5,i0)*zb(i6,i0) + 1/(zb(i1,i0))/(zb(i2,i0))*za(
     &    i1,i3)*za(i2,i7)*zb(i3,i0)*zb(i4,i0)*mt**2 + 1/(zb(i1,i0))/(
     &    zb(i2,i0))*za(i1,i5)*za(i2,i7)*zb(i4,i0)*zb(i5,i0)*mt**2 - 
     &    1/(zb(i2,i0))*za(i1,i2)*za(i1,i3)*za(i6,i7)*zb(i3,i4)*zb(i6,
     &    i0) )
      ttbggmm = ttbggmm + 1/(zb(i1,i0))/(zb(i1,i2))*za(i2,i3)*za(i6,i7)
     & *zb(i3,i4)*zb(i6,i0) + 1/(zb(i1,i0))/(zb(i1,i2))*za(i2,i7)*zb(i4
     &    ,i0)*mt**2 + 1/(zb(i1,i2))/(zb(i2,i0))*za(i1,i3)*za(i6,i7)*
     &    zb(i3,i4)*zb(i6,i0) + 1/(zb(i1,i2))/(zb(i2,i0))*za(i1,i7)*zb(
     &    i4,i0)*mt**2

      return
      end
