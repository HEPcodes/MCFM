      double complex function aqqb_wgg(i1,i2,i4,i5,i6,i7)
      implicit none
C-----Results taken from Trocsanyi and Nagy, hep/ph9806317
C-----A37-A40
c---not finished yet
      include 'constants.f'
      include 'masses.f'
      include 'prods.f'
      integer i1,i2,i4,i3,i5,i6,i7
      double complex t2,Ap_pp,Ap_mm,Ap_pm,Ap_mp
      double precision s123,s234
      aqqb_wgg=czip
      s123=s(i1,i2)+s(i1,i3)+s(i2,i3)
      s234=s(i2,i3)+s(i2,i4)+s(i3,i4)

C A(1+,2+,3+,4-) (A-37)
      Ap_pp=-za(i4,i5)**2*zb(i5,i6)/(za(i1,i2)*za(i2,i3)*za(i3,i4))
C A(1+,2-,3-,4-) (A-38)
      Ap_mm=-zb(i1,i6)**2*za(i5,i6)/(zb(i1,i2)*zb(i2,i3)*zb(i3,i4))
C A(1+,2+,3-,4-) (A-39)
      Ap_pm=
     . -za(i3,i1)*zb(i1,i2)*za(i4,i5)*t2(i3,i1,i2,i6)
     . /(za(i1,i2)*s(i2,i3)*s123)
     . +za(i3,i4)*zb(i4,i2)*za(i1,i6)*t2(i5,i3,i4,i2)
     . /(zb(i3,i4)*s(i2,i3)*s234)
     . +t2(i5,i3,i4,i2)*t2(i3,i1,i2,i6)
     . /(za(i1,i2)*zb(i3,i4)*s(i2,i3))
C A(1+,2-,3+,4-) (A-40)
      Ap_mp=
     . +zb(i1,i3)**2*za(i4,i5)*t2(i2,i1,i3,i6)
     . /(zb(i1,i2)*s(i2,i3)*s123)
     . -za(i2,i4)**2*zb(i1,i6)*t2(i5,i2,i4,i3)
     . /(za(i3,i4)*s(i2,i3)*s234)
     . -zb(i1,i3)*za(i2,i4)*zb(i1,i6)*za(i4,i5)
     . /(zb(i1,i2)*za(i3,i4)*s(i2,i3))
      return 
      end
