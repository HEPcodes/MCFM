      subroutine trocs_4q()
C--Results taken from Nagy and Trocsanyi-hep-ph/9806317
C--Appendix A
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'prods.f'
      integer i1,i2,i3,i4,i5,i6,i7
      double complex t2,tx,A_1,A_2,A_3,A_4

      double precision t123,t134,t234,t345,t567,t167,t267

      t123=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t167=s(i1,i6)+s(i6,i7)+s(i7,i1)
      t267=s(i2,i6)+s(i6,i7)+s(i7,i2)
      t567=s(i5,i6)+s(i6,i7)+s(i7,i5)
      t345=s(i3,i4)+s(i4,i5)+s(i5,i3)
      t234=s(i2,i3)+s(i3,i4)+s(i4,i2)
      t134=s(i1,i3)+s(i3,i4)+s(i4,i1)
C The four-quark one-gluon amplitudes have the form:

C--A50
      A_1= 
     . -zb(i1,i5)*t2(i4,i1,i5,i3)*t2(i4,i2,i6,i7)*za(i6,i2)
     . /(za(i4,i5)*s(i1,i5)*s(i3,i4)*t267)
     . -zb(i1,i7)*t2(i6,i1,i7,i5)*za(i4,i2)**2*zb(i2,i3)
     . /(za(i4,i5)*s(i3,i4)*t234*t167)
     . -t2(i4,i1,i5,i7)*t2(i6,i2,i4,i3)*za(i4,i2)
     . /(za(i1,i5)*za(i5,i4)*s(i3,i4)*t234)
     . +zb(i5,i3)*t2(i4,i3,i5,i1)*t2(i4,i2,i6,i7)*za(i6,i2)
     . /(za(i4,i5)*s(i3,i4)*t345*t267)
     . +zb(i1,i7)*tx(i6,i1,i7,i3,i5,i4)*zb(i3,i5)*za(i4,i2)
     . /(za(i4,i5)*s(i3,i4)*t345*t167)

C--A51
      A_1=
     . (zb(i1,i3)**2*za(i5,i1)*t2(i4,i2,i6,i7)*za(i6,i2))
     . /(zb(i3,i5)*s(i1,i5)*s(i3,i4)*t267)
     . +(zb(i1,i7)*t2(i6,i1,i7,i3)*t2(i5,i2,i4,i3)*za(i4,i2))
     . /(zb(i3,i5)*s(i3,i4)*t234*t167)
     . -zb(i1,i3)*zb(i1,i7)*t2(i6,i2,i4,i3)*za(i4,i2)
     . /(zb(i1,i5)*zb(i5,i3)*s(i3,i4)*t234)
     . +zb(i1,i3)*za(i5,i4)*tx(i3,i4,i5,i2,i6,i7)*za(i6,i2)
     . /(zb(i3,i5)*s(i3,i4)*t345*t267)
     . -zb(i1,i7)*t2(i6,i1,i7,i3)*za(i5,i4)*t2(i2,i4,i5,i3)
     . /(zb(i3,i5)*s(i3,i4)*t345*t167)
 
 
C--A53
      A_2=
     .  -zb(i5,i3)*t2(i4,i3,i5,i1)*t2(i4,i2,i6,i7)*za(i6,i2)
     . /(za(i4,i5)*s(i3,i4)*t345*t267)
     . -zb(i1,i7)*t2(i6,i1,i7,i3)*zb(i2,i5)*za(i4,i2)**2
     . /(za(i4,i5)*s(i2,i5)*s(i3,i4)*t167)
     . -zb(i1,i3)*t2(i4,i1,i3,i5)*t2(i4,i2,i6,i7)*za(i6,i2)
     . /(za(i4,i5)*s(i3,i4)*t134*t267)
     . -zb(i1,i3)*t2(i4,i1,i3,i7)*za(i6,i2)**2
     . /(za(i4,i5)*za(i5,i2)*s(i3,i4)*t134)
     . -zb(i1,i7)*tx(i6,i1,i7,i3,i5,i4)*zb(i3,i5)*za(i4,i2)
     . /(za(i4,i5)*s(i3,i4)*t345*t167)
 
C--A54
      A_2=
     . -zb(i1,i3)*za(i5,i4)*tx(i3,i4,i5,i2,i6,i7)*za(i6,i2)
     . /(zb(i3,i5)*s(i3,i4)*t345*t267)
     . +zb(i1,i3)**2*za(i4,i1)*t2(i5,i2,i6,i7)*za(i6,i2)
     . /(zb(i3,i5)*s(i3,i4)*t134*t267)
     . +zb(i1,i7)*t2(i6,i1,i7,i3)*za(i5,i4)*t2(i2,i4,i5,i3)
     . /(zb(i3,i5)*s(i3,i4)*t345*t167)
     . -zb(i1,i3)*t2(i4,i1,i3,i7)*t2(i6,i2,i5,i3)
     . /(zb(i3,i5)*zb(i5,i2)*s(i3,i4)*t134)
     . +zb(i1,i7)*t2(i6,i1,i7,i3)*za(i5,i2)*t2(i4,i2,i5,i3)
     . /(zb(i3,i5)*s(i2,i5)*s(i3,i4)*t167)
 
C--A56
      A_3=
     . -zb(i5,i3)*t2(i4,i3,i5,i1)*t2(i4,i2,i6,i7)*za(i6,i2)
     . /(za(i4,i5)*s(i3,i5)*t345*t267)
     . -zb(i1,i7)*tx(i6,i1,i7,i3,i5,i4)*zb(i3,i5)*za(i4,i2)
     . /(za(i4,i5)*s(i3,i5)*t345*t167)
 
C--A57
      A_3=
     . -zb(i1,i4)*za(i5,i3)*tx(i4,i3,i5,i2,i6,i7)*za(i6,i2)
     . /(zb(i4,i5)*s(i3,i5)*t345*t267)
     . +zb(i1,i7)*t2(i6,i1,i7,i4)*za(i5,i3)*t2(i2,i3,i5,i4)
     . /(zb(i4,i5)*s(i3,i5)*t345*t167)
 
C--A58
      A_4=
     . zb(i1,i3)*za(i5,i4)*tx(i3,i4,i5,i2,i6,i7)*za(i6,i2)
     . /(zb(i3,i5)*s(i4,i5)*t345*t267)
     . -zb(i1,i7)*t2(i6,i1,i7,i3)*za(i5,i4)*t2(i2,i4,i5,i3)
     . /(zb(i3,i5)*s(i4,i5)*t345*t167)
  
C--A59
      A_4=
     . zb(i5,i4)*t2(i3,i4,i5,i1)*t2(i3,i2,i6,i7)*za(i6,i2)
     . /(za(i3,i5)*s(i4,i5)*t345*t267)
     . +zb(i1,i7)*tx(i6,i1,i7,i4,i5,i3)*zb(i4,i5)*za(i3,i2)
     . /(za(i3,i5)*s(i4,i5)*t345*t167)
      return
      end
