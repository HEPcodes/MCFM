 
      gg_a(2,2,1,2,2)=
     . +za(p2,t6)*za(p2,t6)*za(p3,t5)*zb(p1,t6)*zb(p1,t6)*zb(p4,t6)
     . * (1d0/s34/s126/s26/s12)
 
     . +za(p2,t6)*za(p3,t5)*zb(p1,p4)*zb(p1,t6)*zb(p1,t6)/zb(p1,p2)
     . * (-1d0/s34/s126/s26)
 
      gg_a(2,2,1,1,1)=+ za(p2,t6)*za(p2,t6)*za(p2,p3)*zb(p1,p4)
     . *zb(p1,t6)*zb(p2,t5)
     . * (- mb**2/s34/s126/s26/s25/s12)
 
      gg_a(2,2,1,2,1)=+ za(p2,t6)*za(p2,t6)*za(p3,t5)*zb(p1,p4)
     . *zb(p1,t6)
     . * (- mb/s34/s126/s26/s12)
 
      gg_a(2,2,1,1,2)=
     . +za(p2,t6)*za(p2,p3)*zb(p1,p4)*zb(p1,t6)*zb(p1,t6)*zb(p2,t5)
     . /zb(p1,p2)
     . * (- mb/s34/s126/s26/s25)
 
     . +za(p2,t6)*za(p2,t6)*za(p2,p3)*zb(p1,t6)*zb(p1,t6)*zb(p2,t5)
     . *zb(p4,t6)
     . * (mb/s34/s126/s26/s25/s12)
 
      gg_a(2,1,2,2,2)=
     . +za(p1,t6)*za(p1,t6)*za(p3,t5)*zb(p2,t6)*zb(p2,t6)*zb(p4,t6)
     . * (1d0/s34/s126/s26/s12)
 
     . +za(p1,t6)*za(p3,t5)*zb(p2,p4)*zb(p2,t6)*zb(p2,t6)/zb(p1,p2)
     . * (1d0/s34/s126/s26)
 
      gg_a(2,1,2,1,1)=
     . +za(p1,t6)*za(p1,t6)*za(p1,t6)*za(p2,p3)*zb(p2,t5)*zb(p2,t6)
     . *zb(p4,t6)
     . /za(p1,p2)
     . * (- mb**2/s34/s126/s16/s26/s25)
 
     . +za(p1,t6)*za(p1,t6)*za(p2,p3)*zb(p2,p4)*zb(p2,t5)*zb(p2,t6)
     . * (mb**2/s34/s126/s16/s26/s25 - mb**2/s34/s126
     . /s26/s25/s12)
 
      gg_a(2,1,2,2,1)=
     . +za(p1,t6)*za(p1,t6)*za(p1,t6)*za(p3,t5)*zb(p2,t6)
     . *zb(p4,t6)/za(p1,p2)
     . * (- mb/s34/s126/s16/s26)
 
     . +za(p1,t6)*za(p1,t6)*za(p3,t5)*zb(p2,p4)*zb(p2,t6)
     . * (mb/s34/s126/s16/s26 - mb/s34/s126/s26/s12)
 
      gg_a(2,1,2,1,2)=
     . +za(p1,t6)*za(p1,t6)*za(p2,p3)*zb(p2,t5)*zb(p2,t6)*zb(p2,t6)
     . *zb(p4,t6)
     . * (mb/s34/s126/s26/s25/s12)
 
     . +za(p1,t6)*za(p2,p3)*zb(p2,p4)*zb(p2,t5)*zb(p2,t6)*zb(p2,t6)
     . /zb(p1,p2)
     . * (mb/s34/s126/s26/s25)
 
      gg_a(1,2,1,2,2)=
     . +za(p2,t6)*za(p2,p3)*za(p2,t5)*zb(p1,t6)*zb(p1,t6)*zb(p2,p4)
     . * (mb**2/s34/s126/s16/s26/s25 - mb**2/s34/s126
     . /s26/s25/s12)
 
     . +za(p2,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p1,t6)*zb(p1,t6)
     . *zb(p2,p4)
     . /zb(p1,p2)
     . * (- mb**2/s34/s126/s16/s26/s25)
 
      gg_a(1,2,1,1,1)=
     . +za(p2,t6)*za(p2,t6)*za(p2,p3)*zb(p1,t6)*zb(p4,t5)/za(p1,p2)
     . * (1d0/s34/s126/s26)
 
     . +za(p2,t6)*za(p2,t6)*za(p3,t6)*zb(p1,t6)*zb(p1,t6)*zb(p4,t5)
     . * (1d0/s34/s126/s26/s12)
 
      gg_a(1,2,1,2,1)=
     . +za(p2,t6)*za(p2,t6)*za(p2,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     . /za(p1,p2)
     . * (mb/s34/s126/s26/s25)
 
     . +za(p2,t6)*za(p2,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p1,t6)
     . *zb(p2,p4)
     . * (mb/s34/s126/s26/s25/s12)
 
      gg_a(1,2,1,1,2)=
     . +za(p2,t6)*za(p2,p3)*zb(p1,t6)*zb(p1,t6)*zb(p4,t5)
     . * (mb/s34/s126/s16/s26 - mb/s34/s126/s26/s12)
 
     . +za(p2,t6)*za(p3,t6)*zb(p1,t6)*zb(p1,t6)*zb(p1,t6)*zb(p4,t5)
     . /zb(p1,p2)
     . * (- mb/s34/s126/s16/s26)
 
      gg_a(1,1,2,2,2)=+ za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p2,p4)
     . *zb(p2,t6)*zb(p2,t6)
     . * (- mb**2/s34/s126/s26/s25/s12)
 
      gg_a(1,1,2,1,1)=
     . +za(p1,t6)*za(p1,t6)*za(p1,p3)*zb(p2,t6)*zb(p4,t5)/za(p1,p2)
     . * (-1d0/s34/s126/s26)
 
     . +za(p1,t6)*za(p1,t6)*za(p3,t6)*zb(p2,t6)*zb(p2,t6)*zb(p4,t5)
     . * (1d0/s34/s126/s26/s12)
 
      gg_a(1,1,2,2,1)=
     . +za(p1,t6)*za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p2,p4)*zb(p2,t6)
     . /za(p1,p2)
     . * (- mb/s34/s126/s26/s25)
 
     . +za(p1,t6)*za(p1,t6)*za(p3,t6)*za(p2,t5)*zb(p2,p4)*zb(p2,t6)
     . *zb(p2,t6)
     . * (mb/s34/s126/s26/s25/s12)
 
      gg_a(1,1,2,1,2)=+ za(p1,t6)*za(p1,p3)*zb(p2,t6)*zb(p2,t6)
     . *zb(p4,t5)
     . * (- mb/s34/s126/s26/s12)
 
      gg_a(2,1,1,2,2)=
     . +za(p1,p2)*za(p1,p2)*za(p3,t5)*zb(p1,t6)*zb(p2,p4)/zb(p1,p2)
     . * (- mb**2/s34/s126/s16/s26)
 
     . +za(p2,t6)*za(p1,p2)*za(p3,t5)*zb(p1,t6)*zb(p2,p4)*zb(p2,t6)
     . /zb(p1,p2)/zb(p1,p2)
     . * (-1d0/s34/s126/s26)
 
     . +za(p2,t6)*za(p1,t6)*za(p3,t5)*zb(p1,t6)*zb(p2,t6)*zb(p4,t6)
     . /zb(p1,p2)
     . /zb(p1,p2)
     . * (1d0/s34/s126/s26)
 
      gg_a(2,1,1,1,1)=
     . +za(p2,t6)*za(p1,t6)*za(p2,p3)*zb(p1,t6)*zb(p2,p4)*zb(p2,t5)
     . /zb(p1,p2)/zb(p1,p2)
     . * (- mb**2/s34/s126/s26/s25)
 
     . +za(p2,t6)*za(p1,t6)*za(p2,p3)*zb(p2,t5)*zb(p4,t6)/zb(p1,p2)
     . * (- mb**2/s34/s126/s26/s25)
 
      gg_a(2,1,1,2,1)=
     . +za(p2,t6)*za(p1,t6)*za(p3,t5)*zb(p1,t6)*zb(p2,p4)
     . /zb(p1,p2)/zb(p1,p2)
     . * (- mb/s34/s126/s26)
 
     . +za(p2,t6)*za(p1,t6)*za(p3,t5)*zb(p4,t6)/zb(p1,p2)
     . * (- mb/s34/s126/s26)
 
      gg_a(2,1,1,1,2)=
     . +za(p1,p2)*za(p1,p2)*za(p2,p3)*zb(p1,t6)*zb(p2,p4)*zb(p2,t5)
     . /zb(p1,p2)
     . * (- mb**3/s34/s126/s16/s26/s25)
 
     . +za(p2,t6)*za(p1,p2)*za(p2,p3)*zb(p1,t6)*zb(p2,p4)
     . *zb(p2,t5)*zb(p2,t6)
     . /zb(p1,p2)/zb(p1,p2)
     . * (- mb/s34/s126/s26/s25)
 
     . +za(p2,t6)*za(p1,t6)*za(p2,p3)*zb(p1,t6)*zb(p2,t5)
     . *zb(p2,t6)*zb(p4,t6)
     . /zb(p1,p2)/zb(p1,p2)
     . * (mb/s34/s126/s26/s25)
 
      gg_a(2,2,2,2,2)=
     . +za(p1,t6)*za(p2,t6)*za(p3,t5)*zb(p1,t6)*zb(p2,t6)
     . *zb(p4,t6)/za(p1,p2)
     . /za(p1,p2)
     . * (1d0/s34/s126/s26)
 
     . +za(p1,t6)*za(p3,t5)*zb(p1,p4)*zb(p1,t6)*zb(p2,t6)/za(p1,p2)
     . * (1d0/s34/s126/s26)
 
     . +za(p2,t6)*za(p3,t5)*zb(p1,p2)*zb(p2,t6)*zb(p4,t6)/za(p1,p2)
     . * (1d0/s34/s126/s26)
 
     . +za(p3,t5)*zb(p1,p2)*zb(p1,p4)*zb(p2,t6)
     . * (mb**2/s34/s126/s16/s26 +1d0/s34/s126/s26)
 
      gg_a(2,2,2,1,1)=
     . +za(p1,t6)*za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p1,p4)
     . *zb(p1,t6)*zb(p2,t5)
     . /za(p1,p2)
     . * (mb**4/s34/s126/s16**2/s26/s25)
 
     . +za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p1,p2)*zb(p1,p4)*zb(p2,t5)
     . * (mb**2/s34/s126/s16/s26/s25 + mb**4/s34/s126
     . /s16**2/s26/s25)
 
     . +za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p1,p4)*zb(p2,t5)/za(p1,p2)
     . * (mb**4/s34/s126/s16/s26/s25)
 
     . +za(p1,t6)*za(p2,t6)*za(p2,p3)*zb(p1,p2)
     . *zb(p1,p2)*zb(p2,t5)*zb(p4,t6)/za(p1,p2)
     . * (mb**2/s34/s126/s16/s26/s25)
 
     . +za(p1,t6)*za(p2,t6)*za(p2,p3)*zb(p1,p4)*zb(p2,t5)*zb(p2,t6)
     . /za(p1,p2)/za(p1,p2)
     . * (- mb**2/s34/s126/s26/s25)
 
      gg_a(2,2,2,2,1)=
     . +za(p1,t6)*za(p1,t6)*za(p3,t5)*zb(p1,p2)*zb(p1,p4)*zb(p1,t6)
     . /za(p1,p2)
     . * (mb**3/s34/s126/s16**2/s26)
 
     . +za(p1,t6)*za(p2,t6)*za(p3,t5)*zb(p1,p2)*zb(p1,p2)*zb(p4,t6)
     . /za(p1,p2)
     . * (mb/s34/s126/s16/s26)
 
     . +za(p1,t6)*za(p2,t6)*za(p3,t5)*zb(p1,p4)*zb(p2,t6)
     . /za(p1,p2)/za(p1,p2)
     . * (- mb/s34/s126/s26)
 
     . +za(p1,t6)*za(p3,t5)*zb(p1,p2)*zb(p1,p2)*zb(p1,p4)
     . * (mb/s34/s126/s16/s26 + mb**3/s34/s126/s16**2
     . /s26)
 
     . +za(p1,t6)*za(p3,t5)*zb(p1,p2)*zb(p1,p4)/za(p1,p2)
     . * (mb**3/s34/s126/s16/s26)
 
      gg_a(2,2,2,1,2)=
     . +za(p1,t6)*za(p2,p3)*zb(p1,p4)*zb(p1,t6)*zb(p2,t5)*zb(p2,t6)
     . /za(p1,p2)
     . * (mb/s34/s126/s26/s25)
 
     . +za(p1,t6)*za(p2,t6)*za(p2,p3)*zb(p1,t6)*zb(p2,t5)*zb(p2,t6)
     . *zb(p4,t6)
     . /za(p1,p2)/za(p1,p2)
     . * (mb/s34/s126/s26/s25)
 
     . +za(p2,p3)*zb(p1,p2)*zb(p1,p4)*zb(p2,t5)*zb(p2,t6)
     . * (mb/s34/s126/s26/s25 + mb**3/s34/s126/s16
     . /s26/s25)
 
     . +za(p2,t6)*za(p2,p3)*zb(p1,p2)*zb(p2,t5)*zb(p2,t6)*zb(p4,t6)
     . /za(p1,p2)
     . * (mb/s34/s126/s26/s25)
 
      gg_a(1,1,1,2,2)=
     . +za(p1,p2)*za(p1,p2)*za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     . * (mb**2/s34/s126/s16/s26/s25 + mb**4/s34/s126
     . /s16**2/s26/s25)
 
     . +za(p1,p2)*za(p1,p2)*za(p3,t6)*za(p2,t5)*zb(p1,t6)
     . *zb(p2,p4)*zb(p2,t6)
     . /zb(p1,p2)
     . * (mb**2/s34/s126/s16/s26/s25)
 
     . +za(p1,p2)*za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)/zb(p1,p2)
     . * (mb**4/s34/s126/s16/s26/s25)
 
     . +za(p1,p2)*za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p1,t6)
     . *zb(p2,p4)
     . /zb(p1,p2)
     . * (mb**4/s34/s126/s16**2/s26/s25)
 
     . +za(p2,t6)*za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)*zb(p2,t6)
     . /zb(p1,p2)
     . /zb(p1,p2)
     . * (- mb**2/s34/s126/s26/s25)
 
      gg_a(1,1,1,1,1)=
     . +za(p2,t6)*za(p1,p2)*za(p1,p3)*zb(p4,t5)
     . * (mb**2/s34/s126/s16/s26 +1d0/s34/s126/s26)
 
     . +za(p2,t6)*za(p1,p2)*za(p3,t6)*zb(p2,t6)*zb(p4,t5)/zb(p1,p2)
     . * (1d0/s34/s126/s26)
 
     . +za(p2,t6)*za(p1,t6)*za(p1,p3)*zb(p1,t6)*zb(p4,t5)/zb(p1,p2)
     . * (1d0/s34/s126/s26)
 
     . +za(p2,t6)*za(p1,t6)*za(p3,t6)*zb(p1,t6)*zb(p2,t6)*zb(p4,t5)
     . /zb(p1,p2)
     . /zb(p1,p2)
     . * (1d0/s34/s126/s26)
 
      gg_a(1,1,1,2,1)=
     . +za(p2,t6)*za(p1,p2)*za(p1,p3)*za(p2,t5)*zb(p2,p4)
     . * (mb/s34/s126/s26/s25 + mb**3/s34/s126/s16
     . /s26/s25)
 
     . +za(p2,t6)*za(p1,p2)*za(p3,t6)*za(p2,t5)*zb(p2,p4)*zb(p2,t6)
     . /zb(p1,p2)
     . * (mb/s34/s126/s26/s25)
 
     . +za(p2,t6)*za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     . /zb(p1,p2)
     . * (mb/s34/s126/s26/s25)
 
     . +za(p2,t6)*za(p1,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     . *zb(p2,t6)
     . /zb(p1,p2)/zb(p1,p2)
     . * (mb/s34/s126/s26/s25)
 
      gg_a(1,1,1,1,2)=
     . +za(p1,p2)*za(p1,p2)*za(p1,p3)*zb(p1,t6)*zb(p4,t5)
     . * (mb/s34/s126/s16/s26 + mb**3/s34/s126/s16**2
     . /s26)
 
     . +za(p1,p2)*za(p1,p2)*za(p3,t6)*zb(p1,t6)*zb(p2,t6)*zb(p4,t5)
     .  /zb(p1,p2)
     . * (mb/s34/s126/s16/s26)
 
     . +za(p1,p2)*za(p1,p3)*zb(p1,t6)*zb(p4,t5)/zb(p1,p2)
     . * (mb**3/s34/s126/s16/s26)
 
     . +za(p1,p2)*za(p1,t6)*za(p1,p3)*zb(p1,t6)*zb(p1,t6)*zb(p4,t5)
     . /zb(p1,p2)
     . * (mb**3/s34/s126/s16**2/s26)
 
     . +za(p2,t6)*za(p1,p3)*zb(p1,t6)*zb(p2,t6)*zb(p4,t5)
     . /zb(p1,p2)/zb(p1,p2)
     . * (- mb/s34/s126/s26)
 
      gg_a(1,2,2,2,2)=
     . +za(p1,t6)*za(p2,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)*zb(p2,t6)
     . /za(p1,p2)/za(p1,p2)
     . * (- mb**2/s34/s126/s26/s25)
 
     . +za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)*zb(p2,t6)/za(p1,p2)
     . * (- mb**2/s34/s126/s26/s25)
 
      gg_a(1,2,2,1,1)=
     . +za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p1,p2)*zb(p4,t5)/za(p1,p2)
     . * (- mb**2/s34/s126/s16/s26)
 
     . +za(p1,t6)*za(p2,t6)*za(p2,p3)*zb(p1,p2)*zb(p2,t6)*zb(p4,t5)
     . /za(p1,p2)/za(p1,p2)
     . * (-1d0/s34/s126/s26)
 
     . +za(p1,t6)*za(p2,t6)*za(p3,t6)*zb(p1,t6)*zb(p2,t6)*zb(p4,t5)
     . /za(p1,p2)/za(p1,p2)
     . * (1d0/s34/s126/s26)
 
      gg_a(1,2,2,2,1)=
     . +za(p1,t6)*za(p2,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,p2)*zb(p2,p4)
     . /za(p1,p2)
     . * (- mb**3/s34/s126/s16/s26/s25)
 
     . +za(p1,t6)*za(p2,t6)*za(p2,p3)*za(p2,t5)*zb(p1,p2)*zb(p2,p4)
     . *zb(p2,t6)
     . /za(p1,p2)/za(p1,p2)
     . * (- mb/s34/s126/s26/s25)
 
     . +za(p1,t6)*za(p2,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     . *zb(p2,t6)
     . /za(p1,p2)/za(p1,p2)
     . * (mb/s34/s126/s26/s25)
 
      gg_a(1,2,2,1,2)=
     . +za(p1,t6)*za(p2,p3)*zb(p1,t6)*zb(p2,t6)*zb(p4,t5)
     . /za(p1,p2)/za(p1,p2)
     . * (- mb/s34/s126/s26)
 
     . +za(p3,t6)*zb(p1,t6)*zb(p2,t6)*zb(p4,t5)/za(p1,p2)
     . * (- mb/s34/s126/s26)
 

