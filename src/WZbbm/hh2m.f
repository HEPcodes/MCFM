 
      gg_h(2,2,1,2,2)=czip
 
      gg_h(2,2,1,1,1)=czip
 
      gg_h(2,2,1,2,1)=czip
 
      gg_h(2,2,1,1,2)=czip
 
      gg_h(2,1,2,2,2)=czip
 
      gg_h(2,1,2,1,1)=czip
 
      gg_h(2,1,2,2,1)=czip
 
      gg_h(2,1,2,1,2)=czip
 
      gg_h(1,2,1,2,2)=czip
 
      gg_h(1,2,1,1,1)=czip
 
      gg_h(1,2,1,2,1)=czip
 
      gg_h(1,2,1,1,2)=czip
 
      gg_h(1,1,2,2,2)=czip
 
      gg_h(1,1,2,1,1)=czip
 
      gg_h(1,1,2,2,1)=czip
 
      gg_h(1,1,2,1,2)=czip
 
      gg_h(2,1,1,2,2)=
     . +za(p1,p3)*za(p2,t5)*zb(p4,t6)/zb(p1,p2)
     . * (mb**2/s34/s125/s25)
 
     . +za(p2,p3)*za(p1,t5)*zb(p4,t6)/zb(p1,p2)
     . * (- mb**2/s34/s125/s25 -1d0/s34/s125)
 
     . +za(p3,t5)*za(p1,t5)*zb(p1,t5)*zb(p4,t6)/zb(p1,p2)/zb(p1,p2)
     . * (1d0/s34/s125)
 
      gg_h(2,1,1,1,1)=
     . +za(p1,t6)*za(p1,p3)*zb(p1,p4)*zb(p1,t5)/zb(p1,p2)/zb(p1,p2)
     . * (- mb**2/s34/s125/s16)
 
     . +za(p1,t6)*za(p2,p3)*za(p1,p2)*zb(p1,p4)*zb(p2,t5)/zb(p1,p2)
     . * (mb**2/s34/s125/s16/s25 + mb**4/s34/s125/s16
     . /s25**2)
 
     . +za(p1,t6)*za(p3,t5)*za(p1,p2)*zb(p1,p4)*zb(p1,t5)*zb(p2,t5)
     . /zb(p1,p2)/zb(p1,p2)
     . * (- mb**2/s34/s125/s16/s25)
 
      gg_h(2,1,1,2,1)=
     . +za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p1,p4)/zb(p1,p2)
     . * (mb**3/s34/s125/s16/s25)
 
     . +za(p1,t6)*za(p2,p3)*za(p1,t5)*zb(p1,p4)/zb(p1,p2)
     . * (- mb/s34/s125/s16 - mb**3/s34/s125/s16/s25)
 
     . +za(p1,t6)*za(p3,t5)*za(p1,t5)*zb(p1,p4)*zb(p1,t5)
     . /zb(p1,p2)/zb(p1,p2)
     . * (mb/s34/s125/s16)
 
      gg_h(2,1,1,1,2)=
     . +za(p1,p3)*zb(p1,t5)*zb(p4,t6)/zb(p1,p2)/zb(p1,p2)
     . * (- mb/s34/s125)
 
     . +za(p2,p3)*za(p1,p2)*zb(p2,t5)*zb(p4,t6)/zb(p1,p2)
     . * (mb/s34/s125/s25 + mb**3/s34/s125/s25**2)
 
     . +za(p3,t5)*za(p1,p2)*zb(p1,t5)*zb(p2,t5)*zb(p4,t6)
     . /zb(p1,p2)/zb(p1,p2)
     . * (- mb/s34/s125/s25)
 
      gg_h(2,2,2,2,2)=
     . +za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p4,t6)/za(p1,p2)/za(p1,p2)
     . * (mb**2/s34/s125/s25)
 
     . +za(p2,p3)*za(p1,t5)*zb(p1,p2)*zb(p4,t6)/za(p1,p2)/za(p1,p2)
     . * (- mb**2/s34/s125/s25 -1d0/s34/s125)
 
     . +za(p3,t5)*za(p1,t5)*zb(p1,t5)*zb(p4,t6)/za(p1,p2)/za(p1,p2)
     . * (1d0/s34/s125)
 
      gg_h(2,2,2,1,1)=
     . +za(p1,t6)*za(p1,p3)*zb(p1,p4)*zb(p1,t5)/za(p1,p2)/za(p1,p2)
     . * (- mb**2/s34/s125/s16)
 
     . +za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p1,p4)*zb(p2,t5)/za(p1,p2)
     . * (mb**2/s34/s125/s16/s25 + mb**4/s34/s125/s16
     . /s25**2)
 
     . +za(p1,t6)*za(p3,t5)*zb(p1,p4)*zb(p1,t5)*zb(p2,t5)/za(p1,p2)
     . * (- mb**2/s34/s125/s16/s25)
 
      gg_h(2,2,2,2,1)=
     . +za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,p4)
     . /za(p1,p2)/za(p1,p2)
     . * (mb**3/s34/s125/s16/s25)
 
     . +za(p1,t6)*za(p2,p3)*za(p1,t5)*zb(p1,p2)*zb(p1,p4)
     . /za(p1,p2)/za(p1,p2)
     . * (- mb/s34/s125/s16 - mb**3/s34/s125/s16/s25)
 
     . +za(p1,t6)*za(p3,t5)*za(p1,t5)*zb(p1,p4)*zb(p1,t5)
     . /za(p1,p2)/za(p1,p2)
     . * (mb/s34/s125/s16)
 
      gg_h(2,2,2,1,2)=
     . +za(p1,p3)*zb(p1,t5)*zb(p4,t6)/za(p1,p2)/za(p1,p2)
     . * (- mb/s34/s125)
 
     . +za(p2,p3)*zb(p1,p2)*zb(p2,t5)*zb(p4,t6)/za(p1,p2)
     . * (mb/s34/s125/s25 + mb**3/s34/s125/s25**2)
 
     . +za(p3,t5)*zb(p1,t5)*zb(p2,t5)*zb(p4,t6)/za(p1,p2)
     . * (- mb/s34/s125/s25)
 
      gg_h(1,1,1,2,2)=
     . +za(p1,p3)*za(p1,t5)*za(p2,t5)*zb(p1,t6)*zb(p4,t5)/zb(p1,p2)
     . * (- mb**2/s34/s125/s16/s25)
 
     . +za(p1,p3)*za(p1,t5)*zb(p1,p4)*zb(p1,t6)/zb(p1,p2)/zb(p1,p2)
     . * (- mb**2/s34/s125/s16)
 
     . +za(p1,p3)*za(p2,t5)*za(p1,p2)*zb(p1,t6)*zb(p2,p4)/zb(p1,p2)
     . * (mb**2/s34/s125/s16/s25 + mb**4/s34/s125/s16
     . /s25**2)
 
      gg_h(1,1,1,1,1)=
     . +za(p3,t6)*za(p1,p2)*zb(p1,p4)*zb(p2,t5)/zb(p1,p2)/zb(p1,p2)
     . * (mb**2/s34/s125/s25)
 
     . +za(p3,t6)*za(p1,p2)*zb(p1,t5)*zb(p2,p4)/zb(p1,p2)/zb(p1,p2)
     . * (- mb**2/s34/s125/s25 -1d0/s34/s125)
 
     . +za(p3,t6)*za(p1,t5)*zb(p1,t5)*zb(p4,t5)/zb(p1,p2)/zb(p1,p2)
     . * (1d0/s34/s125)
 
      gg_h(1,1,1,2,1)=
     . +za(p3,t6)*za(p1,t5)*za(p2,t5)*zb(p4,t5)/zb(p1,p2)
     . * (- mb/s34/s125/s25)
 
     . +za(p3,t6)*za(p1,t5)*zb(p1,p4)/zb(p1,p2)/zb(p1,p2)
     . * (- mb/s34/s125)
 
     . +za(p3,t6)*za(p2,t5)*za(p1,p2)*zb(p2,p4)/zb(p1,p2)
     . * (mb/s34/s125/s25 + mb**3/s34/s125/s25**2)
 
      gg_h(1,1,1,1,2)=
     . +za(p1,p3)*za(p1,p2)*zb(p1,p4)*zb(p1,t6)*zb(p2,t5)
     . /zb(p1,p2)/zb(p1,p2)
     . * (mb**3/s34/s125/s16/s25)
 
     . +za(p1,p3)*za(p1,p2)*zb(p1,t5)*zb(p1,t6)*zb(p2,p4)
     . /zb(p1,p2)/zb(p1,p2)
     . * (- mb/s34/s125/s16 - mb**3/s34/s125/s16/s25)
 
     . +za(p1,p3)*za(p1,t5)*zb(p1,t5)*zb(p1,t6)*zb(p4,t5)
     . /zb(p1,p2)/zb(p1,p2)
     . * (mb/s34/s125/s16)
 
      gg_h(1,2,2,2,2)=
     . +za(p1,p3)*za(p1,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,t6)*zb(p4,t5)
     . /za(p1,p2)/za(p1,p2)
     . * (- mb**2/s34/s125/s16/s25)
 
     . +za(p1,p3)*za(p1,t5)*zb(p1,p4)*zb(p1,t6)/za(p1,p2)/za(p1,p2)
     . * (- mb**2/s34/s125/s16)
 
     . +za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,t6)*zb(p2,p4)/za(p1,p2)
     . * (mb**2/s34/s125/s16/s25 + mb**4/s34/s125/s16
     . /s25**2)
 
      gg_h(1,2,2,1,1)=
     . +za(p3,t6)*za(p1,t5)*zb(p1,t5)*zb(p4,t5)/za(p1,p2)/za(p1,p2)
     . * (1d0/s34/s125)
 
     . +za(p3,t6)*zb(p1,p4)*zb(p2,t5)/za(p1,p2)
     . * (mb**2/s34/s125/s25)
 
     . +za(p3,t6)*zb(p1,t5)*zb(p2,p4)/za(p1,p2)
     . * (- mb**2/s34/s125/s25 -1d0/s34/s125)
 
      gg_h(1,2,2,2,1)=
     . +za(p3,t6)*za(p1,t5)*za(p2,t5)*zb(p1,p2)*zb(p4,t5)
     . /za(p1,p2)/za(p1,p2)
     . * (- mb/s34/s125/s25)
 
     . +za(p3,t6)*za(p1,t5)*zb(p1,p4)/za(p1,p2)/za(p1,p2)
     . * (- mb/s34/s125)
 
     . +za(p3,t6)*za(p2,t5)*zb(p1,p2)*zb(p2,p4)/za(p1,p2)
     . * (mb/s34/s125/s25 + mb**3/s34/s125/s25**2)
 
      gg_h(1,2,2,1,2)=
     . +za(p1,p3)*za(p1,t5)*zb(p1,t5)*zb(p1,t6)*zb(p4,t5)
     . /za(p1,p2)/za(p1,p2)
     . * (mb/s34/s125/s16)
 
     . +za(p1,p3)*zb(p1,p4)*zb(p1,t6)*zb(p2,t5)/za(p1,p2)
     . * (mb**3/s34/s125/s16/s25)
 
     . +za(p1,p3)*zb(p1,t5)*zb(p1,t6)*zb(p2,p4)/za(p1,p2)
     . * (- mb/s34/s125/s16 - mb**3/s34/s125/s16/s25)
 

