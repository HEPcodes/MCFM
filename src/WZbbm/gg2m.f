 
      gg_g(2,2,1,2,2)=czip
 
      gg_g(2,2,1,1,1)=czip
 
      gg_g(2,2,1,2,1)=czip
 
      gg_g(2,2,1,1,2)=czip
 
      gg_g(2,1,2,2,2)=czip
 
      gg_g(2,1,2,1,1)=czip
 
      gg_g(2,1,2,2,1)=czip
 
      gg_g(2,1,2,1,2)=czip
 
      gg_g(1,2,1,2,2)=czip
 
      gg_g(1,2,1,1,1)=czip
 
      gg_g(1,2,1,2,1)=czip
 
      gg_g(1,2,1,1,2)=czip
 
      gg_g(1,1,2,2,2)=czip
 
      gg_g(1,1,2,1,1)=czip
 
      gg_g(1,1,2,2,1)=czip
 
      gg_g(1,1,2,1,2)=czip
 
      gg_g(2,1,1,2,2)=
     . +za(p1,t6)*za(p3,t5)*zb(p1,t6)*zb(p4,t6)/zb(p1,p2)/zb(p1,p2)
     . * (-1d0/s34/s126)
 
     . +za(p3,t5)*za(p1,p2)*zb(p1,t6)*zb(p2,p4)/zb(p1,p2)/zb(p1,p2)
     . * (1d0/s34/s126)
 
      gg_g(2,1,1,1,1)=+ za(p1,t6)*za(p2,p3)*zb(p1,p4)*zb(p2,t5)
     . /zb(p1,p2)/zb(p1,p2)
     . * (mb**2/s34/s126/s25)
 
      gg_g(2,1,1,2,1)=+za(p1,t6)*za(p3,t5)*zb(p1,p4)/zb(p1,p2)/zb(p1,p2)
     . * (mb/s34/s126)
 
      gg_g(2,1,1,1,2)=
     . +za(p1,t6)*za(p2,p3)*zb(p1,t6)*zb(p2,t5)*zb(p4,t6)
     . /zb(p1,p2)/zb(p1,p2)
     . * (- mb/s34/s126/s25)
 
     . +za(p2,p3)*za(p1,p2)*zb(p1,t6)*zb(p2,p4)*zb(p2,t5)
     . /zb(p1,p2)/zb(p1,p2)
     . * (mb/s34/s126/s25)
 
      gg_g(2,2,2,2,2)=+ za(p1,t6)*za(p3,t5)*zb(p1,t6)*zb(p4,t6)
     . /za(p1,p2)/za(p1,p2)
     . * (-1d0/s34/s126)
 
     . +za(p3,t5)*zb(p1,t6)*zb(p2,p4)/za(p1,p2)
     . * (1d0/s34/s126)
 
      gg_g(2,2,2,1,1)=+ za(p1,t6)*za(p2,p3)*zb(p1,p4)*zb(p2,t5)
     . /za(p1,p2)/za(p1,p2)
     . * (mb**2/s34/s126/s25)
 
      gg_g(2,2,2,2,1)=+za(p1,t6)*za(p3,t5)*zb(p1,p4)/za(p1,p2)/za(p1,p2)
     . * (mb/s34/s126)
 
      gg_g(2,2,2,1,2)=
     . +za(p1,t6)*za(p2,p3)*zb(p1,t6)*zb(p2,t5)*zb(p4,t6)
     . /za(p1,p2)/za(p1,p2)
     . * (- mb/s34/s126/s25)
 
     . +za(p2,p3)*zb(p1,t6)*zb(p2,p4)*zb(p2,t5)/za(p1,p2)
     . * (mb/s34/s126/s25)
 
      gg_g(1,1,1,2,2)=+ za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     . /zb(p1,p2)/zb(p1,p2)
     . * (mb**2/s34/s126/s25)
 
      gg_g(1,1,1,1,1)=+ za(p1,t6)*za(p2,p3)*zb(p4,t5)/zb(p1,p2)
     . * (1d0/s34/s126)
 
     . +za(p1,t6)*za(p3,t6)*zb(p1,t6)*zb(p4,t5)/zb(p1,p2)/zb(p1,p2)
     . * (-1d0/s34/s126)
 
      gg_g(1,1,1,2,1)=
     . +za(p1,t6)*za(p2,p3)*za(p2,t5)*zb(p2,p4)/zb(p1,p2)
     . * (mb/s34/s126/s25)
 
     . +za(p1,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     . /zb(p1,p2)/zb(p1,p2)
     . * (- mb/s34/s126/s25)
 
      gg_g(1,1,1,1,2)=+ za(p1,p3)*zb(p1,t6)*zb(p4,t5)
     . /zb(p1,p2)/zb(p1,p2)
     . * (mb/s34/s126)
 
      gg_g(1,2,2,2,2)=+ za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     . /za(p1,p2)/za(p1,p2)
     . * (mb**2/s34/s126/s25)
 
      gg_g(1,2,2,1,1)=
     . +za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p4,t5)/za(p1,p2)/za(p1,p2)
     . * (1d0/s34/s126)
 
     . +za(p1,t6)*za(p3,t6)*zb(p1,t6)*zb(p4,t5)/za(p1,p2)/za(p1,p2)
     . * (-1d0/s34/s126)
 
      gg_g(1,2,2,2,1)=
     . +za(p1,t6)*za(p2,p3)*za(p2,t5)*zb(p1,p2)*zb(p2,p4)
     . /za(p1,p2)/za(p1,p2)
     . * (mb/s34/s126/s25)
 
     . +za(p1,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     . /za(p1,p2)/za(p1,p2)
     . * (- mb/s34/s126/s25)
 
      gg_g(1,2,2,1,2)=+ za(p1,p3)*zb(p1,t6)*zb(p4,t5)
     . /za(p1,p2)/za(p1,p2)
     . * (mb/s34/s126)
 

