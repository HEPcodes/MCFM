      subroutine h4qg(P1,P2,P3,P4,P5,za,zb,msq)
      implicit none
c---Matrix element squared for
c     q(-p1)+q(-p2) --> +H+q(p3)+q(p4)+g(p5) 
C     
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'

      integer P1,P2,P3,P4,P5,h1,h2,h5
      double precision s123,s124,s135,s134,s245,s234,S1234,msq
      double precision s12,s14,s15,s23,s25,s34,s35,s45
      double complex XTOTAL,HDPART,TLPART
      double complex a31(2,2,2),a32(2,2,2),a41(2,2,2),a42(2,2,2)

      S123=s(p1,p2)+s(p1,p3)+s(p2,p3)
      S124=s(p1,p2)+s(p1,p4)+s(p2,p4)
      S134=s(p1,p3)+s(p3,p4)+s(p1,p4)
      S135=s(p1,p3)+s(p3,p5)+s(p1,p5)
      S234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      S245=s(p2,p4)+s(p2,p5)+s(p4,p5)

      S1234=s(p1,p2)+s(p1,p3)+s(p1,p4)+s(p2,p3)+s(p2,p4)+s(p3,p4)
      s12=s(p1,p2)
      s14=s(p1,p4)
      s15=s(p1,p5)
      s23=s(p2,p3)
      s25=s(p2,p5)
      s34=s(p3,p4)
      s35=s(p3,p5)
      s45=s(p4,p5)
C Expression ppp31
        XTOTAL=0.
      HDPART=za(P1,P2)*za(P1,P2)*za(P1,P3)*zb(P1,P4)*zb(P2,P3)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*za(P2,P5)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P5)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P5)*zb(P3,P5)
     1 /za(P2,P4)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P3,P5)/za(P1,P5)/za(P2,P4)
     1 *S124**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P1,P3)*za(P1,P4)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P3)*za(P2,P3)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P3)*za(P2,P5)*zb(P2,P3)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P3)*za(P2,P5)*zb(P3,P5)*zb(P3,P5)
     1 /za(P1,P5)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P3)*zb(P1,P4)*zb(P3,P5)*zb(P3,P5)
     1 /za(P2,P4)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P3)*zb(P3,P4)/za(P1,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P2,P5)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P2,P4)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*zb(P2,P5)*zb(P3,P4)*zb(P3,P5)
     1 /za(P2,P4)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*zb(P2,P5)*zb(P3,P4)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P2,P5)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*zb(P2,P5)*zb(P3,P4)*zb(P4,P5)/zb(P2,P4)
     1 *S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*zb(P4,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)
     1 *S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P3)*za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P3)*za(P1,P4)*za(P2,P5)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P3)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)*zb(P3,P5)
     1 /za(P2,P4)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P3)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)/zb(P2,P4)
     1 *S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P2,P4)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*zb(P4,P5)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*zb(P3,P4)*zb(P4,P5)*zb(P4,P5)/zb(P2,P4)
     1 *S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
       a31(2,2,2)=XTOTAL
C Punched 40 terms out of 40.
C Expression ppm31
        XTOTAL=0.
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P3)*zb(P1,P4)*zb(P2,P3)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P3,P5)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P3)*zb(P1,P4)*zb(P3,P4)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P3,P5)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P1,P4)
     1 *zb(P3,P5)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P3,P5)*zb(P1,P3)*zb(P1,P4)
     1 *zb(P3,P5)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P4,P5)*zb(P1,P4)*zb(P1,P4)
     1 *zb(P3,P5)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P3)*zb(P1,P4)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S124**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P1,P3)*zb(P2,P3)*zb(P3,P4)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P3,P5)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*za(P2,P5)*zb(P1,P2)*zb(P3,P5)
     1 /za(P2,P4)/zb(P1,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*za(P3,P5)*zb(P1,P3)*zb(P3,P5)
     1 /za(P2,P4)/zb(P1,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*za(P4,P5)*zb(P1,P4)*zb(P3,P5)
     1 /za(P2,P4)/zb(P1,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*zb(P1,P3)*zb(P2,P3)*zb(P4,P5)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P3,P5)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P1,P3)/za(P2,P4)/zb(P1,P5)
     1 *S124**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P2,P5)*zb(P1,P4)*zb(P2,P3)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P4,P5)*zb(P1,P4)*zb(P3,P4)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*zb(P1,P3)*zb(P3,P4)/za(P2,P4)/zb(P1,P5)
     1 /zb(P2,P4)/zb(P3,P5)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P4)*za(P2,P3)*zb(P1,P3)*zb(P3,P4)*zb(P3,P4)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P3,P5)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*zb(P1,P3)*zb(P3,P4)*zb(P4,P5)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P3,P5)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P3,P5)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*za(P3,P5)*zb(P1,P3)*zb(P3,P4)
     1 *zb(P3,P5)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*za(P4,P5)*zb(P1,P4)*zb(P3,P4)
     1 *zb(P3,P5)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*zb(P1,P3)*zb(P3,P4)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)*zb(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P3,P5)*zb(P1,P3)*zb(P3,P4)*zb(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P4,P5)*zb(P1,P4)*zb(P3,P4)*zb(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*zb(P1,P4)*zb(P3,P4)/zb(P1,P5)/zb(P2,P4)
     1 *S234**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*za(P2,P5)*zb(P2,P3)*zb(P3,P4)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P4,P5)*zb(P3,P4)*zb(P3,P4)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P5)*za(P2,P5)*zb(P2,P3)*zb(P4,P5)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*za(P4,P5)*zb(P3,P4)*zb(P4,P5)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*zb(P3,P4)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)
     1 *S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      a31(2,2,1)=XTOTAL
C Punched 52 terms out of 52.
C Expression ppp32
        XTOTAL=0.
      HDPART=za(P1,P2)*za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P3)
     1 *zb(P2,P5)/za(P1,P3)/za(P1,P5)/za(P2,P4)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P3)
     1 /za(P1,P3)/za(P1,P5)/za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*za(P1,P3)*zb(P1,P4)*zb(P2,P3)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P2,P3)
     1 *zb(P4,P5)/za(P1,P3)/za(P1,P5)/za(P2,P4)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P2,P5)
     1 *zb(P3,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P3)/za(P1,P5)/za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P2,P5)
     1 *zb(P3,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P3)/za(P1,P5)/za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*za(P2,P5)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P5)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P3)*zb(P3,P5)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P5)*zb(P3,P5)
     1 /za(P2,P4)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P2,P3)*zb(P4,P5)/za(P1,P3)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P2,P3)*zb(P4,P5)/za(P1,P3)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= -S15 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P2,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P2,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P3,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P3,P5)/za(P1,P5)/za(P2,P4)
     1 *S124**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P1,P3)*za(P1,P4)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P3)*za(P2,P3)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P3)*za(P2,P5)*zb(P2,P3)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P3)*za(P2,P5)*zb(P3,P5)*zb(P3,P5)
     1 /za(P1,P5)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P3)*zb(P1,P4)*zb(P3,P5)*zb(P3,P5)
     1 /za(P2,P4)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P3)*zb(P3,P4)/za(P1,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P1,P4)*zb(P1,P4)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P1,P5)/za(P2,P4)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P2,P3)*zb(P2,P3)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P1,P5)/za(P2,P4)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*za(P2,P3)*zb(P2,P5)*zb(P3,P4)
     1 *zb(P3,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P4)
     1 /za(P1,P3)/za(P1,P5)/za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P2,P5)*zb(P2,P5)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P1,P5)/zb(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*za(P2,P5)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P2,P4)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P3,P4)*zb(P4,P5)/za(P1,P3)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P3,P4)*zb(P4,P5)/za(P1,P3)
     1 /za(P1,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P2,P5)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P2,P5)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P2,P5)*zb(P3,P4)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*zb(P2,P3)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*zb(P2,P3)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*zb(P2,P5)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*zb(P2,P5)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*zb(P2,P5)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*zb(P3,P4)*zb(P3,P5)/za(P1,P5)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*zb(P3,P4)*zb(P3,P5)/za(P1,P5)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*zb(P3,P4)*zb(P3,P5)/za(P1,P5)/zb(P1,P3)
     1 *S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*zb(P3,P5)*zb(P4,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*zb(P3,P5)*zb(P4,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*zb(P3,P5)*zb(P4,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*zb(P4,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)
     1 *S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P3)*za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*za(P1,P4)*za(P2,P5)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P3)*za(P2,P5)*zb(P3,P4)*zb(P3,P5)*zb(P3,P5)
     1 /za(P1,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P1,P5)/za(P2,P4)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P4)*za(P2,P5)*zb(P3,P4)*zb(P4,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P1,P5)/zb(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-2*za(P1,P4)*za(P2,P5)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P5)*zb(P3,P5)*zb(P4,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P5)*zb(P4,P5)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*zb(P3,P5)*zb(P3,P5)*zb(P4,P5)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
       a32(2,2,2)=XTOTAL
C Punched 118 terms out of 118.
C Expression ppm32
        XTOTAL=0.
      HDPART=-za(P1,P2)*za(P1,P2)*za(P2,P5)*zb(P1,P2)*zb(P1,P4)
     1 *zb(P2,P3)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*za(P3,P5)*zb(P1,P4)*zb(P2,P3)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*za(P4,P5)*zb(P1,P4)*zb(P1,P4)
     1 *zb(P2,P3)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P2)*zb(P1,P4)*zb(P2,P3)
     1 /za(P1,P3)/zb(P1,P3)/zb(P1,P5)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P1,P3)*zb(P1,P4)*zb(P2,P3)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P3,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*za(P2,P5)*zb(P1,P2)*zb(P1,P4)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P3,P5)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*za(P4,P5)*zb(P1,P4)*zb(P1,P4)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P2)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P3)/zb(P1,P3)/zb(P1,P5)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P3)*zb(P1,P4)*zb(P3,P4)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P3,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P1,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P1,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P1,P4)
     1 *zb(P3,P5)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P3,P5)*zb(P1,P3)*zb(P1,P4)
     1 *zb(P3,P5)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P4,P5)*zb(P1,P4)*zb(P1,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P4,P5)*zb(P1,P4)*zb(P1,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P4,P5)*zb(P1,P4)*zb(P1,P4)
     1 *zb(P3,P5)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P2)*zb(P1,P4)*zb(P3,P5)
     1 /za(P1,P3)/zb(P1,P3)/zb(P1,P5)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P3)*zb(P1,P4)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S124**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P1,P5)*zb(P1,P4)*zb(P2,P3)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P4)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P4)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
      TLPART= S25 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P2,P3)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*za(P3,P5)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*za(P4,P5)*zb(P1,P4)*zb(P2,P3)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P1,P2)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P3)/zb(P1,P3)/zb(P1,P5)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*zb(P1,P3)*zb(P2,P3)*zb(P3,P4)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P3,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*za(P2,P5)*zb(P1,P2)*zb(P3,P5)
     1 /za(P2,P4)/zb(P1,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*za(P3,P5)*zb(P1,P3)*zb(P3,P5)
     1 /za(P2,P4)/zb(P1,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*za(P3,P5)*zb(P2,P3)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*za(P3,P5)*zb(P2,P3)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*za(P4,P5)*zb(P1,P4)*zb(P3,P5)
     1 /za(P2,P4)/zb(P1,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P1,P3)*zb(P2,P3)*zb(P4,P5)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P3,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P1,P3)/za(P2,P4)/zb(P1,P5)
     1 *S124**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P1,P4)*zb(P2,P3)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P1,P4)*zb(P2,P3)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)*S245**(-1)
      TLPART= S15 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P3,P5)*zb(P2,P3)*zb(P3,P4)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P5)*zb(P3,P4)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P3,P5)*zb(P3,P4)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P4,P5)*zb(P1,P4)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P4,P5)*zb(P1,P4)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*zb(P1,P2)*zb(P3,P4)/za(P1,P3)/zb(P1,P3)
     1 /zb(P1,P5)/zb(P2,P5)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*zb(P1,P3)*zb(P3,P4)/za(P2,P4)/zb(P1,P5)
     1 /zb(P2,P4)/zb(P3,P5)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P1,P5)*zb(P1,P4)*zb(P3,P4)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*za(P3,P5)*zb(P3,P4)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P4,P5)*zb(P1,P4)*zb(P3,P4)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*zb(P1,P2)*zb(P3,P4)*zb(P3,P4)
     1 /za(P1,P3)/zb(P1,P3)/zb(P1,P5)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*zb(P1,P3)*zb(P3,P4)*zb(P3,P4)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P3,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/zb(P1,P3)/zb(P1,P5)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P5)*za(P3,P5)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P5)*za(P3,P5)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*za(P3,P5)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P3)/zb(P1,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*za(P4,P5)*zb(P1,P4)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/zb(P1,P3)/zb(P1,P5)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P5)*zb(P1,P3)*zb(P3,P4)*zb(P4,P5)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P3,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P5)*zb(P1,P4)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*zb(P1,P4)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)*S245**(-1)
      TLPART= -S15 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P2,P5)*zb(P1,P4)*zb(P3,P4)/za(P1,P3)
     1 /zb(P1,P3)/zb(P1,P5)*S134**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P3,P5)*zb(P3,P4)*zb(P3,P4)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P1,P5)*zb(P1,P4)*zb(P3,P5)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*za(P4,P5)*zb(P1,P4)*zb(P3,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*za(P4,P5)*zb(P1,P4)*zb(P3,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*zb(P1,P2)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P3)/zb(P1,P3)/zb(P1,P5)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*zb(P3,P4)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*zb(P3,P4)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
      TLPART= -S25 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P2,P5)*za(P2,P5)*zb(P1,P2)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P5)*za(P3,P5)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P5)*za(P4,P5)*zb(P1,P4)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P5)*zb(P1,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P5)*zb(P1,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P2,P5)*zb(P1,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P5)*za(P2,P5)*zb(P4,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P5)*zb(P4,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P2,P5)*zb(P4,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P3,P5)*zb(P3,P4)*zb(P3,P5)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*zb(P3,P4)/za(P1,P3)/zb(P1,P3)/zb(P2,P5)
     1 *S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P5)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)*zb(P3,P5)
     1 /zb(P1,P3)/zb(P1,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*za(P3,P5)*zb(P3,P4)*zb(P3,P5)/zb(P1,P5)
     1 *S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*za(P4,P5)*zb(P1,P4)*zb(P3,P4)*zb(P3,P5)
     1 /zb(P1,P3)/zb(P1,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*zb(P3,P4)/zb(P1,P5)*S134**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      a32(2,2,1)=XTOTAL
C Punched 146 terms out of 146.

*end
C Expression pmp31
        XTOTAL=0.
      HDPART=za(P1,P2)*za(P1,P3)*za(P1,P4)*zb(P1,P2)*zb(P2,P3)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P3)*za(P3,P4)*zb(P2,P3)*zb(P2,P3)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P3)*za(P4,P5)*zb(P2,P3)*zb(P2,P5)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P2)*zb(P2,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P4,P5)*zb(P2,P5)*zb(P2,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P3)*za(P1,P4)*za(P1,P4)*zb(P1,P2)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*za(P1,P4)*za(P1,P5)*zb(P1,P2)*zb(P1,P5)
     1 *zb(P3,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*za(P1,P4)*za(P3,P4)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P3)*za(P1,P4)*za(P4,P5)*zb(P1,P5)*zb(P3,P5)
     1 /za(P2,P4)/za(P3,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P3)*za(P1,P4)*za(P4,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*za(P1,P4)*zb(P1,P2)*zb(P3,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S124**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P3)*za(P1,P4)*zb(P2,P3)/za(P1,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P3)*za(P1,P5)*za(P3,P4)*zb(P1,P5)*zb(P2,P3)
     1 *zb(P3,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*za(P1,P5)*zb(P1,P5)*zb(P2,P3)*zb(P2,P5)
     1 /za(P3,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*za(P3,P4)*zb(P2,P3)*zb(P3,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P3)*zb(P2,P3)*zb(P2,P5)/za(P3,P5)/zb(P2,P4)
     1 *S234**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P1,P4)*zb(P1,P2)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P5)*za(P2,P3)*zb(P1,P2)*zb(P2,P5)
     1 *zb(P3,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*za(P3,P4)*zb(P1,P2)*zb(P3,P5)
     1 *zb(P4,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P4,P5)*zb(P2,P5)*zb(P3,P5)
     1 /za(P2,P4)/za(P3,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P3,P4)*za(P4,P5)*zb(P3,P5)*zb(P4,P5)
     1 /za(P2,P4)/za(P3,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P3,P4)*zb(P2,P3)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P3,P5)/za(P2,P4)/za(P3,P5)
     1 *S124**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P4,P5)*zb(P2,P5)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*zb(P2,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)
     1 *S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P2,P3)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)
     1 *zb(P3,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*zb(P2,P3)*zb(P2,P5)*zb(P2,P5)
     1 /za(P3,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P3,P4)*za(P3,P4)*zb(P2,P3)*zb(P3,P5)
     1 *zb(P4,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)*zb(P4,P5)
     1 /za(P3,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      a31(2,1,2)=XTOTAL
C Punched 52 terms out of 52.
C Expression mpp31
        XTOTAL=0.
      HDPART=za(P1,P2)*za(P1,P3)*za(P2,P3)*zb(P1,P2)*zb(P1,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P3)*za(P3,P4)*zb(P1,P4)*zb(P1,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P3)*zb(P1,P4)*zb(P1,P5)*zb(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P1,P4)*zb(P1,P5)*zb(P2,P5)
     1 /za(P2,P4)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P1,P4)*zb(P2,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)*zb(P4,P5)
     1 /za(P2,P4)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P1,P4)*zb(P4,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P3)*za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*za(P2,P3)*za(P2,P5)*zb(P1,P5)*zb(P1,P5)
     1 /za(P2,P4)/za(P3,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*za(P2,P3)*zb(P1,P4)/za(P1,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P3)*za(P2,P3)*zb(P1,P5)*zb(P1,P5)*zb(P3,P4)
     1 /za(P2,P4)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P3)*za(P2,P5)*za(P3,P4)*zb(P1,P4)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P3)*zb(P1,P4)*zb(P1,P5)*zb(P4,P5)/zb(P2,P4)
     1 *S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P3)*za(P2,P5)*zb(P1,P5)*zb(P2,P5)
     1 /za(P2,P4)/za(P3,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P3)*zb(P1,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P2,P4)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P3)*zb(P1,P5)/za(P2,P4)/za(P3,P5)
     1 *S234**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P3)*za(P2,P3)*zb(P2,P5)*zb(P3,P4)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P5)*za(P3,P4)*zb(P1,P5)*zb(P4,P5)
     1 /za(P2,P4)/za(P3,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P5)*zb(P2,P5)*zb(P4,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P3,P4)*zb(P1,P5)*zb(P3,P4)*zb(P4,P5)
     1 /za(P2,P4)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P3,P4)*zb(P4,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*zb(P1,P4)*zb(P2,P5)*zb(P4,P5)/zb(P2,P4)
     1 *S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*zb(P4,P5)/za(P1,P5)/za(P2,P4)/zb(P2,P4)
     1 *S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P5)*za(P3,P4)*zb(P4,P5)*zb(P4,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P3,P4)*zb(P1,P4)*zb(P4,P5)*zb(P4,P5)/zb(P2,P4)
     1 *S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      a31(1,2,2)=XTOTAL
C Punched 40 terms out of 40.
C Expression pmp32
        XTOTAL=0.
      HDPART=-za(P1,P2)*za(P1,P4)*za(P2,P3)*zb(P1,P2)*zb(P2,P3)
     1 *zb(P2,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P2,P3)*zb(P1,P2)*zb(P2,P3)
     1 /za(P1,P3)/za(P2,P5)/za(P3,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*za(P3,P4)*zb(P1,P2)*zb(P2,P3)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P2)*zb(P1,P5)*zb(P2,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P2)*zb(P2,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P5)*zb(P2,P3)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*za(P3,P4)*zb(P2,P3)*zb(P2,P3)
     1 *zb(P2,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*za(P3,P4)*zb(P2,P3)*zb(P2,P3)
     1 /za(P1,P3)/za(P2,P5)/za(P3,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*za(P3,P4)*zb(P2,P3)*zb(P2,P3)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P1,P5)*zb(P2,P3)*zb(P2,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= -S15 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P3,P4)*zb(P2,P3)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P4,P5)*zb(P1,P5)*zb(P2,P3)*zb(P2,P5)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P4,P5)*zb(P1,P5)*zb(P2,P3)*zb(P2,P5)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P4,P5)*zb(P2,P5)*zb(P2,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P3)*za(P1,P4)*za(P1,P5)*zb(P1,P2)*zb(P1,P5)
     1 *zb(P3,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*za(P1,P4)*za(P4,P5)*zb(P1,P5)*zb(P3,P5)
     1 /za(P2,P4)/za(P3,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*za(P1,P4)*zb(P1,P2)*zb(P3,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S124**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P1,P4)*za(P2,P3)*zb(P1,P2)*zb(P2,P5)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P4)*za(P2,P3)*zb(P1,P2)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P3,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P4)*za(P3,P4)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P4)*zb(P1,P2)*zb(P1,P5)*zb(P3,P4)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P4)*zb(P1,P2)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P4)*zb(P1,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P5)*za(P2,P3)*zb(P1,P2)*zb(P2,P5)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P5)*za(P2,P3)*zb(P1,P2)*zb(P2,P5)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*za(P2,P3)*zb(P1,P2)*zb(P2,P5)
     1 *zb(P3,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P5)*za(P2,P3)*zb(P1,P2)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P5)/za(P3,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*za(P2,P3)*zb(P2,P5)*zb(P2,P5)
     1 /za(P1,P3)/za(P3,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*za(P3,P4)*zb(P1,P2)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*za(P3,P4)*zb(P1,P2)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P5)*za(P3,P4)*zb(P1,P2)*zb(P3,P5)
     1 *zb(P4,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P5)*za(P3,P4)*zb(P2,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P3,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*zb(P1,P5)*zb(P2,P5)/za(P3,P5)
     1 *S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*zb(P1,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P3,P4)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P3,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*za(P4,P5)*zb(P2,P5)*zb(P2,P5)
     1 *zb(P3,P4)/za(P1,P3)/za(P3,P5)/zb(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*za(P4,P5)*zb(P2,P5)*zb(P3,P5)
     1 /za(P2,P4)/za(P3,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*zb(P2,P3)*zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*zb(P2,P3)*zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P2,P3)*zb(P2,P3)/za(P1,P3)/za(P2,P5)
     1 /za(P3,P5)/zb(P1,P3)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P3,P4)*za(P3,P4)*zb(P2,P3)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P3,P4)*za(P4,P5)*zb(P2,P5)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P3,P5)/zb(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P3,P4)*za(P4,P5)*zb(P3,P5)*zb(P4,P5)
     1 /za(P2,P4)/za(P3,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P1,P5)*zb(P2,P3)*zb(P3,P4)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P2,P3)*zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P2,P3)*zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P2,P3)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P2,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P2,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= S15 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P2,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P3,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P4)*za(P3,P4)*zb(P3,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P3,P5)/za(P2,P4)/za(P3,P5)
     1 *S124**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P4)*za(P4,P5)*zb(P1,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P4,P5)*zb(P1,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P4,P5)*zb(P1,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P3,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P4,P5)*zb(P2,P5)*zb(P4,P5)/za(P2,P4)
     1 /za(P3,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*zb(P1,P2)*zb(P3,P5)/za(P2,P4)/za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*zb(P1,P2)*zb(P3,P5)/za(P2,P4)/za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= S25 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*zb(P1,P5)*zb(P2,P3)/za(P2,P4)/za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*zb(P1,P5)*zb(P2,P3)/za(P2,P4)/za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*zb(P2,P5)/za(P2,P4)/za(P3,P5)/zb(P2,P4)
     1 *S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*zb(P2,P5)/za(P3,P5)*S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*zb(P3,P5)/za(P1,P3)/za(P2,P5)/zb(P1,P3)
     1 *S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P5)*za(P2,P3)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*za(P3,P4)*zb(P2,P3)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P5)/za(P3,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*za(P4,P5)*zb(P2,P5)*zb(P2,P5)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P3,P4)*za(P3,P4)*zb(P2,P3)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P3,P4)*za(P3,P4)*zb(P2,P3)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P3,P4)*za(P4,P5)*zb(P2,P5)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P3,P4)*zb(P2,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P3,P4)*zb(P2,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P3,P4)*zb(P2,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P5)*za(P3,P4)*zb(P3,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P4,P5)*zb(P1,P5)*zb(P2,P5)*zb(P3,P5)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P3,P4)*zb(P2,P3)*zb(P3,P5)/za(P2,P4)/za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P3,P4)*zb(P2,P3)*zb(P3,P5)/za(P2,P4)/za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S25 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P4,P5)*zb(P2,P5)*zb(P3,P5)/za(P2,P4)/za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P4,P5)*zb(P2,P5)*zb(P3,P5)/za(P2,P4)/za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P4,P5)*zb(P2,P5)*zb(P3,P5)/za(P2,P4)/za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
        a32(2,1,2)=XTOTAL
C Punched 146 terms out of 146.
C Expression mpp32
        XTOTAL=0.
      HDPART=-za(P1,P2)*za(P1,P3)*zb(P1,P4)*zb(P1,P5)*zb(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P1,P4)
     1 *zb(P2,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P1,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P3,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*za(P3,P4)*zb(P1,P2)*zb(P1,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P1,P4)
     1 *zb(P2,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P1,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P3,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P1,P2)*zb(P1,P4)*zb(P1,P5)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P1,P4)*zb(P1,P5)*zb(P2,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P1,P4)*zb(P1,P5)*zb(P2,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*zb(P1,P4)*zb(P1,P5)*zb(P2,P5)
     1 /za(P2,P4)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P1,P4)*zb(P1,P5)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*za(P3,P4)*zb(P1,P4)*zb(P1,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P1,P4)*zb(P1,P4)*zb(P1,P5)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)*zb(P4,P5)
     1 /za(P2,P4)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P3)*za(P2,P5)*zb(P1,P4)*zb(P1,P5)*zb(P1,P5)
     1 /za(P3,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P3)*zb(P1,P4)*zb(P1,P5)*zb(P4,P5)/zb(P2,P4)
     1 *S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P2,P5)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P3,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P3)*za(P3,P4)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P2,P5)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P3,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P1,P5)*zb(P3,P4)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= S15 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P3)*za(P2,P3)*zb(P1,P4)*zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P3)*zb(P1,P4)*zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*za(P2,P3)*zb(P1,P4)/za(P1,P3)/za(P2,P5)
     1 /za(P3,P5)/zb(P1,P3)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P3)*za(P2,P3)*zb(P1,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P3)*zb(P1,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P3)*zb(P1,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P5)*za(P3,P4)*zb(P1,P4)*zb(P2,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P3,P5)/zb(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P1,P5)*zb(P4,P5)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P1,P5)*zb(P4,P5)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P5)*zb(P1,P4)*zb(P1,P5)*zb(P2,P5)
     1 /za(P3,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P5)*zb(P1,P5)*zb(P2,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*za(P3,P4)*zb(P1,P4)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)*zb(P3,P4)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P4,P5)/za(P1,P3)
     1 /za(P3,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P1,P5)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P1,P5)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*zb(P1,P4)*zb(P1,P5)/za(P2,P4)/za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*zb(P1,P4)*zb(P1,P5)/za(P2,P4)/za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*zb(P1,P4)*zb(P1,P5)/za(P3,P5)/zb(P1,P3)
     1 *S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*zb(P1,P4)*zb(P2,P5)*zb(P4,P5)/zb(P2,P4)
     1 *S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*zb(P1,P5)*zb(P4,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*zb(P1,P5)*zb(P4,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*zb(P1,P5)*zb(P4,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P5)*za(P3,P4)*za(P3,P4)*zb(P1,P4)*zb(P4,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P3,P5)/zb(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P5)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)*zb(P4,P5)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P5)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)*zb(P4,P5)
     1 /za(P2,P4)/za(P3,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=2*za(P2,P5)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)*zb(P4,P5)
     1 /za(P3,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P5)*za(P3,P4)*zb(P1,P5)*zb(P4,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P5)*zb(P1,P5)*zb(P1,P5)*zb(P4,P5)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P3,P4)*zb(P1,P4)*zb(P4,P5)*zb(P4,P5)/zb(P2,P4)
     1 *S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
       a32(1,2,2)=XTOTAL
C Punched 90 terms out of 90.

*end
C Expression ppp42
        XTOTAL=0.
      HDPART=-za(P1,P2)*za(P1,P2)*za(P1,P5)*zb(P1,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*za(P2,P4)*zb(P1,P4)*zb(P2,P3)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P5)*zb(P2,P3)*zb(P4,P5)
     1 /za(P1,P3)/zb(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P1,P5)*zb(P2,P3)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P4,P5)/za(P1,P3)/za(P2,P5)
     1 *S123**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P2,P4)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P5)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P3)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P2,P3)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P2,P4)*zb(P1,P4)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P2,P4)*zb(P4,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*zb(P1,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*za(P2,P4)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/zb(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P4)*zb(P2,P3)*zb(P4,P5)*zb(P4,P5)
     1 /za(P1,P3)/zb(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P4)*zb(P3,P4)/za(P1,P3)/za(P2,P5)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*zb(P1,P5)*zb(P3,P4)*zb(P3,P5)/zb(P1,P3)
     1 *S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*zb(P3,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)
     1 *S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P2,P3)*za(P2,P4)*zb(P3,P4)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P4)*zb(P3,P4)*zb(P4,P5)*zb(P4,P5)
     1 /za(P1,P3)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*za(P2,P4)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*zb(P3,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*zb(P3,P4)*zb(P3,P5)*zb(P3,P5)/zb(P1,P3)
     1 *S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P4)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)/zb(P1,P3)
     1 *S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      a42(2,2,2)=XTOTAL
C Punched 40 terms out of 40.
C Expression ppm42
        XTOTAL=0.
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P3)*zb(P2,P4)
     1 /za(P1,P3)/zb(P1,P3)/zb(P2,P5)/zb(P4,P5)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P2,P4)*zb(P3,P4)
     1 /za(P1,P3)/zb(P1,P3)/zb(P2,P5)/zb(P4,P5)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P1,P5)*zb(P1,P2)*zb(P4,P5)
     1 /za(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P2,P3)
     1 *zb(P4,P5)/za(P1,P3)/zb(P1,P3)/zb(P2,P5)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P3,P5)*zb(P2,P3)*zb(P4,P5)
     1 /za(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P4,P5)*zb(P2,P4)*zb(P4,P5)
     1 /za(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P4)*zb(P2,P3)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P4)*zb(P2,P4)*zb(P3,P5)
     1 /za(P1,P3)/zb(P1,P3)/zb(P2,P5)/zb(P4,P5)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P2,P4)/za(P1,P3)/zb(P2,P5)
     1 *S123**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P2,P4)*zb(P3,P4)
     1 /za(P1,P3)/zb(P1,P3)/zb(P2,P5)/zb(P4,P5)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*za(P3,P5)*zb(P2,P3)*zb(P2,P3)
     1 *zb(P4,P5)/za(P1,P3)/zb(P1,P3)/zb(P2,P5)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*za(P4,P5)*zb(P2,P3)*zb(P2,P4)
     1 *zb(P4,P5)/za(P1,P3)/zb(P1,P3)/zb(P2,P5)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P2,P3)*zb(P2,P4)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P3,P5)*zb(P2,P3)*zb(P3,P4)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*zb(P2,P4)*zb(P3,P4)/za(P1,P3)/zb(P1,P3)
     1 /zb(P2,P5)/zb(P4,P5)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/zb(P1,P3)/zb(P2,P5)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*zb(P1,P4)*zb(P3,P4)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*zb(P2,P4)*zb(P3,P4)*zb(P3,P4)
     1 /za(P1,P3)/zb(P1,P3)/zb(P2,P5)/zb(P4,P5)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P5)*za(P3,P5)*zb(P2,P3)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/zb(P1,P3)/zb(P2,P5)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P5)*za(P4,P5)*zb(P2,P4)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/zb(P1,P3)/zb(P2,P5)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*zb(P2,P4)*zb(P3,P4)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P4)*za(P3,P5)*zb(P3,P4)*zb(P3,P4)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P1,P5)*zb(P1,P4)*zb(P3,P5)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*zb(P2,P4)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P3)/zb(P1,P3)/zb(P2,P5)/zb(P4,P5)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)*zb(P3,P5)
     1 /zb(P1,P3)/zb(P2,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P3,P5)*zb(P3,P4)*zb(P3,P5)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*zb(P3,P4)/za(P1,P3)/zb(P1,P3)/zb(P2,P5)
     1 *S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P5)*za(P3,P5)*zb(P2,P3)*zb(P3,P4)*zb(P3,P5)
     1 /zb(P1,P3)/zb(P2,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P5)*za(P4,P5)*zb(P2,P4)*zb(P3,P4)*zb(P3,P5)
     1 /zb(P1,P3)/zb(P2,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*zb(P2,P3)*zb(P3,P4)/zb(P1,P3)/zb(P2,P5)
     1 *S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      a42(2,2,1)=XTOTAL
C Punched 52 terms out of 52.
C Expression ppp41
        XTOTAL=0.
      HDPART=-za(P1,P2)*za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P1,P5)
     1 *zb(P2,P3)/za(P1,P3)/za(P2,P4)/za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P3)
     1 /za(P1,P5)/za(P2,P4)/za(P2,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P1,P5)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P2,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*za(P1,P5)*zb(P1,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*za(P2,P3)*zb(P1,P4)*zb(P2,P3)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*za(P2,P3)*zb(P1,P5)*zb(P2,P3)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P2,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*za(P2,P4)*zb(P1,P4)*zb(P2,P3)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P3)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S25 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P5)*zb(P2,P3)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P5)*zb(P2,P3)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P1,P5)*zb(P2,P3)*zb(P4,P5)
     1 /za(P1,P3)/zb(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P1,P5)*zb(P2,P3)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P1,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P1,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P1,P2)*zb(P2,P3)*zb(P4,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P3,P4)/za(P1,P5)/za(P2,P4)
     1 /za(P2,P5)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P4,P5)/za(P1,P3)/za(P2,P5)
     1 *S123**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P2,P3)*zb(P1,P4)*zb(P3,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P2,P3)*zb(P1,P5)*zb(P3,P4)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P2,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*za(P2,P4)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P5)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P5)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P3,P4)*zb(P4,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P2,P3)*zb(P1,P5)*zb(P3,P4)
     1 *zb(P3,P5)/za(P2,P4)/za(P2,P5)/zb(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P2,P3)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P2,P4)*zb(P1,P4)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P2,P4)*zb(P4,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*zb(P1,P5)*zb(P3,P4)*zb(P4,P5)
     1 /za(P2,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P5)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*za(P2,P3)*zb(P2,P3)*zb(P3,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*za(P2,P4)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/zb(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)/za(P2,P4)
     1 /za(P2,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P2,P4)*zb(P2,P3)*zb(P4,P5)*zb(P4,P5)
     1 /za(P1,P3)/zb(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P4)*zb(P3,P4)/za(P1,P3)/za(P2,P5)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*zb(P3,P4)*zb(P4,P5)/za(P1,P3)/za(P2,P5)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*zb(P3,P4)*zb(P4,P5)/za(P1,P3)/za(P2,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*zb(P3,P4)*zb(P4,P5)/za(P2,P5)/zb(P2,P4)
     1 *S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*zb(P3,P5)*zb(P4,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*zb(P3,P5)*zb(P4,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*zb(P3,P5)*zb(P4,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*zb(P3,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)
     1 *S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P2,P3)*za(P2,P3)*zb(P3,P4)*zb(P3,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P2,P4)*zb(P3,P4)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P4)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)
     1 *zb(P3,P5)/za(P2,P4)/za(P2,P5)/zb(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*za(P2,P4)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=2*za(P1,P5)*za(P2,P3)*zb(P3,P4)*zb(P3,P5)*zb(P4,P5)
     1 /za(P2,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*zb(P3,P5)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*zb(P3,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P4)*zb(P3,P4)*zb(P4,P5)*zb(P4,P5)
     1 /za(P2,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*zb(P3,P5)*zb(P4,P5)*zb(P4,P5)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      a41(2,2,2)=XTOTAL
C Punched 118 terms out of 118.
C Expression ppm41
        XTOTAL=0.
      HDPART=za(P1,P2)*za(P1,P2)*za(P1,P5)*zb(P1,P2)*zb(P1,P4)
     1 *zb(P2,P3)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P2)*za(P3,P5)*zb(P1,P4)*zb(P2,P3)
     1 *zb(P2,P3)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*za(P4,P5)*zb(P1,P4)*zb(P2,P3)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P1,P2)*zb(P1,P4)*zb(P2,P3)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P2)*zb(P1,P4)*zb(P2,P3)*zb(P2,P4)
     1 /za(P1,P3)/zb(P1,P3)/zb(P2,P5)/zb(P4,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P1,P5)*zb(P1,P2)*zb(P1,P4)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*za(P3,P5)*zb(P1,P4)*zb(P2,P3)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P4,P5)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P2)*zb(P1,P4)*zb(P3,P4)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P4)*zb(P2,P4)*zb(P3,P4)
     1 /za(P1,P3)/zb(P1,P3)/zb(P2,P5)/zb(P4,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P1,P5)*zb(P1,P2)*zb(P4,P5)
     1 /za(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P2,P3)*zb(P1,P2)*zb(P2,P3)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P2,P3)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P2,P3)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P2,P3)
     1 *zb(P4,P5)/za(P1,P3)/zb(P1,P3)/zb(P2,P5)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P3,P5)*zb(P2,P3)*zb(P4,P5)
     1 /za(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P4,P5)*zb(P1,P4)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*za(P4,P5)*zb(P1,P4)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*za(P4,P5)*zb(P2,P4)*zb(P4,P5)
     1 /za(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P5)*zb(P1,P2)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P2)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P1,P5)*zb(P1,P4)*zb(P2,P3)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P1,P4)*zb(P2,P3)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)*S135**(-1)
      TLPART= -S25 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P1,P5)*zb(P1,P4)*zb(P2,P4)*zb(P3,P5)
     1 /za(P1,P3)/zb(P1,P3)/zb(P2,P5)/zb(P4,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P5)*zb(P2,P4)/za(P1,P3)/zb(P2,P5)
     1 *S123**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P2,P3)*za(P3,P5)*zb(P2,P3)*zb(P2,P3)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*za(P4,P5)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*zb(P1,P2)*zb(P2,P3)*zb(P3,P4)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*zb(P2,P3)*zb(P2,P4)*zb(P3,P4)
     1 /za(P1,P3)/zb(P1,P3)/zb(P2,P5)/zb(P4,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*za(P3,P5)*zb(P2,P3)*zb(P2,P3)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P5)*za(P3,P5)*zb(P2,P3)*zb(P2,P3)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*za(P3,P5)*zb(P2,P3)*zb(P2,P3)
     1 *zb(P4,P5)/za(P1,P3)/zb(P1,P3)/zb(P2,P5)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*za(P4,P5)*zb(P2,P3)*zb(P2,P4)
     1 *zb(P4,P5)/za(P1,P3)/zb(P1,P3)/zb(P2,P5)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P1,P2)*zb(P2,P3)*zb(P4,P5)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P1,P4)*zb(P2,P3)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P2,P3)*zb(P2,P4)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P2,P5)*zb(P2,P3)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P5)*zb(P2,P3)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
      TLPART= -S15 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P3,P5)*zb(P2,P3)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P3,P5)*zb(P2,P3)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P2)*za(P4,P5)*zb(P1,P4)*zb(P3,P4)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P4,P5)*zb(P3,P4)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P4,P5)*zb(P3,P4)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*zb(P1,P2)*zb(P3,P4)/za(P2,P4)/zb(P1,P5)
     1 /zb(P2,P4)/zb(P2,P5)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*zb(P2,P4)*zb(P3,P4)/za(P1,P3)/zb(P1,P3)
     1 /zb(P2,P5)/zb(P4,P5)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P1,P5)*za(P2,P3)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P3,P5)*zb(P2,P3)*zb(P3,P4)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*za(P4,P5)*zb(P3,P4)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*zb(P1,P2)*zb(P3,P4)*zb(P3,P4)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*zb(P2,P4)*zb(P3,P4)*zb(P3,P4)
     1 /za(P1,P3)/zb(P1,P3)/zb(P2,P5)/zb(P4,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*za(P3,P5)*zb(P2,P3)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*za(P3,P5)*zb(P2,P3)*zb(P3,P4)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P5)*zb(P1,P2)*zb(P3,P4)*zb(P4,P5)
     1 /za(P2,P4)/zb(P1,P5)/zb(P2,P4)/zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*zb(P3,P4)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*zb(P3,P4)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
      TLPART= S15 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P1,P5)*za(P2,P3)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P3,P5)/za(P2,P4)/zb(P2,P4)/zb(P2,P5)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P1,P5)*za(P2,P5)*zb(P1,P2)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P1,P5)*zb(P1,P2)*zb(P3,P4)*zb(P4,P5)
     1 /zb(P2,P4)/zb(P2,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*za(P3,P5)*zb(P2,P3)*zb(P3,P4)
     1 *zb(P3,P5)/za(P2,P4)/zb(P2,P4)/zb(P2,P5)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*za(P4,P5)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*za(P4,P5)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P3)*za(P4,P5)*zb(P3,P4)*zb(P3,P5)
     1 /za(P2,P4)/zb(P2,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*zb(P2,P3)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P3)*zb(P2,P3)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)*S135**(-1)
      TLPART= S25 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P2,P3)*zb(P2,P3)*zb(P3,P4)/za(P2,P4)
     1 /zb(P2,P4)/zb(P2,P5)*S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P2,P3)*zb(P2,P4)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P3)/zb(P1,P3)/zb(P2,P5)/zb(P4,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P5)*za(P3,P5)*zb(P2,P3)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)
     1 /zb(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P2,P5)*za(P4,P5)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P5)*zb(P2,P3)*zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P5)*zb(P2,P3)*zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P2,P5)*zb(P2,P3)*zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)*S135**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P5)*za(P2,P5)*zb(P3,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P2,P5)*zb(P3,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P5)*S245**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P2,P5)*zb(P3,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P5)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P5)*za(P3,P5)*zb(P2,P3)*zb(P3,P4)*zb(P4,P5)
     1 /zb(P2,P4)/zb(P2,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P4,P5)*zb(P3,P4)*zb(P4,P5)/zb(P2,P5)
     1 *S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*zb(P3,P4)/zb(P2,P5)*S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P3)*za(P2,P5)*zb(P2,P3)*zb(P3,P4)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P4,P5)*zb(P3,P4)*zb(P3,P4)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*za(P2,P5)*zb(P2,P3)*zb(P4,P5)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P5)*za(P4,P5)*zb(P3,P4)*zb(P4,P5)/za(P2,P4)
     1 /zb(P1,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*zb(P3,P4)/za(P2,P4)/zb(P1,P5)/zb(P2,P4)
     1 *S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      a41(2,2,1)=XTOTAL
C Punched 146 terms out of 146.

*end
C Expression pmp42
        XTOTAL=0.
      HDPART=za(P1,P2)*za(P1,P4)*za(P2,P4)*zb(P1,P2)*zb(P2,P3)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P5)*zb(P2,P3)*zb(P2,P5)
     1 /za(P1,P3)/zb(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P5)*zb(P2,P3)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P4)*za(P3,P4)*zb(P2,P3)*zb(P2,P3)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P4)*zb(P2,P3)*zb(P2,P5)*zb(P2,P5)
     1 /za(P1,P3)/zb(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P3)/zb(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P2,P3)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P4)*za(P1,P5)*zb(P1,P5)*zb(P2,P5)
     1 /za(P1,P3)/za(P4,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P4)*za(P2,P4)*zb(P1,P2)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P4)*zb(P1,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P3)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P4)*zb(P1,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P4)*zb(P2,P5)/za(P1,P3)/za(P4,P5)
     1 *S134**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P4)*za(P1,P5)*za(P2,P4)*zb(P1,P2)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P5)*za(P2,P4)*zb(P2,P5)*zb(P2,P5)
     1 /za(P1,P3)/za(P4,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P5)*za(P3,P4)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P3)/za(P4,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P5)*zb(P1,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P4)*za(P3,P4)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P4)*zb(P2,P3)/za(P1,P3)/za(P2,P5)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P4)*za(P2,P4)*zb(P2,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P3)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P3,P4)*zb(P2,P5)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P3)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P3,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*zb(P1,P5)*zb(P2,P3)*zb(P3,P5)/zb(P1,P3)
     1 *S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*zb(P3,P5)/za(P1,P3)/za(P2,P5)/zb(P1,P3)
     1 *S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P5)*za(P2,P4)*za(P3,P4)*zb(P2,P3)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P3,P4)*zb(P3,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P4)*zb(P2,P3)*zb(P2,P5)*zb(P3,P5)/zb(P1,P3)
     1 *S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P3,P4)*zb(P2,P3)*zb(P3,P5)*zb(P3,P5)/zb(P1,P3)
     1 *S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      a42(2,1,2)=XTOTAL
C Punched 40 terms out of 40.
C Expression mpp42
        XTOTAL=0.
      HDPART=za(P1,P2)*za(P2,P3)*za(P2,P4)*zb(P1,P2)*zb(P1,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P1,P2)*zb(P1,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P4)*za(P3,P4)*zb(P1,P4)*zb(P1,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P4)*za(P3,P5)*zb(P1,P4)*zb(P1,P5)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P5)*zb(P1,P5)*zb(P1,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P1,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P3,P5)*zb(P1,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P4,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P5)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*zb(P1,P4)*zb(P1,P5)*zb(P1,P5)
     1 /za(P4,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P3)*za(P2,P4)*zb(P1,P2)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P4)*za(P2,P5)*zb(P1,P2)*zb(P2,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P4)*za(P3,P4)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P4)*za(P3,P5)*zb(P1,P5)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P5)/za(P4,P5)/zb(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P4)*za(P3,P5)*zb(P2,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P4,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P4)*zb(P1,P2)*zb(P4,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S123**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*za(P2,P4)*zb(P1,P4)/za(P1,P3)/za(P2,P5)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P3)*za(P2,P5)*za(P3,P4)*zb(P1,P2)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P3,P4)*za(P3,P5)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P4,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P4,P5)/za(P1,P3)/za(P4,P5)
     1 *S123**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*za(P3,P5)*zb(P1,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*zb(P1,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)
     1 *S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P4)*za(P2,P5)*za(P3,P4)*zb(P1,P4)*zb(P2,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P4)*za(P2,P5)*zb(P1,P4)*zb(P1,P5)*zb(P2,P5)
     1 /za(P4,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P4)*za(P3,P4)*zb(P1,P4)*zb(P4,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P4)*zb(P1,P4)*zb(P1,P5)/za(P4,P5)/zb(P1,P3)
     1 *S134**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P5)*za(P3,P4)*za(P3,P4)*zb(P1,P4)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)*zb(P3,P5)
     1 /za(P4,P5)/zb(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      a42(1,2,2)=XTOTAL
C Punched 52 terms out of 52.
C Expression pmp41
        XTOTAL=0.
      HDPART=-za(P1,P2)*za(P1,P4)*za(P1,P4)*zb(P1,P2)*zb(P1,P5)
     1 *zb(P2,P3)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P1,P4)*zb(P1,P2)*zb(P2,P3)
     1 /za(P1,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P3,P4)*zb(P1,P2)*zb(P2,P3)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*za(P3,P4)*zb(P1,P5)*zb(P2,P3)
     1 *zb(P2,P3)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*za(P3,P4)*zb(P2,P3)*zb(P2,P3)
     1 /za(P1,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P2)*zb(P2,P3)*zb(P2,P5)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P5)*zb(P2,P3)*zb(P2,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P1,P5)*zb(P2,P3)*zb(P2,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*zb(P1,P5)*zb(P2,P3)*zb(P2,P5)
     1 /za(P1,P3)/zb(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*zb(P2,P3)*zb(P2,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P4)*zb(P2,P3)*zb(P2,P5)*zb(P2,P5)
     1 /za(P1,P3)/zb(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P3,P4)*za(P3,P4)*zb(P2,P3)*zb(P2,P3)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P2,P3)*zb(P2,P3)*zb(P2,P5)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P3)/zb(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P4)*za(P1,P4)*zb(P1,P2)*zb(P1,P5)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P4)*za(P1,P4)*zb(P1,P2)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P4)*za(P3,P4)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P4)*za(P3,P4)*zb(P1,P5)*zb(P2,P3)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P4)*za(P3,P4)*zb(P2,P3)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P4)*zb(P1,P2)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P4)*zb(P1,P2)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P4)*zb(P1,P2)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S25 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P4)*za(P1,P4)*zb(P1,P5)*zb(P2,P3)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P4)*zb(P1,P5)*zb(P2,P3)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P1,P4)*zb(P1,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P4)*zb(P1,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P4)*zb(P2,P3)/za(P1,P5)/za(P2,P4)
     1 /za(P4,P5)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P1,P4)*zb(P2,P5)*zb(P3,P4)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*za(P3,P4)*zb(P1,P5)*zb(P2,P3)
     1 *zb(P3,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*zb(P1,P2)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*zb(P1,P2)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P1,P5)*zb(P1,P5)*zb(P2,P3)*zb(P2,P5)
     1 /za(P4,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P1,P5)*zb(P1,P5)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P3,P4)*za(P3,P4)*zb(P2,P3)*zb(P3,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P2,P3)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P2,P3)*zb(P3,P5)/za(P2,P4)
     1 /za(P4,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P2,P5)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P3,P4)*zb(P2,P5)*zb(P3,P4)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*zb(P1,P5)*zb(P2,P3)*zb(P3,P5)/zb(P1,P3)
     1 *S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*zb(P2,P3)*zb(P2,P5)/za(P1,P3)/za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*zb(P2,P3)*zb(P2,P5)/za(P1,P3)/za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*zb(P2,P3)*zb(P2,P5)/za(P4,P5)/zb(P2,P4)
     1 *S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*zb(P2,P5)*zb(P3,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*zb(P2,P5)*zb(P3,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*zb(P2,P5)*zb(P3,P5)/za(P1,P3)/za(P2,P4)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P1,P5)*za(P2,P4)*zb(P2,P3)*zb(P2,P5)*zb(P2,P5)
     1 /za(P4,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*za(P3,P4)*za(P3,P4)*zb(P2,P3)*zb(P3,P5)
     1 *zb(P3,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)*zb(P3,P5)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-2*za(P1,P5)*za(P3,P4)*zb(P2,P3)*zb(P2,P5)*zb(P3,P5)
     1 /za(P4,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P5)*za(P3,P4)*zb(P2,P5)*zb(P3,P5)*zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P5)*zb(P2,P5)*zb(P2,P5)*zb(P3,P5)/za(P1,P3)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P4)*zb(P2,P3)*zb(P2,P5)*zb(P3,P5)/zb(P1,P3)
     1 *S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P3,P4)*zb(P2,P3)*zb(P3,P5)*zb(P3,P5)/zb(P1,P3)
     1 *S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      a41(2,1,2)=XTOTAL
C Punched 90 terms out of 90.
C Expression mpp41
        XTOTAL=0.
      HDPART=-za(P1,P2)*za(P1,P4)*za(P2,P3)*zb(P1,P2)*zb(P1,P4)
     1 *zb(P1,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P1,P4)*za(P2,P3)*zb(P1,P2)*zb(P1,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*za(P3,P4)*zb(P1,P4)*zb(P1,P4)
     1 *zb(P1,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P1,P4)*za(P3,P4)*zb(P1,P4)*zb(P1,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*za(P3,P4)*zb(P1,P2)*zb(P1,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P2,P3)*zb(P1,P2)*zb(P1,P4)*zb(P2,P5)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*zb(P1,P2)*zb(P1,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P2,P3)*zb(P1,P4)*zb(P2,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P3,P4)*za(P3,P4)*zb(P1,P4)*zb(P1,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P1,P4)*zb(P1,P4)*zb(P2,P5)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S25 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P2)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P3,P4)*zb(P1,P4)*zb(P4,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P5)*zb(P1,P4)*zb(P1,P5)*zb(P2,P5)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P2)*za(P3,P5)*zb(P1,P4)*zb(P1,P5)*zb(P2,P5)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P2)*za(P3,P5)*zb(P1,P5)*zb(P1,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P1,P5)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P1,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P1,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P1,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P2,P5)*zb(P1,P2)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P2,P5)*zb(P1,P5)*zb(P1,P5)
     1 /za(P2,P4)/za(P4,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)
     1 *zb(P3,P4)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P3,P4)
     1 /za(P1,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*za(P3,P5)*zb(P1,P5)*zb(P1,P5)
     1 *zb(P3,P4)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*za(P3,P5)*zb(P1,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P4,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P1,P4)*za(P2,P3)*zb(P1,P4)*zb(P1,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P3)*zb(P1,P4)*zb(P1,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P2,P3)*zb(P1,P4)/za(P1,P5)/za(P2,P4)
     1 /za(P4,P5)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P1,P4)*za(P2,P5)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*za(P3,P4)*zb(P1,P4)*zb(P1,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*za(P3,P4)*zb(P1,P4)*zb(P4,P5)
     1 /za(P1,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P1,P4)*za(P2,P5)*za(P3,P5)*zb(P1,P5)*zb(P1,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P3)*za(P3,P4)*zb(P1,P2)*zb(P3,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P3)*zb(P1,P2)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P3)*zb(P2,P5)*zb(P3,P4)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P4)*za(P2,P5)*zb(P1,P2)*zb(P2,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P4)*za(P3,P5)*zb(P2,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P4,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P4)*zb(P1,P2)*zb(P4,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S123**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P3)*za(P2,P5)*za(P3,P4)*zb(P1,P2)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P5)*za(P3,P4)*zb(P1,P2)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P2,P5)*za(P3,P4)*zb(P1,P2)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P5)*za(P3,P4)*zb(P1,P5)*zb(P3,P5)
     1 /za(P2,P4)/za(P4,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P5)*zb(P1,P5)*zb(P2,P5)/za(P4,P5)
     1 *S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P2,P5)*zb(P2,P5)*zb(P4,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P3,P4)*za(P3,P4)*zb(P1,P4)*zb(P3,P4)
     1 *zb(P3,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P3,P4)*za(P3,P5)*zb(P1,P5)*zb(P3,P4)
     1 *zb(P3,P5)/za(P2,P4)/za(P4,P5)/zb(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*za(P3,P5)*zb(P3,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P4,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P3)*za(P3,P4)*zb(P1,P4)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P1,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P1,P5)*zb(P3,P4)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= S25 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P1,P5)*zb(P3,P4)/za(P2,P4)
     1 /za(P4,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P2,P3)*za(P3,P4)*zb(P3,P4)*zb(P4,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P4)*zb(P4,P5)/za(P1,P3)/za(P4,P5)
     1 *S123**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*za(P3,P5)*zb(P1,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*za(P3,P5)*zb(P1,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P3,P5)*zb(P1,P5)*zb(P2,P5)*zb(P3,P4)
     1 /za(P4,P5)/zb(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*za(P3,P5)*zb(P1,P5)*zb(P3,P5)/za(P1,P3)
     1 /za(P4,P5)/zb(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P3)*zb(P1,P2)*zb(P4,P5)/za(P1,P3)/za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*zb(P1,P2)*zb(P4,P5)/za(P1,P3)/za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= -S15 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*zb(P1,P4)*zb(P2,P5)/za(P1,P3)/za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P3)*zb(P1,P4)*zb(P2,P5)/za(P1,P3)/za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*zb(P1,P5)/za(P1,P3)/za(P4,P5)/zb(P1,P3)
     1 *S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*zb(P1,P5)/za(P4,P5)*S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P3)*zb(P4,P5)/za(P1,P5)/za(P2,P4)/zb(P2,P4)
     1 *S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P5)*za(P3,P4)*za(P3,P4)*zb(P1,P4)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*za(P3,P4)*za(P3,P4)*zb(P1,P4)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*za(P3,P4)*za(P3,P5)*zb(P1,P5)*zb(P3,P5)
     1 *zb(P4,P5)/za(P1,P3)/za(P2,P4)/za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*za(P3,P4)*zb(P1,P5)*zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P2,P5)*za(P3,P4)*zb(P1,P5)*zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P5)*za(P3,P4)*zb(P1,P5)*zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P2,P5)*za(P3,P4)*zb(P4,P5)*zb(P4,P5)/za(P1,P5)
     1 /za(P2,P4)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-za(P2,P5)*za(P3,P5)*zb(P1,P5)*zb(P2,P5)*zb(P4,P5)
     1 /za(P1,P3)/za(P4,P5)/zb(P1,P3)/zb(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P3,P4)*zb(P1,P4)*zb(P4,P5)/za(P1,P3)/za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P3,P4)*zb(P1,P4)*zb(P4,P5)/za(P1,P3)/za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= S15 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-za(P3,P5)*zb(P1,P5)*zb(P4,P5)/za(P1,P3)/za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=za(P3,P5)*zb(P1,P5)*zb(P4,P5)/za(P1,P3)/za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=za(P3,P5)*zb(P1,P5)*zb(P4,P5)/za(P1,P3)/za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      a41(1,2,2)=XTOTAL
C Punched 146 terms out of 146.

*end

      write(6,*) 
      write(6,*) 'ppp31',a31(2,2,2)
      write(6,*) 'ppm31',a31(2,2,1)
      write(6,*) 'ppp32',a32(2,2,2)
      write(6,*) 'ppm32',a32(2,2,1)

      write(6,*) 
      write(6,*) 'pmp31',a31(2,1,2)
      write(6,*) 'mpp31',a31(1,2,2)
      write(6,*) 'pmp32',a32(2,1,2)
      write(6,*) 'mpp32',a32(1,2,2)

      write(6,*) 
      write(6,*) 'ppp42',a42(2,2,2)
      write(6,*) 'ppm42',a42(2,2,1)
      write(6,*) 'ppp41',a41(2,2,2)
      write(6,*) 'ppm41',a41(2,2,1)

      write(6,*) 
      write(6,*) 'pmp42',a42(2,1,2)
      write(6,*) 'mpp42',a42(1,2,2)
      write(6,*) 'pmp41',a41(2,1,2)
      write(6,*) 'mpp41',a41(1,2,2)

C Expression ppp31
        XTOTAL=0.
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P1,P3)*za(P1,P4)*za(P2,P3)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P2,P5)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P5)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P5)*za(P3,P5)
     1 /zb(P2,P4)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P3,P5)/zb(P1,P5)/zb(P2,P4)
     1 *S124**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P1,P3)*zb(P1,P4)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P3)*zb(P2,P3)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P3)*zb(P2,P5)*za(P2,P3)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P3)*zb(P2,P5)*za(P3,P5)*za(P3,P5)
     1 /zb(P1,P5)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P3)*za(P1,P4)*za(P3,P5)*za(P3,P5)
     1 /zb(P2,P4)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P3)*za(P3,P4)/zb(P1,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P2,P5)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P2,P4)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*za(P2,P5)*za(P3,P4)*za(P3,P5)
     1 /zb(P2,P4)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*za(P2,P5)*za(P3,P4)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P2,P5)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*za(P2,P5)*za(P3,P4)*za(P4,P5)/za(P2,P4)
     1 *S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*za(P4,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)
     1 *S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P3)*zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P3)*zb(P1,P4)*zb(P2,P5)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P3)*zb(P2,P3)*za(P3,P4)*za(P3,P5)*za(P3,P5)
     1 /zb(P2,P4)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P3)*za(P3,P4)*za(P3,P5)*za(P4,P5)/za(P2,P4)
     1 *S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P2,P4)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*za(P4,P5)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*za(P3,P4)*za(P4,P5)*za(P4,P5)/za(P2,P4)
     1 *S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
       a31(1,1,1)=XTOTAL
C Punched 40 terms out of 40.
C Expression ppm31
        XTOTAL=0.
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P3)*za(P1,P4)*za(P2,P3)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P3,P5)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P3)*za(P1,P4)*za(P3,P4)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P3,P5)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P1,P4)
     1 *za(P3,P5)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P3,P5)*za(P1,P3)*za(P1,P4)
     1 *za(P3,P5)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P4,P5)*za(P1,P4)*za(P1,P4)
     1 *za(P3,P5)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P3)*za(P1,P4)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S124**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P1,P3)*za(P2,P3)*za(P3,P4)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P3,P5)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*zb(P2,P5)*za(P1,P2)*za(P3,P5)
     1 /zb(P2,P4)/za(P1,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*zb(P3,P5)*za(P1,P3)*za(P3,P5)
     1 /zb(P2,P4)/za(P1,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*zb(P4,P5)*za(P1,P4)*za(P3,P5)
     1 /zb(P2,P4)/za(P1,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*za(P1,P3)*za(P2,P3)*za(P4,P5)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P3,P5)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P1,P3)/zb(P2,P4)/za(P1,P5)
     1 *S124**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*za(P1,P4)*za(P2,P3)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P4,P5)*za(P1,P4)*za(P3,P4)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*za(P1,P3)*za(P3,P4)/zb(P2,P4)/za(P1,P5)
     1 /za(P2,P4)/za(P3,P5)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*za(P1,P3)*za(P3,P4)*za(P3,P4)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P3,P5)*S135**(-1)
     
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*za(P1,P3)*za(P3,P4)*za(P4,P5)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P3,P5)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P3,P4)
     1 *za(P3,P5)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*zb(P3,P5)*za(P1,P3)*za(P3,P4)
     1 *za(P3,P5)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*zb(P4,P5)*za(P1,P4)*za(P3,P4)
     1 *za(P3,P5)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*za(P1,P3)*za(P3,P4)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P3,P4)*za(P4,P5)
     1 /za(P1,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P3,P5)*za(P1,P3)*za(P3,P4)*za(P4,P5)
     1 /za(P1,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P4,P5)*za(P1,P4)*za(P3,P4)*za(P4,P5)
     1 /za(P1,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*za(P1,P4)*za(P3,P4)/za(P1,P5)/za(P2,P4)
     1 *S234**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*zb(P2,P5)*za(P2,P3)*za(P3,P4)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P4,P5)*za(P3,P4)*za(P3,P4)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P5)*zb(P2,P5)*za(P2,P3)*za(P4,P5)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*zb(P4,P5)*za(P3,P4)*za(P4,P5)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*za(P3,P4)/zb(P2,P4)/za(P1,P5)/za(P2,P4)
     1 *S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      a31(1,1,2)=XTOTAL
C Punched 52 terms out of 52.
C Expression ppp32
        XTOTAL=0.
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P3)
     1 *za(P2,P5)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P3)
     1 /zb(P1,P3)/zb(P1,P5)/zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P1,P3)*za(P1,P4)*za(P2,P3)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P2,P3)
     1 *za(P4,P5)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P2,P5)
     1 *za(P3,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P3)/zb(P1,P5)/zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P2,P5)
     1 *za(P3,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P3)/zb(P1,P5)/zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P2,P5)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P5)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P3)*za(P3,P5)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P5)*za(P3,P5)
     1 /zb(P2,P4)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P2,P3)*za(P4,P5)/zb(P1,P3)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P2,P3)*za(P4,P5)/zb(P1,P3)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= -S15 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P2,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P2,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P3,P4)/zb(P1,P3)/zb(P1,P5)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P3,P5)/zb(P1,P5)/zb(P2,P4)
     1 *S124**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P1,P3)*zb(P1,P4)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P3)*zb(P2,P3)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P3)*zb(P2,P5)*za(P2,P3)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P3)*zb(P2,P5)*za(P3,P5)*za(P3,P5)
     1 /zb(P1,P5)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P3)*za(P1,P4)*za(P3,P5)*za(P3,P5)
     1 /zb(P2,P4)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P3)*za(P3,P4)/zb(P1,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P1,P4)*za(P1,P4)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P2,P3)*za(P2,P3)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P2,P3)*za(P2,P5)*za(P3,P4)
     1 *za(P3,P4)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P4)
     1 /zb(P1,P3)/zb(P1,P5)/zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P2,P5)*za(P2,P5)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P1,P5)/za(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P2,P5)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P2,P4)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P3,P4)*za(P4,P5)/zb(P1,P3)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P3,P4)*za(P4,P5)/zb(P1,P3)
     1 /zb(P1,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P2,P5)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P2,P5)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P2,P5)*za(P3,P4)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P3,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*za(P2,P3)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*za(P2,P3)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*za(P2,P5)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*za(P2,P5)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*za(P2,P5)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*za(P3,P4)*za(P3,P5)/zb(P1,P5)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*za(P3,P4)*za(P3,P5)/zb(P1,P5)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*za(P3,P4)*za(P3,P5)/zb(P1,P5)/za(P1,P3)
     1 *S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*za(P3,P5)*za(P4,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*za(P3,P5)*za(P4,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*za(P3,P5)*za(P4,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*za(P4,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)
     1 *S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P3)*zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*zb(P1,P4)*zb(P2,P5)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P3)*zb(P2,P5)*za(P3,P4)*za(P3,P5)*za(P3,P5)
     1 /zb(P1,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P1,P5)/zb(P2,P4)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P4)*zb(P2,P5)*za(P3,P4)*za(P4,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P1,P5)/za(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*za(P3,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*za(P3,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-2*zb(P1,P4)*zb(P2,P5)*za(P3,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P5)*za(P3,P5)*za(P4,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P5)*za(P4,P5)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*za(P3,P5)*za(P3,P5)*za(P4,P5)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
       a32(1,1,1)=XTOTAL
C Punched 118 terms out of 118.
C Expression ppm32
        XTOTAL=0.
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P2,P5)*za(P1,P2)*za(P1,P4)
     1 *za(P2,P3)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P3,P5)*za(P1,P4)*za(P2,P3)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P4,P5)*za(P1,P4)*za(P1,P4)
     1 *za(P2,P3)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P2)*za(P1,P4)*za(P2,P3)
     1 /zb(P1,P3)/za(P1,P3)/za(P1,P5)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P1,P3)*za(P1,P4)*za(P2,P3)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P3,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P2,P5)*za(P1,P2)*za(P1,P4)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P3,P5)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P4,P5)*za(P1,P4)*za(P1,P4)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P2)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P3)/za(P1,P3)/za(P1,P5)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P3)*za(P1,P4)*za(P3,P4)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P3,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P1,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P1,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P1,P4)
     1 *za(P3,P5)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P3,P5)*za(P1,P3)*za(P1,P4)
     1 *za(P3,P5)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P4,P5)*za(P1,P4)*za(P1,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P4,P5)*za(P1,P4)*za(P1,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P4,P5)*za(P1,P4)*za(P1,P4)
     1 *za(P3,P5)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P2)*za(P1,P4)*za(P3,P5)
     1 /zb(P1,P3)/za(P1,P3)/za(P1,P5)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P3)*za(P1,P4)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S124**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*za(P1,P4)*za(P2,P3)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P4)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P4)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
      TLPART= S25 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P2,P3)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*zb(P3,P5)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*zb(P4,P5)*za(P1,P4)*za(P2,P3)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P1,P2)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P3)/za(P1,P3)/za(P1,P5)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*za(P1,P3)*za(P2,P3)*za(P3,P4)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P3,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*zb(P2,P5)*za(P1,P2)*za(P3,P5)
     1 /zb(P2,P4)/za(P1,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*zb(P3,P5)*za(P1,P3)*za(P3,P5)
     1 /zb(P2,P4)/za(P1,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*zb(P3,P5)*za(P2,P3)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*zb(P3,P5)*za(P2,P3)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*zb(P4,P5)*za(P1,P4)*za(P3,P5)
     1 /zb(P2,P4)/za(P1,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P1,P2)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P1,P2)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P1,P3)*za(P2,P3)*za(P4,P5)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P3,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P1,P3)/zb(P2,P4)/za(P1,P5)
     1 *S124**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P1,P4)*za(P2,P3)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P1,P4)*za(P2,P3)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)*S245**(-1)
      TLPART= S15 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P3,P5)*za(P2,P3)*za(P3,P4)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P5)*za(P3,P4)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P3,P5)*za(P3,P4)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P4,P5)*za(P1,P4)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P4,P5)*za(P1,P4)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*za(P1,P2)*za(P3,P4)/zb(P1,P3)/za(P1,P3)
     1 /za(P1,P5)/za(P2,P5)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*za(P1,P3)*za(P3,P4)/zb(P2,P4)/za(P1,P5)
     1 /za(P2,P4)/za(P3,P5)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P1,P5)*za(P1,P4)*za(P3,P4)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P3,P4)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*zb(P3,P5)*za(P3,P4)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P4,P5)*za(P1,P4)*za(P3,P4)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*za(P1,P2)*za(P3,P4)*za(P3,P4)
     1 /zb(P1,P3)/za(P1,P3)/za(P1,P5)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*za(P1,P3)*za(P3,P4)*za(P3,P4)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P3,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*zb(P2,P5)*za(P1,P2)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/za(P1,P3)/za(P1,P5)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P5)*zb(P3,P5)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P5)*zb(P3,P5)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*zb(P3,P5)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P3)/za(P1,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*zb(P4,P5)*za(P1,P4)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/za(P1,P3)/za(P1,P5)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P5)*za(P1,P3)*za(P3,P4)*za(P4,P5)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P3,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P5)*za(P1,P4)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*za(P1,P4)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)*S245**(-1)
      TLPART= -S15 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P2,P5)*za(P1,P4)*za(P3,P4)/zb(P1,P3)
     1 /za(P1,P3)/za(P1,P5)*S134**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P3,P5)*za(P3,P4)*za(P3,P4)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P1,P5)*za(P1,P4)*za(P3,P5)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P3,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P3,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*zb(P4,P5)*za(P1,P4)*za(P3,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*zb(P4,P5)*za(P1,P4)*za(P3,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*za(P1,P2)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P3)/za(P1,P3)/za(P1,P5)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*za(P3,P4)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*za(P3,P4)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
      TLPART= -S25 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P2,P5)*zb(P2,P5)*za(P1,P2)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P5)*zb(P3,P5)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P5)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P5)*zb(P4,P5)*za(P1,P4)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P1,P5)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P5)*za(P1,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P5)*za(P1,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P2,P5)*za(P1,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P1,P5)/za(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P5)*zb(P2,P5)*za(P4,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P5)*za(P4,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P2,P5)*za(P4,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P3,P5)*za(P3,P4)*za(P3,P5)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*za(P3,P4)/zb(P1,P3)/za(P1,P3)/za(P2,P5)
     1 *S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P5)*zb(P2,P5)*za(P1,P2)*za(P3,P4)*za(P3,P5)
     1 /za(P1,P3)/za(P1,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*zb(P3,P5)*za(P3,P4)*za(P3,P5)/za(P1,P5)
     1 *S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*zb(P4,P5)*za(P1,P4)*za(P3,P4)*za(P3,P5)
     1 /za(P1,P3)/za(P1,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*za(P3,P4)/za(P1,P5)*S134**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      a32(1,1,2)=XTOTAL
C Punched 146 terms out of 146.

*end
C Expression pmp31
        XTOTAL=0.
      HDPART=zb(P1,P2)*zb(P1,P3)*zb(P1,P4)*za(P1,P2)*za(P2,P3)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P3)*zb(P3,P4)*za(P2,P3)*za(P2,P3)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P3)*zb(P4,P5)*za(P2,P3)*za(P2,P5)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P2)*za(P2,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P2,P3)*za(P2,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P4,P5)*za(P2,P5)*za(P2,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P3)*zb(P1,P4)*zb(P1,P4)*za(P1,P2)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*zb(P1,P4)*zb(P1,P5)*za(P1,P2)*za(P1,P5)
     1 *za(P3,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*zb(P1,P4)*zb(P3,P4)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P3)*zb(P1,P4)*zb(P4,P5)*za(P1,P5)*za(P3,P5)
     1 /zb(P2,P4)/zb(P3,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P3)*zb(P1,P4)*zb(P4,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*zb(P1,P4)*za(P1,P2)*za(P3,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S124**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P3)*zb(P1,P4)*za(P2,P3)/zb(P1,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P3)*zb(P1,P5)*zb(P3,P4)*za(P1,P5)*za(P2,P3)
     1 *za(P3,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*zb(P1,P5)*za(P1,P5)*za(P2,P3)*za(P2,P5)
     1 /zb(P3,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*zb(P3,P4)*za(P2,P3)*za(P3,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P3)*za(P2,P3)*za(P2,P5)/zb(P3,P5)/za(P2,P4)
     1 *S234**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P1,P4)*za(P1,P2)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P5)*zb(P2,P3)*za(P1,P2)*za(P2,P5)
     1 *za(P3,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*zb(P3,P4)*za(P1,P2)*za(P3,P5)
     1 *za(P4,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P4,P5)*za(P2,P5)*za(P3,P5)
     1 /zb(P2,P4)/zb(P3,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P3,P4)*zb(P4,P5)*za(P3,P5)*za(P4,P5)
     1 /zb(P2,P4)/zb(P3,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P3,P4)*za(P2,P3)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P3,P5)/zb(P2,P4)/zb(P3,P5)
     1 *S124**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P4,P5)*za(P2,P5)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*za(P2,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)
     1 *S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P2,P3)*zb(P3,P4)*za(P2,P3)*za(P2,P5)
     1 *za(P3,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*za(P2,P3)*za(P2,P5)*za(P2,P5)
     1 /zb(P3,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P3,P4)*zb(P3,P4)*za(P2,P3)*za(P3,P5)
     1 *za(P4,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P3,P4)*za(P2,P3)*za(P2,P5)*za(P4,P5)
     1 /zb(P3,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      a31(1,2,1)=XTOTAL
C Punched 52 terms out of 52.
C Expression mpp31
        XTOTAL=0.
      HDPART=zb(P1,P2)*zb(P1,P3)*zb(P2,P3)*za(P1,P2)*za(P1,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P3)*zb(P3,P4)*za(P1,P4)*za(P1,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P3)*za(P1,P4)*za(P1,P5)*za(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P1,P4)*za(P1,P5)*za(P2,P5)
     1 /zb(P2,P4)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P1,P4)*za(P2,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P1,P4)*za(P1,P5)*za(P4,P5)
     1 /zb(P2,P4)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P1,P4)*za(P4,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P3)*zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*zb(P2,P3)*zb(P2,P5)*za(P1,P5)*za(P1,P5)
     1 /zb(P2,P4)/zb(P3,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*zb(P2,P3)*za(P1,P4)/zb(P1,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P3)*zb(P2,P3)*za(P1,P5)*za(P1,P5)*za(P3,P4)
     1 /zb(P2,P4)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P3)*zb(P2,P5)*zb(P3,P4)*za(P1,P4)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S135**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P3)*za(P1,P4)*za(P1,P5)*za(P4,P5)/za(P2,P4)
     1 *S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P3)*zb(P2,P5)*za(P1,P5)*za(P2,P5)
     1 /zb(P2,P4)/zb(P3,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P3)*za(P1,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P2,P4)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P3)*za(P1,P5)/zb(P2,P4)/zb(P3,P5)
     1 *S234**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P3)*zb(P2,P3)*za(P2,P5)*za(P3,P4)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P5)*zb(P3,P4)*za(P1,P5)*za(P4,P5)
     1 /zb(P2,P4)/zb(P3,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P5)*za(P2,P5)*za(P4,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P3,P4)*za(P1,P5)*za(P3,P4)*za(P4,P5)
     1 /zb(P2,P4)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P3,P4)*za(P4,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*za(P1,P4)*za(P2,P5)*za(P4,P5)/za(P2,P4)
     1 *S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*za(P4,P5)/zb(P1,P5)/zb(P2,P4)/za(P2,P4)
     1 *S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P5)*zb(P3,P4)*za(P4,P5)*za(P4,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P3,P4)*za(P1,P4)*za(P4,P5)*za(P4,P5)/za(P2,P4)
     1 *S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      a31(2,1,1)=XTOTAL
C Punched 40 terms out of 40.
C Expression pmp32
        XTOTAL=0.
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P2,P3)*za(P1,P2)*za(P2,P3)
     1 *za(P2,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P2,P3)*za(P1,P2)*za(P2,P3)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P3,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P3,P4)*za(P1,P2)*za(P2,P3)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P2)*za(P1,P5)*za(P2,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P2)*za(P2,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P5)*za(P2,P3)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*zb(P3,P4)*za(P2,P3)*za(P2,P3)
     1 *za(P2,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*zb(P3,P4)*za(P2,P3)*za(P2,P3)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P3,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*zb(P3,P4)*za(P2,P3)*za(P2,P3)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P1,P5)*za(P2,P3)*za(P2,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P2,P3)*za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P3,P4)*za(P2,P3)*za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= -S15 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P3,P4)*za(P2,P3)*za(P2,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P3,P4)*za(P2,P3)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P4,P5)*za(P1,P5)*za(P2,P3)*za(P2,P5)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P4,P5)*za(P1,P5)*za(P2,P3)*za(P2,P5)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P4,P5)*za(P2,P5)*za(P2,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P3)*zb(P1,P4)*zb(P1,P5)*za(P1,P2)*za(P1,P5)
     1 *za(P3,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*zb(P1,P4)*zb(P4,P5)*za(P1,P5)*za(P3,P5)
     1 /zb(P2,P4)/zb(P3,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*zb(P1,P4)*za(P1,P2)*za(P3,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S124**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P1,P4)*zb(P2,P3)*za(P1,P2)*za(P2,P5)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P4)*zb(P2,P3)*za(P1,P2)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P3,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P4)*zb(P3,P4)*za(P1,P2)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P4)*za(P1,P2)*za(P1,P5)*za(P3,P4)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P4)*za(P1,P2)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P4)*za(P1,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P5)*zb(P2,P3)*za(P1,P2)*za(P2,P5)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P5)*zb(P2,P3)*za(P1,P2)*za(P2,P5)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*zb(P2,P3)*za(P1,P2)*za(P2,P5)
     1 *za(P3,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P5)*zb(P2,P3)*za(P1,P2)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P3,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*zb(P2,P3)*za(P2,P5)*za(P2,P5)
     1 /zb(P1,P3)/zb(P3,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*zb(P3,P4)*za(P1,P2)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*zb(P3,P4)*za(P1,P2)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P5)*zb(P3,P4)*za(P1,P2)*za(P3,P5)
     1 *za(P4,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)*S124**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P5)*zb(P3,P4)*za(P2,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P3,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*za(P1,P5)*za(P2,P5)/zb(P3,P5)
     1 *S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*za(P1,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P3,P4)*za(P2,P3)*za(P2,P5)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P3,P4)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P3,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*zb(P4,P5)*za(P2,P5)*za(P2,P5)
     1 *za(P3,P4)/zb(P1,P3)/zb(P3,P5)/za(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*zb(P4,P5)*za(P2,P5)*za(P3,P5)
     1 /zb(P2,P4)/zb(P3,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*za(P2,P3)*za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*za(P2,P3)*za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P2,P3)*za(P2,P3)/zb(P1,P3)/zb(P2,P5)
     1 /zb(P3,P5)/za(P1,P3)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P3,P4)*zb(P3,P4)*za(P2,P3)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P3,P4)*zb(P4,P5)*za(P2,P5)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P3,P5)/za(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P3,P4)*zb(P4,P5)*za(P3,P5)*za(P4,P5)
     1 /zb(P2,P4)/zb(P3,P5)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P1,P5)*za(P2,P3)*za(P3,P4)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P2,P3)*za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P2,P3)*za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P2,P3)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P2,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P2,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= S15 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P2,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P3,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P4)*zb(P3,P4)*za(P3,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P3,P5)/zb(P2,P4)/zb(P3,P5)
     1 *S124**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P4)*zb(P4,P5)*za(P1,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P4,P5)*za(P1,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P4,P5)*za(P1,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P3,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P4,P5)*za(P2,P5)*za(P4,P5)/zb(P2,P4)
     1 /zb(P3,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*za(P1,P2)*za(P3,P5)/zb(P2,P4)/zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*za(P1,P2)*za(P3,P5)/zb(P2,P4)/zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= S25 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*za(P1,P5)*za(P2,P3)/zb(P2,P4)/zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*za(P1,P5)*za(P2,P3)/zb(P2,P4)/zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*za(P2,P5)/zb(P2,P4)/zb(P3,P5)/za(P2,P4)
     1 *S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*za(P2,P5)/zb(P3,P5)*S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*za(P3,P5)/zb(P1,P3)/zb(P2,P5)/za(P1,P3)
     1 *S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*zb(P3,P4)*za(P2,P3)*za(P2,P5)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*zb(P3,P4)*za(P2,P3)*za(P2,P5)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*zb(P3,P4)*za(P2,P3)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P3,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*zb(P4,P5)*za(P2,P5)*za(P2,P5)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P3,P4)*zb(P3,P4)*za(P2,P3)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P3,P4)*zb(P3,P4)*za(P2,P3)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P3,P4)*zb(P4,P5)*za(P2,P5)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P3,P4)*za(P2,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P3,P4)*za(P2,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P3,P4)*za(P2,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P5)*zb(P3,P4)*za(P3,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P4,P5)*za(P1,P5)*za(P2,P5)*za(P3,P5)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P3,P4)*za(P2,P3)*za(P3,P5)/zb(P2,P4)/zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P3,P4)*za(P2,P3)*za(P3,P5)/zb(P2,P4)/zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S25 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P4,P5)*za(P2,P5)*za(P3,P5)/zb(P2,P4)/zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P4,P5)*za(P2,P5)*za(P3,P5)/zb(P2,P4)/zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P4,P5)*za(P2,P5)*za(P3,P5)/zb(P2,P4)/zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
        a32(1,2,1)=XTOTAL
C Punched 146 terms out of 146.
C Expression mpp32
        XTOTAL=0.
      HDPART=-zb(P1,P2)*zb(P1,P3)*za(P1,P4)*za(P1,P5)*za(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P1,P4)
     1 *za(P2,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P1,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P3,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*zb(P3,P4)*za(P1,P2)*za(P1,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P1,P4)
     1 *za(P2,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P1,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P3,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P1,P2)*za(P1,P4)*za(P1,P5)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P1,P4)*za(P1,P5)*za(P2,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P1,P4)*za(P1,P5)*za(P2,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*za(P1,P4)*za(P1,P5)*za(P2,P5)
     1 /zb(P2,P4)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P1,P4)*za(P1,P5)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*zb(P3,P4)*za(P1,P4)*za(P1,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P1,P4)*za(P1,P4)*za(P1,P5)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P1,P4)*za(P1,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P1,P4)*za(P1,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P3,P4)*za(P1,P4)*za(P1,P5)*za(P4,P5)
     1 /zb(P2,P4)/za(P2,P4)*S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P3)*zb(P2,P5)*za(P1,P4)*za(P1,P5)*za(P1,P5)
     1 /zb(P3,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P3)*za(P1,P4)*za(P1,P5)*za(P4,P5)/za(P2,P4)
     1 *S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P2,P5)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P3,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P3)*zb(P3,P4)*za(P1,P2)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P2,P5)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P3,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P1,P5)*za(P3,P4)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= S15 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P3)*zb(P2,P3)*za(P1,P4)*za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P3)*za(P1,P4)*za(P2,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*zb(P2,P3)*za(P1,P4)/zb(P1,P3)/zb(P2,P5)
     1 /zb(P3,P5)/za(P1,P3)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P3)*zb(P2,P3)*za(P1,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P3)*za(P1,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P3)*za(P1,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P5)*zb(P3,P4)*za(P1,P4)*za(P2,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P3,P5)/za(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P1,P5)*za(P4,P5)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P1,P5)*za(P4,P5)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P5)*za(P1,P4)*za(P1,P5)*za(P2,P5)
     1 /zb(P3,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P5)*za(P1,P5)*za(P2,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*zb(P3,P4)*za(P1,P4)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P3,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P1,P5)*za(P3,P4)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P4,P5)/zb(P1,P3)
     1 /zb(P3,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P1,P5)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P1,P5)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*za(P1,P4)*za(P1,P5)/zb(P2,P4)/zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*za(P1,P4)*za(P1,P5)/zb(P2,P4)/zb(P3,P5)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*za(P1,P4)*za(P1,P5)/zb(P3,P5)/za(P1,P3)
     1 *S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*za(P1,P4)*za(P2,P5)*za(P4,P5)/za(P2,P4)
     1 *S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*za(P1,P5)*za(P4,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*za(P1,P5)*za(P4,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*za(P1,P5)*za(P4,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P5)*zb(P3,P4)*zb(P3,P4)*za(P1,P4)*za(P4,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P3,P5)/za(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P5)*zb(P3,P4)*za(P1,P4)*za(P1,P5)*za(P4,P5)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P5)*zb(P3,P4)*za(P1,P4)*za(P1,P5)*za(P4,P5)
     1 /zb(P2,P4)/zb(P3,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=2*zb(P2,P5)*zb(P3,P4)*za(P1,P4)*za(P1,P5)*za(P4,P5)
     1 /zb(P3,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P5)*zb(P3,P4)*za(P1,P5)*za(P4,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P5)*za(P1,P5)*za(P1,P5)*za(P4,P5)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P3,P4)*za(P1,P4)*za(P4,P5)*za(P4,P5)/za(P2,P4)
     1 *S124**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
       a32(2,1,1)=XTOTAL
C Punched 90 terms out of 90.

*end
C Expression ppp42
        XTOTAL=0.
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P1,P5)*za(P1,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P2,P4)*za(P1,P4)*za(P2,P3)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P5)*za(P2,P3)*za(P4,P5)
     1 /zb(P1,P3)/za(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P1,P5)*za(P2,P3)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P4,P5)/zb(P1,P3)/zb(P2,P5)
     1 *S123**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P2,P4)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P5)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P3)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P2,P3)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P2,P4)*za(P1,P4)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P2,P4)*za(P4,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*za(P1,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*zb(P2,P4)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/za(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P4)*za(P2,P3)*za(P4,P5)*za(P4,P5)
     1 /zb(P1,P3)/za(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P4)*za(P3,P4)/zb(P1,P3)/zb(P2,P5)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*za(P1,P5)*za(P3,P4)*za(P3,P5)/za(P1,P3)
     1 *S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*za(P3,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)
     1 *S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P2,P3)*zb(P2,P4)*za(P3,P4)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P4)*za(P3,P4)*za(P4,P5)*za(P4,P5)
     1 /zb(P1,P3)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*zb(P2,P4)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*za(P3,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*za(P3,P4)*za(P3,P5)*za(P3,P5)/za(P1,P3)
     1 *S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P4)*za(P3,P4)*za(P3,P5)*za(P4,P5)/za(P1,P3)
     1 *S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      a42(1,1,1)=XTOTAL
C Punched 40 terms out of 40.
C Expression ppm42
        XTOTAL=0.
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P3)*za(P2,P4)
     1 /zb(P1,P3)/za(P1,P3)/za(P2,P5)/za(P4,P5)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P2,P4)*za(P3,P4)
     1 /zb(P1,P3)/za(P1,P3)/za(P2,P5)/za(P4,P5)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P1,P5)*za(P1,P2)*za(P4,P5)
     1 /zb(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P2,P3)
     1 *za(P4,P5)/zb(P1,P3)/za(P1,P3)/za(P2,P5)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P3,P5)*za(P2,P3)*za(P4,P5)
     1 /zb(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P4,P5)*za(P2,P4)*za(P4,P5)
     1 /zb(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P4)*za(P2,P3)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P4)*za(P2,P4)*za(P3,P5)
     1 /zb(P1,P3)/za(P1,P3)/za(P2,P5)/za(P4,P5)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P2,P4)/zb(P1,P3)/za(P2,P5)
     1 *S123**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P2,P4)*za(P3,P4)
     1 /zb(P1,P3)/za(P1,P3)/za(P2,P5)/za(P4,P5)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*zb(P3,P5)*za(P2,P3)*za(P2,P3)
     1 *za(P4,P5)/zb(P1,P3)/za(P1,P3)/za(P2,P5)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*zb(P4,P5)*za(P2,P3)*za(P2,P4)
     1 *za(P4,P5)/zb(P1,P3)/za(P1,P3)/za(P2,P5)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P2,P3)*za(P2,P4)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P3,P5)*za(P2,P3)*za(P3,P4)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*za(P2,P4)*za(P3,P4)/zb(P1,P3)/za(P1,P3)
     1 /za(P2,P5)/za(P4,P5)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/za(P1,P3)/za(P2,P5)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*za(P1,P4)*za(P3,P4)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*za(P2,P4)*za(P3,P4)*za(P3,P4)
     1 /zb(P1,P3)/za(P1,P3)/za(P2,P5)/za(P4,P5)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P5)*zb(P3,P5)*za(P2,P3)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/za(P1,P3)/za(P2,P5)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P5)*zb(P4,P5)*za(P2,P4)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/za(P1,P3)/za(P2,P5)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*za(P2,P4)*za(P3,P4)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P4)*zb(P3,P5)*za(P3,P4)*za(P3,P4)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P1,P5)*za(P1,P4)*za(P3,P5)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*za(P2,P4)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P3)/za(P1,P3)/za(P2,P5)/za(P4,P5)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P3,P4)*za(P3,P5)
     1 /za(P1,P3)/za(P2,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P3,P5)*za(P3,P4)*za(P3,P5)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*za(P3,P4)/zb(P1,P3)/za(P1,P3)/za(P2,P5)
     1 *S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P5)*zb(P3,P5)*za(P2,P3)*za(P3,P4)*za(P3,P5)
     1 /za(P1,P3)/za(P2,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P5)*zb(P4,P5)*za(P2,P4)*za(P3,P4)*za(P3,P5)
     1 /za(P1,P3)/za(P2,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*za(P2,P3)*za(P3,P4)/za(P1,P3)/za(P2,P5)
     1 *S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      a42(1,1,2)=XTOTAL
C Punched 52 terms out of 52.
C Expression ppp41
        XTOTAL=0.
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P1,P5)
     1 *za(P2,P3)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P3)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P2,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P1,P5)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P2,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P1,P5)*za(P1,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P2,P3)*za(P1,P4)*za(P2,P3)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P2,P3)*za(P1,P5)*za(P2,P3)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P2,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P2,P4)*za(P1,P4)*za(P2,P3)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P3)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P2,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P2,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S25 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P5)*za(P2,P3)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P5)*za(P2,P3)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P1,P5)*za(P2,P3)*za(P4,P5)
     1 /zb(P1,P3)/za(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P1,P5)*za(P2,P3)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P1,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P2,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P1,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P2,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*za(P2,P3)*za(P4,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P3,P4)/zb(P1,P5)/zb(P2,P4)
     1 /zb(P2,P5)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P4,P5)/zb(P1,P3)/zb(P2,P5)
     1 *S123**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P2,P3)*za(P1,P4)*za(P3,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P2,P3)*za(P1,P5)*za(P3,P4)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P2,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P2,P4)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P5)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P5)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P3,P4)*za(P4,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P2,P3)*za(P1,P5)*za(P3,P4)
     1 *za(P3,P5)/zb(P2,P4)/zb(P2,P5)/za(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P2,P3)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P2,P4)*za(P1,P4)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P2,P4)*za(P4,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*za(P1,P5)*za(P3,P4)*za(P4,P5)
     1 /zb(P2,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P5)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*zb(P2,P3)*za(P2,P3)*za(P3,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*zb(P2,P4)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/za(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P3,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P2,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P3,P4)*za(P3,P5)/zb(P2,P4)
     1 /zb(P2,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P2,P4)*za(P2,P3)*za(P4,P5)*za(P4,P5)
     1 /zb(P1,P3)/za(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P4)*za(P3,P4)/zb(P1,P3)/zb(P2,P5)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*za(P3,P4)*za(P4,P5)/zb(P1,P3)/zb(P2,P5)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*za(P3,P4)*za(P4,P5)/zb(P1,P3)/zb(P2,P5)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*za(P3,P4)*za(P4,P5)/zb(P2,P5)/za(P2,P4)
     1 *S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*za(P3,P5)*za(P4,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*za(P3,P5)*za(P4,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*za(P3,P5)*za(P4,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*za(P3,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)
     1 *S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P2,P3)*zb(P2,P3)*za(P3,P4)*za(P3,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P2,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P2,P4)*za(P3,P4)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P4)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*za(P3,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*zb(P2,P3)*za(P3,P4)*za(P3,P5)
     1 *za(P3,P5)/zb(P2,P4)/zb(P2,P5)/za(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*zb(P2,P4)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*za(P3,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*za(P3,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=2*zb(P1,P5)*zb(P2,P3)*za(P3,P4)*za(P3,P5)*za(P4,P5)
     1 /zb(P2,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*za(P3,P5)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*za(P3,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P4)*za(P3,P4)*za(P4,P5)*za(P4,P5)
     1 /zb(P2,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*za(P3,P5)*za(P4,P5)*za(P4,P5)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      a41(1,1,1)=XTOTAL
C Punched 118 terms out of 118.
C Expression ppm41
        XTOTAL=0.
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P1,P5)*za(P1,P2)*za(P1,P4)
     1 *za(P2,P3)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P2)*zb(P3,P5)*za(P1,P4)*za(P2,P3)
     1 *za(P2,P3)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*zb(P4,P5)*za(P1,P4)*za(P2,P3)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P1,P2)*za(P1,P4)*za(P2,P3)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P2)*za(P1,P4)*za(P2,P3)*za(P2,P4)
     1 /zb(P1,P3)/za(P1,P3)/za(P2,P5)/za(P4,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P1,P5)*za(P1,P2)*za(P1,P4)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P3,P5)*za(P1,P4)*za(P2,P3)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P4,P5)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P2)*za(P1,P4)*za(P3,P4)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P4)*za(P2,P4)*za(P3,P4)
     1 /zb(P1,P3)/za(P1,P3)/za(P2,P5)/za(P4,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P1,P5)*za(P1,P2)*za(P4,P5)
     1 /zb(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P2,P3)*za(P1,P2)*za(P2,P3)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P2,P3)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P2,P3)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P2,P3)
     1 *za(P4,P5)/zb(P1,P3)/za(P1,P3)/za(P2,P5)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P3,P5)*za(P2,P3)*za(P4,P5)
     1 /zb(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P4,P5)*za(P1,P4)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*zb(P4,P5)*za(P1,P4)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*zb(P4,P5)*za(P2,P4)*za(P4,P5)
     1 /zb(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*za(P1,P2)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P2)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*za(P1,P4)*za(P2,P3)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P1,P4)*za(P2,P3)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)*S135**(-1)
      TLPART= -S25 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P1,P5)*za(P1,P4)*za(P2,P4)*za(P3,P5)
     1 /zb(P1,P3)/za(P1,P3)/za(P2,P5)/za(P4,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P5)*za(P2,P4)/zb(P1,P3)/za(P2,P5)
     1 *S123**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P2,P3)*zb(P3,P5)*za(P2,P3)*za(P2,P3)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*zb(P4,P5)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*za(P1,P2)*za(P2,P3)*za(P3,P4)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*za(P2,P3)*za(P2,P4)*za(P3,P4)
     1 /zb(P1,P3)/za(P1,P3)/za(P2,P5)/za(P4,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*zb(P3,P5)*za(P2,P3)*za(P2,P3)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*zb(P3,P5)*za(P2,P3)*za(P2,P3)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*zb(P3,P5)*za(P2,P3)*za(P2,P3)
     1 *za(P4,P5)/zb(P1,P3)/za(P1,P3)/za(P2,P5)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*zb(P4,P5)*za(P2,P3)*za(P2,P4)
     1 *za(P4,P5)/zb(P1,P3)/za(P1,P3)/za(P2,P5)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P1,P2)*za(P2,P3)*za(P4,P5)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P1,P4)*za(P2,P3)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P2,P3)*za(P2,P4)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P5)*S123**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P2,P5)*za(P2,P3)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P5)*za(P2,P3)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
      TLPART= -S15 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P3,P5)*za(P2,P3)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P3,P5)*za(P2,P3)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)*S135**(-1)
      TLPART= S12 +S14 +S23 +S25 +S34 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P2)*zb(P4,P5)*za(P1,P4)*za(P3,P4)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P4,P5)*za(P3,P4)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P4,P5)*za(P3,P4)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*za(P1,P2)*za(P3,P4)/zb(P2,P4)/za(P1,P5)
     1 /za(P2,P4)/za(P2,P5)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*za(P2,P4)*za(P3,P4)/zb(P1,P3)/za(P1,P3)
     1 /za(P2,P5)/za(P4,P5)*S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P1,P5)*zb(P2,P3)*za(P1,P2)*za(P3,P4)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P3,P5)*za(P2,P3)*za(P3,P4)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*zb(P4,P5)*za(P3,P4)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*za(P1,P2)*za(P3,P4)*za(P3,P4)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*za(P2,P4)*za(P3,P4)*za(P3,P4)
     1 /zb(P1,P3)/za(P1,P3)/za(P2,P5)/za(P4,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*zb(P3,P5)*za(P2,P3)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*zb(P3,P5)*za(P2,P3)*za(P3,P4)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P5)*za(P1,P2)*za(P3,P4)*za(P4,P5)
     1 /zb(P2,P4)/za(P1,P5)/za(P2,P4)/za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*za(P3,P4)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*za(P3,P4)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
      TLPART= S15 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P1,P5)*zb(P2,P3)*za(P1,P2)*za(P3,P4)
     1 *za(P3,P5)/zb(P2,P4)/za(P2,P4)/za(P2,P5)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P1,P5)*zb(P2,P5)*za(P1,P2)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P1,P5)*za(P1,P2)*za(P3,P4)*za(P4,P5)
     1 /za(P2,P4)/za(P2,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*zb(P3,P5)*za(P2,P3)*za(P3,P4)
     1 *za(P3,P5)/zb(P2,P4)/za(P2,P4)/za(P2,P5)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*zb(P4,P5)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*zb(P4,P5)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P3)*zb(P4,P5)*za(P3,P4)*za(P3,P5)
     1 /zb(P2,P4)/za(P2,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*za(P2,P3)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P3)*za(P2,P3)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)*S135**(-1)
      TLPART= S25 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P2,P3)*za(P2,P3)*za(P3,P4)/zb(P2,P4)
     1 /za(P2,P4)/za(P2,P5)*S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P2,P3)*za(P2,P4)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P3)/za(P1,P3)/za(P2,P5)/za(P4,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P5)*zb(P3,P5)*za(P2,P3)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)
     1 /za(P2,P5)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P2,P5)*zb(P4,P5)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P5)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P5)*za(P2,P3)*za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P5)*za(P2,P3)*za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P2,P5)*za(P2,P3)*za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)/za(P1,P3)/za(P2,P4)/za(P2,P5)*S135**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P5)*zb(P2,P5)*za(P3,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P5)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P2,P5)*za(P3,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P5)*S245**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P2,P5)*za(P3,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P5)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P5)*zb(P3,P5)*za(P2,P3)*za(P3,P4)*za(P4,P5)
     1 /za(P2,P4)/za(P2,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P4,P5)*za(P3,P4)*za(P4,P5)/za(P2,P5)
     1 *S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*za(P3,P4)/za(P2,P5)*S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P3)*zb(P2,P5)*za(P2,P3)*za(P3,P4)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P4,P5)*za(P3,P4)*za(P3,P4)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*zb(P2,P5)*za(P2,P3)*za(P4,P5)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P5)*zb(P4,P5)*za(P3,P4)*za(P4,P5)/zb(P2,P4)
     1 /za(P1,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*za(P3,P4)/zb(P2,P4)/za(P1,P5)/za(P2,P4)
     1 *S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      a41(1,1,2)=XTOTAL
C Punched 146 terms out of 146.

*end
C Expression pmp42
        XTOTAL=0.
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P2,P4)*za(P1,P2)*za(P2,P3)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P5)*za(P2,P3)*za(P2,P5)
     1 /zb(P1,P3)/za(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P5)*za(P2,P3)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P4)*zb(P3,P4)*za(P2,P3)*za(P2,P3)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P4)*za(P2,P3)*za(P2,P5)*za(P2,P5)
     1 /zb(P1,P3)/za(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P2,P3)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P3)/za(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P2,P3)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P4)*zb(P1,P5)*za(P1,P5)*za(P2,P5)
     1 /zb(P1,P3)/zb(P4,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P4)*zb(P2,P4)*za(P1,P2)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P4)*za(P1,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P3)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P4)*za(P1,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P4)*za(P2,P5)/zb(P1,P3)/zb(P4,P5)
     1 *S134**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*zb(P2,P4)*za(P1,P2)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P5)*zb(P2,P4)*za(P2,P5)*za(P2,P5)
     1 /zb(P1,P3)/zb(P4,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P5)*zb(P3,P4)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P3)/zb(P4,P5)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P5)*za(P1,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P4)*zb(P3,P4)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P4)*za(P2,P3)/zb(P1,P3)/zb(P2,P5)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P4)*zb(P2,P4)*za(P2,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P3)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P3,P4)*za(P2,P5)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P3)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P3,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*za(P1,P5)*za(P2,P3)*za(P3,P5)/za(P1,P3)
     1 *S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*za(P3,P5)/zb(P1,P3)/zb(P2,P5)/za(P1,P3)
     1 *S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P5)*zb(P2,P4)*zb(P3,P4)*za(P2,P3)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P3,P4)*za(P3,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P4)*za(P2,P3)*za(P2,P5)*za(P3,P5)/za(P1,P3)
     1 *S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P3,P4)*za(P2,P3)*za(P3,P5)*za(P3,P5)/za(P1,P3)
     1 *S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      a42(1,2,1)=XTOTAL
C Punched 40 terms out of 40.
C Expression mpp42
        XTOTAL=0.
      HDPART=zb(P1,P2)*zb(P2,P3)*zb(P2,P4)*za(P1,P2)*za(P1,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P1,P2)*za(P1,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P4)*zb(P3,P4)*za(P1,P4)*za(P1,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P4)*zb(P3,P5)*za(P1,P4)*za(P1,P5)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P1,P4)*za(P1,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P5)*za(P1,P5)*za(P1,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P1,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P3,P5)*za(P1,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P4,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P5)*zb(P3,P4)*za(P1,P4)*za(P1,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*za(P1,P4)*za(P1,P5)*za(P1,P5)
     1 /zb(P4,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P3)*zb(P2,P4)*za(P1,P2)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P4)*zb(P2,P5)*za(P1,P2)*za(P2,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P4)*zb(P3,P4)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P4)*zb(P3,P5)*za(P1,P5)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P5)/zb(P4,P5)/za(P1,P3)*S245**(-1)
     1 
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P4)*zb(P3,P5)*za(P2,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P4,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P4)*za(P1,P2)*za(P4,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S123**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*zb(P2,P4)*za(P1,P4)/zb(P1,P3)/zb(P2,P5)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P3)*zb(P2,P5)*zb(P3,P4)*za(P1,P2)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P3,P4)*zb(P3,P5)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P4,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P4,P5)/zb(P1,P3)/zb(P4,P5)
     1 *S123**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*zb(P3,P5)*za(P1,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*za(P1,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)
     1 *S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P4)*zb(P2,P5)*zb(P3,P4)*za(P1,P4)*za(P2,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P4)*zb(P2,P5)*za(P1,P4)*za(P1,P5)*za(P2,P5)
     1 /zb(P4,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P4)*zb(P3,P4)*za(P1,P4)*za(P4,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P4)*za(P1,P4)*za(P1,P5)/zb(P4,P5)/za(P1,P3)
     1 *S134**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P5)*zb(P3,P4)*zb(P3,P4)*za(P1,P4)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)*S134**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*zb(P3,P4)*za(P1,P4)*za(P1,P5)*za(P3,P5)
     1 /zb(P4,P5)/za(P1,P3)*S134**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      a42(2,1,1)=XTOTAL
C Punched 52 terms out of 52.
C Expression pmp41
        XTOTAL=0.
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P1,P4)*za(P1,P2)*za(P1,P5)
     1 *za(P2,P3)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P1,P4)*za(P1,P2)*za(P2,P3)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P3,P4)*za(P1,P2)*za(P2,P3)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P3,P4)*za(P1,P5)*za(P2,P3)
     1 *za(P2,P3)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P3,P4)*za(P2,P3)*za(P2,P3)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P2)*za(P2,P3)*za(P2,P5)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P5)*za(P2,P3)*za(P2,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P1,P5)*za(P2,P3)*za(P2,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*za(P1,P5)*za(P2,P3)*za(P2,P5)
     1 /zb(P1,P3)/za(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*za(P2,P3)*za(P2,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P4)*za(P2,P3)*za(P2,P5)*za(P2,P5)
     1 /zb(P1,P3)/za(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P3,P4)*zb(P3,P4)*za(P2,P3)*za(P2,P3)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P2,P3)*za(P2,P3)*za(P2,P5)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P2,P3)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P2,P3)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P3,P4)*za(P2,P3)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P3)/za(P1,P3)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P4)*zb(P1,P4)*za(P1,P2)*za(P1,P5)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P4)*zb(P1,P4)*za(P1,P2)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P4)*zb(P3,P4)*za(P1,P2)*za(P3,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P4)*zb(P3,P4)*za(P1,P5)*za(P2,P3)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P4)*zb(P3,P4)*za(P2,P3)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P4)*za(P1,P2)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P4)*za(P1,P2)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P4)*za(P1,P2)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S25 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P4)*zb(P1,P4)*za(P1,P5)*za(P2,P3)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P4)*za(P1,P5)*za(P2,P3)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P1,P4)*za(P1,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P4)*za(P1,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P4)*za(P2,P3)/zb(P1,P5)/zb(P2,P4)
     1 /zb(P4,P5)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P1,P4)*za(P2,P5)*za(P3,P4)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*zb(P3,P4)*za(P1,P5)*za(P2,P3)
     1 *za(P3,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*za(P1,P2)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*za(P1,P2)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P1,P5)*za(P1,P5)*za(P2,P3)*za(P2,P5)
     1 /zb(P4,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P1,P5)*za(P1,P5)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P3,P4)*zb(P3,P4)*za(P2,P3)*za(P3,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P2,P3)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P2,P3)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P2,P3)*za(P3,P5)/zb(P2,P4)
     1 /zb(P4,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P2,P5)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P3,P4)*za(P2,P5)*za(P3,P4)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*za(P1,P5)*za(P2,P3)*za(P3,P5)/za(P1,P3)
     1 *S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*za(P2,P3)*za(P2,P5)/zb(P1,P3)/zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*za(P2,P3)*za(P2,P5)/zb(P1,P3)/zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*za(P2,P3)*za(P2,P5)/zb(P4,P5)/za(P2,P4)
     1 *S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*za(P2,P5)*za(P3,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*za(P2,P5)*za(P3,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*za(P2,P5)*za(P3,P5)/zb(P1,P3)/zb(P2,P4)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P1,P5)*zb(P2,P4)*za(P2,P3)*za(P2,P5)*za(P2,P5)
     1 /zb(P4,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*zb(P3,P4)*zb(P3,P4)*za(P2,P3)*za(P3,P5)
     1 *za(P3,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P3,P4)*za(P2,P3)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P3,P4)*za(P2,P3)*za(P2,P5)*za(P3,P5)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-2*zb(P1,P5)*zb(P3,P4)*za(P2,P3)*za(P2,P5)*za(P3,P5)
     1 /zb(P4,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P5)*zb(P3,P4)*za(P2,P5)*za(P3,P5)*za(P3,P5)
     1 /zb(P1,P3)/zb(P2,P4)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P5)*za(P2,P5)*za(P2,P5)*za(P3,P5)/zb(P1,P3)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P4)*za(P2,P3)*za(P2,P5)*za(P3,P5)/za(P1,P3)
     1 *S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P3,P4)*za(P2,P3)*za(P3,P5)*za(P3,P5)/za(P1,P3)
     1 *S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      a41(1,2,1)=XTOTAL
C Punched 90 terms out of 90.
C Expression mpp41
        XTOTAL=0.
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P2,P3)*za(P1,P2)*za(P1,P4)
     1 *za(P1,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P1,P4)*zb(P2,P3)*za(P1,P2)*za(P1,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P3,P4)*za(P1,P4)*za(P1,P4)
     1 *za(P1,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P1,P4)*zb(P3,P4)*za(P1,P4)*za(P1,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*zb(P3,P4)*za(P1,P2)*za(P1,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P2,P3)*za(P1,P2)*za(P1,P4)*za(P2,P5)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*za(P1,P2)*za(P1,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P2,P3)*za(P1,P4)*za(P2,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P3,P4)*zb(P3,P4)*za(P1,P4)*za(P1,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P1,P4)*za(P1,P4)*za(P2,P5)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P4)*za(P1,P4)*za(P1,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P3,P4)*za(P1,P4)*za(P1,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S25 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P2)*zb(P3,P4)*za(P1,P4)*za(P1,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P3,P4)*za(P1,P4)*za(P4,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P5)*za(P1,P4)*za(P1,P5)*za(P2,P5)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P2)*zb(P3,P5)*za(P1,P4)*za(P1,P5)*za(P2,P5)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P2)*zb(P3,P5)*za(P1,P5)*za(P1,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P1,P5)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P1,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P1,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P1,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P2,P5)*za(P1,P2)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P2,P5)*za(P1,P5)*za(P1,P5)
     1 /zb(P2,P4)/zb(P4,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P1,P5)
     1 *za(P3,P4)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P3,P4)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*zb(P3,P5)*za(P1,P5)*za(P1,P5)
     1 *za(P3,P4)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*zb(P3,P5)*za(P1,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P4,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P1,P4)*zb(P2,P3)*za(P1,P4)*za(P1,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P3)*za(P1,P4)*za(P1,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P2,P3)*za(P1,P4)/zb(P1,P5)/zb(P2,P4)
     1 /zb(P4,P5)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P1,P4)*zb(P2,P5)*zb(P3,P4)*za(P1,P4)*za(P1,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*zb(P3,P4)*za(P1,P4)*za(P1,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*zb(P3,P4)*za(P1,P4)*za(P4,P5)
     1 /zb(P1,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P1,P4)*zb(P2,P5)*zb(P3,P5)*za(P1,P5)*za(P1,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P3)*zb(P3,P4)*za(P1,P2)*za(P3,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P3)*za(P1,P2)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P3)*za(P2,P5)*za(P3,P4)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P4)*zb(P2,P5)*za(P1,P2)*za(P2,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P4)*zb(P3,P5)*za(P2,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P4,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P4)*za(P1,P2)*za(P4,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S123**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P3)*zb(P2,P5)*zb(P3,P4)*za(P1,P2)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P5)*zb(P3,P4)*za(P1,P2)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P2,P5)*zb(P3,P4)*za(P1,P2)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)*S123**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P5)*zb(P3,P4)*za(P1,P5)*za(P3,P5)
     1 /zb(P2,P4)/zb(P4,P5)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P5)*za(P1,P5)*za(P2,P5)/zb(P4,P5)
     1 *S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P2,P5)*za(P2,P5)*za(P4,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P3,P4)*zb(P3,P4)*za(P1,P4)*za(P3,P4)
     1 *za(P3,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P3,P4)*zb(P3,P5)*za(P1,P5)*za(P3,P4)
     1 *za(P3,P5)/zb(P2,P4)/zb(P4,P5)/za(P2,P4)*S234**(-1)
     1 *S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*zb(P3,P5)*za(P3,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P4,P5)*S123**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P3)*zb(P3,P4)*za(P1,P4)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P1,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P1,P5)*za(P3,P4)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= S25 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P1,P5)*za(P3,P4)/zb(P2,P4)
     1 /zb(P4,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P2,P3)*zb(P3,P4)*za(P3,P4)*za(P4,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P4)*za(P4,P5)/zb(P1,P3)/zb(P4,P5)
     1 *S123**(-1)*S1234**(-1)
      TLPART= S15 +S25 +S35 +S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*zb(P3,P5)*za(P1,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*zb(P3,P5)*za(P1,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P3,P5)*za(P1,P5)*za(P2,P5)*za(P3,P4)
     1 /zb(P4,P5)/za(P2,P4)*S234**(-1)*S1234**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*zb(P3,P5)*za(P1,P5)*za(P3,P5)/zb(P1,P3)
     1 /zb(P4,P5)/za(P1,P3)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P3)*za(P1,P2)*za(P4,P5)/zb(P1,P3)/zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*za(P1,P2)*za(P4,P5)/zb(P1,P3)/zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= -S15 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*za(P1,P4)*za(P2,P5)/zb(P1,P3)/zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P3)*za(P1,P4)*za(P2,P5)/zb(P1,P3)/zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*za(P1,P5)/zb(P1,P3)/zb(P4,P5)/za(P1,P3)
     1 *S245**(-1)
      TLPART= -S12 -S14 -S15 -S23 -S34 -S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*za(P1,P5)/zb(P4,P5)*S234**(-1)*S1234**(-1)
      TLPART= -S15 -S25 -S35 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P3)*za(P4,P5)/zb(P1,P5)/zb(P2,P4)/za(P2,P4)
     1 *S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P5)*zb(P3,P4)*zb(P3,P4)*za(P1,P4)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*zb(P3,P4)*zb(P3,P4)*za(P1,P4)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*zb(P3,P4)*zb(P3,P5)*za(P1,P5)*za(P3,P5)
     1 *za(P4,P5)/zb(P1,P3)/zb(P2,P4)/zb(P4,P5)/za(P1,P3)
     1 /za(P2,P4)*S245**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*zb(P3,P4)*za(P1,P5)*za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P2,P5)*zb(P3,P4)*za(P1,P5)*za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= S12 +S14 +S15 +S23 +S34 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P5)*zb(P3,P4)*za(P1,P5)*za(P4,P5)/zb(P1,P3)
     1 /zb(P2,P4)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= S12 +S14 +S23 +S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P2,P5)*zb(P3,P4)*za(P4,P5)*za(P4,P5)/zb(P1,P5)
     1 /zb(P2,P4)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=-zb(P2,P5)*zb(P3,P5)*za(P1,P5)*za(P2,P5)*za(P4,P5)
     1 /zb(P1,P3)/zb(P4,P5)/za(P1,P3)/za(P2,P4)*S135**(-1)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P3,P4)*za(P1,P4)*za(P4,P5)/zb(P1,P3)/zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P3,P4)*za(P1,P4)*za(P4,P5)/zb(P1,P3)/zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= S15 +S35
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=-zb(P3,P5)*za(P1,P5)*za(P4,P5)/zb(P1,P3)/zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)
        XTOTAL=XTOTAL+HDPART
      HDPART=zb(P3,P5)*za(P1,P5)*za(P4,P5)/zb(P1,P3)/zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)*S245**(-1)
      TLPART= -S12 -S14 -S23 -S34
        XTOTAL=XTOTAL+HDPART*TLPART
      HDPART=zb(P3,P5)*za(P1,P5)*za(P4,P5)/zb(P1,P3)/zb(P4,P5)
     1 /za(P1,P3)/za(P2,P4)*S135**(-1)
      TLPART= -S12 -S14 -S23 -S25 -S34 -S45
        XTOTAL=XTOTAL+HDPART*TLPART
      a41(2,1,1)=XTOTAL
C Punched 146 terms out of 146.
*end

      write(6,*) 
      write(6,*) 'mmm31',a31(1,1,1)
      write(6,*) 'mmp31',a31(1,1,2)
      write(6,*) 'mmm32',a32(1,1,1)
      write(6,*) 'mmp32',a32(1,1,2)

      write(6,*) 
      write(6,*) 'mpm31',a31(1,2,1)
      write(6,*) 'pmm31',a31(2,1,1)
      write(6,*) 'mpm32',a32(1,2,1)
      write(6,*) 'pmm32',a32(2,1,1)

      write(6,*) 
      write(6,*) 'mmm42',a42(1,1,1)
      write(6,*) 'mmp42',a42(1,1,2)
      write(6,*) 'mmm41',a41(1,1,1)
      write(6,*) 'mmp41',a41(1,1,2)

      write(6,*) 
      write(6,*) 'mpm42',a42(1,2,1)
      write(6,*) 'pmm42',a42(2,1,1)
      write(6,*) 'mpm41',a41(1,2,1)
      write(6,*) 'pmm41',a41(2,1,1)
      msq=0d0
      do h1=1,2
      do h2=1,2
      do h5=1,2
      msq=msq+0.5d0*V*(
     . +(abs(a32(h1,h2,h5))**2+abs(a41(h1,h2,h5))**2)*xn
     . +(abs(a31(h1,h2,h5))**2+abs(a42(h1,h2,h5))**2)/xn
     . +two*(
     . +dble(a31(h1,h2,h5)*Dconjg(a32(h1,h2,h5)))
     . +dble(a31(h1,h2,h5)*Dconjg(a41(h1,h2,h5)))
     . +dble(a32(h1,h2,h5)*Dconjg(a42(h1,h2,h5)))
     . +dble(a41(h1,h2,h5)*Dconjg(a42(h1,h2,h5))))/xn)
      enddo
      enddo
      enddo
C full matrix element is given by the above multiplied by
C      g^6*A^2/2d0
      return
      end
