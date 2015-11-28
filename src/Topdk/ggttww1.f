      subroutine ggttww1(s1t,s2t,s12,loab,loba)
c--- Helicity amplitudes for the sub-process gg -> t+tbar,
c--- including decays (see qqb_QQbdk.f for calling details).
      implicit none
      include 'constants.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer q1,q2,q3,q4,q5,q6,q7,q8
      double complex loab(2,2),loba(2,2)
      double complex x1mm,x1pp,x1pm,x1mp
      double complex x2mm,x2pp,x2pm,x2mp
      double complex x3mm,x3pp
      double precision mt2,s1t,s2t,s12
      parameter(q1=1,q2=2,q3=3,q4=4,q5=5,q6=6,q7=7,q8=8)

      mt2=mt**2

      x1pm=-zb(q1,q6)*za(q2,q6)/(s12*s1t)
     &*(za(q7,q5)*za(q2,q3)*zb(q1,q5)*zb(q3,q4)-za(q7,q2)*zb(q1,q4)*mt2)
      x2pm=-za(q2,q8)*zb(q1,q8)/(s12*s2t)
     &*(za(q7,q5)*za(q2,q3)*zb(q1,q5)*zb(q3,q4)-za(q7,q2)*zb(q1,q4)*mt2)

      x1mp=-za(q1,q6)*zb(q2,q6)/(s12*s1t)
     &*(za(q7,q5)*za(q1,q3)*zb(q2,q5)*zb(q3,q4)-za(q7,q1)*zb(q2,q4)*mt2)
      x2mp=-za(q1,q8)*zb(q2,q8)/(s12*s2t)
     &*(za(q7,q5)*za(q1,q3)*zb(q2,q5)*zb(q3,q4)-za(q7,q1)*zb(q2,q4)*mt2)

      x1mm=
     & za(q7,q5)*za(q1,q3)*zb(q3,q4)*zb(q4,q5)/(zb(q1,q4)*zb(q2,q4)*s1t)
     & *(za(q1,q2)*zb(q1,q4)+za(q2,q8)*zb(q4,q8))
      x2mm=
     & za(q7,q5)*za(q2,q3)*zb(q3,q4)*zb(q4,q5)/(zb(q1,q4)*zb(q2,q4)*s2t)
     & *(-za(q1,q2)*zb(q2,q4)+za(q1,q6)*zb(q4,q6))
      x3mm=
     & za(q7,q5)*za(q1,q2)*zb(q3,q4)*zb(q4,q5)/(zb(q1,q4)*zb(q2,q4)*s12)
     & *(za(q1,q3)*zb(q1,q4)+za(q2,q3)*zb(q2,q4))

      x1pp=za(q7,q3)*za(q7,q5)*za(q7,q6)*zb(q1,q6)*zb(q2,q5)*zb(q3,q4)
     & /(za(q7,q1)*za(q7,q2)*s1t)
      x2pp=za(q7,q3)*za(q7,q5)*za(q7,q8)*zb(q1,q5)*zb(q2,q8)*zb(q3,q4)
     & /(za(q7,q1)*za(q7,q2)*s2t)
      x3pp=-za(q7,q3)*za(q7,q5)*zb(q1,q2)*zb(q3,q4)
     &    *(za(q7,q1)*zb(q1,q5)+za(q7,q2)*zb(q2,q5))
     & /(za(q7,q1)*za(q7,q2)*s12)

      loab(1,1)=x1mm+x3mm
      loab(1,2)=x1mp
      loab(2,1)=x1pm
      loab(2,2)=x1pp+x3pp

      loba(1,1)=x2mm-x3mm
      loba(1,2)=x2mp
      loba(2,1)=x2pm
      loba(2,2)=x2pp-x3pp
 
    
      return
      end
