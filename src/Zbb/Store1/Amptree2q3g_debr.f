      double complex Function mPPPPM(p,k1,k2,k3,q,l,l_)
      implicit none
      include 'constants.f'
      include 'debr.f'
      double complex zU, zV
      call spassign(p,k1,k2,k3,q,l,l_)
      zV=HL_PH
      zU=(HPQH*TQLT+HP1H*T1LT+HP2H*T2LT+HP3H*T3LT)/(HQ1H*H12H*H23H*H3PH)
      mPPPPM=-zV*zU
      return
      end

      double complex Function mPMMMM(p,k1,k2,k3,q,l,l_)
      implicit none
      include 'constants.f'
      include 'debr.f'
      double complex zU, zV
      call spassign(p,k1,k2,k3,q,l,l_)
      zU=TQLT
      zV=(HL_PH*TPQT+HL_1H*T1QT+HL_2H*T2QT+HL_3H*T3QT)
     . /(TQ1T*T12T*T23T*T3PT)
      mPMMMM=zU*zV
      return
      end

      double complex Function mPPPMM(p,k1,k2,k3,q,l,l_) 
      implicit none
      include 'constants.f'
      include 'debr.f'
      include 'sprods_com.f'
      double complex m1, m2, m3, m41, m42, m43
      double precision s3,s4
      integer a,b,c,d
      S3(a,b,c)=(s(a,b)+s(a,c)+s(b,c))
      S4(a,b,c,d)=(s(a,b)+s(a,c)+s(a,d)+s(b,c)+s(b,d)+s(c,d))
      call spassign(p,k1,k2,k3,q,l,l_)

      m1=-HL_PH*(H3QH*TQLT+H31H*T1LT+H32H*T2LT)
     . /(H12H*H23H*T23T*S4(q,k1,k2,k3))*
     . (T12T*(TQ1T*H13H+TQ2T*H23H)/S3(k1,k2,k3)
     . +(H3QH*TQ2T+H31H*T12T)/HQ1H)
     
      m2=-(H3QH*TQLT+H31H*T1LT+H32H*T2LT)*(HL_PH*TP2T+HL_3H*T32T)/
     . (HQ1H*H12H*H23H*T23T*T3PT)

      m3=-(H3QH*TQLT+H31H*T1LT)*(HL_PH*TP2T+HL_3H*T32T)*H3PH*TP2T/
     . (HQ1H*H13H*T3PT*S23*S3(p,k2,k3))
  
      m41=TQLT*(HL_PH*TP1T+HL_2H*T21T+HL_3H*T31T)*H3PH*H3PH*TP2T*TP2T/
     . (T3PT*H31H*S23*S3(p,k2,k3)*S4(p,k1,k2,k3))
  
      m42=-TQLT*H3PH*TP2T/(H12H*T23T*T3PT*S4(p,k1,k2,k3))*
     . ((HL_PH*TP1T+HL_2H*T21T+HL_3H*T31T)/H23H+
     . (HL_PH*TP2T+HL_1H*T12T+HL_3H*T32T)/H13H)
  
      m43=TQLT*H3PH*T12T/(H12H*T23T*S3(k1,k2,k3)*S4(p,k1,k2,k3))*
     . ((HL_PH*TP1T+HL_2H*T21T+HL_3H*T31T)*H13H/H23H+
     . (HL_PH*TP2T+HL_1H*T12T+HL_3H*T32T))
  
      mPPPMM=(m1+m2+m3+m41+m42+m43)
      return
      end

      double complex Function mPPMPM(p,k1,k2,k3,q,l,l_) 
      implicit none
      include 'constants.f'
      include 'debr.f'
      include 'sprods_com.f'
      double complex m11, m12, m13, m2, m3, m41, m43
      double precision S3,S4
      integer a,b,c,d
      S3(a,b,c)=(s(a,b)+s(a,c)+s(b,c))
      S4(a,b,c,d)=(s(a,b)+s(a,c)+s(a,d)+s(b,c)+s(b,d)+s(c,d))
      call spassign(p,k1,k2,k3,q,l,l_)

      m11=-T31T*T31T*(TQ1T*H12H+TQ3T*H32H)
     . *(H2QH*TQLT+H21H*T1LT+H23H*T3LT)*HL_PH/
     . (S12*S23*S3(k1,k2,k3)*S4(q,k1,k2,k3))
  
      m12=-T31T*(H2QH*TQ3T+H21H*T13T)
     . *(H2QH*TQLT+H21H*T1LT+H23H*T3LT)*HL_PH/
     . (T12T*HQ1H*H12H*S23*S4(q,k1,k2,k3))
  
      m13=-H2QH*TQ1T*(H2QH*TQ3T+H21H*T13T)
     . *(H2QH*TQLT+H21H*T1LT+H23H*T3LT)*HL_PH/
     . (HQ1H*H23H*S12*S3(q,k1,k2)*S4(q,k1,k2,k3))

      m2=H2QH*TQ1T*(H2QH*TQLT+H21H*T1LT)*HL_PH*H2PH/
     . (HQ1H*H23H*H3PH*S12*S3(q,k1,k2))
  
      m3=(H2QH*TQLT+H21H*T1LT)
     . *(HL_PH*TP1T+HL_2H*T21T+HL_3H*T31T)*H2PH*H2PH/
     . (HQ1H*H12H*H23H*H3PH*T12T*S3(p,k2,k3))+
     . (H2QH*TQLT+H21H*T1LT)*(HL_PH*TP3T+HL_2H*T23T)*H2PH*T31T/
     . (T12T*HQ1H*H12H*S23*S3(p,k3,k2))

      m41=-TQLT*(HL_PH*TP1T+HL_2H*T21T+HL_3H*T31T)*H2PH*H2PH/
     . (S12*S3(p,k2,k3)*S4(p,k1,k2,k3))*
     . ((H2PH*TP1T+H23H*T31T)/(H23H*H3PH)+T31T*TP3T/S23)
  
      m43=T31T*T31T*TQLT/(S12*S23*S3(k1,k2,k3)*S4(p,k1,k2,k3))*
     . ((HL_PH*TP1T+HL_2H*T21T+HL_3H*T31T)*H12H*H2PH+
     . (HL_PH*TP3T+HL_1H*T13T+HL_2H*T23T)*H32H*H2PH)
  
      mPPMPM=(m11+m12+m13+m2+m3+m41+m43)
      return
      end

      double complex Function mPMPPM(p,k1,k2,k3,q,l,l_) 
      implicit none
      include 'constants.f'
      include 'debr.f'
      include 'sprods_com.f'
      double complex m11, m12, m13, m2, m3, m41, m43
      double precision S3,S4
      integer a,b,c,d
      S3(a,b,c)=(s(a,b)+s(a,c)+s(b,c))
      S4(a,b,c,d)=(s(a,b)+s(a,c)+s(a,d)+s(b,c)+s(b,d)+s(c,d))  
      call spassign(p,k1,k2,k3,q,l,l_)
  
      m11=-T32T*T32T*(TQ2T*H21H+TQ3T*H31H)
     . *(H1QH*TQLT+H12H*T2LT+H13H*T3LT)*HL_PH/
     . (S12*S23*S3(k1,k2,k3)*S4(q,k1,k2,k3))

      m12=-TQ2T*(TQ2T*H21H+TQ3T*H31H)
     . *(H1QH*TQLT+H12H*T2LT+H13H*T3LT)*HL_PH/
     . (TQ1T*H23H*H31H*S12*S4(q,k1,k2,k3))

      m13=TQ2T*TQ2T*(H1QH*TQ3T+H12H*T23T)
     . *(H1QH*TQLT+H12H*T2LT+H13H*T3LT)*HL_PH/
     . (TQ1T*H13H*S12*S3(q,k1,k2)*S4(q,k1,k2,k3))
  
      m2=-TQ2T*TQ2T*(H1QH*TQLT+H12H*T2LT)*HL_PH*H1PH/
     . (TQ1T*H13H*H3PH*S12*S3(q,k1,k2))

      m3=TQ2T*TQLT*H1PH/(TQ1T*S12*H13H*S3(p,k2,k3))*
     . ((HL_PH*TP2T+HL_3H*T32T)*H1PH/H3PH-
     . (HL_PH*TP2T+HL_3H*T32T)*H21H/H23H-
     . (HL_PH*TP3T+HL_2H*T23T)*H31H/H23H)
  
      m41=-TQLT*(HL_PH*TP2T+HL_1H*T12T+HL_3H*T32T)*H1PH/
     . (H13H*S12*S3(p,k2,k3)*S4(p,k1,k2,k3))*
     . ((H1PH*TP2T+H13H*T32T)*H1PH/H3PH-
     . (H1PH*TP2T+H13H*T32T)*H21H/H23H-
     . (H1PH*TP3T+H12H*T23T)*H31H/H23H)

      m43=TQLT*H1PH*T32T*T32T/(S12*S23*S3(k1,k2,k3)*S4(p,k1,k2,k3))*
     . ((HL_PH*TP2T+HL_1H*T12T+HL_3H*T32T)*H21H+
     . (HL_PH*TP3T+HL_1H*T13T+HL_2H*T23T)*H31H)
  
      mPMPPM=(m11+m12+m13+m2+m3+m41+m43)
      return
      end

      double complex Function mPPMMM(p,k1,k2,k3,q,l,l_) 
      implicit none
      include 'constants.f'
      include 'debr.f'
      include 'sprods_com.f'
      double complex m11, m13, m2, m3, m41, m43, TEMP
      double precision S3,S4
      integer a,b,c,d
      S3(a,b,c)=(s(a,b)+s(a,c)+s(b,c))
      S4(a,b,c,d)=(s(a,b)+s(a,c)+s(a,d)+s(b,c)+s(b,d)+s(c,d))  
      call spassign(p,k1,k2,k3,q,l,l_)
  
      TEMP=T12T*(H2QH*TQLT+H21H*T1LT+H23H*T3LT)
     . +T13T*(H3QH*TQLT+H31H*T1LT+H32H*T2LT)
    
      m11=-TQ1T*TEMP*HL_PH/(S12*S4(q,k1,k2,k3))*
     . (H32H*H32H/(S23*S3(k1,k2,k3))-H2QH/(T23T*T31T*HQ1H))

  
      m13=H2QH*H2QH*TQ1T*TQ1T*(H3QH*TQLT+H31H*T1LT+H32H*T2LT)*HL_PH/
     . (HQ1H*T13T*S12*S3(q,k1,k2)*S4(q,k1,k2,k3))

      m2=H2QH*TQ1T*(H2QH*TQLT+H21H*T1LT)*(HL_PH*TP1T+HL_3H*T31T)/
     . (HQ1H*T13T*T3PT*S12*S3(q,k1,k2))
  
      m3=-(H2QH*TQLT+H21H*T1LT)*(HL_PH*TP1T+HL_2H*T21T+HL_3H*T31T)/
     . (HQ1H*T23T*T3PT*S12)

      m41=-TQLT*(HL_PH*TP1T+HL_2H*T21T+HL_3H*T31T)
     . *(H2PH*TP1T+H23H*T31T)/
     . (T23T*T3PT*S12*S4(p,k1,k2,k3))

      m43=TQLT*(HL_PH*TP1T+HL_2H*T21T+HL_3H*T31T)
     . *(HP2H*T21T+HP3H*T31T)*H32H*H32H/
     . (S12*S23*S3(k1,k2,k3)*S4(p,k1,k2,k3))
  
       mPPMMM=(m11+m13+m2+m3+m41+m43)
      return
      end

      double complex Function mPMPMM(p,k1,k2,k3,q,l,l_) 
      implicit none
      include 'constants.f'
      include 'debr.f'
      include 'sprods_com.f'
      double complex m11, m12, m13, m2, m3, m41, m43
      double precision S3,S4
      integer a,b,c,d
      S3(a,b,c)=(s(a,b)+s(a,c)+s(b,c))
      S4(a,b,c,d)=(s(a,b)+s(a,c)+s(a,d)+s(b,c)+s(b,d)+s(c,d))
      call spassign(p,k1,k2,k3,q,l,l_)

      m11=-H31H*H31H*TQ2T*(T21T*(H1QH*TQLT+H12H*T2LT+H13H*T3LT)+
     . T23T*(H3QH*TQLT+H32H*T2LT+H31H*T1LT))*HL_PH/
     . (S12*S23*S3(k1,k2,k3)*S4(q,k1,k2,k3))


      m12=H31H*TQ2T*TQ2T*(H3QH*TQLT+H31H*T1LT+H32H*T2LT)*HL_PH/
     . (TQ1T*S12*S23*S4(q,k1,k2,k3))
  
      m13=-TQ2T*TQ2T*TQ2T*H1QH*(H3QH*TQLT+H31H*T1LT+H32H*T2LT)*HL_PH/
     . (TQ1T*T23T*S12*S3(q,k1,k2)*S4(q,k1,k2,k3))

      m2=-TQ2T*TQ2T*(H1QH*TQLT+H12H*T2LT)*(HL_PH*TP2T+HL_3H*T32T)/
     . (TQ1T*T23T*T3PT*S12*S3(k1,k2,q))

      m3=-TQ2T*TQLT*(HL_PH*TP2T+HL_3H*T32T)/(TQ1T*S12*S3(k3,k2,p))*
     . (H3PH*H31H/S23-(H1PH*TP2T+H13H*T32T)/(T23T*T3PT))

      m41=TQLT*(HL_PH*TP2T+HL_1H*T12T+HL_3H*T32T)*(H1PH*TP2T+H13H*T32T)/
     . (S12*S3(k3,k2,p)*S4(k1,k3,k2,p))*
     . (H3PH*H31H/S23-(H1PH*TP2T+H13H*T32T)/(T23T*T3PT))

      m43=TQLT*(HL_PH*TP2T+HL_1H*T12T+HL_3H*T32T)
     . *H31H*H31H*(HP1H*T12T+HP3H*T32T)/
     . (S12*S23*S3(k1,k2,k3)*S4(p,k1,k2,k3))

      mPMPMM=(m11+m12+m13+m2+m3+m41+m43)
      return
      end

      double complex Function mPMMPM(p,k1,k2,k3,q,l,l_) 
      implicit none
      include 'constants.f'
      include 'debr.f'
      include 'sprods_com.f'
      double complex m11, m13, m2, m3, m41, m42, m43
      double precision S3,S4
      integer a,b,c,d
      S3(a,b,c)=(s(a,b)+s(a,c)+s(b,c))
      S4(a,b,c,d)=(s(a,b)+s(a,c)+s(a,d)+s(b,c)+s(b,d)+s(c,d))
      call spassign(p,k1,k2,k3,q,l,l_)

      m11=-H12H*H12H*TQ3T*(T31T*(H1QH*TQLT+H12H*T2LT+H13H*T3LT)+
     . T32T*(H2QH*TQLT+H21H*T1LT+H23H*T3LT))*HL_PH/
     . (S12*S23*S3(k1,k2,k3)*S4(q,k1,k2,k3))

      m13=TQ3T*(H2QH*TQLT+H21H*T1LT+H23H*T3LT)*HL_PH/
     . (T31T*S23*S3(q,k1,k2)*S4(q,k1,k2,k3))*
     . ((T31T*(H1QH*TQ3T+H12H*T23T)+T32T*(H2QH*TQ3T+H21H*T13T))/T12T-
     . TQ3T*(H2QH*TQ3T+H21H*T13T)/TQ1T)

      m2=-TQ3T*HL_PH*H2PH/(T31T*H3PH*S23*S3(q,k1,k2))*
     . ((T31T*(H1QH*TQLT+H12H*T2LT)+T32T*(H2QH*TQLT+H21H*T1LT))/T12T-
     . TQ3T*(H2QH*TQLT+H21H*T1LT)/TQ1T)

      m3=TQ3T*TQLT*(HL_PH*TP3T+HL_2H*T23T)*H2PH*H2PH/
     . (TQ1T*T13T*H3PH*S23*S3(p,k2,k3))

      m41=TQLT*(HL_PH*TP3T+HL_1H*T13T+HL_2H*T23T)
     . *(H1PH*TP3T+H12H*T23T)*H2PH*H2PH/
     . (T31T*H3PH*S23*S3(p,k2,k3)*S4(p,k1,k2,k3))

      m42=TQLT*(HL_PH*TP3T+HL_1H*T13T+HL_2H*T23T)
     . *(T31T*H1PH+T32T*H2PH)*H2PH/
     . (H3PH*T31T*T12T*S23*S4(p,k1,k2,k3))
  
      m43=TQLT*(HL_PH*TP3T+HL_1H*T13T+HL_2H*T23T)
     . *(T31T*H1PH+T32T*H2PH)*H12H*H12H/
     . (S12*S23*S3(k1,k2,k3)*S4(p,k1,k2,k3))
  
      mPMMPM=(m11+m13+m2+m3+m41+m42+m43)
      return
      end
