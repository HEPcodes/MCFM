      double complex function Fsc4(j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      integer j1,j2,j3,j4,j5,j6
      double complex L0,L1,Lsm1_2mh,I3m,Lnrat
      double precision t  
      Fsc4= 
     .-(Lnrat(-s(j2,j3),-s(j5,j6))*zb(j1,j2)*
     .((za(j1,j5)*za(j3,j5)*zb(j1,j3)+za(j1,j5)*za(j4,j5)*zb(j1,j4))/
     .za(j5,j6)-(-(za(j1,j4)*zb(j1,j6)*zb(j4,j6))-
     .za(j2,j4)*zb(j2,j6)*zb(j4,j6))/zb(j5,j6)))/
     .(2d0*zb(j1,j4)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4)))+
     .(za(j1,j4)*(-(za(j2,j4)*((za(j4,j5)*
     .(-((za(j3,j5)*
     .((s(j1,j4)-s(j2,j3)-s(j5,j6))*za(j1,j4)*
     .zb(j1,j3)-
     .(-s(j1,j4)+s(j2,j3)-s(j5,j6))*za(j2,j4)*
     .zb(j2,j3)))/za(j5,j6))+
     .((-s(j1,j4)-s(j2,j3)+s(j5,j6))*za(j2,j4)+
     .2d0*za(j1,j4)*za(j2,j3)*zb(j1,j3))*zb(j2,j6)))/
     .za(j1,j4)-(zb(j1,j6)*
     .(-(za(j3,j5)*
     .(-((-s(j1,j4)-s(j2,j3)+s(j5,j6))*zb(j1,j3))-
     .2d0*za(j2,j4)*zb(j1,j4)*zb(j2,j3)))-
     .((-((-s(j1,j4)+s(j2,j3)-s(j5,j6))*za(j2,j3)*
     .zb(j1,j3))+
     .(s(j1,j4)-s(j2,j3)-s(j5,j6))*za(j2,j4)*
     .zb(j1,j4))*zb(j2,j6))/zb(j5,j6)))/zb(j1,j4)))/
     .(2d0*(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)
     .*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4)))+
     .(za(j2,j4)*zb(j1,j2)*(-(za(j1,j2)*zb(j1,j6))+
     .za(j2,j4)*zb(j4,j6))*
     .((L0(-t(j1,j2,j4),-s(j1,j4))*
     .((-2d0*za(j2,j4)*zb(j3,j6))/
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))+
     .(-(za(j1,j4)*zb(j1,j6))-za(j2,j4)*zb(j2,j6))/
     .t(j1,j2,j4)))/s(j1,j4)-
     .(L1(-s(j1,j4),-t(j1,j2,j4))*za(j2,j4)*zb(j2,j6))/t(j1,j2,j4)**2))/
     .(2d0*za(j1,j2)*(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*
     .zb(j5,j6))+(za(j2,j4)*za(j3,j4)*zb(j3,j6)**2*
     .((-2d0*L0(-s(j5,j6),-t(j1,j2,j4))*za(j2,j4))/
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))+
     .(L1(-s(j5,j6),-t(j1,j2,j4))*za(j3,j4))/t(j1,j2,j4)))/
     .(2d0*za(j1,j2)*za(j1,j4)*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*zb(j5,j6))-
     .(Lsm1_2mh(s(j2,j3),t(j1,j2,j4),s(j1,j4),s(j5,j6))*za(j2,j4)**3*
     .zb(j3,j6)**2*t(j1,j2,j4))/
     .(za(j1,j2)*za(j1,j4)*(-(za(j1,j2)*zb(j1,j3))-
     .za(j2,j4)*zb(j3,j4))**3*zb(j5,j6))-
     .(Lnrat(-s(j2,j3),-s(j5,j6))*za(j2,j4)*
     .(-((za(j2,j5)*za(j4,j5)*zb(j1,j2))/za(j5,j6))+
     .(2d0*(-s(j1,j4)-s(j2,j3)+s(j5,j6))*
     .(za(j2,j5)*za(j3,j4)+za(j2,j4)*za(j3,j5))*zb(j1,j6)*
     .zb(j2,j3))/
     .(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)+
     .((-s(j1,j4)+s(j2,j3)-s(j5,j6))*zb(j2,j3)*
     .(za(j2,j5)*za(j3,j4)*
     .(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))-
     .za(j2,j3)*za(j4,j5)*
     .(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))+
     .za(j2,j4)*zb(j1,j3)*
     .(za(j3,j5)**2+
     .(za(j2,j5)**2*
     .(-(za(j1,j3)*zb(j1,j2))-za(j3,j4)*zb(j2,j4)))/
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4)))))/
     .((s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)*
     .za(j5,j6))+(za(j2,j4)*zb(j1,j3)*
     .(-((za(j2,j5)*za(j3,j5)*zb(j2,j3))/za(j5,j6))-
     .(za(j2,j3)*zb(j2,j6)*zb(j3,j6))/zb(j5,j6)))/
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))+
     .((-s(j1,j4)+s(j2,j3)-s(j5,j6))*za(j2,j3)*
     .(-((za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .zb(j2,j6)*zb(j3,j6))-
     .(za(j2,j4)*zb(j1,j3)*
     .(-(za(j1,j3)*zb(j1,j2))-za(j3,j4)*zb(j2,j4))*
     .zb(j3,j6)**2)/
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))-
     .za(j2,j4)*zb(j2,j6)*
     .(-(zb(j1,j6)*zb(j2,j3))+zb(j1,j2)*zb(j3,j6))))/
     .((s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)*
     .zb(j5,j6))-(s(j2,j3)*
     .(-((za(j4,j5)**2*zb(j1,j4))/za(j5,j6))-
     .(za(j1,j4)*zb(j1,j6)**2)/zb(j5,j6))*t(j1,j2,j4))/
     .(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)+
     .(6d0*s(j1,j4)*s(j2,j3)*
     .(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(-(za(j2,j5)*zb(j2,j6))-za(j3,j5)*zb(j3,j6))*
     .(-t(j1,j2,j4)+t(j1,j3,j4)))/
     .(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)**2))/
     .(2d0*s(j1,j4)*za(j1,j2)*(-(za(j1,j2)*zb(j1,j3))-
     .za(j2,j4)*zb(j3,j4)))-
     .(I3m(s(j2,j3),s(j1,j4),s(j5,j6))*za(j2,j4)*
     .((za(j2,j4)*zb(j1,j3)*
     .(-(za(j1,j3)*zb(j1,j2))-za(j3,j4)*zb(j2,j4))*
     .(za(j1,j2)*(za(j2,j5)*zb(j1,j2)+
     .za(j3,j5)*zb(j1,j3))-
     .za(j2,j4)*(-(za(j2,j5)*zb(j2,j4))-
     .za(j3,j5)*zb(j3,j4))+
     .za(j2,j3)*(-(za(j1,j5)*zb(j1,j3))+
     .za(j4,j5)*zb(j3,j4)))*zb(j3,j6))/
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))-
     .zb(j1,j2)*(za(j2,j4)*
     .(za(j1,j3)*(za(j2,j5)*zb(j1,j2)+
     .za(j3,j5)*zb(j1,j3))-
     .za(j2,j3)*
     .(-(za(j1,j5)*zb(j1,j2))+za(j4,j5)*zb(j2,j4))-
     .za(j3,j4)*
     .(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4)))+
     .za(j2,j3)*(-(za(j2,j4)*
     .(-(za(j1,j5)*zb(j1,j2))+za(j4,j5)*zb(j2,j4)))-
     .za(j3,j4)*
     .(-(za(j1,j5)*zb(j1,j3))+za(j4,j5)*zb(j3,j4))))*
     .zb(j3,j6)+(s(j2,j3)*
     .(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(-(za(j2,j5)*zb(j2,j6))-za(j3,j5)*zb(j3,j6)))/2+
     .za(j3,j4)*zb(j2,j3)*
     .(za(j2,j3)*za(j4,j5)*zb(j1,j4)*zb(j3,j6)+
     .(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j2,j5)*zb(j2,j6))-za(j3,j5)*zb(j3,j6)))+
     .(3d0*s(j2,j3)*(-s(j1,j4)+s(j2,j3)-s(j5,j6))*
     .(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(-(za(j2,j5)*zb(j2,j6))-za(j3,j5)*zb(j3,j6))*
     .(-t(j1,j2,j4)+t(j1,j3,j4)))/
     .(2d0*(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2))))/
     .((s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)
     .*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4)))-
     .(Lnrat(-s(j1,j4),-s(j5,j6))*za(j2,j4)*
     .(4d0*(za(j2,j5)*za(j3,j4)*zb(j1,j6)*zb(j2,j3)+
     .za(j2,j3)*za(j4,j5)*zb(j1,j3)*zb(j2,j6))-
     .(8d0*s(j2,j3)*za(j2,j4)*za(j5,j6)*zb(j1,j6)*zb(j3,j6))/
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))-
     .2d0*(-(za(j2,j4)*zb(j1,j2))+za(j3,j4)*zb(j1,j3))*
     .(-(za(j2,j5)*zb(j2,j6))-za(j3,j5)*zb(j3,j6))+
     .(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3)-
     .(4d0*s(j5,j6)*za(j2,j4)*zb(j1,j3))/
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4)))*
     .((za(j2,j5)*za(j3,j5)*zb(j2,j3))/za(j5,j6)-
     .(za(j2,j3)*zb(j2,j6)*zb(j3,j6))/zb(j5,j6))-
     .(s(j1,j4)-s(j2,j3)-s(j5,j6))*
     .(-((za(j2,j5)*za(j4,j5)*zb(j1,j2))/za(j5,j6))-
     .(za(j3,j4)*zb(j1,j6)*zb(j3,j6))/zb(j5,j6)-
     .(2d0*za(j2,j4)*zb(j3,j6)*
     .(za(j2,j3)*zb(j1,j2)*zb(j3,j6)-
     .za(j3,j4)*zb(j1,j4)*zb(j3,j6)-
     .za(j2,j5)*zb(j1,j2)*zb(j5,j6)))/
     .((-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*
     .zb(j5,j6)))+
     .(3d0*(-s(j1,j4)-s(j2,j3)+s(j5,j6))*
     .(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(-(za(j2,j5)*zb(j2,j6))-za(j3,j5)*zb(j3,j6))*
     .(-t(j1,j2,j4)+t(j1,j3,j4)))/
     .(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)))/
     .(2d0*(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)
     .*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4)))))/za(j2,j4)
      Fsc4=Fsc4-
     .(za(j1,j5)*zb(j1,j2)*((L0(-t(j2,j3,j4),-s(j5,j6))*
     .(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4)))/(s(j5,j6)*za(j5,j6))
     .+(L1(-s(j5,j6),-t(j2,j3,j4))*(-(za(j1,j3)*zb(j1,j6)*zb(j2,j3))-
     .za(j1,j4)*zb(j1,j6)*zb(j2,j4)))/t(j2,j3,j4)**2))/
     .(2d0*zb(j3,j4)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))+
     .(za(j3,j4)*zb(j2,j3)*(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4))*
     .((L1(-t(j2,j3,j4),-s(j2,j3))*za(j4,j5)*zb(j2,j4))/s(j2,j3)**2+
     .(L0(-t(j2,j3,j4),-s(j2,j3))*
     .(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4)))/(s(j2,j3)*t(j2,j3,j4))
     .))/(2d0*za(j5,j6)*zb(j3,j4)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j
     .4)))+
     .(-((za(j1,j3)*za(j3,j5)*za(j4,j5))/
     .(za(j1,j2)*za(j2,j3)*za(j5,j6)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))))-
     .(za(j1,j5)*za(j3,j5)*zb(j1,j2))/
     .(za(j2,j3)*za(j5,j6)*zb(j3,j4)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))-
     .(za(j3,j5)**2*zb(j1,j3)**2)/
     .(za(j2,j3)*za(j5,j6)*zb(j1,j4)*zb(j3,j4)*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4)))-
     .(s(j3,j4)*za(j1,j5)*za(j4,j5))/
     .(za(j1,j2)*za(j5,j6)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4)))+
     .(za(j1,j5)*za(j3,j4)*za(j3,j5)*zb(j1,j3))/
     .(za(j2,j3)*za(j5,j6)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4)))-
     .(za(j1,j3)*zb(j1,j6)**2)/
     .(za(j1,j2)*za(j2,j3)*zb(j1,j4)*zb(j3,j4)*zb(j5,j6))-
     .(za(j1,j3)*zb(j1,j2)*zb(j1,j6)*zb(j4,j6))/
     .(za(j2,j3)*zb(j1,j4)*zb(j3,j4)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*zb(j5,j6))+
     .(za(j1,j4)*za(j3,j4)*zb(j3,j6)*zb(j4,j6))/
     .(za(j1,j2)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*zb(j5,j6))-
     .(s(j1,j3)*za(j3,j4)*zb(j3,j6)*zb(j4,j6))/
     .(za(j2,j3)*zb(j3,j4)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*zb(j5,j6))+
     .(za(j3,j4)*(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6))*zb(j4,j6))/
     .(za(j1,j2)*za(j2,j3)*zb(j3,j4)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*zb(j5,j6))-
     .(za(j2,j4)*zb(j2,j6)*(-(za(j1,j4)*zb(j1,j6))-za(j2,j4)*zb(j2,j6)))
     ./
     .(za(j1,j2)*(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*zb(j5,j6)*
     .t(j1,j2,j4))-(za(j1,j5)*zb(j1,j2)*
     .(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4)))/
     .(za(j5,j6)*zb(j3,j4)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*
     .t(j2,j3,j4)))/2+(za(j1,j3)*
     .((za(j4,j5)*zb(j1,j6)*(-((-s(j1,j4)-s(j2,j3)+s(j5,j6))*za(j1,j3))-
     .2d0*za(j1,j4)*za(j2,j3)*zb(j2,j4)))/
     .((s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-2d0*s(j1,j4)*s(j5,
     .j6)-
     .2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))+
     .(((s(j1,j4)-s(j2,j3)-s(j5,j6))*za(j1,j3)*zb(j1,j4)-
     .(-s(j1,j4)+s(j2,j3)-s(j5,j6))*za(j2,j3)*zb(j2,j4))*
     .(-((za(j4,j5)**2*zb(j1,j4))/za(j5,j6))-
     .(za(j1,j4)*zb(j1,j6)**2)/zb(j5,j6)))/
     .(2d0*(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)
     .*
     .zb(j1,j4)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))+
     .(Lsm1_2mh(s(j1,j4),t(j1,j2,j3),s(j2,j3),s(j5,j6))*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6))**2)/
     .(za(j1,j2)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))**3*zb(j5,j6))
     .-(L0(-t(j1,j2,j3),-s(j2,j3))*za(j2,j3)*zb(j1,j2)*
     .(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6))**2)/
     .(s(j2,j3)*za(j1,j2)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))**2*
     .zb(j5,j6))+(Lnrat(-s(j2,j3),-s(j5,j6))*za(j2,j3)*
     .((-6d0*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*zb(j2,j3)*
     .(-((-s(j1,j4)-s(j2,j3)+s(j5,j6))*za(j1,j3))-
     .2d0*za(j1,j4)*za(j2,j3)*zb(j2,j4))*
     .(-(za(j2,j5)*zb(j2,j6))-za(j3,j5)*zb(j3,j6)))/
     .(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)**2-
     .((za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6))*
     .(zb(j1,j4)*zb(j2,j6)+zb(j1,j2)*zb(j4,j6)))/
     .(zb(j1,j4)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*
     .zb(j5,j6))+(-((za(j4,j5)*
     .(za(j4,j5)*zb(j1,j4)*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))
     .+2d0*(-((s(j1,j4)-s(j2,j3)-s(j5,j6))*za(j1,j5)*zb(j1,j2))+
     .(-s(j1,j4)-s(j2,j3)+s(j5,j6))*za(j5,j6)*
     .zb(j2,j6))))/za(j5,j6))+
     .((-s(j1,j4)+s(j2,j3)-s(j5,j6))*za(j1,j3)*zb(j2,j3)*
     .(-(za(j1,j5)*zb(j1,j6))+za(j4,j5)*zb(j4,j6)))/
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))-
     .(za(j3,j4)*zb(j1,j6)*zb(j2,j3)*
     .(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6)))/zb(j5,j6)
     .-(zb(j2,j3)*zb(j4,j6)*(-((2d0*(-s(j1,j4)+s(j2,j3)-s(j5,j6))*
     .za(j3,j4)-za(j1,j4)*za(j2,j3)*zb(j1,j2))*
     .zb(j1,j6))-
     .((-s(j1,j4)+s(j2,j3)-s(j5,j6))*za(j1,j3)*
     .(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .zb(j4,j6))/
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))))/
     .(zb(j1,j4)*zb(j5,j6))+
     .((-s(j1,j4)+s(j2,j3)-s(j5,j6))*zb(j1,j6)*zb(j2,j4)*
     .(-(za(j2,j4)*zb(j2,j6))-za(j3,j4)*zb(j3,j6)-
     .za(j4,j5)*zb(j5,j6)))/(zb(j1,j4)*zb(j5,j6))+
     .za(j4,j5)*(za(j1,j3)*zb(j1,j6)*zb(j2,j3)+
     .za(j4,j5)*zb(j2,j4)*zb(j5,j6))+
     .4d0*za(j1,j5)*zb(j2,j3)*
     .(za(j3,j4)*zb(j1,j6)+
     .(za(j1,j3)*za(j4,j5)*zb(j1,j4)*zb(j5,j6))/
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))))/
     .(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)))/
     .(2d0*za(j1,j2)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))-
     .(L1(-s(j5,j6),-t(j1,j2,j3))*za(j3,j4)*za(j4,j5)*zb(j4,j6))/
     .(za(j1,j2)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*t(j1,j2,j3))
     .-(L0(-t(j1,j2,j3),-s(j5,j6))*za(j3,j4)*
     .(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6))*zb(j4,j6)*t(j1,j2,j3)
     .)/(s(j5,j6)*za(j1,j2)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))**2
     .*
     .zb(j5,j6))+(Lnrat(-s(j1,j4),-s(j5,j6))*
     .((-6d0*za(j1,j4)*za(j2,j3)*
     .(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(2d0*za(j1,j3)*zb(j1,j4)*zb(j2,j3)+
     .(-s(j1,j4)-s(j2,j3)+s(j5,j6))*zb(j2,j4))*
     .(-(za(j2,j5)*zb(j2,j6))-za(j3,j5)*zb(j3,j6)))/
     .(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)+
     .(za(j4,j5)*(-3d0*(s(j1,j4)-s(j2,j3)-s(j5,j6))*za(j1,j3)*
     .za(j4,j5)*zb(j1,j4)+
     .za(j1,j3)*za(j1,j4)*
     .(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))*zb(j1,j4)
     .+za(j1,j4)*za(j2,j3)*(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))*
     .zb(j2,j4)+
     .2d0*za(j2,j3)*zb(j1,j2)*
     .((-s(j1,j4)+s(j2,j3)-s(j5,j6))*za(j1,j5)+
     .2d0*za(j1,j4)*za(j5,j6)*zb(j4,j6))))/za(j5,j6)+
     .(za(j1,j4)*zb(j1,j6)*
     .((-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(-(za(j2,j4)*zb(j2,j6))-za(j3,j4)*zb(j3,j6))+
     .2d0*((s(j1,j4)-s(j2,j3)-s(j5,j6))*za(j3,j4)*zb(j4,j6)-
     .(-s(j1,j4)-s(j2,j3)+s(j5,j6))*za(j3,j5)*zb(j5,j6)))
     .)/zb(j5,j6)-(za(j1,j3)*((za(j1,j5)*za(j4,j5)*zb(j1,j4))/za(j5,j6)+
     .(za(j1,j4)*zb(j1,j6)*zb(j4,j6))/zb(j5,j6))*
     .(-s(j1,j4)**2+2d0*s(j1,j4)*s(j2,j3)-s(j2,j3)**2+
     .2d0*s(j1,j4)*s(j5,j6)+2d0*s(j2,j3)*s(j5,j6)-s(j5,j6)**2+
     .(s(j1,j4)-s(j2,j3)-s(j5,j6))*(t(j1,j2,j3)-t(j2,j3,j4)))
     .)/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))))/
     .(2d0*(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)
     .*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))+
     .(I3m(s(j2,j3),s(j1,j4),s(j5,j6))*za(j2,j3)*
     .((-3d0*za(j1,j4)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .zb(j2,j3)*((s(j1,j4)-s(j2,j3)-s(j5,j6))*za(j1,j3)*
     .zb(j1,j4)-
     .(-s(j1,j4)+s(j2,j3)-s(j5,j6))*za(j2,j3)*zb(j2,j4))*
     .(-(za(j2,j5)*zb(j2,j6))-za(j3,j5)*zb(j3,j6)))/
     .(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)**2-
     .(za(j4,j5)*zb(j1,j2)*
     .(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6)))/t(j1,j2,j3)-
     .(-(za(j3,j4)*zb(j1,j6)*zb(j2,j3)*
     .((-s(j1,j4)+s(j2,j3)-s(j5,j6))*za(j1,j5)-
     .za(j4,j5)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))+
     .2d0*za(j1,j4)*za(j5,j6)*zb(j4,j6)))+
     .(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(-(za(j1,j4)*zb(j2,j4)*
     .(-(za(j2,j5)*zb(j2,j6))-za(j3,j5)*zb(j3,j6)))+
     .3d0*za(j1,j3)*za(j4,j5)*zb(j2,j3)*zb(j4,j6))-
     .za(j4,j5)*zb(j1,j2)*
     .((s(j1,j4)-s(j2,j3)-s(j5,j6))*za(j1,j4)*zb(j4,j6)-
     .(-s(j1,j4)-s(j2,j3)+s(j5,j6))*za(j1,j5)*zb(j5,j6))-
     .(za(j1,j3)*za(j1,j5)*zb(j2,j3)*
     .(za(j2,j4)*zb(j1,j2)*zb(j4,j6)+
     .za(j3,j4)*zb(j1,j3)*zb(j4,j6))*
     .(-t(j1,j2,j3)+t(j2,j3,j4)))/
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))/
     .(s(j1,j4)**2-2d0*s(j1,j4)*s(j2,j3)+s(j2,j3)**2-
     .2d0*s(j1,j4)*s(j5,j6)-2d0*s(j2,j3)*s(j5,j6)+s(j5,j6)**2)))/
     .(za(j1,j2)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))))/za(j2,j3)
      return
      end
