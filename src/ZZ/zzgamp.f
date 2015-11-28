      subroutine zzgamp(p1,p2,p3,p4,p5,p6,p7,za,zb,zaa,z34,z56)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'srdiags.f'
      double precision tinv,s12,s34,s56
      double complex zaa(2,2,2,2),z34(2,2,2,2),z56(2,2,2,2),zab2
      double complex iza(7,7),izb(7,7)
      integer p1,p2,p3,p4,p5,p6,p7,i,j
C--   order of indices is gluon helicity 
      tinv(p1,p2,p3)=1d0/(s(p1,p2)+s(p2,p3)+s(p1,p3))
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      do i=1,7
      do j=i+1,7
         iza(i,j)=cone/za(i,j)
         izb(i,j)=cone/zb(i,j)
         iza(j,i)=-iza(i,j)
         izb(j,i)=-izb(i,j)
      enddo
      enddo
      zaa(2,2,2,2)= + tinv(p1,p3,p4)*tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * (  - za(p2,p6)**2*zb(p1,p3)*zb(p5,p6)*iza(p2,p7)*zab2(p4,p1,
     &    p3,p7) )
      zaa(2,2,2,2) = zaa(2,2,2,2) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * (  - za(p2,p4)**2*zb(p1,p5)*zb(p3,p4)*iza(
     &    p2,p7)*zab2(p6,p1,p5,p7) )
      zaa(2,2,2,2) = zaa(2,2,2,2) + tinv(p2,p3,p4)*s34**(-1)*s56**(-1)
     &  * (  - za(p2,p4)*iza(p1,p7)*iza(p2,p7)*zab2(p2,p1,p7,p5)*zab2(
     &    p6,p2,p4,p3) )
      zaa(2,2,2,2) = zaa(2,2,2,2) + tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * (  - za(p2,p6)*iza(p1,p7)*iza(p2,p7)*zab2(p2,p1,p7,p3)*zab2(
     &    p4,p2,p6,p5) )

      zaa(1,1,1,1)= + tinv(p1,p3,p4)*tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * ( za(p1,p3)*za(p5,p6)*zb(p2,p6)**2*izb(p2,p7)*zab2(p7,p1,p3,
     &    p4) )
      zaa(1,1,1,1) = zaa(1,1,1,1) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * ( za(p1,p5)*za(p3,p4)*zb(p2,p4)**2*izb(p2,
     &    p7)*zab2(p7,p1,p5,p6) )
      zaa(1,1,1,1) = zaa(1,1,1,1) + tinv(p2,p3,p4)*s34**(-1)*s56**(-1)
     &  * ( zb(p2,p4)*izb(p1,p7)*izb(p2,p7)*zab2(p3,p2,p4,p6)*zab2(p5,
     &    p1,p7,p2) )
      zaa(1,1,1,1) = zaa(1,1,1,1) + tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * ( zb(p2,p6)*izb(p1,p7)*izb(p2,p7)*zab2(p3,p1,p7,p2)*zab2(p5,
     &    p2,p6,p4) )

      zaa(1,2,2,2)= + tinv(p1,p3,p4)*s34**(-1)*s56**(-1) * (  - za(p1,
     &    p4)*iza(p1,p7)*iza(p2,p7)*zab2(p1,p2,p7,p5)*zab2(p6,p1,p4,p3)
     &     )
      zaa(1,2,2,2) = zaa(1,2,2,2) + tinv(p1,p3,p4)*tinv(p2,p5,p6)*
     & s34**(-1)*s56**(-1) * (  - za(p1,p4)**2*zb(p2,p5)*zb(p3,p4)*iza(
     &    p1,p7)*zab2(p6,p2,p5,p7) )
      zaa(1,2,2,2) = zaa(1,2,2,2) + tinv(p1,p5,p6)*s34**(-1)*s56**(-1)
     &  * (  - za(p1,p6)*iza(p1,p7)*iza(p2,p7)*zab2(p1,p2,p7,p3)*zab2(
     &    p4,p1,p6,p5) )
      zaa(1,2,2,2) = zaa(1,2,2,2) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * (  - za(p1,p6)**2*zb(p2,p3)*zb(p5,p6)*iza(
     &    p1,p7)*zab2(p4,p2,p3,p7) )

      zaa(2,1,2,2)= + tinv(p1,p3,p4)*tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * (  - za(p2,p6)**2*zb(p1,p4)*zb(p5,p6)*iza(p2,p7)*zab2(p3,p1,
     &    p4,p7) )
      zaa(2,1,2,2) = zaa(2,1,2,2) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * ( za(p2,p3)**2*zb(p1,p5)*zb(p3,p4)*iza(p2,
     &    p7)*zab2(p6,p1,p5,p7) )
      zaa(2,1,2,2) = zaa(2,1,2,2) + tinv(p2,p3,p4)*s34**(-1)*s56**(-1)
     &  * (  - za(p2,p3)*iza(p1,p7)*iza(p2,p7)*zab2(p2,p1,p7,p5)*zab2(
     &    p6,p2,p3,p4) )
      zaa(2,1,2,2) = zaa(2,1,2,2) + tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * (  - za(p2,p6)*iza(p1,p7)*iza(p2,p7)*zab2(p2,p1,p7,p4)*zab2(
     &    p3,p2,p6,p5) )

      zaa(2,2,1,2)= + tinv(p1,p3,p4)*tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * ( za(p2,p5)**2*zb(p1,p3)*zb(p5,p6)*iza(p2,p7)*zab2(p4,p1,p3,
     &    p7) )
      zaa(2,2,1,2) = zaa(2,2,1,2) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * (  - za(p2,p4)**2*zb(p1,p6)*zb(p3,p4)*iza(
     &    p2,p7)*zab2(p5,p1,p6,p7) )
      zaa(2,2,1,2) = zaa(2,2,1,2) + tinv(p2,p3,p4)*s34**(-1)*s56**(-1)
     &  * (  - za(p2,p4)*iza(p1,p7)*iza(p2,p7)*zab2(p2,p1,p7,p6)*zab2(
     &    p5,p2,p4,p3) )
      zaa(2,2,1,2) = zaa(2,2,1,2) + tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * (  - za(p2,p5)*iza(p1,p7)*iza(p2,p7)*zab2(p2,p1,p7,p3)*zab2(
     &    p4,p2,p5,p6) )

      zaa(2,2,2,1)= + tinv(p1,p3,p4)*s34**(-1)*s56**(-1) * ( zb(p1,p3)*
     &    izb(p1,p7)*izb(p2,p7)*zab2(p4,p1,p3,p5)*zab2(p6,p2,p7,p1) )
      zaa(2,2,2,1) = zaa(2,2,2,1) + tinv(p1,p3,p4)*tinv(p2,p5,p6)*
     & s34**(-1)*s56**(-1) * (  - za(p2,p6)*za(p3,p4)*zb(p1,p3)**2*izb(
     &    p1,p7)*zab2(p7,p2,p6,p5) )
      zaa(2,2,2,1) = zaa(2,2,2,1) + tinv(p1,p5,p6)*s34**(-1)*s56**(-1)
     &  * ( zb(p1,p5)*izb(p1,p7)*izb(p2,p7)*zab2(p4,p2,p7,p1)*zab2(p6,
     &    p1,p5,p3) )
      zaa(2,2,2,1) = zaa(2,2,2,1) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * (  - za(p2,p4)*za(p5,p6)*zb(p1,p5)**2*izb(
     &    p1,p7)*zab2(p7,p2,p4,p3) )

      zaa(1,1,2,2)= + tinv(p1,p3,p4)*s34**(-1)*s56**(-1) * (  - za(p1,
     &    p3)*iza(p1,p7)*iza(p2,p7)*zab2(p1,p2,p7,p5)*zab2(p6,p1,p3,p4)
     &     )
      zaa(1,1,2,2) = zaa(1,1,2,2) + tinv(p1,p3,p4)*tinv(p2,p5,p6)*
     & s34**(-1)*s56**(-1) * ( za(p1,p3)**2*zb(p2,p5)*zb(p3,p4)*iza(p1,
     &    p7)*zab2(p6,p2,p5,p7) )
      zaa(1,1,2,2) = zaa(1,1,2,2) + tinv(p1,p5,p6)*s34**(-1)*s56**(-1)
     &  * (  - za(p1,p6)*iza(p1,p7)*iza(p2,p7)*zab2(p1,p2,p7,p4)*zab2(
     &    p3,p1,p6,p5) )
      zaa(1,1,2,2) = zaa(1,1,2,2) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * (  - za(p1,p6)**2*zb(p2,p4)*zb(p5,p6)*iza(
     &    p1,p7)*zab2(p3,p2,p4,p7) )

      zaa(1,2,1,2)= + tinv(p1,p3,p4)*s34**(-1)*s56**(-1) * (  - za(p1,
     &    p4)*iza(p1,p7)*iza(p2,p7)*zab2(p1,p2,p7,p6)*zab2(p5,p1,p4,p3)
     &     )
      zaa(1,2,1,2) = zaa(1,2,1,2) + tinv(p1,p3,p4)*tinv(p2,p5,p6)*
     & s34**(-1)*s56**(-1) * (  - za(p1,p4)**2*zb(p2,p6)*zb(p3,p4)*iza(
     &    p1,p7)*zab2(p5,p2,p6,p7) )
      zaa(1,2,1,2) = zaa(1,2,1,2) + tinv(p1,p5,p6)*s34**(-1)*s56**(-1)
     &  * (  - za(p1,p5)*iza(p1,p7)*iza(p2,p7)*zab2(p1,p2,p7,p3)*zab2(
     &    p4,p1,p5,p6) )
      zaa(1,2,1,2) = zaa(1,2,1,2) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * ( za(p1,p5)**2*zb(p2,p3)*zb(p5,p6)*iza(p1,
     &    p7)*zab2(p4,p2,p3,p7) )

      zaa(1,2,2,1)= + tinv(p1,p3,p4)*tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * (  - za(p1,p4)*za(p5,p6)*zb(p2,p5)**2*izb(p2,p7)*zab2(p7,p1,
     &    p4,p3) )
      zaa(1,2,2,1) = zaa(1,2,2,1) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * (  - za(p1,p6)*za(p3,p4)*zb(p2,p3)**2*izb(
     &    p2,p7)*zab2(p7,p1,p6,p5) )
      zaa(1,2,2,1) = zaa(1,2,2,1) + tinv(p2,p3,p4)*s34**(-1)*s56**(-1)
     &  * ( zb(p2,p3)*izb(p1,p7)*izb(p2,p7)*zab2(p4,p2,p3,p5)*zab2(p6,
     &    p1,p7,p2) )
      zaa(1,2,2,1) = zaa(1,2,2,1) + tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * ( zb(p2,p5)*izb(p1,p7)*izb(p2,p7)*zab2(p4,p1,p7,p2)*zab2(p6,
     &    p2,p5,p3) )

      zaa(2,1,1,2)= + tinv(p1,p3,p4)*tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * ( za(p2,p5)**2*zb(p1,p4)*zb(p5,p6)*iza(p2,p7)*zab2(p3,p1,p4,
     &    p7) )
      zaa(2,1,1,2) = zaa(2,1,1,2) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * ( za(p2,p3)**2*zb(p1,p6)*zb(p3,p4)*iza(p2,
     &    p7)*zab2(p5,p1,p6,p7) )
      zaa(2,1,1,2) = zaa(2,1,1,2) + tinv(p2,p3,p4)*s34**(-1)*s56**(-1)
     &  * (  - za(p2,p3)*iza(p1,p7)*iza(p2,p7)*zab2(p2,p1,p7,p6)*zab2(
     &    p5,p2,p3,p4) )
      zaa(2,1,1,2) = zaa(2,1,1,2) + tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * (  - za(p2,p5)*iza(p1,p7)*iza(p2,p7)*zab2(p2,p1,p7,p4)*zab2(
     &    p3,p2,p5,p6) )

      zaa(2,1,2,1)= + tinv(p1,p3,p4)*s34**(-1)*s56**(-1) * ( zb(p1,p4)*
     &    izb(p1,p7)*izb(p2,p7)*zab2(p3,p1,p4,p5)*zab2(p6,p2,p7,p1) )
      zaa(2,1,2,1) = zaa(2,1,2,1) + tinv(p1,p3,p4)*tinv(p2,p5,p6)*
     & s34**(-1)*s56**(-1) * ( za(p2,p6)*za(p3,p4)*zb(p1,p4)**2*izb(p1,
     &    p7)*zab2(p7,p2,p6,p5) )
      zaa(2,1,2,1) = zaa(2,1,2,1) + tinv(p1,p5,p6)*s34**(-1)*s56**(-1)
     &  * ( zb(p1,p5)*izb(p1,p7)*izb(p2,p7)*zab2(p3,p2,p7,p1)*zab2(p6,
     &    p1,p5,p4) )
      zaa(2,1,2,1) = zaa(2,1,2,1) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * (  - za(p2,p3)*za(p5,p6)*zb(p1,p5)**2*izb(
     &    p1,p7)*zab2(p7,p2,p3,p4) )

      zaa(2,2,1,1)= + tinv(p1,p3,p4)*s34**(-1)*s56**(-1) * ( zb(p1,p3)*
     &    izb(p1,p7)*izb(p2,p7)*zab2(p4,p1,p3,p6)*zab2(p5,p2,p7,p1) )
      zaa(2,2,1,1) = zaa(2,2,1,1) + tinv(p1,p3,p4)*tinv(p2,p5,p6)*
     & s34**(-1)*s56**(-1) * (  - za(p2,p5)*za(p3,p4)*zb(p1,p3)**2*izb(
     &    p1,p7)*zab2(p7,p2,p5,p6) )
      zaa(2,2,1,1) = zaa(2,2,1,1) + tinv(p1,p5,p6)*s34**(-1)*s56**(-1)
     &  * ( zb(p1,p6)*izb(p1,p7)*izb(p2,p7)*zab2(p4,p2,p7,p1)*zab2(p5,
     &    p1,p6,p3) )
      zaa(2,2,1,1) = zaa(2,2,1,1) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * ( za(p2,p4)*za(p5,p6)*zb(p1,p6)**2*izb(p1,
     &    p7)*zab2(p7,p2,p4,p3) )

      zaa(2,1,1,1)= + tinv(p1,p3,p4)*s34**(-1)*s56**(-1) * ( zb(p1,p4)*
     &    izb(p1,p7)*izb(p2,p7)*zab2(p3,p1,p4,p6)*zab2(p5,p2,p7,p1) )
      zaa(2,1,1,1) = zaa(2,1,1,1) + tinv(p1,p3,p4)*tinv(p2,p5,p6)*
     & s34**(-1)*s56**(-1) * ( za(p2,p5)*za(p3,p4)*zb(p1,p4)**2*izb(p1,
     &    p7)*zab2(p7,p2,p5,p6) )
      zaa(2,1,1,1) = zaa(2,1,1,1) + tinv(p1,p5,p6)*s34**(-1)*s56**(-1)
     &  * ( zb(p1,p6)*izb(p1,p7)*izb(p2,p7)*zab2(p3,p2,p7,p1)*zab2(p5,
     &    p1,p6,p4) )
      zaa(2,1,1,1) = zaa(2,1,1,1) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * ( za(p2,p3)*za(p5,p6)*zb(p1,p6)**2*izb(p1,
     &    p7)*zab2(p7,p2,p3,p4) )

      zaa(1,2,1,1)= + tinv(p1,p3,p4)*tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * ( za(p1,p4)*za(p5,p6)*zb(p2,p6)**2*izb(p2,p7)*zab2(p7,p1,p4,
     &    p3) )
      zaa(1,2,1,1) = zaa(1,2,1,1) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * (  - za(p1,p5)*za(p3,p4)*zb(p2,p3)**2*izb(
     &    p2,p7)*zab2(p7,p1,p5,p6) )
      zaa(1,2,1,1) = zaa(1,2,1,1) + tinv(p2,p3,p4)*s34**(-1)*s56**(-1)
     &  * ( zb(p2,p3)*izb(p1,p7)*izb(p2,p7)*zab2(p4,p2,p3,p6)*zab2(p5,
     &    p1,p7,p2) )
      zaa(1,2,1,1) = zaa(1,2,1,1) + tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * ( zb(p2,p6)*izb(p1,p7)*izb(p2,p7)*zab2(p4,p1,p7,p2)*zab2(p5,
     &    p2,p6,p3) )

      zaa(1,1,2,1)= + tinv(p1,p3,p4)*tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * (  - za(p1,p3)*za(p5,p6)*zb(p2,p5)**2*izb(p2,p7)*zab2(p7,p1,
     &    p3,p4) )
      zaa(1,1,2,1) = zaa(1,1,2,1) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * ( za(p1,p6)*za(p3,p4)*zb(p2,p4)**2*izb(p2,
     &    p7)*zab2(p7,p1,p6,p5) )
      zaa(1,1,2,1) = zaa(1,1,2,1) + tinv(p2,p3,p4)*s34**(-1)*s56**(-1)
     &  * ( zb(p2,p4)*izb(p1,p7)*izb(p2,p7)*zab2(p3,p2,p4,p5)*zab2(p6,
     &    p1,p7,p2) )
      zaa(1,1,2,1) = zaa(1,1,2,1) + tinv(p2,p5,p6)*s34**(-1)*s56**(-1)
     &  * ( zb(p2,p5)*izb(p1,p7)*izb(p2,p7)*zab2(p3,p1,p7,p2)*zab2(p6,
     &    p2,p5,p4) )

      zaa(1,1,1,2)= + tinv(p1,p3,p4)*s34**(-1)*s56**(-1) * (  - za(p1,
     &    p3)*iza(p1,p7)*iza(p2,p7)*zab2(p1,p2,p7,p6)*zab2(p5,p1,p3,p4)
     &     )
      zaa(1,1,1,2) = zaa(1,1,1,2) + tinv(p1,p3,p4)*tinv(p2,p5,p6)*
     & s34**(-1)*s56**(-1) * ( za(p1,p3)**2*zb(p2,p6)*zb(p3,p4)*iza(p1,
     &    p7)*zab2(p5,p2,p6,p7) )
      zaa(1,1,1,2) = zaa(1,1,1,2) + tinv(p1,p5,p6)*s34**(-1)*s56**(-1)
     &  * (  - za(p1,p5)*iza(p1,p7)*iza(p2,p7)*zab2(p1,p2,p7,p4)*zab2(
     &    p3,p1,p5,p6) )
      zaa(1,1,1,2) = zaa(1,1,1,2) + tinv(p1,p5,p6)*tinv(p2,p3,p4)*
     & s34**(-1)*s56**(-1) * ( za(p1,p5)**2*zb(p2,p4)*zb(p5,p6)*iza(p1,
     &    p7)*zab2(p3,p2,p4,p7) )


c--- do not calculate single resonant diagrams unnecessarily
      if (srdiags .eqv. .false.) return

      z34(2,2,2,2)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * (  - 
     &    za(p1,p2)*za(p2,p6)*za(p3,p4)*zb(p1,p3)*zb(p3,p5)*iza(p1,p7)*
     &    iza(p2,p7) + za(p1,p2)*za(p2,p6)*za(p4,p5)*zb(p1,p5)*zb(p3,p5
     &    )*iza(p1,p7)*iza(p2,p7) - za(p2,p6)*za(p3,p4)*zb(p3,p5)*zb(p3
     &    ,p7)*iza(p1,p7) + za(p2,p6)*za(p4,p5)*zb(p3,p5)*zb(p5,p7)*
     &    iza(p1,p7) )
      z34(2,2,2,2) = z34(2,2,2,2) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * (  - za(p1,p2)*za(p2,p4)*za(p4,p6)*zb(p1,p5)*zb(p3,
     &    p4)*iza(p1,p7)*iza(p2,p7) - za(p1,p2)*za(p2,p6)*za(p4,p6)*zb(
     &    p1,p5)*zb(p3,p6)*iza(p1,p7)*iza(p2,p7) - za(p2,p4)*za(p4,p6)*
     &    zb(p3,p4)*zb(p5,p7)*iza(p1,p7) - za(p2,p6)*za(p4,p6)*zb(p3,p6
     &    )*zb(p5,p7)*iza(p1,p7) )

      z34(1,1,1,1)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * ( za(p1
     &    ,p3)*za(p3,p5)*zb(p1,p2)*zb(p2,p6)*zb(p3,p4)*izb(p1,p7)*izb(
     &    p2,p7) - za(p1,p5)*za(p3,p5)*zb(p1,p2)*zb(p2,p6)*zb(p4,p5)*
     &    izb(p1,p7)*izb(p2,p7) + za(p3,p5)*za(p3,p7)*zb(p2,p6)*zb(p3,
     &    p4)*izb(p1,p7) - za(p3,p5)*za(p5,p7)*zb(p2,p6)*zb(p4,p5)*izb(
     &    p1,p7) )
      z34(1,1,1,1) = z34(1,1,1,1) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * ( za(p1,p5)*za(p3,p4)*zb(p1,p2)*zb(p2,p4)*zb(p4,p6)*
     &    izb(p1,p7)*izb(p2,p7) + za(p1,p5)*za(p3,p6)*zb(p1,p2)*zb(p2,
     &    p6)*zb(p4,p6)*izb(p1,p7)*izb(p2,p7) + za(p3,p4)*za(p5,p7)*zb(
     &    p2,p4)*zb(p4,p6)*izb(p1,p7) + za(p3,p6)*za(p5,p7)*zb(p2,p6)*
     &    zb(p4,p6)*izb(p1,p7) )

      z34(1,2,2,2)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * (  - 
     &    za(p1,p2)*za(p1,p6)*za(p3,p4)*zb(p2,p3)*zb(p3,p5)*iza(p1,p7)*
     &    iza(p2,p7) + za(p1,p2)*za(p1,p6)*za(p4,p5)*zb(p2,p5)*zb(p3,p5
     &    )*iza(p1,p7)*iza(p2,p7) + za(p1,p6)*za(p3,p4)*zb(p3,p5)*zb(p3
     &    ,p7)*iza(p2,p7) - za(p1,p6)*za(p4,p5)*zb(p3,p5)*zb(p5,p7)*
     &    iza(p2,p7) )
      z34(1,2,2,2) = z34(1,2,2,2) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * (  - za(p1,p2)*za(p1,p4)*za(p4,p6)*zb(p2,p5)*zb(p3,
     &    p4)*iza(p1,p7)*iza(p2,p7) - za(p1,p2)*za(p1,p6)*za(p4,p6)*zb(
     &    p2,p5)*zb(p3,p6)*iza(p1,p7)*iza(p2,p7) + za(p1,p4)*za(p4,p6)*
     &    zb(p3,p4)*zb(p5,p7)*iza(p2,p7) + za(p1,p6)*za(p4,p6)*zb(p3,p6
     &    )*zb(p5,p7)*iza(p2,p7) )

      z34(2,1,2,2)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * ( za(p1
     &    ,p2)*za(p2,p6)*za(p3,p4)*zb(p1,p4)*zb(p4,p5)*iza(p1,p7)*iza(
     &    p2,p7) + za(p1,p2)*za(p2,p6)*za(p3,p5)*zb(p1,p5)*zb(p4,p5)*
     &    iza(p1,p7)*iza(p2,p7) + za(p2,p6)*za(p3,p4)*zb(p4,p5)*zb(p4,
     &    p7)*iza(p1,p7) + za(p2,p6)*za(p3,p5)*zb(p4,p5)*zb(p5,p7)*iza(
     &    p1,p7) )
      z34(2,1,2,2) = z34(2,1,2,2) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * ( za(p1,p2)*za(p2,p3)*za(p3,p6)*zb(p1,p5)*zb(p3,p4)*
     &    iza(p1,p7)*iza(p2,p7) - za(p1,p2)*za(p2,p6)*za(p3,p6)*zb(p1,
     &    p5)*zb(p4,p6)*iza(p1,p7)*iza(p2,p7) + za(p2,p3)*za(p3,p6)*zb(
     &    p3,p4)*zb(p5,p7)*iza(p1,p7) - za(p2,p6)*za(p3,p6)*zb(p4,p6)*
     &    zb(p5,p7)*iza(p1,p7) )

      z34(2,2,1,2)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * ( za(p1
     &    ,p2)*za(p2,p4)*za(p4,p5)*zb(p1,p6)*zb(p3,p4)*iza(p1,p7)*iza(
     &    p2,p7) + za(p1,p2)*za(p2,p5)*za(p4,p5)*zb(p1,p6)*zb(p3,p5)*
     &    iza(p1,p7)*iza(p2,p7) + za(p2,p4)*za(p4,p5)*zb(p3,p4)*zb(p6,
     &    p7)*iza(p1,p7) + za(p2,p5)*za(p4,p5)*zb(p3,p5)*zb(p6,p7)*iza(
     &    p1,p7) )
      z34(2,2,1,2) = z34(2,2,1,2) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * ( za(p1,p2)*za(p2,p5)*za(p3,p4)*zb(p1,p3)*zb(p3,p6)*
     &    iza(p1,p7)*iza(p2,p7) - za(p1,p2)*za(p2,p5)*za(p4,p6)*zb(p1,
     &    p6)*zb(p3,p6)*iza(p1,p7)*iza(p2,p7) + za(p2,p5)*za(p3,p4)*zb(
     &    p3,p6)*zb(p3,p7)*iza(p1,p7) - za(p2,p5)*za(p4,p6)*zb(p3,p6)*
     &    zb(p6,p7)*iza(p1,p7) )

      z34(2,2,2,1)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * ( za(p2
     &    ,p6)*za(p3,p4)*zb(p1,p2)*zb(p1,p3)*zb(p3,p5)*izb(p1,p7)*izb(
     &    p2,p7) - za(p2,p6)*za(p4,p5)*zb(p1,p2)*zb(p1,p5)*zb(p3,p5)*
     &    izb(p1,p7)*izb(p2,p7) - za(p3,p4)*za(p6,p7)*zb(p1,p3)*zb(p3,
     &    p5)*izb(p2,p7) + za(p4,p5)*za(p6,p7)*zb(p1,p5)*zb(p3,p5)*izb(
     &    p2,p7) )
      z34(2,2,2,1) = z34(2,2,2,1) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * ( za(p2,p4)*za(p4,p6)*zb(p1,p2)*zb(p1,p5)*zb(p3,p4)*
     &    izb(p1,p7)*izb(p2,p7) + za(p2,p6)*za(p4,p6)*zb(p1,p2)*zb(p1,
     &    p5)*zb(p3,p6)*izb(p1,p7)*izb(p2,p7) - za(p4,p6)*za(p4,p7)*zb(
     &    p1,p5)*zb(p3,p4)*izb(p2,p7) - za(p4,p6)*za(p6,p7)*zb(p1,p5)*
     &    zb(p3,p6)*izb(p2,p7) )

      z34(1,1,2,2)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * ( za(p1
     &    ,p2)*za(p1,p6)*za(p3,p4)*zb(p2,p4)*zb(p4,p5)*iza(p1,p7)*iza(
     &    p2,p7) + za(p1,p2)*za(p1,p6)*za(p3,p5)*zb(p2,p5)*zb(p4,p5)*
     &    iza(p1,p7)*iza(p2,p7) - za(p1,p6)*za(p3,p4)*zb(p4,p5)*zb(p4,
     &    p7)*iza(p2,p7) - za(p1,p6)*za(p3,p5)*zb(p4,p5)*zb(p5,p7)*iza(
     &    p2,p7) )
      z34(1,1,2,2) = z34(1,1,2,2) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * ( za(p1,p2)*za(p1,p3)*za(p3,p6)*zb(p2,p5)*zb(p3,p4)*
     &    iza(p1,p7)*iza(p2,p7) - za(p1,p2)*za(p1,p6)*za(p3,p6)*zb(p2,
     &    p5)*zb(p4,p6)*iza(p1,p7)*iza(p2,p7) - za(p1,p3)*za(p3,p6)*zb(
     &    p3,p4)*zb(p5,p7)*iza(p2,p7) + za(p1,p6)*za(p3,p6)*zb(p4,p6)*
     &    zb(p5,p7)*iza(p2,p7) )

      z34(1,2,1,2)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * ( za(p1
     &    ,p2)*za(p1,p4)*za(p4,p5)*zb(p2,p6)*zb(p3,p4)*iza(p1,p7)*iza(
     &    p2,p7) + za(p1,p2)*za(p1,p5)*za(p4,p5)*zb(p2,p6)*zb(p3,p5)*
     &    iza(p1,p7)*iza(p2,p7) - za(p1,p4)*za(p4,p5)*zb(p3,p4)*zb(p6,
     &    p7)*iza(p2,p7) - za(p1,p5)*za(p4,p5)*zb(p3,p5)*zb(p6,p7)*iza(
     &    p2,p7) )
      z34(1,2,1,2) = z34(1,2,1,2) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * ( za(p1,p2)*za(p1,p5)*za(p3,p4)*zb(p2,p3)*zb(p3,p6)*
     &    iza(p1,p7)*iza(p2,p7) - za(p1,p2)*za(p1,p5)*za(p4,p6)*zb(p2,
     &    p6)*zb(p3,p6)*iza(p1,p7)*iza(p2,p7) - za(p1,p5)*za(p3,p4)*zb(
     &    p3,p6)*zb(p3,p7)*iza(p2,p7) + za(p1,p5)*za(p4,p6)*zb(p3,p6)*
     &    zb(p6,p7)*iza(p2,p7) )

      z34(1,2,2,1)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * ( za(p1
     &    ,p6)*za(p3,p4)*zb(p1,p2)*zb(p2,p3)*zb(p3,p5)*izb(p1,p7)*izb(
     &    p2,p7) - za(p1,p6)*za(p4,p5)*zb(p1,p2)*zb(p2,p5)*zb(p3,p5)*
     &    izb(p1,p7)*izb(p2,p7) + za(p3,p4)*za(p6,p7)*zb(p2,p3)*zb(p3,
     &    p5)*izb(p1,p7) - za(p4,p5)*za(p6,p7)*zb(p2,p5)*zb(p3,p5)*izb(
     &    p1,p7) )
      z34(1,2,2,1) = z34(1,2,2,1) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * ( za(p1,p4)*za(p4,p6)*zb(p1,p2)*zb(p2,p5)*zb(p3,p4)*
     &    izb(p1,p7)*izb(p2,p7) + za(p1,p6)*za(p4,p6)*zb(p1,p2)*zb(p2,
     &    p5)*zb(p3,p6)*izb(p1,p7)*izb(p2,p7) + za(p4,p6)*za(p4,p7)*zb(
     &    p2,p5)*zb(p3,p4)*izb(p1,p7) + za(p4,p6)*za(p6,p7)*zb(p2,p5)*
     &    zb(p3,p6)*izb(p1,p7) )

      z34(2,1,1,2)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * (  - 
     &    za(p1,p2)*za(p2,p3)*za(p3,p5)*zb(p1,p6)*zb(p3,p4)*iza(p1,p7)*
     &    iza(p2,p7) + za(p1,p2)*za(p2,p5)*za(p3,p5)*zb(p1,p6)*zb(p4,p5
     &    )*iza(p1,p7)*iza(p2,p7) - za(p2,p3)*za(p3,p5)*zb(p3,p4)*zb(p6
     &    ,p7)*iza(p1,p7) + za(p2,p5)*za(p3,p5)*zb(p4,p5)*zb(p6,p7)*
     &    iza(p1,p7) )
      z34(2,1,1,2) = z34(2,1,1,2) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * (  - za(p1,p2)*za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(p4,
     &    p6)*iza(p1,p7)*iza(p2,p7) - za(p1,p2)*za(p2,p5)*za(p3,p6)*zb(
     &    p1,p6)*zb(p4,p6)*iza(p1,p7)*iza(p2,p7) - za(p2,p5)*za(p3,p4)*
     &    zb(p4,p6)*zb(p4,p7)*iza(p1,p7) - za(p2,p5)*za(p3,p6)*zb(p4,p6
     &    )*zb(p6,p7)*iza(p1,p7) )

      z34(2,1,2,1)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * (  - 
     &    za(p2,p6)*za(p3,p4)*zb(p1,p2)*zb(p1,p4)*zb(p4,p5)*izb(p1,p7)*
     &    izb(p2,p7) - za(p2,p6)*za(p3,p5)*zb(p1,p2)*zb(p1,p5)*zb(p4,p5
     &    )*izb(p1,p7)*izb(p2,p7) + za(p3,p4)*za(p6,p7)*zb(p1,p4)*zb(p4
     &    ,p5)*izb(p2,p7) + za(p3,p5)*za(p6,p7)*zb(p1,p5)*zb(p4,p5)*
     &    izb(p2,p7) )
      z34(2,1,2,1) = z34(2,1,2,1) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * (  - za(p2,p3)*za(p3,p6)*zb(p1,p2)*zb(p1,p5)*zb(p3,
     &    p4)*izb(p1,p7)*izb(p2,p7) + za(p2,p6)*za(p3,p6)*zb(p1,p2)*zb(
     &    p1,p5)*zb(p4,p6)*izb(p1,p7)*izb(p2,p7) + za(p3,p6)*za(p3,p7)*
     &    zb(p1,p5)*zb(p3,p4)*izb(p2,p7) - za(p3,p6)*za(p6,p7)*zb(p1,p5
     &    )*zb(p4,p6)*izb(p2,p7) )

      z34(2,2,1,1)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * (  - 
     &    za(p2,p4)*za(p4,p5)*zb(p1,p2)*zb(p1,p6)*zb(p3,p4)*izb(p1,p7)*
     &    izb(p2,p7) - za(p2,p5)*za(p4,p5)*zb(p1,p2)*zb(p1,p6)*zb(p3,p5
     &    )*izb(p1,p7)*izb(p2,p7) + za(p4,p5)*za(p4,p7)*zb(p1,p6)*zb(p3
     &    ,p4)*izb(p2,p7) + za(p4,p5)*za(p5,p7)*zb(p1,p6)*zb(p3,p5)*
     &    izb(p2,p7) )
      z34(2,2,1,1) = z34(2,2,1,1) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * (  - za(p2,p5)*za(p3,p4)*zb(p1,p2)*zb(p1,p3)*zb(p3,
     &    p6)*izb(p1,p7)*izb(p2,p7) + za(p2,p5)*za(p4,p6)*zb(p1,p2)*zb(
     &    p1,p6)*zb(p3,p6)*izb(p1,p7)*izb(p2,p7) + za(p3,p4)*za(p5,p7)*
     &    zb(p1,p3)*zb(p3,p6)*izb(p2,p7) - za(p4,p6)*za(p5,p7)*zb(p1,p6
     &    )*zb(p3,p6)*izb(p2,p7) )

      z34(2,1,1,1)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * ( za(p2
     &    ,p3)*za(p3,p5)*zb(p1,p2)*zb(p1,p6)*zb(p3,p4)*izb(p1,p7)*izb(
     &    p2,p7) - za(p2,p5)*za(p3,p5)*zb(p1,p2)*zb(p1,p6)*zb(p4,p5)*
     &    izb(p1,p7)*izb(p2,p7) - za(p3,p5)*za(p3,p7)*zb(p1,p6)*zb(p3,
     &    p4)*izb(p2,p7) + za(p3,p5)*za(p5,p7)*zb(p1,p6)*zb(p4,p5)*izb(
     &    p2,p7) )
      z34(2,1,1,1) = z34(2,1,1,1) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * ( za(p2,p5)*za(p3,p4)*zb(p1,p2)*zb(p1,p4)*zb(p4,p6)*
     &    izb(p1,p7)*izb(p2,p7) + za(p2,p5)*za(p3,p6)*zb(p1,p2)*zb(p1,
     &    p6)*zb(p4,p6)*izb(p1,p7)*izb(p2,p7) - za(p3,p4)*za(p5,p7)*zb(
     &    p1,p4)*zb(p4,p6)*izb(p2,p7) - za(p3,p6)*za(p5,p7)*zb(p1,p6)*
     &    zb(p4,p6)*izb(p2,p7) )

      z34(1,2,1,1)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p4,p5)*zb(p1,p2)*zb(p2,p6)*zb(p3,p4)*izb(p1,p7)*
     &    izb(p2,p7) - za(p1,p5)*za(p4,p5)*zb(p1,p2)*zb(p2,p6)*zb(p3,p5
     &    )*izb(p1,p7)*izb(p2,p7) - za(p4,p5)*za(p4,p7)*zb(p2,p6)*zb(p3
     &    ,p4)*izb(p1,p7) - za(p4,p5)*za(p5,p7)*zb(p2,p6)*zb(p3,p5)*
     &    izb(p1,p7) )
      z34(1,2,1,1) = z34(1,2,1,1) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * (  - za(p1,p5)*za(p3,p4)*zb(p1,p2)*zb(p2,p3)*zb(p3,
     &    p6)*izb(p1,p7)*izb(p2,p7) + za(p1,p5)*za(p4,p6)*zb(p1,p2)*zb(
     &    p2,p6)*zb(p3,p6)*izb(p1,p7)*izb(p2,p7) - za(p3,p4)*za(p5,p7)*
     &    zb(p2,p3)*zb(p3,p6)*izb(p1,p7) + za(p4,p6)*za(p5,p7)*zb(p2,p6
     &    )*zb(p3,p6)*izb(p1,p7) )

      z34(1,1,2,1)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * (  - 
     &    za(p1,p6)*za(p3,p4)*zb(p1,p2)*zb(p2,p4)*zb(p4,p5)*izb(p1,p7)*
     &    izb(p2,p7) - za(p1,p6)*za(p3,p5)*zb(p1,p2)*zb(p2,p5)*zb(p4,p5
     &    )*izb(p1,p7)*izb(p2,p7) - za(p3,p4)*za(p6,p7)*zb(p2,p4)*zb(p4
     &    ,p5)*izb(p1,p7) - za(p3,p5)*za(p6,p7)*zb(p2,p5)*zb(p4,p5)*
     &    izb(p1,p7) )
      z34(1,1,2,1) = z34(1,1,2,1) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * (  - za(p1,p3)*za(p3,p6)*zb(p1,p2)*zb(p2,p5)*zb(p3,
     &    p4)*izb(p1,p7)*izb(p2,p7) + za(p1,p6)*za(p3,p6)*zb(p1,p2)*zb(
     &    p2,p5)*zb(p4,p6)*izb(p1,p7)*izb(p2,p7) - za(p3,p6)*za(p3,p7)*
     &    zb(p2,p5)*zb(p3,p4)*izb(p1,p7) + za(p3,p6)*za(p6,p7)*zb(p2,p5
     &    )*zb(p4,p6)*izb(p1,p7) )

      z34(1,1,1,2)= + tinv(p1,p2,p7)*tinv(p3,p4,p5)*s34**(-1) * (  - 
     &    za(p1,p2)*za(p1,p3)*za(p3,p5)*zb(p2,p6)*zb(p3,p4)*iza(p1,p7)*
     &    iza(p2,p7) + za(p1,p2)*za(p1,p5)*za(p3,p5)*zb(p2,p6)*zb(p4,p5
     &    )*iza(p1,p7)*iza(p2,p7) + za(p1,p3)*za(p3,p5)*zb(p3,p4)*zb(p6
     &    ,p7)*iza(p2,p7) - za(p1,p5)*za(p3,p5)*zb(p4,p5)*zb(p6,p7)*
     &    iza(p2,p7) )
      z34(1,1,1,2) = z34(1,1,1,2) + tinv(p1,p2,p7)*tinv(p3,p4,p6)*
     & s34**(-1) * (  - za(p1,p2)*za(p1,p5)*za(p3,p4)*zb(p2,p4)*zb(p4,
     &    p6)*iza(p1,p7)*iza(p2,p7) - za(p1,p2)*za(p1,p5)*za(p3,p6)*zb(
     &    p2,p6)*zb(p4,p6)*iza(p1,p7)*iza(p2,p7) + za(p1,p5)*za(p3,p4)*
     &    zb(p4,p6)*zb(p4,p7)*iza(p2,p7) + za(p1,p5)*za(p3,p6)*zb(p4,p6
     &    )*zb(p6,p7)*iza(p2,p7) )

      z56(2,2,2,2)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * ( za(p1
     &    ,p2)*za(p2,p4)*za(p3,p6)*zb(p1,p3)*zb(p3,p5)*iza(p1,p7)*iza(
     &    p2,p7) + za(p1,p2)*za(p2,p4)*za(p5,p6)*zb(p1,p5)*zb(p3,p5)*
     &    iza(p1,p7)*iza(p2,p7) + za(p2,p4)*za(p3,p6)*zb(p3,p5)*zb(p3,
     &    p7)*iza(p1,p7) + za(p2,p4)*za(p5,p6)*zb(p3,p5)*zb(p5,p7)*iza(
     &    p1,p7) )
      z56(2,2,2,2) = z56(2,2,2,2) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * (  - za(p1,p2)*za(p2,p4)*za(p4,p6)*zb(p1,p3)*zb(p4,
     &    p5)*iza(p1,p7)*iza(p2,p7) + za(p1,p2)*za(p2,p6)*za(p4,p6)*zb(
     &    p1,p3)*zb(p5,p6)*iza(p1,p7)*iza(p2,p7) - za(p2,p4)*za(p4,p6)*
     &    zb(p3,p7)*zb(p4,p5)*iza(p1,p7) + za(p2,p6)*za(p4,p6)*zb(p3,p7
     &    )*zb(p5,p6)*iza(p1,p7) )

      z56(1,1,1,1)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * (  - 
     &    za(p1,p3)*za(p3,p5)*zb(p1,p2)*zb(p2,p4)*zb(p3,p6)*izb(p1,p7)*
     &    izb(p2,p7) - za(p1,p5)*za(p3,p5)*zb(p1,p2)*zb(p2,p4)*zb(p5,p6
     &    )*izb(p1,p7)*izb(p2,p7) - za(p3,p5)*za(p3,p7)*zb(p2,p4)*zb(p3
     &    ,p6)*izb(p1,p7) - za(p3,p5)*za(p5,p7)*zb(p2,p4)*zb(p5,p6)*
     &    izb(p1,p7) )
      z56(1,1,1,1) = z56(1,1,1,1) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * ( za(p1,p3)*za(p4,p5)*zb(p1,p2)*zb(p2,p4)*zb(p4,p6)*
     &    izb(p1,p7)*izb(p2,p7) - za(p1,p3)*za(p5,p6)*zb(p1,p2)*zb(p2,
     &    p6)*zb(p4,p6)*izb(p1,p7)*izb(p2,p7) + za(p3,p7)*za(p4,p5)*zb(
     &    p2,p4)*zb(p4,p6)*izb(p1,p7) - za(p3,p7)*za(p5,p6)*zb(p2,p6)*
     &    zb(p4,p6)*izb(p1,p7) )

      z56(1,2,2,2)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * ( za(p1
     &    ,p2)*za(p1,p4)*za(p3,p6)*zb(p2,p3)*zb(p3,p5)*iza(p1,p7)*iza(
     &    p2,p7) + za(p1,p2)*za(p1,p4)*za(p5,p6)*zb(p2,p5)*zb(p3,p5)*
     &    iza(p1,p7)*iza(p2,p7) - za(p1,p4)*za(p3,p6)*zb(p3,p5)*zb(p3,
     &    p7)*iza(p2,p7) - za(p1,p4)*za(p5,p6)*zb(p3,p5)*zb(p5,p7)*iza(
     &    p2,p7) )
      z56(1,2,2,2) = z56(1,2,2,2) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * (  - za(p1,p2)*za(p1,p4)*za(p4,p6)*zb(p2,p3)*zb(p4,
     &    p5)*iza(p1,p7)*iza(p2,p7) + za(p1,p2)*za(p1,p6)*za(p4,p6)*zb(
     &    p2,p3)*zb(p5,p6)*iza(p1,p7)*iza(p2,p7) + za(p1,p4)*za(p4,p6)*
     &    zb(p3,p7)*zb(p4,p5)*iza(p2,p7) - za(p1,p6)*za(p4,p6)*zb(p3,p7
     &    )*zb(p5,p6)*iza(p2,p7) )

      z56(2,1,2,2)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * ( za(p1
     &    ,p2)*za(p2,p3)*za(p3,p6)*zb(p1,p4)*zb(p3,p5)*iza(p1,p7)*iza(
     &    p2,p7) - za(p1,p2)*za(p2,p6)*za(p3,p6)*zb(p1,p4)*zb(p5,p6)*
     &    iza(p1,p7)*iza(p2,p7) + za(p2,p3)*za(p3,p6)*zb(p3,p5)*zb(p4,
     &    p7)*iza(p1,p7) - za(p2,p6)*za(p3,p6)*zb(p4,p7)*zb(p5,p6)*iza(
     &    p1,p7) )
      z56(2,1,2,2) = z56(2,1,2,2) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * (  - za(p1,p2)*za(p2,p3)*za(p4,p6)*zb(p1,p4)*zb(p4,
     &    p5)*iza(p1,p7)*iza(p2,p7) - za(p1,p2)*za(p2,p3)*za(p5,p6)*zb(
     &    p1,p5)*zb(p4,p5)*iza(p1,p7)*iza(p2,p7) - za(p2,p3)*za(p4,p6)*
     &    zb(p4,p5)*zb(p4,p7)*iza(p1,p7) - za(p2,p3)*za(p5,p6)*zb(p4,p5
     &    )*zb(p5,p7)*iza(p1,p7) )

      z56(2,2,1,2)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * ( za(p1
     &    ,p2)*za(p2,p4)*za(p3,p5)*zb(p1,p3)*zb(p3,p6)*iza(p1,p7)*iza(
     &    p2,p7) - za(p1,p2)*za(p2,p4)*za(p5,p6)*zb(p1,p6)*zb(p3,p6)*
     &    iza(p1,p7)*iza(p2,p7) + za(p2,p4)*za(p3,p5)*zb(p3,p6)*zb(p3,
     &    p7)*iza(p1,p7) - za(p2,p4)*za(p5,p6)*zb(p3,p6)*zb(p6,p7)*iza(
     &    p1,p7) )
      z56(2,2,1,2) = z56(2,2,1,2) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * (  - za(p1,p2)*za(p2,p4)*za(p4,p5)*zb(p1,p3)*zb(p4,
     &    p6)*iza(p1,p7)*iza(p2,p7) - za(p1,p2)*za(p2,p5)*za(p4,p5)*zb(
     &    p1,p3)*zb(p5,p6)*iza(p1,p7)*iza(p2,p7) - za(p2,p4)*za(p4,p5)*
     &    zb(p3,p7)*zb(p4,p6)*iza(p1,p7) - za(p2,p5)*za(p4,p5)*zb(p3,p7
     &    )*zb(p5,p6)*iza(p1,p7) )

      z56(2,2,2,1)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * (  - 
     &    za(p2,p4)*za(p3,p6)*zb(p1,p2)*zb(p1,p3)*zb(p3,p5)*izb(p1,p7)*
     &    izb(p2,p7) - za(p2,p4)*za(p5,p6)*zb(p1,p2)*zb(p1,p5)*zb(p3,p5
     &    )*izb(p1,p7)*izb(p2,p7) + za(p3,p6)*za(p4,p7)*zb(p1,p3)*zb(p3
     &    ,p5)*izb(p2,p7) + za(p4,p7)*za(p5,p6)*zb(p1,p5)*zb(p3,p5)*
     &    izb(p2,p7) )
      z56(2,2,2,1) = z56(2,2,2,1) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * ( za(p2,p4)*za(p4,p6)*zb(p1,p2)*zb(p1,p3)*zb(p4,p5)*
     &    izb(p1,p7)*izb(p2,p7) - za(p2,p6)*za(p4,p6)*zb(p1,p2)*zb(p1,
     &    p3)*zb(p5,p6)*izb(p1,p7)*izb(p2,p7) - za(p4,p6)*za(p4,p7)*zb(
     &    p1,p3)*zb(p4,p5)*izb(p2,p7) + za(p4,p6)*za(p6,p7)*zb(p1,p3)*
     &    zb(p5,p6)*izb(p2,p7) )

      z56(1,1,2,2)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * ( za(p1
     &    ,p2)*za(p1,p3)*za(p3,p6)*zb(p2,p4)*zb(p3,p5)*iza(p1,p7)*iza(
     &    p2,p7) - za(p1,p2)*za(p1,p6)*za(p3,p6)*zb(p2,p4)*zb(p5,p6)*
     &    iza(p1,p7)*iza(p2,p7) - za(p1,p3)*za(p3,p6)*zb(p3,p5)*zb(p4,
     &    p7)*iza(p2,p7) + za(p1,p6)*za(p3,p6)*zb(p4,p7)*zb(p5,p6)*iza(
     &    p2,p7) )
      z56(1,1,2,2) = z56(1,1,2,2) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * (  - za(p1,p2)*za(p1,p3)*za(p4,p6)*zb(p2,p4)*zb(p4,
     &    p5)*iza(p1,p7)*iza(p2,p7) - za(p1,p2)*za(p1,p3)*za(p5,p6)*zb(
     &    p2,p5)*zb(p4,p5)*iza(p1,p7)*iza(p2,p7) + za(p1,p3)*za(p4,p6)*
     &    zb(p4,p5)*zb(p4,p7)*iza(p2,p7) + za(p1,p3)*za(p5,p6)*zb(p4,p5
     &    )*zb(p5,p7)*iza(p2,p7) )

      z56(1,2,1,2)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * ( za(p1
     &    ,p2)*za(p1,p4)*za(p3,p5)*zb(p2,p3)*zb(p3,p6)*iza(p1,p7)*iza(
     &    p2,p7) - za(p1,p2)*za(p1,p4)*za(p5,p6)*zb(p2,p6)*zb(p3,p6)*
     &    iza(p1,p7)*iza(p2,p7) - za(p1,p4)*za(p3,p5)*zb(p3,p6)*zb(p3,
     &    p7)*iza(p2,p7) + za(p1,p4)*za(p5,p6)*zb(p3,p6)*zb(p6,p7)*iza(
     &    p2,p7) )
      z56(1,2,1,2) = z56(1,2,1,2) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * (  - za(p1,p2)*za(p1,p4)*za(p4,p5)*zb(p2,p3)*zb(p4,
     &    p6)*iza(p1,p7)*iza(p2,p7) - za(p1,p2)*za(p1,p5)*za(p4,p5)*zb(
     &    p2,p3)*zb(p5,p6)*iza(p1,p7)*iza(p2,p7) + za(p1,p4)*za(p4,p5)*
     &    zb(p3,p7)*zb(p4,p6)*iza(p2,p7) + za(p1,p5)*za(p4,p5)*zb(p3,p7
     &    )*zb(p5,p6)*iza(p2,p7) )

      z56(1,2,2,1)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * (  - 
     &    za(p1,p4)*za(p3,p6)*zb(p1,p2)*zb(p2,p3)*zb(p3,p5)*izb(p1,p7)*
     &    izb(p2,p7) - za(p1,p4)*za(p5,p6)*zb(p1,p2)*zb(p2,p5)*zb(p3,p5
     &    )*izb(p1,p7)*izb(p2,p7) - za(p3,p6)*za(p4,p7)*zb(p2,p3)*zb(p3
     &    ,p5)*izb(p1,p7) - za(p4,p7)*za(p5,p6)*zb(p2,p5)*zb(p3,p5)*
     &    izb(p1,p7) )
      z56(1,2,2,1) = z56(1,2,2,1) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * ( za(p1,p4)*za(p4,p6)*zb(p1,p2)*zb(p2,p3)*zb(p4,p5)*
     &    izb(p1,p7)*izb(p2,p7) - za(p1,p6)*za(p4,p6)*zb(p1,p2)*zb(p2,
     &    p3)*zb(p5,p6)*izb(p1,p7)*izb(p2,p7) + za(p4,p6)*za(p4,p7)*zb(
     &    p2,p3)*zb(p4,p5)*izb(p1,p7) - za(p4,p6)*za(p6,p7)*zb(p2,p3)*
     &    zb(p5,p6)*izb(p1,p7) )

      z56(2,1,1,2)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * ( za(p1
     &    ,p2)*za(p2,p3)*za(p3,p5)*zb(p1,p4)*zb(p3,p6)*iza(p1,p7)*iza(
     &    p2,p7) + za(p1,p2)*za(p2,p5)*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*
     &    iza(p1,p7)*iza(p2,p7) + za(p2,p3)*za(p3,p5)*zb(p3,p6)*zb(p4,
     &    p7)*iza(p1,p7) + za(p2,p5)*za(p3,p5)*zb(p4,p7)*zb(p5,p6)*iza(
     &    p1,p7) )
      z56(2,1,1,2) = z56(2,1,1,2) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * (  - za(p1,p2)*za(p2,p3)*za(p4,p5)*zb(p1,p4)*zb(p4,
     &    p6)*iza(p1,p7)*iza(p2,p7) + za(p1,p2)*za(p2,p3)*za(p5,p6)*zb(
     &    p1,p6)*zb(p4,p6)*iza(p1,p7)*iza(p2,p7) - za(p2,p3)*za(p4,p5)*
     &    zb(p4,p6)*zb(p4,p7)*iza(p1,p7) + za(p2,p3)*za(p5,p6)*zb(p4,p6
     &    )*zb(p6,p7)*iza(p1,p7) )

      z56(2,1,2,1)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * (  - 
     &    za(p2,p3)*za(p3,p6)*zb(p1,p2)*zb(p1,p4)*zb(p3,p5)*izb(p1,p7)*
     &    izb(p2,p7) + za(p2,p6)*za(p3,p6)*zb(p1,p2)*zb(p1,p4)*zb(p5,p6
     &    )*izb(p1,p7)*izb(p2,p7) + za(p3,p6)*za(p3,p7)*zb(p1,p4)*zb(p3
     &    ,p5)*izb(p2,p7) - za(p3,p6)*za(p6,p7)*zb(p1,p4)*zb(p5,p6)*
     &    izb(p2,p7) )
      z56(2,1,2,1) = z56(2,1,2,1) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * ( za(p2,p3)*za(p4,p6)*zb(p1,p2)*zb(p1,p4)*zb(p4,p5)*
     &    izb(p1,p7)*izb(p2,p7) + za(p2,p3)*za(p5,p6)*zb(p1,p2)*zb(p1,
     &    p5)*zb(p4,p5)*izb(p1,p7)*izb(p2,p7) - za(p3,p7)*za(p4,p6)*zb(
     &    p1,p4)*zb(p4,p5)*izb(p2,p7) - za(p3,p7)*za(p5,p6)*zb(p1,p5)*
     &    zb(p4,p5)*izb(p2,p7) )

      z56(2,2,1,1)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * (  - 
     &    za(p2,p4)*za(p3,p5)*zb(p1,p2)*zb(p1,p3)*zb(p3,p6)*izb(p1,p7)*
     &    izb(p2,p7) + za(p2,p4)*za(p5,p6)*zb(p1,p2)*zb(p1,p6)*zb(p3,p6
     &    )*izb(p1,p7)*izb(p2,p7) + za(p3,p5)*za(p4,p7)*zb(p1,p3)*zb(p3
     &    ,p6)*izb(p2,p7) - za(p4,p7)*za(p5,p6)*zb(p1,p6)*zb(p3,p6)*
     &    izb(p2,p7) )
      z56(2,2,1,1) = z56(2,2,1,1) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * ( za(p2,p4)*za(p4,p5)*zb(p1,p2)*zb(p1,p3)*zb(p4,p6)*
     &    izb(p1,p7)*izb(p2,p7) + za(p2,p5)*za(p4,p5)*zb(p1,p2)*zb(p1,
     &    p3)*zb(p5,p6)*izb(p1,p7)*izb(p2,p7) - za(p4,p5)*za(p4,p7)*zb(
     &    p1,p3)*zb(p4,p6)*izb(p2,p7) - za(p4,p5)*za(p5,p7)*zb(p1,p3)*
     &    zb(p5,p6)*izb(p2,p7) )

      z56(2,1,1,1)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * (  - 
     &    za(p2,p3)*za(p3,p5)*zb(p1,p2)*zb(p1,p4)*zb(p3,p6)*izb(p1,p7)*
     &    izb(p2,p7) - za(p2,p5)*za(p3,p5)*zb(p1,p2)*zb(p1,p4)*zb(p5,p6
     &    )*izb(p1,p7)*izb(p2,p7) + za(p3,p5)*za(p3,p7)*zb(p1,p4)*zb(p3
     &    ,p6)*izb(p2,p7) + za(p3,p5)*za(p5,p7)*zb(p1,p4)*zb(p5,p6)*
     &    izb(p2,p7) )
      z56(2,1,1,1) = z56(2,1,1,1) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * ( za(p2,p3)*za(p4,p5)*zb(p1,p2)*zb(p1,p4)*zb(p4,p6)*
     &    izb(p1,p7)*izb(p2,p7) - za(p2,p3)*za(p5,p6)*zb(p1,p2)*zb(p1,
     &    p6)*zb(p4,p6)*izb(p1,p7)*izb(p2,p7) - za(p3,p7)*za(p4,p5)*zb(
     &    p1,p4)*zb(p4,p6)*izb(p2,p7) + za(p3,p7)*za(p5,p6)*zb(p1,p6)*
     &    zb(p4,p6)*izb(p2,p7) )

      z56(1,2,1,1)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * (  - 
     &    za(p1,p4)*za(p3,p5)*zb(p1,p2)*zb(p2,p3)*zb(p3,p6)*izb(p1,p7)*
     &    izb(p2,p7) + za(p1,p4)*za(p5,p6)*zb(p1,p2)*zb(p2,p6)*zb(p3,p6
     &    )*izb(p1,p7)*izb(p2,p7) - za(p3,p5)*za(p4,p7)*zb(p2,p3)*zb(p3
     &    ,p6)*izb(p1,p7) + za(p4,p7)*za(p5,p6)*zb(p2,p6)*zb(p3,p6)*
     &    izb(p1,p7) )
      z56(1,2,1,1) = z56(1,2,1,1) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * ( za(p1,p4)*za(p4,p5)*zb(p1,p2)*zb(p2,p3)*zb(p4,p6)*
     &    izb(p1,p7)*izb(p2,p7) + za(p1,p5)*za(p4,p5)*zb(p1,p2)*zb(p2,
     &    p3)*zb(p5,p6)*izb(p1,p7)*izb(p2,p7) + za(p4,p5)*za(p4,p7)*zb(
     &    p2,p3)*zb(p4,p6)*izb(p1,p7) + za(p4,p5)*za(p5,p7)*zb(p2,p3)*
     &    zb(p5,p6)*izb(p1,p7) )

      z56(1,1,2,1)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * (  - 
     &    za(p1,p3)*za(p3,p6)*zb(p1,p2)*zb(p2,p4)*zb(p3,p5)*izb(p1,p7)*
     &    izb(p2,p7) + za(p1,p6)*za(p3,p6)*zb(p1,p2)*zb(p2,p4)*zb(p5,p6
     &    )*izb(p1,p7)*izb(p2,p7) - za(p3,p6)*za(p3,p7)*zb(p2,p4)*zb(p3
     &    ,p5)*izb(p1,p7) + za(p3,p6)*za(p6,p7)*zb(p2,p4)*zb(p5,p6)*
     &    izb(p1,p7) )
      z56(1,1,2,1) = z56(1,1,2,1) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * ( za(p1,p3)*za(p4,p6)*zb(p1,p2)*zb(p2,p4)*zb(p4,p5)*
     &    izb(p1,p7)*izb(p2,p7) + za(p1,p3)*za(p5,p6)*zb(p1,p2)*zb(p2,
     &    p5)*zb(p4,p5)*izb(p1,p7)*izb(p2,p7) + za(p3,p7)*za(p4,p6)*zb(
     &    p2,p4)*zb(p4,p5)*izb(p1,p7) + za(p3,p7)*za(p5,p6)*zb(p2,p5)*
     &    zb(p4,p5)*izb(p1,p7) )

      z56(1,1,1,2)= + tinv(p1,p2,p7)*tinv(p3,p5,p6)*s56**(-1) * ( za(p1
     &    ,p2)*za(p1,p3)*za(p3,p5)*zb(p2,p4)*zb(p3,p6)*iza(p1,p7)*iza(
     &    p2,p7) + za(p1,p2)*za(p1,p5)*za(p3,p5)*zb(p2,p4)*zb(p5,p6)*
     &    iza(p1,p7)*iza(p2,p7) - za(p1,p3)*za(p3,p5)*zb(p3,p6)*zb(p4,
     &    p7)*iza(p2,p7) - za(p1,p5)*za(p3,p5)*zb(p4,p7)*zb(p5,p6)*iza(
     &    p2,p7) )
      z56(1,1,1,2) = z56(1,1,1,2) + tinv(p1,p2,p7)*tinv(p4,p5,p6)*
     & s56**(-1) * (  - za(p1,p2)*za(p1,p3)*za(p4,p5)*zb(p2,p4)*zb(p4,
     &    p6)*iza(p1,p7)*iza(p2,p7) + za(p1,p2)*za(p1,p3)*za(p5,p6)*zb(
     &    p2,p6)*zb(p4,p6)*iza(p1,p7)*iza(p2,p7) + za(p1,p3)*za(p4,p5)*
     &    zb(p4,p6)*zb(p4,p7)*iza(p2,p7) - za(p1,p3)*za(p5,p6)*zb(p4,p6
     &    )*zb(p6,p7)*iza(p2,p7) )

      return
      end
