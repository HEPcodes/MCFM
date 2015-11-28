      subroutine wbbgmassf(ju,jd,jg,jb,jc,jn,je,m,qedf)
      implicit none
      include 'constants.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      double precision m,sned,sneu,sneud
      double complex qedf(2,2,2)
      double complex iza(7,7),izb(7,7)
      integer i,j,ju,jd,jg,jb,jc,jn,je
C--   order of indices is gluon helicity hg,sb,sc
      sned=s(jn,je)+s(jn,jd)+s(je,jd)
      sneu=s(jn,je)+s(jn,ju)+s(je,ju)
      sneud=sneu+s(jn,jd)+s(je,jd)+s(ju,jd)
      do i=1,7
      do j=i+1,7
         iza(i,j)=cone/za(i,j)
         izb(i,j)=cone/zb(i,j)
         iza(j,i)=-iza(i,j)
         izb(j,i)=-izb(i,j)
      enddo
      enddo
      qedf(2,2,2)= + iza(jb,jg)*sned**(-1)*sneud**(-1)*fourrt2 * (  - 
     &    za(jd,jc)*za(jd,jn)*zb(ju,jg)*zb(jd,je) - za(jd,jn)*za(jn,jc)
     &    *zb(ju,jg)*zb(jn,je) )
      qedf(2,2,2) = qedf(2,2,2) + iza(jb,jg)*sneu**(-1)*sneud**(-1)*
     & fourrt2 * ( za(ju,jn)*za(jd,jc)*zb(ju,jg)*zb(ju,je) - za(jd,jc)*
     &    za(jn,je)*zb(ju,je)*zb(je,jg) )
      qedf(2,2,2) = qedf(2,2,2) + iza(jb,jg)**2*iza(jc,jg)*izb(jc,jg)*
     & sned**(-1)*sneud**(-1)*fourrt2*m**2 * ( za(jd,jg)*za(jd,jn)*za(
     &    jb,jc)*zb(ju,jg)*zb(jd,je) + za(jd,jn)*za(jb,jc)*za(jn,jg)*
     &    zb(ju,jg)*zb(jn,je) )
      qedf(2,2,2) = qedf(2,2,2) + iza(jb,jg)**2*iza(jc,jg)*izb(jc,jg)*
     & sneu**(-1)*sneud**(-1)*fourrt2*m**2 * (  - za(ju,jn)*za(jd,jg)*
     &    za(jb,jc)*zb(ju,jg)*zb(ju,je) + za(jd,jg)*za(jb,jc)*za(jn,je)
     &    *zb(ju,je)*zb(je,jg) )
      qedf(2,2,2) = qedf(2,2,2) + iza(jb,jg)*iza(jc,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2 * ( za(jd,jc)*za(jd,jn)*za(jb,jc)*zb(ju,jb)*
     &    zb(jd,je) + za(jd,jn)*za(jb,jc)*za(jn,jc)*zb(ju,jb)*zb(jn,je)
     &     )
      qedf(2,2,2) = qedf(2,2,2) + iza(jb,jg)*iza(jc,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2 * (  - za(ju,jn)*za(jd,jc)*za(jb,jc)*zb(ju,
     &    jb)*zb(ju,je) + za(jd,jc)*za(jb,jc)*za(jn,je)*zb(ju,je)*zb(je
     &    ,jb) )

      qedf(1,2,2)= + iza(jb,jg)*izb(jb,jg)*izb(jc,jg)**2*sned**(-1)*
     & sneud**(-1)*fourrt2*m**2 * (  - za(jd,jg)*za(jd,jn)*zb(ju,jg)*
     &    zb(jd,je)*zb(jb,jc) - za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(jb,jc)
     &    *zb(jn,je) )
      qedf(1,2,2) = qedf(1,2,2) + iza(jb,jg)*izb(jb,jg)*izb(jc,jg)**2*
     & sneu**(-1)*sneud**(-1)*fourrt2*m**2 * ( za(ju,jn)*za(jd,jg)*zb(
     &    ju,jg)*zb(ju,je)*zb(jb,jc) - za(jd,jg)*za(jn,je)*zb(ju,je)*
     &    zb(jb,jc)*zb(je,jg) )
      qedf(1,2,2) = qedf(1,2,2) + izb(jb,jg)*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2 * (  - za(jd,jc)*za(jd,jn)*zb(ju,jb)*zb(jd,
     &    je)*zb(jb,jc) - za(jd,jn)*za(jn,jc)*zb(ju,jb)*zb(jb,jc)*zb(jn
     &    ,je) )
      qedf(1,2,2) = qedf(1,2,2) + izb(jb,jg)*izb(jc,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2 * ( za(ju,jn)*za(jd,jc)*zb(ju,jb)*zb(ju,je)*
     &    zb(jb,jc) - za(jd,jc)*za(jn,je)*zb(ju,je)*zb(jb,jc)*zb(je,jb)
     &     )
      qedf(1,2,2) = qedf(1,2,2) + izb(jc,jg)*sned**(-1)*sneud**(-1)*
     & fourrt2 * (  - za(jd,jg)*za(jd,jn)*zb(ju,jb)*zb(jd,je) - za(jd,
     &    jn)*za(jn,jg)*zb(ju,jb)*zb(jn,je) )
      qedf(1,2,2) = qedf(1,2,2) + izb(jc,jg)*sneu**(-1)*sneud**(-1)*
     & fourrt2 * ( za(ju,jn)*za(jd,jg)*zb(ju,jb)*zb(ju,je) - za(jd,jg)*
     &    za(jn,je)*zb(ju,je)*zb(je,jb) )

      qedf(2,1,2)= + iza(jb,jg)*iza(jc,jg)*izb(jb,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2*m * (  - za(jd,jc)*za(jd,jn)*za(jb,jc)*zb(ju
     &    ,jg)*zb(jd,je) - za(jd,jn)*za(jb,jc)*za(jn,jc)*zb(ju,jg)*zb(
     &    jn,je) )
      qedf(2,1,2) = qedf(2,1,2) + iza(jb,jg)*iza(jc,jg)*izb(jb,jg)*
     & sneu**(-1)*sneud**(-1)*fourrt2*m * ( za(ju,jn)*za(jd,jc)*za(jb,
     &    jc)*zb(ju,jg)*zb(ju,je) - za(jd,jc)*za(jb,jc)*za(jn,je)*zb(ju
     &    ,je)*zb(je,jg) )
      qedf(2,1,2) = qedf(2,1,2) + iza(jb,jg)*iza(jc,jg)*izb(jc,jg)*
     & sned**(-1)*sneud**(-1)*fourrt2*m * (  - za(jd,jb)*za(jd,jn)*za(
     &    jb,jc)*zb(ju,jg)*zb(jd,je) - za(jd,jn)*za(jb,jc)*za(jn,jb)*
     &    zb(ju,jg)*zb(jn,je) )
      qedf(2,1,2) = qedf(2,1,2) + iza(jb,jg)*iza(jc,jg)*izb(jc,jg)*
     & sneu**(-1)*sneud**(-1)*fourrt2*m * ( za(ju,jn)*za(jd,jb)*za(jb,
     &    jc)*zb(ju,jg)*zb(ju,je) - za(jd,jb)*za(jb,jc)*za(jn,je)*zb(ju
     &    ,je)*zb(je,jg) )

      qedf(2,2,1)= + iza(jb,jg)**2*iza(jc,jg)*sned**(-1)*sneud**(-1)*
     & fourrt2*m * (  - za(jd,jg)*za(jd,jn)*za(jb,jc)*zb(ju,jc)*zb(jd,
     &    je) - za(jd,jn)*za(jb,jc)*za(jn,jg)*zb(ju,jc)*zb(jn,je) )
      qedf(2,2,1) = qedf(2,2,1) + iza(jb,jg)**2*iza(jc,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2*m * ( za(ju,jn)*za(jd,jg)*za(jb,jc)*zb(ju,jc
     &    )*zb(ju,je) - za(jd,jg)*za(jb,jc)*za(jn,je)*zb(ju,je)*zb(je,
     &    jc) )
      qedf(2,2,1) = qedf(2,2,1) + iza(jb,jg)*iza(jc,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2*m * ( za(jd,jb)*za(jd,jn)*zb(ju,jb)*zb(jd,je
     &    ) + za(jd,jn)*za(jn,jb)*zb(ju,jb)*zb(jn,je) )
      qedf(2,2,1) = qedf(2,2,1) + iza(jb,jg)*iza(jc,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2*m * (  - za(ju,jn)*za(jd,jb)*zb(ju,jb)*zb(ju
     &    ,je) + za(jd,jb)*za(jn,je)*zb(ju,je)*zb(je,jb) )
      qedf(2,2,1) = qedf(2,2,1) + iza(jc,jg)**2*sned**(-1)*sneud**(-1)*
     & fourrt2*m * (  - za(jd,jc)*za(jd,jn)*zb(ju,jb)*zb(jd,je) - za(jd
     &    ,jn)*za(jn,jc)*zb(ju,jb)*zb(jn,je) )
      qedf(2,2,1) = qedf(2,2,1) + iza(jc,jg)**2*sneu**(-1)*sneud**(-1)*
     & fourrt2*m * ( za(ju,jn)*za(jd,jc)*zb(ju,jb)*zb(ju,je) - za(jd,jc
     &    )*za(jn,je)*zb(ju,je)*zb(je,jb) )

      qedf(1,1,1)= + iza(jc,jg)*izb(jb,jg)**2*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2*m**2 * (  - za(jd,jg)*za(jd,jn)*zb(ju,jg)*
     &    zb(jd,je)*zb(jb,jc) - za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(jb,jc)
     &    *zb(jn,je) )
      qedf(1,1,1) = qedf(1,1,1) + iza(jc,jg)*izb(jb,jg)**2*izb(jc,jg)*
     & sneu**(-1)*sneud**(-1)*fourrt2*m**2 * ( za(ju,jn)*za(jd,jg)*zb(
     &    ju,jg)*zb(ju,je)*zb(jb,jc) - za(jd,jg)*za(jn,je)*zb(ju,je)*
     &    zb(jb,jc)*zb(je,jg) )
      qedf(1,1,1) = qedf(1,1,1) + izb(jb,jg)*sned**(-1)*sneud**(-1)*
     & fourrt2 * ( za(jd,jg)*za(jd,jn)*zb(ju,jc)*zb(jd,je) + za(jd,jn)*
     &    za(jn,jg)*zb(ju,jc)*zb(jn,je) )
      qedf(1,1,1) = qedf(1,1,1) + izb(jb,jg)*sneu**(-1)*sneud**(-1)*
     & fourrt2 * (  - za(ju,jn)*za(jd,jg)*zb(ju,jc)*zb(ju,je) + za(jd,
     &    jg)*za(jn,je)*zb(ju,je)*zb(je,jc) )
      qedf(1,1,1) = qedf(1,1,1) + izb(jb,jg)*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2 * (  - za(jd,jb)*za(jd,jn)*zb(ju,jc)*zb(jd,
     &    je)*zb(jb,jc) - za(jd,jn)*za(jn,jb)*zb(ju,jc)*zb(jb,jc)*zb(jn
     &    ,je) )
      qedf(1,1,1) = qedf(1,1,1) + izb(jb,jg)*izb(jc,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2 * ( za(ju,jn)*za(jd,jb)*zb(ju,jc)*zb(ju,je)*
     &    zb(jb,jc) - za(jd,jb)*za(jn,je)*zb(ju,je)*zb(jb,jc)*zb(je,jc)
     &     )

      qedf(2,1,1)= + iza(jb,jg)*iza(jc,jg)*sned**(-1)*sneud**(-1)*
     & fourrt2 * ( za(jd,jb)*za(jd,jn)*za(jb,jc)*zb(ju,jc)*zb(jd,je) + 
     &    za(jd,jn)*za(jb,jc)*za(jn,jb)*zb(ju,jc)*zb(jn,je) )
      qedf(2,1,1) = qedf(2,1,1) + iza(jb,jg)*iza(jc,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2 * (  - za(ju,jn)*za(jd,jb)*za(jb,jc)*zb(ju,
     &    jc)*zb(ju,je) + za(jd,jb)*za(jb,jc)*za(jn,je)*zb(ju,je)*zb(je
     &    ,jc) )
      qedf(2,1,1) = qedf(2,1,1) + iza(jb,jg)*iza(jc,jg)**2*izb(jb,jg)*
     & sned**(-1)*sneud**(-1)*fourrt2*m**2 * ( za(jd,jg)*za(jd,jn)*za(
     &    jb,jc)*zb(ju,jg)*zb(jd,je) + za(jd,jn)*za(jb,jc)*za(jn,jg)*
     &    zb(ju,jg)*zb(jn,je) )
      qedf(2,1,1) = qedf(2,1,1) + iza(jb,jg)*iza(jc,jg)**2*izb(jb,jg)*
     & sneu**(-1)*sneud**(-1)*fourrt2*m**2 * (  - za(ju,jn)*za(jd,jg)*
     &    za(jb,jc)*zb(ju,jg)*zb(ju,je) + za(jd,jg)*za(jb,jc)*za(jn,je)
     &    *zb(ju,je)*zb(je,jg) )
      qedf(2,1,1) = qedf(2,1,1) + iza(jc,jg)*sned**(-1)*sneud**(-1)*
     & fourrt2 * ( za(jd,jb)*za(jd,jn)*zb(ju,jg)*zb(jd,je) + za(jd,jn)*
     &    za(jn,jb)*zb(ju,jg)*zb(jn,je) )
      qedf(2,1,1) = qedf(2,1,1) + iza(jc,jg)*sneu**(-1)*sneud**(-1)*
     & fourrt2 * (  - za(ju,jn)*za(jd,jb)*zb(ju,jg)*zb(ju,je) + za(jd,
     &    jb)*za(jn,je)*zb(ju,je)*zb(je,jg) )

      qedf(1,2,1)= + iza(jb,jg)*izb(jb,jg)*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2*m * ( za(jd,jg)*za(jd,jn)*zb(ju,jc)*zb(jd,je
     &    )*zb(jb,jc) + za(jd,jn)*za(jn,jg)*zb(ju,jc)*zb(jb,jc)*zb(jn,
     &    je) )
      qedf(1,2,1) = qedf(1,2,1) + iza(jb,jg)*izb(jb,jg)*izb(jc,jg)*
     & sneu**(-1)*sneud**(-1)*fourrt2*m * (  - za(ju,jn)*za(jd,jg)*zb(
     &    ju,jc)*zb(ju,je)*zb(jb,jc) + za(jd,jg)*za(jn,je)*zb(ju,je)*
     &    zb(jb,jc)*zb(je,jc) )
      qedf(1,2,1) = qedf(1,2,1) + iza(jc,jg)*izb(jb,jg)*izb(jc,jg)*
     & sned**(-1)*sneud**(-1)*fourrt2*m * ( za(jd,jg)*za(jd,jn)*zb(ju,
     &    jb)*zb(jd,je)*zb(jb,jc) + za(jd,jn)*za(jn,jg)*zb(ju,jb)*zb(jb
     &    ,jc)*zb(jn,je) )
      qedf(1,2,1) = qedf(1,2,1) + iza(jc,jg)*izb(jb,jg)*izb(jc,jg)*
     & sneu**(-1)*sneud**(-1)*fourrt2*m * (  - za(ju,jn)*za(jd,jg)*zb(
     &    ju,jb)*zb(ju,je)*zb(jb,jc) + za(jd,jg)*za(jn,je)*zb(ju,je)*
     &    zb(jb,jc)*zb(je,jb) )

      qedf(1,1,2)= + izb(jb,jg)**2*izb(jc,jg)*sned**(-1)*sneud**(-1)*
     & fourrt2*m * ( za(jd,jc)*za(jd,jn)*zb(ju,jg)*zb(jd,je)*zb(jb,jc)
     &     + za(jd,jn)*za(jn,jc)*zb(ju,jg)*zb(jb,jc)*zb(jn,je) )
      qedf(1,1,2) = qedf(1,1,2) + izb(jb,jg)**2*izb(jc,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2*m * (  - za(ju,jn)*za(jd,jc)*zb(ju,jg)*zb(ju
     &    ,je)*zb(jb,jc) + za(jd,jc)*za(jn,je)*zb(ju,je)*zb(jb,jc)*zb(
     &    je,jg) )
      qedf(1,1,2) = qedf(1,1,2) + izb(jb,jg)*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2*m * (  - za(jd,jb)*za(jd,jn)*zb(ju,jb)*zb(jd
     &    ,je) - za(jd,jn)*za(jn,jb)*zb(ju,jb)*zb(jn,je) )
      qedf(1,1,2) = qedf(1,1,2) + izb(jb,jg)*izb(jc,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2*m * ( za(ju,jn)*za(jd,jb)*zb(ju,jb)*zb(ju,je
     &    ) - za(jd,jb)*za(jn,je)*zb(ju,je)*zb(je,jb) )
      qedf(1,1,2) = qedf(1,1,2) + izb(jc,jg)**2*sned**(-1)*sneud**(-1)*
     & fourrt2*m * ( za(jd,jb)*za(jd,jn)*zb(ju,jc)*zb(jd,je) + za(jd,jn
     &    )*za(jn,jb)*zb(ju,jc)*zb(jn,je) )
      qedf(1,1,2) = qedf(1,1,2) + izb(jc,jg)**2*sneu**(-1)*sneud**(-1)*
     & fourrt2*m * (  - za(ju,jn)*za(jd,jb)*zb(ju,jc)*zb(ju,je) + za(jd
     &    ,jb)*za(jn,je)*zb(ju,je)*zb(je,jc) )

      return
      end
