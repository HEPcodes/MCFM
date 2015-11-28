      subroutine wbbgmassa(ju,jd,jg,jb,jc,jn,je,m,qcda)
      implicit none
      include 'constants.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      double precision m,sned,sneu,sneud,snedg,sneug,sbc
      double complex qcda(2,2,2)
      double complex iza(7,7),izb(7,7)
      integer i,j,ju,jd,jg,jb,jc,jn,je
C--   order of indices is gluon helicity hg,sb,sc
      sned=s(jn,je)+s(jn,jd)+s(je,jd)
      sneu=s(jn,je)+s(jn,ju)+s(je,ju)
      sneud=sneu+s(jn,jd)+s(je,jd)+s(ju,jd)
      snedg=sned+s(jn,jg)+s(je,jg)+s(jd,jg)
      sneug=sneu+s(jn,jg)+s(je,jg)+s(ju,jg)
      sbc=sneug+s(jn,jd)+s(je,jd)+s(ju,jd)+s(jd,jg)
      do i=1,7
      do j=i+1,7
         iza(i,j)=cone/za(i,j)
         izb(i,j)=cone/zb(i,j)
         iza(j,i)=-iza(i,j)
         izb(j,i)=-izb(i,j)
      enddo
      enddo
      qcda(2,2,2)= + sned**(-1)*sneud**(-1)*sbc**(-1)*fourrt2 * (  - 
     &    za(jd,jc)*za(jd,jn)*zb(ju,jg)*zb(jd,je)*zb(jb,jg) )
      qcda(2,2,2) = qcda(2,2,2) + sneu**(-1)*sneud**(-1)*sbc**(-1)*
     & fourrt2 * ( za(ju,jn)*za(jd,jc)*zb(ju,jg)*zb(ju,je)*zb(jb,jg) - 
     &    za(jd,jc)*za(jn,je)*zb(ju,je)*zb(jb,jg)*zb(je,jg) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jd,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jc)*za(jd,jn)**2*zb(ju,jb)*zb(jd,jg)
     &    *zb(jn,je) + za(jd,jn)**2*za(jn,jc)*zb(ju,jb)*zb(jn,jg)*zb(jn
     &    ,je) + za(jd,jn)**2*za(je,jc)*zb(ju,jb)*zb(jn,je)*zb(je,jg) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jd,jg)*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(ju,jd)*za(jd,jc)*za(jd,jn)*zb(ju,jg)
     &    *zb(ju,jb)*zb(jd,je) - za(ju,jd)*za(jd,jn)*za(jn,jc)*zb(ju,jg
     &    )*zb(ju,jb)*zb(jn,je) - 1.D0/2.D0*za(ju,jc)*za(jd,jn)**2*zb(
     &    ju,jg)*zb(ju,jb)*zb(jn,je) + 1.D0/2.D0*za(jd,jc)*za(jd,jn)**2
     &    *zb(ju,jd)*zb(jb,jg)*zb(jn,je) - 1.D0/2.D0*za(jd,jc)*za(jd,jn
     &    )**2*zb(ju,jg)*zb(jd,jb)*zb(jn,je) + za(jd,jc)*za(jd,jn)**2*
     &    zb(ju,jb)*zb(jd,je)*zb(jn,jg) - 1.D0/2.D0*za(jd,jc)*za(jd,jn)
     &    **2*zb(ju,jn)*zb(jd,je)*zb(jb,jg) + za(jd,jc)*za(jd,jn)*za(jd
     &    ,je)*zb(ju,jb)*zb(jd,je)*zb(je,jg) - 1.D0/2.D0*za(jd,jc)*za(
     &    jd,jn)*za(jd,je)*zb(ju,je)*zb(jd,je)*zb(jb,jg) - za(jd,jc)*
     &    za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(jb,jg)*zb(jn,je) - 1.D0/2.D0
     &    *za(jd,jc)*za(jd,jn)*za(jn,je)*zb(ju,je)*zb(jb,jg)*zb(jn,je)
     &     + 1.D0/2.D0*za(jd,jn)**2*za(jc,jg)*zb(ju,jg)*zb(jb,jg)*zb(jn
     &    ,je) - 1.D0/2.D0*za(jd,jn)**2*za(jn,jc)*zb(ju,jg)*zb(jn,jb)*
     &    zb(jn,je) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jd,jg)*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jn)**2*za(jn,jc)*zb(ju,jb)*zb(jn,jg)
     &    *zb(jn,je) - 1.D0/2.D0*za(jd,jn)**2*za(je,jc)*zb(ju,jg)*zb(jn
     &    ,je)*zb(je,jb) + za(jd,jn)*za(jd,je)*za(jn,jc)*zb(ju,jb)*zb(
     &    jn,je)*zb(je,jg) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jd,jg)*sneu**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * ( za(ju,jd)*za(ju,jn)*za(jd,jc)*zb(ju,jg)*
     &    zb(ju,jb)*zb(ju,je) - za(ju,jd)*za(jd,jc)*za(jn,je)*zb(ju,jg)
     &    *zb(ju,je)*zb(je,jb) - 1.D0/2.D0*za(ju,jd)*za(jd,jc)*za(jn,je
     &    )*zb(ju,je)**2*zb(jb,jg) - za(ju,jn)*za(jd,jc)*za(jd,jn)*zb(
     &    ju,jb)*zb(ju,je)*zb(jn,jg) + 1.D0/2.D0*za(ju,jn)*za(jd,jc)*
     &    za(jd,jn)*zb(ju,jn)*zb(ju,je)*zb(jb,jg) - za(ju,jn)*za(jd,jc)
     &    *za(jd,je)*zb(ju,jb)*zb(ju,je)*zb(je,jg) + 1.D0/2.D0*za(ju,jn
     &    )*za(jd,jc)*za(jd,je)*zb(ju,je)**2*zb(jb,jg) + 1.D0/2.D0*za(
     &    jd,jc)*za(jd,jn)*za(jn,je)*zb(ju,je)*zb(jb,jg)*zb(jn,je) + 
     &    za(jd,jc)*za(jd,jn)*za(jn,je)*zb(ju,je)*zb(jn,jg)*zb(je,jb)
     &     + za(jd,jc)*za(jd,je)*za(jn,je)*zb(ju,je)*zb(je,jg)*zb(je,jb
     &    ) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jd,jg)*iza(jb,jg)*iza(jc,jg)*izb(
     & jc,jg)*sned**(-1)*sneud**(-1)*fourrt2*m**2 * ( za(jd,jc)*za(jd,
     &    jn)*za(jn,jg)*zb(ju,jg)*zb(jn,je) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jd,jg)*iza(jb,jg)*izb(jc,jg)*
     & sned**(-1)*snedg**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jn)**2*
     &    za(jn,jg)*zb(ju,jg)*zb(jn,jg)*zb(jn,je) + za(jd,jn)**2*za(je,
     &    jg)*zb(ju,jg)*zb(jn,je)*zb(je,jg) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jd,jg)*iza(jb,jg)*izb(jc,jg)*
     & sned**(-1)*sneud**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(ju,jd)*
     &    za(jd,jn)*za(jn,jg)*zb(ju,jg)**2*zb(jn,je) - 1.D0/2.D0*za(ju,
     &    jg)*za(jd,jn)**2*zb(ju,jg)**2*zb(jn,je) + 1.D0/2.D0*za(jd,jn)
     &    **2*za(jn,jg)*zb(ju,jg)*zb(jn,jg)*zb(jn,je) - 1.D0/2.D0*za(jd
     &    ,jn)**2*za(je,jg)*zb(ju,jg)*zb(jn,je)*zb(je,jg) + za(jd,jn)*
     &    za(jd,je)*za(jn,jg)*zb(ju,jg)*zb(jn,je)*zb(je,jg) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jd,jg)*iza(jc,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2 * ( za(jd,jc)**2*za(jd,jn)*zb(ju,jb)*zb(jd,
     &    je) + za(jd,jc)*za(jd,jn)*za(jn,jc)*zb(ju,jb)*zb(jn,je) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jd,jg)*iza(jc,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2 * (  - za(ju,jn)*za(jd,jc)**2*zb(ju,jb)*zb(
     &    ju,je) + za(jd,jc)**2*za(jn,je)*zb(ju,je)*zb(je,jb) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jb,jg)*iza(jc,jg)*izb(jc,jg)*
     & sned**(-1)*sneud**(-1)*fourrt2*m**2 * ( za(jd,jc)*za(jd,jn)*zb(
     &    ju,jg)*zb(jd,je) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jb,jg)*iza(jc,jg)*izb(jc,jg)*
     & sneu**(-1)*sneud**(-1)*fourrt2*m**2 * (  - za(ju,jn)*za(jd,jc)*
     &    zb(ju,jg)*zb(ju,je) + za(jd,jc)*za(jn,je)*zb(ju,je)*zb(je,jg)
     &     )
      qcda(2,2,2) = qcda(2,2,2) + iza(jb,jg)*izb(jc,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jn)**2*zb(ju,jg)*
     &    zb(jd,jg)*zb(jn,je) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jb,jg)*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(ju,jd)*za(jd,jn)*
     &    zb(ju,jg)**2*zb(jd,je) - 1.D0/2.D0*za(jd,jn)**2*zb(ju,jg)*zb(
     &    jd,jg)*zb(jn,je) + za(jd,jn)**2*zb(ju,jg)*zb(jd,je)*zb(jn,jg)
     &     + za(jd,jn)*za(jd,je)*zb(ju,jg)*zb(jd,je)*zb(je,jg) )
      qcda(2,2,2) = qcda(2,2,2) + iza(jb,jg)*izb(jc,jg)*sneu**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(ju,jd)*za(ju,jn)*zb(ju
     &    ,jg)**2*zb(ju,je) - za(ju,jd)*za(jn,je)*zb(ju,jg)*zb(ju,je)*
     &    zb(je,jg) - za(ju,jn)*za(jd,jn)*zb(ju,jg)*zb(ju,je)*zb(jn,jg)
     &     - za(ju,jn)*za(jd,je)*zb(ju,jg)*zb(ju,je)*zb(je,jg) + za(jd,
     &    jn)*za(jn,je)*zb(ju,je)*zb(jn,jg)*zb(je,jg) + za(jd,je)*za(jn
     &    ,je)*zb(ju,je)*zb(je,jg)**2 )

      qcda(1,2,2)= + sned**(-1)*snedg**(-1)*sbc**(-1)*fourrt2 * ( za(jd
     &    ,jg)*za(jd,jn)*za(jc,jg)*zb(ju,jb)*zb(jd,je) + za(jd,jn)*za(
     &    jc,jg)*za(jn,jg)*zb(ju,jb)*zb(jn,je) )
      qcda(1,2,2) = qcda(1,2,2) + iza(jb,jg)*izb(jd,jg)*izb(jc,jg)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(jd,jg)*za(jn,jg)*
     &    zb(ju,jg)*zb(jd,je) - za(jn,jg)**2*zb(ju,jg)*zb(jn,je) )
      qcda(1,2,2) = qcda(1,2,2) + iza(jb,jg)*izb(jd,jg)*izb(jc,jg)**2*
     & sneu**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(ju,jn)*za(jd,jg)*zb(ju,
     &    jg)*zb(ju,je)*zb(jd,jc) - za(jd,jg)*za(jn,je)*zb(ju,je)*zb(jd
     &    ,jc)*zb(je,jg) )
      qcda(1,2,2) = qcda(1,2,2) + iza(jb,jg)*izb(jd,jg)*izb(jc,jg)**2*
     & snedg**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(jd,jg)*za(jd,jn)*
     &    zb(ju,jg)*zb(jd,jc)*zb(jd,je) - za(jd,jn)*za(jn,jg)*zb(ju,jg)
     &    *zb(jd,jc)*zb(jn,je) )
      qcda(1,2,2) = qcda(1,2,2) + iza(jb,jg)*izb(jc,jg)**2*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jg)**2*za(jd,jn)*
     &    zb(ju,jg)*zb(jd,jc)*zb(jd,je) + za(jd,jg)*za(jd,jn)*za(jn,jg)
     &    *zb(ju,jg)*zb(jd,jc)*zb(jn,je) + za(jd,jg)*za(jd,jn)*za(jn,jg
     &    )*zb(ju,jg)*zb(jd,je)*zb(jn,jc) + za(jd,jg)*za(jd,jn)*za(je,
     &    jg)*zb(ju,jg)*zb(jd,je)*zb(je,jc) + za(jd,jn)*za(jn,jg)**2*
     &    zb(ju,jg)*zb(jn,jc)*zb(jn,je) + za(jd,jn)*za(jn,jg)*za(je,jg)
     &    *zb(ju,jg)*zb(jn,je)*zb(je,jc) )
      qcda(1,2,2) = qcda(1,2,2) + iza(jb,jg)*izb(jc,jg)**2*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m**2 * ( 1.D0/2.D0*za(ju,jg)*za(jd
     &    ,jg)*za(jd,jn)*zb(ju,jg)*zb(ju,jc)*zb(jd,je) + 1.D0/2.D0*za(
     &    ju,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(ju,jc)*zb(jn,je) + 
     &    za(jd,jg)**2*za(jd,jn)*zb(ju,jg)*zb(jd,jc)*zb(jd,je) - 1.D0/2.
     &    D0*za(jd,jg)**2*za(jd,jn)*zb(ju,jc)*zb(jd,jg)*zb(jd,je) + za(
     &    jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(jd,jc)*zb(jn,je) + 
     &    za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(jd,je)*zb(jn,jc)
     &     - 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jc)*zb(jd,jg
     &    )*zb(jn,je) - 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,
     &    jc)*zb(jd,je)*zb(jn,jg) + za(jd,jg)*za(jd,jn)*za(je,jg)*zb(ju
     &    ,jg)*zb(jd,je)*zb(je,jc) - 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(
     &    je,jg)*zb(ju,jc)*zb(jd,je)*zb(je,jg) + za(jd,jn)*za(jn,jg)**2
     &    *zb(ju,jg)*zb(jn,jc)*zb(jn,je) - 1.D0/2.D0*za(jd,jn)*za(jn,jg
     &    )**2*zb(ju,jc)*zb(jn,jg)*zb(jn,je) + za(jd,jn)*za(jn,jg)*za(
     &    je,jg)*zb(ju,jg)*zb(jn,je)*zb(je,jc) )
      qcda(1,2,2) = qcda(1,2,2) + iza(jb,jg)*izb(jc,jg)**2*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m**2 * (  - 1.D0/2.D0*za(jd,jn)*
     &    za(jn,jg)*za(je,jg)*zb(ju,jc)*zb(jn,je)*zb(je,jg) )
      qcda(1,2,2) = qcda(1,2,2) + iza(jb,jg)*izb(jc,jg)**2*sneu**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m**2 * (  - 1.D0/2.D0*za(ju,jg)*
     &    za(ju,jn)*za(jd,jg)*zb(ju,jg)*zb(ju,jc)*zb(ju,je) - 1.D0/2.D0
     &    *za(ju,jg)*za(jd,jg)*za(jn,je)*zb(ju,jg)*zb(ju,je)*zb(je,jc)
     &     + za(ju,jg)*za(jd,jg)*za(jn,je)*zb(ju,jc)*zb(ju,je)*zb(je,jg
     &    ) - za(ju,jn)*za(jd,jg)**2*zb(ju,jg)*zb(ju,je)*zb(jd,jc) + 1.D
     &    0/2.D0*za(ju,jn)*za(jd,jg)**2*zb(ju,jc)*zb(ju,je)*zb(jd,jg)
     &     - za(ju,jn)*za(jd,jg)*za(jn,jg)*zb(ju,jg)*zb(ju,je)*zb(jn,jc
     &    ) + 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(jn,jg)*zb(ju,jc)*zb(ju,
     &    je)*zb(jn,jg) - za(ju,jn)*za(jd,jg)*za(je,jg)*zb(ju,jg)*zb(ju
     &    ,je)*zb(je,jc) + 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(je,jg)*zb(
     &    ju,jc)*zb(ju,je)*zb(je,jg) - 1.D0/2.D0*za(jd,jg)**2*za(jn,je)
     &    *zb(ju,je)*zb(jd,jg)*zb(je,jc) + za(jd,jg)**2*za(jn,je)*zb(ju
     &    ,je)*zb(jd,jc)*zb(je,jg) - 1.D0/2.D0*za(jd,jg)*za(jn,jg)*za(
     &    jn,je)*zb(ju,je)*zb(jn,jg)*zb(je,jc) + za(jd,jg)*za(jn,jg)*
     &    za(jn,je)*zb(ju,je)*zb(jn,jc)*zb(je,jg) )
      qcda(1,2,2) = qcda(1,2,2) + iza(jb,jg)*izb(jc,jg)**2*sneu**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m**2 * ( 1.D0/2.D0*za(jd,jg)*za(jn
     &    ,je)*za(je,jg)*zb(ju,je)*zb(je,jg)*zb(je,jc) )
      qcda(1,2,2) = qcda(1,2,2) + izb(jd,jg)*sneu**(-1)*sbc**(-1)*
     & fourrt2 * ( za(ju,jn)*za(jc,jg)*zb(ju,jb)*zb(ju,je) - za(jc,jg)*
     &    za(jn,je)*zb(ju,je)*zb(je,jb) )
      qcda(1,2,2) = qcda(1,2,2) + izb(jd,jg)*snedg**(-1)*sbc**(-1)*
     & fourrt2 * (  - za(jd,jc)*za(jn,jg)*zb(ju,jb)*zb(jd,je) - za(jc,
     &    jg)*za(jn,jg)*zb(ju,jb)*zb(je,jg) - za(jn,jg)*za(jn,jc)*zb(ju
     &    ,jb)*zb(jn,je) )
      qcda(1,2,2) = qcda(1,2,2) + izb(jd,jg)*izb(jc,jg)*sneu**(-1)*
     & sbc**(-1)*fourrt2 * ( za(ju,jn)*za(jd,jc)*zb(ju,jb)*zb(ju,je)*
     &    zb(jd,jc) - za(jd,jc)*za(jn,je)*zb(ju,je)*zb(jd,jc)*zb(je,jb)
     &     )
      qcda(1,2,2) = qcda(1,2,2) + izb(jd,jg)*izb(jc,jg)*snedg**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(jd,jc)*za(jd,jn)*zb(ju,jb)*zb(jd,jc)
     &    *zb(jd,je) - za(jd,jn)*za(jc,jg)*zb(ju,jb)*zb(jd,jc)*zb(je,jg
     &    ) - za(jd,jn)*za(jn,jc)*zb(ju,jb)*zb(jd,jc)*zb(jn,je) )
      qcda(1,2,2) = qcda(1,2,2) + izb(jc,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jg)*za(jd,jc)*za(jd,jn)*zb(ju,jb)*
     &    zb(jd,jc)*zb(jd,je) + za(jd,jg)*za(jd,jn)*za(jn,jc)*zb(ju,jb)
     &    *zb(jd,je)*zb(jn,jc) + za(jd,jg)*za(jd,jn)*za(je,jc)*zb(ju,jb
     &    )*zb(jd,je)*zb(je,jc) + za(jd,jc)*za(jd,jn)*za(jn,jg)*zb(ju,
     &    jb)*zb(jd,jc)*zb(jn,je) + za(jd,jn)*za(jn,jg)*za(jn,jc)*zb(ju
     &    ,jb)*zb(jn,jc)*zb(jn,je) + za(jd,jn)*za(jn,jg)*za(je,jc)*zb(
     &    ju,jb)*zb(jn,je)*zb(je,jc) )
      qcda(1,2,2) = qcda(1,2,2) + izb(jc,jg)*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * ( za(ju,jg)*za(jd,jc)*za(jd,jn)*zb(ju,jb)*
     &    zb(ju,jc)*zb(jd,je) + za(ju,jg)*za(jd,jn)*za(jn,jc)*zb(ju,jb)
     &    *zb(ju,jc)*zb(jn,je) - 1.D0/2.D0*za(ju,jc)*za(jd,jg)*za(jd,jn
     &    )*zb(ju,jb)*zb(ju,jc)*zb(jd,je) - 1.D0/2.D0*za(ju,jc)*za(jd,
     &    jn)*za(jn,jg)*zb(ju,jb)*zb(ju,jc)*zb(jn,je) + za(jd,jg)*za(jd
     &    ,jc)*za(jd,jn)*zb(ju,jb)*zb(jd,jc)*zb(jd,je) - 1.D0/2.D0*za(
     &    jd,jg)*za(jd,jc)*za(jd,jn)*zb(ju,jc)*zb(jd,jb)*zb(jd,je) - 
     &    za(jd,jg)*za(jd,jn)*za(jc,jg)*zb(ju,jg)*zb(jd,je)*zb(jb,jc)
     &     + 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(jc,jg)*zb(ju,jc)*zb(jd,je
     &    )*zb(jb,jg) + za(jd,jg)*za(jd,jn)*za(jn,jc)*zb(ju,jb)*zb(jd,
     &    jc)*zb(jn,je) - 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(jn,jc)*zb(ju
     &    ,jc)*zb(jd,je)*zb(jn,jb) - 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(
     &    je,jc)*zb(ju,jc)*zb(jd,je)*zb(je,jb) + za(jd,jc)*za(jd,jn)*
     &    za(jn,jg)*zb(ju,jb)*zb(jd,je)*zb(jn,jc) - 1.D0/2.D0*za(jd,jc)
     &    *za(jd,jn)*za(jn,jg)*zb(ju,jc)*zb(jd,jb)*zb(jn,je) )
      qcda(1,2,2) = qcda(1,2,2) + izb(jc,jg)*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jc)*za(jd,jn)*za(je,jg)*zb(ju,jb)*
     &    zb(jd,je)*zb(je,jc) + 1.D0/2.D0*za(jd,jn)**2*za(jc,jg)*zb(ju,
     &    jd)*zb(jb,jc)*zb(jn,je) - 1.D0/2.D0*za(jd,jn)**2*za(jc,jg)*
     &    zb(ju,jn)*zb(jd,je)*zb(jb,jc) - 1.D0/2.D0*za(jd,jn)*za(jd,je)
     &    *za(jc,jg)*zb(ju,je)*zb(jd,je)*zb(jb,jc) - za(jd,jn)*za(jc,jg
     &    )*za(jn,jg)*zb(ju,jg)*zb(jb,jc)*zb(jn,je) + 1.D0/2.D0*za(jd,
     &    jn)*za(jc,jg)*za(jn,jg)*zb(ju,jc)*zb(jb,jg)*zb(jn,je) - 1.D0/
     &    2.D0*za(jd,jn)*za(jc,jg)*za(jn,je)*zb(ju,je)*zb(jb,jc)*zb(jn,
     &    je) + za(jd,jn)*za(jn,jg)*za(jn,jc)*zb(ju,jb)*zb(jn,jc)*zb(jn
     &    ,je) - 1.D0/2.D0*za(jd,jn)*za(jn,jg)*za(jn,jc)*zb(ju,jc)*zb(
     &    jn,jb)*zb(jn,je) - 1.D0/2.D0*za(jd,jn)*za(jn,jg)*za(je,jc)*
     &    zb(ju,jc)*zb(jn,je)*zb(je,jb) + za(jd,jn)*za(jn,jc)*za(je,jg)
     &    *zb(ju,jb)*zb(jn,je)*zb(je,jc) )
      qcda(1,2,2) = qcda(1,2,2) + izb(jc,jg)*sned**(-1)*sneud**(-1)*
     & fourrt2 * (  - za(jd,jg)*za(jd,jn)*zb(ju,jb)*zb(jd,je) - za(jd,
     &    jn)*za(jn,jg)*zb(ju,jb)*zb(jn,je) )
      qcda(1,2,2) = qcda(1,2,2) + izb(jc,jg)*sneu**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * (  - 1.D0/2.D0*za(ju,jd)*za(jc,jg)*za(jn,je)
     &    *zb(ju,je)**2*zb(jb,jc) - za(ju,jg)*za(ju,jn)*za(jd,jc)*zb(ju
     &    ,jb)*zb(ju,jc)*zb(ju,je) + za(ju,jg)*za(jd,jc)*za(jn,je)*zb(
     &    ju,jc)*zb(ju,je)*zb(je,jb) + 1.D0/2.D0*za(ju,jc)*za(ju,jn)*
     &    za(jd,jg)*zb(ju,jb)*zb(ju,jc)*zb(ju,je) - 1.D0/2.D0*za(ju,jc)
     &    *za(jd,jg)*za(jn,je)*zb(ju,jb)*zb(ju,je)*zb(je,jc) - za(ju,jn
     &    )*za(jd,jg)*za(jd,jc)*zb(ju,jb)*zb(ju,je)*zb(jd,jc) + 1.D0/2.D
     &    0*za(ju,jn)*za(jd,jg)*za(jd,jc)*zb(ju,jc)*zb(ju,je)*zb(jd,jb)
     &     + za(ju,jn)*za(jd,jg)*za(jc,jg)*zb(ju,jg)*zb(ju,je)*zb(jb,jc
     &    ) - 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(jc,jg)*zb(ju,jc)*zb(ju,
     &    je)*zb(jb,jg) + 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(jn,jc)*zb(ju
     &    ,jc)*zb(ju,je)*zb(jn,jb) + 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(
     &    je,jc)*zb(ju,jc)*zb(ju,je)*zb(je,jb) - za(ju,jn)*za(jd,jc)*
     &    za(jn,jg)*zb(ju,jb)*zb(ju,je)*zb(jn,jc) - za(ju,jn)*za(jd,jc)
     &    *za(je,jg)*zb(ju,jb)*zb(ju,je)*zb(je,jc) )
      qcda(1,2,2) = qcda(1,2,2) + izb(jc,jg)*sneu**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * ( 1.D0/2.D0*za(ju,jn)*za(jd,jn)*za(jc,jg)*
     &    zb(ju,jn)*zb(ju,je)*zb(jb,jc) + 1.D0/2.D0*za(ju,jn)*za(jd,je)
     &    *za(jc,jg)*zb(ju,je)**2*zb(jb,jc) - 1.D0/2.D0*za(jd,jg)*za(jd
     &    ,jc)*za(jn,je)*zb(ju,je)*zb(jd,jb)*zb(je,jc) + za(jd,jg)*za(
     &    jd,jc)*za(jn,je)*zb(ju,je)*zb(jd,jc)*zb(je,jb) + 1.D0/2.D0*
     &    za(jd,jg)*za(jc,jg)*za(jn,je)*zb(ju,je)*zb(jb,jg)*zb(je,jc)
     &     - za(jd,jg)*za(jc,jg)*za(jn,je)*zb(ju,je)*zb(jb,jc)*zb(je,jg
     &    ) - 1.D0/2.D0*za(jd,jg)*za(jn,jc)*za(jn,je)*zb(ju,je)*zb(jn,
     &    jb)*zb(je,jc) - 1.D0/2.D0*za(jd,jg)*za(jn,je)*za(je,jc)*zb(ju
     &    ,je)*zb(je,jb)*zb(je,jc) + za(jd,jc)*za(jn,jg)*za(jn,je)*zb(
     &    ju,je)*zb(jn,jc)*zb(je,jb) + za(jd,jc)*za(jn,je)*za(je,jg)*
     &    zb(ju,je)*zb(je,jb)*zb(je,jc) + 1.D0/2.D0*za(jd,jn)*za(jc,jg)
     &    *za(jn,je)*zb(ju,je)*zb(jb,jc)*zb(jn,je) )
      qcda(1,2,2) = qcda(1,2,2) + izb(jc,jg)*sneu**(-1)*sneud**(-1)*
     & fourrt2 * ( za(ju,jn)*za(jd,jg)*zb(ju,jb)*zb(ju,je) - za(jd,jg)*
     &    za(jn,je)*zb(ju,je)*zb(je,jb) )

      qcda(2,1,2)= + iza(jd,jg)*iza(jc,jg)*izb(jb,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2*m * (  - za(jd,jc)**2*za(jd,jn)*zb(ju,jg)*
     &    zb(jd,je) - za(jd,jc)*za(jd,jn)*za(jn,jc)*zb(ju,jg)*zb(jn,je)
     &     )
      qcda(2,1,2) = qcda(2,1,2) + iza(jd,jg)*iza(jc,jg)*izb(jb,jg)*
     & sneu**(-1)*sneud**(-1)*fourrt2*m * ( za(ju,jn)*za(jd,jc)**2*zb(
     &    ju,jg)*zb(ju,je) - za(jd,jc)**2*za(jn,je)*zb(ju,je)*zb(je,jg)
     &     )
      qcda(2,1,2) = qcda(2,1,2) + iza(jd,jg)*iza(jc,jg)*izb(jc,jg)*
     & sned**(-1)*sneud**(-1)*fourrt2*m * (  - za(jd,jb)*za(jd,jc)*za(
     &    jd,jn)*zb(ju,jg)*zb(jd,je) - za(jd,jc)*za(jd,jn)*za(jn,jb)*
     &    zb(ju,jg)*zb(jn,je) )
      qcda(2,1,2) = qcda(2,1,2) + iza(jd,jg)*iza(jc,jg)*izb(jc,jg)*
     & sneu**(-1)*sneud**(-1)*fourrt2*m * ( za(ju,jn)*za(jd,jb)*za(jd,
     &    jc)*zb(ju,jg)*zb(ju,je) - za(jd,jb)*za(jd,jc)*za(jn,je)*zb(ju
     &    ,je)*zb(je,jg) )
      qcda(2,1,2) = qcda(2,1,2) + iza(jd,jg)*izb(jb,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jc)*za(jd,jn)**2*
     &    zb(ju,jg)*zb(jd,jg)*zb(jn,je) - za(jd,jn)**2*za(jn,jc)*zb(ju,
     &    jg)*zb(jn,jg)*zb(jn,je) - za(jd,jn)**2*za(je,jc)*zb(ju,jg)*
     &    zb(jn,je)*zb(je,jg) )
      qcda(2,1,2) = qcda(2,1,2) + iza(jd,jg)*izb(jb,jg)*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(jd,jc)*za(jd,jn
     &    )*zb(ju,jg)**2*zb(jd,je) + za(ju,jd)*za(jd,jn)*za(jn,jc)*zb(
     &    ju,jg)**2*zb(jn,je) + 1.D0/2.D0*za(ju,jc)*za(jd,jn)**2*zb(ju,
     &    jg)**2*zb(jn,je) + 1.D0/2.D0*za(jd,jc)*za(jd,jn)**2*zb(ju,jg)
     &    *zb(jd,jg)*zb(jn,je) - za(jd,jc)*za(jd,jn)**2*zb(ju,jg)*zb(jd
     &    ,je)*zb(jn,jg) - za(jd,jc)*za(jd,jn)*za(jd,je)*zb(ju,jg)*zb(
     &    jd,je)*zb(je,jg) - 1.D0/2.D0*za(jd,jn)**2*za(jn,jc)*zb(ju,jg)
     &    *zb(jn,jg)*zb(jn,je) + 1.D0/2.D0*za(jd,jn)**2*za(je,jc)*zb(ju
     &    ,jg)*zb(jn,je)*zb(je,jg) - za(jd,jn)*za(jd,je)*za(jn,jc)*zb(
     &    ju,jg)*zb(jn,je)*zb(je,jg) )
      qcda(2,1,2) = qcda(2,1,2) + iza(jd,jg)*izb(jb,jg)*sneu**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * (  - za(ju,jd)*za(ju,jn)*za(jd
     &    ,jc)*zb(ju,jg)**2*zb(ju,je) + za(ju,jd)*za(jd,jc)*za(jn,je)*
     &    zb(ju,jg)*zb(ju,je)*zb(je,jg) + za(ju,jn)*za(jd,jc)*za(jd,jn)
     &    *zb(ju,jg)*zb(ju,je)*zb(jn,jg) + za(ju,jn)*za(jd,jc)*za(jd,je
     &    )*zb(ju,jg)*zb(ju,je)*zb(je,jg) - za(jd,jc)*za(jd,jn)*za(jn,
     &    je)*zb(ju,je)*zb(jn,jg)*zb(je,jg) - za(jd,jc)*za(jd,je)*za(jn
     &    ,je)*zb(ju,je)*zb(je,jg)**2 )
      qcda(2,1,2) = qcda(2,1,2) + iza(jd,jg)*izb(jc,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jb)*za(jd,jn)**2*
     &    zb(ju,jg)*zb(jd,jg)*zb(jn,je) - za(jd,jn)**2*za(jn,jb)*zb(ju,
     &    jg)*zb(jn,jg)*zb(jn,je) - za(jd,jn)**2*za(je,jb)*zb(ju,jg)*
     &    zb(jn,je)*zb(je,jg) )
      qcda(2,1,2) = qcda(2,1,2) + iza(jd,jg)*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(jd,jb)*za(jd,jn
     &    )*zb(ju,jg)**2*zb(jd,je) + za(ju,jd)*za(jd,jn)*za(jn,jb)*zb(
     &    ju,jg)**2*zb(jn,je) + 1.D0/2.D0*za(ju,jb)*za(jd,jn)**2*zb(ju,
     &    jg)**2*zb(jn,je) + 1.D0/2.D0*za(jd,jb)*za(jd,jn)**2*zb(ju,jg)
     &    *zb(jd,jg)*zb(jn,je) - za(jd,jb)*za(jd,jn)**2*zb(ju,jg)*zb(jd
     &    ,je)*zb(jn,jg) - za(jd,jb)*za(jd,jn)*za(jd,je)*zb(ju,jg)*zb(
     &    jd,je)*zb(je,jg) - 1.D0/2.D0*za(jd,jn)**2*za(jn,jb)*zb(ju,jg)
     &    *zb(jn,jg)*zb(jn,je) + 1.D0/2.D0*za(jd,jn)**2*za(je,jb)*zb(ju
     &    ,jg)*zb(jn,je)*zb(je,jg) - za(jd,jn)*za(jd,je)*za(jn,jb)*zb(
     &    ju,jg)*zb(jn,je)*zb(je,jg) )
      qcda(2,1,2) = qcda(2,1,2) + iza(jd,jg)*izb(jc,jg)*sneu**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * (  - za(ju,jd)*za(ju,jn)*za(jd
     &    ,jb)*zb(ju,jg)**2*zb(ju,je) + za(ju,jd)*za(jd,jb)*za(jn,je)*
     &    zb(ju,jg)*zb(ju,je)*zb(je,jg) + za(ju,jn)*za(jd,jb)*za(jd,jn)
     &    *zb(ju,jg)*zb(ju,je)*zb(jn,jg) + za(ju,jn)*za(jd,jb)*za(jd,je
     &    )*zb(ju,jg)*zb(ju,je)*zb(je,jg) - za(jd,jb)*za(jd,jn)*za(jn,
     &    je)*zb(ju,je)*zb(jn,jg)*zb(je,jg) - za(jd,jb)*za(jd,je)*za(jn
     &    ,je)*zb(ju,je)*zb(je,jg)**2 )

      qcda(2,2,1)= + iza(jd,jg)*iza(jb,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jn)**2*za(jn,jg)*zb(ju,jc)*zb(
     &    jn,jg)*zb(jn,je) - za(jd,jn)**2*za(je,jg)*zb(ju,jc)*zb(jn,je)
     &    *zb(je,jg) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jd,jg)*iza(jb,jg)*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(jd,jn)*za(jn,jg
     &    )*zb(ju,jg)*zb(ju,jc)*zb(jn,je) + 1.D0/2.D0*za(ju,jg)*za(jd,
     &    jn)**2*zb(ju,jg)*zb(ju,jc)*zb(jn,je) + 1.D0/2.D0*za(jd,jn)**2
     &    *za(jn,jg)*zb(ju,jg)*zb(jn,jc)*zb(jn,je) - za(jd,jn)**2*za(jn
     &    ,jg)*zb(ju,jc)*zb(jn,jg)*zb(jn,je) + 1.D0/2.D0*za(jd,jn)**2*
     &    za(je,jg)*zb(ju,jg)*zb(jn,je)*zb(je,jc) - za(jd,jn)*za(jd,je)
     &    *za(jn,jg)*zb(ju,jc)*zb(jn,je)*zb(je,jg) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jd,jg)*iza(jb,jg)*iza(jc,jg)*
     & sned**(-1)*sneud**(-1)*fourrt2*m * (  - za(jd,jc)*za(jd,jn)*za(
     &    jn,jg)*zb(ju,jc)*zb(jn,je) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jd,jg)*iza(jc,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jn)**2*za(jn,jg)*
     &    zb(ju,jb)*zb(jn,jg)*zb(jn,je) - za(jd,jn)**2*za(je,jg)*zb(ju,
     &    jb)*zb(jn,je)*zb(je,jg) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jd,jg)*iza(jc,jg)*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(jd,jn)*za(jn,jg
     &    )*zb(ju,jg)*zb(ju,jb)*zb(jn,je) + 1.D0/2.D0*za(ju,jg)*za(jd,
     &    jn)**2*zb(ju,jg)*zb(ju,jb)*zb(jn,je) + 1.D0/2.D0*za(jd,jn)**2
     &    *za(jn,jg)*zb(ju,jg)*zb(jn,jb)*zb(jn,je) - za(jd,jn)**2*za(jn
     &    ,jg)*zb(ju,jb)*zb(jn,jg)*zb(jn,je) + 1.D0/2.D0*za(jd,jn)**2*
     &    za(je,jg)*zb(ju,jg)*zb(jn,je)*zb(je,jb) - za(jd,jn)*za(jd,je)
     &    *za(jn,jg)*zb(ju,jb)*zb(jn,je)*zb(je,jg) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jd,jg)*iza(jc,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2*m * (  - za(jd,jn)**2*zb(ju,jb)*zb(jn,je) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jb,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jn)**2*zb(ju,jc)*zb(jd,jg)*zb(
     &    jn,je) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jb,jg)*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(jd,jn)*zb(ju,jg)*zb(ju,jc)*
     &    zb(jd,je) + za(jd,jg)*za(jd,jn)*zb(ju,jg)*zb(jd,je)*zb(jc,jg)
     &     - 1.D0/2.D0*za(jd,jn)**2*zb(ju,jd)*zb(jc,jg)*zb(jn,je) + 1.D0
     &    /2.D0*za(jd,jn)**2*zb(ju,jg)*zb(jd,jc)*zb(jn,je) - za(jd,jn)
     &    **2*zb(ju,jc)*zb(jd,je)*zb(jn,jg) + 1.D0/2.D0*za(jd,jn)**2*
     &    zb(ju,jn)*zb(jd,je)*zb(jc,jg) - za(jd,jn)*za(jd,je)*zb(ju,jc)
     &    *zb(jd,je)*zb(je,jg) + 1.D0/2.D0*za(jd,jn)*za(jd,je)*zb(ju,je
     &    )*zb(jd,je)*zb(jc,jg) + za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(jc,
     &    jg)*zb(jn,je) + 1.D0/2.D0*za(jd,jn)*za(jn,je)*zb(ju,je)*zb(jc
     &    ,jg)*zb(jn,je) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jb,jg)*sneu**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(ju,jd)*za(ju,jn)*zb(ju,jg)*zb(ju,
     &    jc)*zb(ju,je) + za(ju,jd)*za(jn,je)*zb(ju,jg)*zb(ju,je)*zb(je
     &    ,jc) + 1.D0/2.D0*za(ju,jd)*za(jn,je)*zb(ju,je)**2*zb(jc,jg)
     &     - za(ju,jn)*za(jd,jg)*zb(ju,jg)*zb(ju,je)*zb(jc,jg) + za(ju,
     &    jn)*za(jd,jn)*zb(ju,jc)*zb(ju,je)*zb(jn,jg) - 1.D0/2.D0*za(ju
     &    ,jn)*za(jd,jn)*zb(ju,jn)*zb(ju,je)*zb(jc,jg) + za(ju,jn)*za(
     &    jd,je)*zb(ju,jc)*zb(ju,je)*zb(je,jg) - 1.D0/2.D0*za(ju,jn)*
     &    za(jd,je)*zb(ju,je)**2*zb(jc,jg) + za(jd,jg)*za(jn,je)*zb(ju,
     &    je)*zb(jc,jg)*zb(je,jg) - 1.D0/2.D0*za(jd,jn)*za(jn,je)*zb(ju
     &    ,je)*zb(jc,jg)*zb(jn,je) - za(jd,jn)*za(jn,je)*zb(ju,je)*zb(
     &    jn,jg)*zb(je,jc) - za(jd,je)*za(jn,je)*zb(ju,je)*zb(je,jg)*
     &    zb(je,jc) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jb,jg)*iza(jc,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2*m * (  - za(jd,jg)*za(jd,jn)*zb(ju,jg)*zb(jd
     &    ,je) - za(jd,jc)*za(jd,jn)*zb(ju,jc)*zb(jd,je) - za(jd,jn)*
     &    za(jn,jg)*zb(ju,jg)*zb(jn,je) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jb,jg)*iza(jc,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2*m * ( za(ju,jn)*za(jd,jg)*zb(ju,jg)*zb(ju,je
     &    ) + za(ju,jn)*za(jd,jc)*zb(ju,jc)*zb(ju,je) - za(jd,jg)*za(jn
     &    ,je)*zb(ju,je)*zb(je,jg) - za(jd,jc)*za(jn,je)*zb(ju,je)*zb(
     &    je,jc) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jc,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jn)**2*zb(ju,jb)*zb(jd,jg)*zb(
     &    jn,je) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jc,jg)*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(jd,jn)*zb(ju,jg)*zb(ju,jb)*
     &    zb(jd,je) + za(jd,jg)*za(jd,jn)*zb(ju,jg)*zb(jd,je)*zb(jb,jg)
     &     - 1.D0/2.D0*za(jd,jn)**2*zb(ju,jd)*zb(jb,jg)*zb(jn,je) + 1.D0
     &    /2.D0*za(jd,jn)**2*zb(ju,jg)*zb(jd,jb)*zb(jn,je) - za(jd,jn)
     &    **2*zb(ju,jb)*zb(jd,je)*zb(jn,jg) + 1.D0/2.D0*za(jd,jn)**2*
     &    zb(ju,jn)*zb(jd,je)*zb(jb,jg) - za(jd,jn)*za(jd,je)*zb(ju,jb)
     &    *zb(jd,je)*zb(je,jg) + 1.D0/2.D0*za(jd,jn)*za(jd,je)*zb(ju,je
     &    )*zb(jd,je)*zb(jb,jg) + za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(jb,
     &    jg)*zb(jn,je) + 1.D0/2.D0*za(jd,jn)*za(jn,je)*zb(ju,je)*zb(jb
     &    ,jg)*zb(jn,je) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jc,jg)*sneu**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(ju,jd)*za(ju,jn)*zb(ju,jg)*zb(ju,
     &    jb)*zb(ju,je) + za(ju,jd)*za(jn,je)*zb(ju,jg)*zb(ju,je)*zb(je
     &    ,jb) + 1.D0/2.D0*za(ju,jd)*za(jn,je)*zb(ju,je)**2*zb(jb,jg)
     &     - za(ju,jn)*za(jd,jg)*zb(ju,jg)*zb(ju,je)*zb(jb,jg) + za(ju,
     &    jn)*za(jd,jn)*zb(ju,jb)*zb(ju,je)*zb(jn,jg) - 1.D0/2.D0*za(ju
     &    ,jn)*za(jd,jn)*zb(ju,jn)*zb(ju,je)*zb(jb,jg) + za(ju,jn)*za(
     &    jd,je)*zb(ju,jb)*zb(ju,je)*zb(je,jg) - 1.D0/2.D0*za(ju,jn)*
     &    za(jd,je)*zb(ju,je)**2*zb(jb,jg) + za(jd,jg)*za(jn,je)*zb(ju,
     &    je)*zb(jb,jg)*zb(je,jg) - 1.D0/2.D0*za(jd,jn)*za(jn,je)*zb(ju
     &    ,je)*zb(jb,jg)*zb(jn,je) - za(jd,jn)*za(jn,je)*zb(ju,je)*zb(
     &    jn,jg)*zb(je,jb) - za(jd,je)*za(jn,je)*zb(ju,je)*zb(je,jg)*
     &    zb(je,jb) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jc,jg)**2*sned**(-1)*sneud**(-1)*
     & fourrt2*m * (  - za(jd,jc)*za(jd,jn)*zb(ju,jb)*zb(jd,je) - za(jd
     &    ,jn)*za(jn,jc)*zb(ju,jb)*zb(jn,je) )
      qcda(2,2,1) = qcda(2,2,1) + iza(jc,jg)**2*sneu**(-1)*sneud**(-1)*
     & fourrt2*m * ( za(ju,jn)*za(jd,jc)*zb(ju,jb)*zb(ju,je) - za(jd,jc
     &    )*za(jn,je)*zb(ju,je)*zb(je,jb) )

      qcda(1,1,1)= + sned**(-1)*snedg**(-1)*sbc**(-1)*fourrt2 * ( za(jd
     &    ,jg)*za(jd,jn)*za(jb,jg)*zb(ju,jc)*zb(jd,je) + za(jd,jn)*za(
     &    jb,jg)*za(jn,jg)*zb(ju,jc)*zb(jn,je) )
      qcda(1,1,1) = qcda(1,1,1) + sned**(-1)*sneud**(-1)*sbc**(-1)*
     & fourrt2 * ( 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(jb,jg)*zb(ju,jc)*
     &    zb(jd,je) + 1.D0/2.D0*za(jd,jn)*za(jb,jg)*za(jn,jg)*zb(ju,jc)
     &    *zb(jn,je) )
      qcda(1,1,1) = qcda(1,1,1) + sneu**(-1)*sneud**(-1)*sbc**(-1)*
     & fourrt2 * (  - 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(jb,jg)*zb(ju,jc)
     &    *zb(ju,je) + 1.D0/2.D0*za(jd,jg)*za(jb,jg)*za(jn,je)*zb(ju,je
     &    )*zb(je,jc) )
      qcda(1,1,1) = qcda(1,1,1) + iza(jc,jg)*izb(jd,jg)*izb(jb,jg)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(jd,jg)*za(jn,jg)*
     &    zb(ju,jg)*zb(jd,je) - za(jn,jg)**2*zb(ju,jg)*zb(jn,je) )
      qcda(1,1,1) = qcda(1,1,1) + iza(jc,jg)*izb(jd,jg)*izb(jb,jg)*izb(
     & jc,jg)*sneu**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(ju,jn)*za(jd,jg)
     &    *zb(ju,jg)*zb(ju,je)*zb(jd,jc) - za(jd,jg)*za(jn,je)*zb(ju,je
     &    )*zb(jd,jc)*zb(je,jg) )
      qcda(1,1,1) = qcda(1,1,1) + iza(jc,jg)*izb(jd,jg)*izb(jb,jg)*izb(
     & jc,jg)*snedg**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(jd,jg)*za(jd
     &    ,jn)*zb(ju,jg)*zb(jd,jc)*zb(jd,je) - za(jd,jn)*za(jn,jg)*zb(
     &    ju,jg)*zb(jd,jc)*zb(jn,je) )
      qcda(1,1,1) = qcda(1,1,1) + iza(jc,jg)*izb(jb,jg)*izb(jc,jg)*
     & sned**(-1)*snedg**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jg)**2*
     &    za(jd,jn)*zb(ju,jg)*zb(jd,jc)*zb(jd,je) + za(jd,jg)*za(jd,jn)
     &    *za(jn,jg)*zb(ju,jg)*zb(jd,jc)*zb(jn,je) + za(jd,jg)*za(jd,jn
     &    )*za(jn,jg)*zb(ju,jg)*zb(jd,je)*zb(jn,jc) + za(jd,jg)*za(jd,
     &    jn)*za(je,jg)*zb(ju,jg)*zb(jd,je)*zb(je,jc) + za(jd,jn)*za(jn
     &    ,jg)**2*zb(ju,jg)*zb(jn,jc)*zb(jn,je) + za(jd,jn)*za(jn,jg)*
     &    za(je,jg)*zb(ju,jg)*zb(jn,je)*zb(je,jc) )
      qcda(1,1,1) = qcda(1,1,1) + iza(jc,jg)*izb(jb,jg)*izb(jc,jg)*
     & sned**(-1)*sneud**(-1)*sbc**(-1)*fourrt2*m**2 * ( 1.D0/2.D0*za(
     &    ju,jg)*za(jd,jg)*za(jd,jn)*zb(ju,jg)*zb(ju,jc)*zb(jd,je) + 1.D
     &    0/2.D0*za(ju,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(ju,jc)*zb(
     &    jn,je) + za(jd,jg)**2*za(jd,jn)*zb(ju,jg)*zb(jd,jc)*zb(jd,je)
     &     - 1.D0/2.D0*za(jd,jg)**2*za(jd,jn)*zb(ju,jc)*zb(jd,jg)*zb(jd
     &    ,je) + za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(jd,jc)*zb(
     &    jn,je) + za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(jd,je)*
     &    zb(jn,jc) - 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jc)
     &    *zb(jd,jg)*zb(jn,je) - 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(jn,jg
     &    )*zb(ju,jc)*zb(jd,je)*zb(jn,jg) + za(jd,jg)*za(jd,jn)*za(je,
     &    jg)*zb(ju,jg)*zb(jd,je)*zb(je,jc) - 1.D0/2.D0*za(jd,jg)*za(jd
     &    ,jn)*za(je,jg)*zb(ju,jc)*zb(jd,je)*zb(je,jg) + za(jd,jn)*za(
     &    jn,jg)**2*zb(ju,jg)*zb(jn,jc)*zb(jn,je) - 1.D0/2.D0*za(jd,jn)
     &    *za(jn,jg)**2*zb(ju,jc)*zb(jn,jg)*zb(jn,je) + za(jd,jn)*za(jn
     &    ,jg)*za(je,jg)*zb(ju,jg)*zb(jn,je)*zb(je,jc) )
      qcda(1,1,1) = qcda(1,1,1) + iza(jc,jg)*izb(jb,jg)*izb(jc,jg)*
     & sned**(-1)*sneud**(-1)*sbc**(-1)*fourrt2*m**2 * (  - 1.D0/2.D0*
     &    za(jd,jn)*za(jn,jg)*za(je,jg)*zb(ju,jc)*zb(jn,je)*zb(je,jg) )
      qcda(1,1,1) = qcda(1,1,1) + iza(jc,jg)*izb(jb,jg)*izb(jc,jg)*
     & sneu**(-1)*sneud**(-1)*sbc**(-1)*fourrt2*m**2 * (  - 1.D0/2.D0*
     &    za(ju,jg)*za(ju,jn)*za(jd,jg)*zb(ju,jg)*zb(ju,jc)*zb(ju,je)
     &     - 1.D0/2.D0*za(ju,jg)*za(jd,jg)*za(jn,je)*zb(ju,jg)*zb(ju,je
     &    )*zb(je,jc) + za(ju,jg)*za(jd,jg)*za(jn,je)*zb(ju,jc)*zb(ju,
     &    je)*zb(je,jg) - za(ju,jn)*za(jd,jg)**2*zb(ju,jg)*zb(ju,je)*
     &    zb(jd,jc) + 1.D0/2.D0*za(ju,jn)*za(jd,jg)**2*zb(ju,jc)*zb(ju,
     &    je)*zb(jd,jg) - za(ju,jn)*za(jd,jg)*za(jn,jg)*zb(ju,jg)*zb(ju
     &    ,je)*zb(jn,jc) + 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(jn,jg)*zb(
     &    ju,jc)*zb(ju,je)*zb(jn,jg) - za(ju,jn)*za(jd,jg)*za(je,jg)*
     &    zb(ju,jg)*zb(ju,je)*zb(je,jc) + 1.D0/2.D0*za(ju,jn)*za(jd,jg)
     &    *za(je,jg)*zb(ju,jc)*zb(ju,je)*zb(je,jg) - 1.D0/2.D0*za(jd,jg
     &    )**2*za(jn,je)*zb(ju,je)*zb(jd,jg)*zb(je,jc) + za(jd,jg)**2*
     &    za(jn,je)*zb(ju,je)*zb(jd,jc)*zb(je,jg) - 1.D0/2.D0*za(jd,jg)
     &    *za(jn,jg)*za(jn,je)*zb(ju,je)*zb(jn,jg)*zb(je,jc) + za(jd,jg
     &    )*za(jn,jg)*za(jn,je)*zb(ju,je)*zb(jn,jc)*zb(je,jg) )
      qcda(1,1,1) = qcda(1,1,1) + iza(jc,jg)*izb(jb,jg)*izb(jc,jg)*
     & sneu**(-1)*sneud**(-1)*sbc**(-1)*fourrt2*m**2 * ( 1.D0/2.D0*za(
     &    jd,jg)*za(jn,je)*za(je,jg)*zb(ju,je)*zb(je,jg)*zb(je,jc) )
      qcda(1,1,1) = qcda(1,1,1) + izb(jd,jg)*sneu**(-1)*sbc**(-1)*
     & fourrt2 * ( za(ju,jn)*za(jb,jg)*zb(ju,jc)*zb(ju,je) - za(jb,jg)*
     &    za(jn,je)*zb(ju,je)*zb(je,jc) )
      qcda(1,1,1) = qcda(1,1,1) + izb(jd,jg)*snedg**(-1)*sbc**(-1)*
     & fourrt2 * (  - za(jd,jb)*za(jn,jg)*zb(ju,jc)*zb(jd,je) - za(jb,
     &    jg)*za(jn,jg)*zb(ju,jc)*zb(je,jg) - za(jn,jg)*za(jn,jb)*zb(ju
     &    ,jc)*zb(jn,je) )
      qcda(1,1,1) = qcda(1,1,1) + izb(jd,jg)*izb(jc,jg)*sneu**(-1)*
     & sbc**(-1)*fourrt2 * ( za(ju,jn)*za(jd,jb)*zb(ju,jc)*zb(ju,je)*
     &    zb(jd,jc) - za(jd,jb)*za(jn,je)*zb(ju,je)*zb(jd,jc)*zb(je,jc)
     &     )
      qcda(1,1,1) = qcda(1,1,1) + izb(jd,jg)*izb(jc,jg)*snedg**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(jd,jb)*za(jd,jn)*zb(ju,jc)*zb(jd,jc)
     &    *zb(jd,je) - za(jd,jn)*za(jb,jg)*zb(ju,jc)*zb(jd,jc)*zb(je,jg
     &    ) - za(jd,jn)*za(jn,jb)*zb(ju,jc)*zb(jd,jc)*zb(jn,je) )
      qcda(1,1,1) = qcda(1,1,1) + izb(jc,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jg)*za(jd,jb)*za(jd,jn)*zb(ju,jc)*
     &    zb(jd,jc)*zb(jd,je) + za(jd,jg)*za(jd,jn)*za(jn,jb)*zb(ju,jc)
     &    *zb(jd,je)*zb(jn,jc) + za(jd,jg)*za(jd,jn)*za(je,jb)*zb(ju,jc
     &    )*zb(jd,je)*zb(je,jc) + za(jd,jb)*za(jd,jn)*za(jn,jg)*zb(ju,
     &    jc)*zb(jd,jc)*zb(jn,je) + za(jd,jn)*za(jn,jg)*za(jn,jb)*zb(ju
     &    ,jc)*zb(jn,jc)*zb(jn,je) + za(jd,jn)*za(jn,jg)*za(je,jb)*zb(
     &    ju,jc)*zb(jn,je)*zb(je,jc) )
      qcda(1,1,1) = qcda(1,1,1) + izb(jc,jg)*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * ( za(ju,jg)*za(jd,jb)*za(jd,jn)*zb(ju,jc)**2
     &    *zb(jd,je) + za(ju,jg)*za(jd,jn)*za(jn,jb)*zb(ju,jc)**2*zb(jn
     &    ,je) - 1.D0/2.D0*za(ju,jb)*za(jd,jg)*za(jd,jn)*zb(ju,jc)**2*
     &    zb(jd,je) - 1.D0/2.D0*za(ju,jb)*za(jd,jn)*za(jn,jg)*zb(ju,jc)
     &    **2*zb(jn,je) + 1.D0/2.D0*za(jd,jg)*za(jd,jb)*za(jd,jn)*zb(ju
     &    ,jc)*zb(jd,jc)*zb(jd,je) + za(jd,jg)*za(jd,jn)*za(jn,jb)*zb(
     &    ju,jc)*zb(jd,jc)*zb(jn,je) - 1.D0/2.D0*za(jd,jg)*za(jd,jn)*
     &    za(jn,jb)*zb(ju,jc)*zb(jd,je)*zb(jn,jc) - 1.D0/2.D0*za(jd,jg)
     &    *za(jd,jn)*za(je,jb)*zb(ju,jc)*zb(jd,je)*zb(je,jc) - 1.D0/2.D0
     &    *za(jd,jb)*za(jd,jn)*za(jn,jg)*zb(ju,jc)*zb(jd,jc)*zb(jn,je)
     &     + za(jd,jb)*za(jd,jn)*za(jn,jg)*zb(ju,jc)*zb(jd,je)*zb(jn,jc
     &    ) + za(jd,jb)*za(jd,jn)*za(je,jg)*zb(ju,jc)*zb(jd,je)*zb(je,
     &    jc) + 1.D0/2.D0*za(jd,jn)*za(jn,jg)*za(jn,jb)*zb(ju,jc)*zb(jn
     &    ,jc)*zb(jn,je) - 1.D0/2.D0*za(jd,jn)*za(jn,jg)*za(je,jb)*zb(
     &    ju,jc)*zb(jn,je)*zb(je,jc) )
      qcda(1,1,1) = qcda(1,1,1) + izb(jc,jg)*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jn)*za(jn,jb)*za(je,jg)*zb(ju,jc)*
     &    zb(jn,je)*zb(je,jc) )
      qcda(1,1,1) = qcda(1,1,1) + izb(jc,jg)*sneu**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(ju,jg)*za(ju,jn)*za(jd,jb)*zb(ju,jc)
     &    **2*zb(ju,je) + za(ju,jg)*za(jd,jb)*za(jn,je)*zb(ju,jc)*zb(ju
     &    ,je)*zb(je,jc) + 1.D0/2.D0*za(ju,jb)*za(ju,jn)*za(jd,jg)*zb(
     &    ju,jc)**2*zb(ju,je) - 1.D0/2.D0*za(ju,jb)*za(jd,jg)*za(jn,je)
     &    *zb(ju,jc)*zb(ju,je)*zb(je,jc) - 1.D0/2.D0*za(ju,jn)*za(jd,jg
     &    )*za(jd,jb)*zb(ju,jc)*zb(ju,je)*zb(jd,jc) + 1.D0/2.D0*za(ju,
     &    jn)*za(jd,jg)*za(jn,jb)*zb(ju,jc)*zb(ju,je)*zb(jn,jc) + 1.D0/
     &    2.D0*za(ju,jn)*za(jd,jg)*za(je,jb)*zb(ju,jc)*zb(ju,je)*zb(je,
     &    jc) - za(ju,jn)*za(jd,jb)*za(jn,jg)*zb(ju,jc)*zb(ju,je)*zb(jn
     &    ,jc) - za(ju,jn)*za(jd,jb)*za(je,jg)*zb(ju,jc)*zb(ju,je)*zb(
     &    je,jc) + 1.D0/2.D0*za(jd,jg)*za(jd,jb)*za(jn,je)*zb(ju,je)*
     &    zb(jd,jc)*zb(je,jc) - 1.D0/2.D0*za(jd,jg)*za(jn,jb)*za(jn,je)
     &    *zb(ju,je)*zb(jn,jc)*zb(je,jc) - 1.D0/2.D0*za(jd,jg)*za(jn,je
     &    )*za(je,jb)*zb(ju,je)*zb(je,jc)**2 + za(jd,jb)*za(jn,jg)*za(
     &    jn,je)*zb(ju,je)*zb(jn,jc)*zb(je,jc) )
      qcda(1,1,1) = qcda(1,1,1) + izb(jc,jg)*sneu**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jb)*za(jn,je)*za(je,jg)*zb(ju,je)*
     &    zb(je,jc)**2 )

      qcda(2,1,1)= + sned**(-1)*sneud**(-1)*sbc**(-1)*fourrt2 * (  - 
     &    za(jd,jb)*za(jd,jn)*zb(ju,jg)*zb(jd,je)*zb(jc,jg) )
      qcda(2,1,1) = qcda(2,1,1) + sneu**(-1)*sneud**(-1)*sbc**(-1)*
     & fourrt2 * ( za(ju,jn)*za(jd,jb)*zb(ju,jg)*zb(ju,je)*zb(jc,jg) - 
     &    za(jd,jb)*za(jn,je)*zb(ju,je)*zb(jc,jg)*zb(je,jg) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jd,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jb)*za(jd,jn)**2*zb(ju,jc)*zb(jd,jg)
     &    *zb(jn,je) + za(jd,jn)**2*za(jn,jb)*zb(ju,jc)*zb(jn,jg)*zb(jn
     &    ,je) + za(jd,jn)**2*za(je,jb)*zb(ju,jc)*zb(jn,je)*zb(je,jg) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jd,jg)*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(ju,jd)*za(jd,jb)*za(jd,jn)*zb(ju,jg)
     &    *zb(ju,jc)*zb(jd,je) - za(ju,jd)*za(jd,jn)*za(jn,jb)*zb(ju,jg
     &    )*zb(ju,jc)*zb(jn,je) - 1.D0/2.D0*za(ju,jb)*za(jd,jn)**2*zb(
     &    ju,jg)*zb(ju,jc)*zb(jn,je) + 1.D0/2.D0*za(jd,jb)*za(jd,jn)**2
     &    *zb(ju,jd)*zb(jc,jg)*zb(jn,je) - 1.D0/2.D0*za(jd,jb)*za(jd,jn
     &    )**2*zb(ju,jg)*zb(jd,jc)*zb(jn,je) + za(jd,jb)*za(jd,jn)**2*
     &    zb(ju,jc)*zb(jd,je)*zb(jn,jg) - 1.D0/2.D0*za(jd,jb)*za(jd,jn)
     &    **2*zb(ju,jn)*zb(jd,je)*zb(jc,jg) + za(jd,jb)*za(jd,jn)*za(jd
     &    ,je)*zb(ju,jc)*zb(jd,je)*zb(je,jg) - 1.D0/2.D0*za(jd,jb)*za(
     &    jd,jn)*za(jd,je)*zb(ju,je)*zb(jd,je)*zb(jc,jg) - za(jd,jb)*
     &    za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(jc,jg)*zb(jn,je) - 1.D0/2.D0
     &    *za(jd,jb)*za(jd,jn)*za(jn,je)*zb(ju,je)*zb(jc,jg)*zb(jn,je)
     &     + 1.D0/2.D0*za(jd,jn)**2*za(jb,jg)*zb(ju,jg)*zb(jc,jg)*zb(jn
     &    ,je) - 1.D0/2.D0*za(jd,jn)**2*za(jn,jb)*zb(ju,jg)*zb(jn,jc)*
     &    zb(jn,je) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jd,jg)*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jn)**2*za(jn,jb)*zb(ju,jc)*zb(jn,jg)
     &    *zb(jn,je) - 1.D0/2.D0*za(jd,jn)**2*za(je,jb)*zb(ju,jg)*zb(jn
     &    ,je)*zb(je,jc) + za(jd,jn)*za(jd,je)*za(jn,jb)*zb(ju,jc)*zb(
     &    jn,je)*zb(je,jg) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jd,jg)*sneu**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2 * ( za(ju,jd)*za(ju,jn)*za(jd,jb)*zb(ju,jg)*
     &    zb(ju,jc)*zb(ju,je) - za(ju,jd)*za(jd,jb)*za(jn,je)*zb(ju,jg)
     &    *zb(ju,je)*zb(je,jc) - 1.D0/2.D0*za(ju,jd)*za(jd,jb)*za(jn,je
     &    )*zb(ju,je)**2*zb(jc,jg) - za(ju,jn)*za(jd,jb)*za(jd,jn)*zb(
     &    ju,jc)*zb(ju,je)*zb(jn,jg) + 1.D0/2.D0*za(ju,jn)*za(jd,jb)*
     &    za(jd,jn)*zb(ju,jn)*zb(ju,je)*zb(jc,jg) - za(ju,jn)*za(jd,jb)
     &    *za(jd,je)*zb(ju,jc)*zb(ju,je)*zb(je,jg) + 1.D0/2.D0*za(ju,jn
     &    )*za(jd,jb)*za(jd,je)*zb(ju,je)**2*zb(jc,jg) + 1.D0/2.D0*za(
     &    jd,jb)*za(jd,jn)*za(jn,je)*zb(ju,je)*zb(jc,jg)*zb(jn,je) + 
     &    za(jd,jb)*za(jd,jn)*za(jn,je)*zb(ju,je)*zb(jn,jg)*zb(je,jc)
     &     + za(jd,jb)*za(jd,je)*za(jn,je)*zb(ju,je)*zb(je,jg)*zb(je,jc
     &    ) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jd,jg)*iza(jc,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2 * ( za(jd,jb)*za(jd,jc)*za(jd,jn)*zb(ju,jc)*
     &    zb(jd,je) + za(jd,jc)*za(jd,jn)*za(jn,jb)*zb(ju,jc)*zb(jn,je)
     &     )
      qcda(2,1,1) = qcda(2,1,1) + iza(jd,jg)*iza(jc,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2 * (  - za(ju,jn)*za(jd,jb)*za(jd,jc)*zb(ju,
     &    jc)*zb(ju,je) + za(jd,jb)*za(jd,jc)*za(jn,je)*zb(ju,je)*zb(je
     &    ,jc) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jd,jg)*iza(jc,jg)*izb(jb,jg)*
     & sned**(-1)*snedg**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jn)**2*
     &    za(jn,jg)*zb(ju,jg)*zb(jn,jg)*zb(jn,je) + za(jd,jn)**2*za(je,
     &    jg)*zb(ju,jg)*zb(jn,je)*zb(je,jg) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jd,jg)*iza(jc,jg)*izb(jb,jg)*
     & sned**(-1)*sneud**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(ju,jd)*
     &    za(jd,jn)*za(jn,jg)*zb(ju,jg)**2*zb(jn,je) - 1.D0/2.D0*za(ju,
     &    jg)*za(jd,jn)**2*zb(ju,jg)**2*zb(jn,je) + 1.D0/2.D0*za(jd,jn)
     &    **2*za(jn,jg)*zb(ju,jg)*zb(jn,jg)*zb(jn,je) - 1.D0/2.D0*za(jd
     &    ,jn)**2*za(je,jg)*zb(ju,jg)*zb(jn,je)*zb(je,jg) + za(jd,jn)*
     &    za(jd,je)*za(jn,jg)*zb(ju,jg)*zb(jn,je)*zb(je,jg) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jd,jg)*iza(jc,jg)*izb(jb,jg)*
     & sned**(-1)*sneud**(-1)*fourrt2*m**2 * ( za(jd,jn)**2*zb(ju,jg)*
     &    zb(jn,je) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jc,jg)*sned**(-1)*sneud**(-1)*
     & fourrt2 * ( za(jd,jb)*za(jd,jn)*zb(ju,jg)*zb(jd,je) + za(jd,jn)*
     &    za(jn,jb)*zb(ju,jg)*zb(jn,je) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jc,jg)*sneu**(-1)*sneud**(-1)*
     & fourrt2 * (  - za(ju,jn)*za(jd,jb)*zb(ju,jg)*zb(ju,je) + za(jd,
     &    jb)*za(jn,je)*zb(ju,je)*zb(je,jg) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jc,jg)**2*izb(jb,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2*m**2 * ( za(jd,jc)*za(jd,jn)*zb(ju,jg)*zb(jd
     &    ,je) + za(jd,jn)*za(jn,jc)*zb(ju,jg)*zb(jn,je) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jc,jg)**2*izb(jb,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2*m**2 * (  - za(ju,jn)*za(jd,jc)*zb(ju,jg)*
     &    zb(ju,je) + za(jd,jc)*za(jn,je)*zb(ju,je)*zb(je,jg) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jc,jg)*izb(jb,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jn)**2*zb(ju,jg)*
     &    zb(jd,jg)*zb(jn,je) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jc,jg)*izb(jb,jg)*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(ju,jd)*za(jd,jn)*
     &    zb(ju,jg)**2*zb(jd,je) - 1.D0/2.D0*za(jd,jn)**2*zb(ju,jg)*zb(
     &    jd,jg)*zb(jn,je) + za(jd,jn)**2*zb(ju,jg)*zb(jd,je)*zb(jn,jg)
     &     + za(jd,jn)*za(jd,je)*zb(ju,jg)*zb(jd,je)*zb(je,jg) )
      qcda(2,1,1) = qcda(2,1,1) + iza(jc,jg)*izb(jb,jg)*sneu**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(ju,jd)*za(ju,jn)*zb(ju
     &    ,jg)**2*zb(ju,je) - za(ju,jd)*za(jn,je)*zb(ju,jg)*zb(ju,je)*
     &    zb(je,jg) - za(ju,jn)*za(jd,jn)*zb(ju,jg)*zb(ju,je)*zb(jn,jg)
     &     - za(ju,jn)*za(jd,je)*zb(ju,jg)*zb(ju,je)*zb(je,jg) + za(jd,
     &    jn)*za(jn,je)*zb(ju,je)*zb(jn,jg)*zb(je,jg) + za(jd,je)*za(jn
     &    ,je)*zb(ju,je)*zb(je,jg)**2 )

      qcda(1,2,1)= + iza(jb,jg)*izb(jd,jg)*snedg**(-1)*sbc**(-1)*
     & fourrt2*m * ( za(jd,jg)*za(jn,jg)*zb(ju,jc)*zb(jd,je) + za(jn,jg
     &    )**2*zb(ju,jc)*zb(jn,je) )
      qcda(1,2,1) = qcda(1,2,1) + iza(jb,jg)*izb(jd,jg)*izb(jc,jg)*
     & sneu**(-1)*sbc**(-1)*fourrt2*m * (  - za(ju,jn)*za(jd,jg)*zb(ju,
     &    jc)*zb(ju,je)*zb(jd,jc) + za(jd,jg)*za(jn,je)*zb(ju,je)*zb(jd
     &    ,jc)*zb(je,jc) )
      qcda(1,2,1) = qcda(1,2,1) + iza(jb,jg)*izb(jd,jg)*izb(jc,jg)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * ( za(jd,jg)*za(jd,jn)*zb(ju,jc
     &    )*zb(jd,jc)*zb(jd,je) + za(jd,jn)*za(jn,jg)*zb(ju,jc)*zb(jd,
     &    jc)*zb(jn,je) )
      qcda(1,2,1) = qcda(1,2,1) + iza(jb,jg)*izb(jc,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jg)**2*za(jd,jn)*
     &    zb(ju,jc)*zb(jd,jc)*zb(jd,je) - za(jd,jg)*za(jd,jn)*za(jn,jg)
     &    *zb(ju,jc)*zb(jd,jc)*zb(jn,je) - za(jd,jg)*za(jd,jn)*za(jn,jg
     &    )*zb(ju,jc)*zb(jd,je)*zb(jn,jc) - za(jd,jg)*za(jd,jn)*za(je,
     &    jg)*zb(ju,jc)*zb(jd,je)*zb(je,jc) - za(jd,jn)*za(jn,jg)**2*
     &    zb(ju,jc)*zb(jn,jc)*zb(jn,je) - za(jd,jn)*za(jn,jg)*za(je,jg)
     &    *zb(ju,jc)*zb(jn,je)*zb(je,jc) )
      qcda(1,2,1) = qcda(1,2,1) + iza(jb,jg)*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * (  - 1.D0/2.D0*za(ju,jg)*za(jd
     &    ,jg)*za(jd,jn)*zb(ju,jc)**2*zb(jd,je) - 1.D0/2.D0*za(ju,jg)*
     &    za(jd,jn)*za(jn,jg)*zb(ju,jc)**2*zb(jn,je) - 1.D0/2.D0*za(jd,
     &    jg)**2*za(jd,jn)*zb(ju,jc)*zb(jd,jc)*zb(jd,je) - 1.D0/2.D0*
     &    za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jc)*zb(jd,jc)*zb(jn,je)
     &     - 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jc)*zb(jd,je
     &    )*zb(jn,jc) - 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(je,jg)*zb(ju,
     &    jc)*zb(jd,je)*zb(je,jc) - 1.D0/2.D0*za(jd,jn)*za(jn,jg)**2*
     &    zb(ju,jc)*zb(jn,jc)*zb(jn,je) - 1.D0/2.D0*za(jd,jn)*za(jn,jg)
     &    *za(je,jg)*zb(ju,jc)*zb(jn,je)*zb(je,jc) )
      qcda(1,2,1) = qcda(1,2,1) + iza(jb,jg)*izb(jc,jg)*sneu**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * ( 1.D0/2.D0*za(ju,jg)*za(ju,jn
     &    )*za(jd,jg)*zb(ju,jc)**2*zb(ju,je) - 1.D0/2.D0*za(ju,jg)*za(
     &    jd,jg)*za(jn,je)*zb(ju,jc)*zb(ju,je)*zb(je,jc) + 1.D0/2.D0*
     &    za(ju,jn)*za(jd,jg)**2*zb(ju,jc)*zb(ju,je)*zb(jd,jc) + 1.D0/2.
     &    D0*za(ju,jn)*za(jd,jg)*za(jn,jg)*zb(ju,jc)*zb(ju,je)*zb(jn,jc
     &    ) + 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(je,jg)*zb(ju,jc)*zb(ju,
     &    je)*zb(je,jc) - 1.D0/2.D0*za(jd,jg)**2*za(jn,je)*zb(ju,je)*
     &    zb(jd,jc)*zb(je,jc) - 1.D0/2.D0*za(jd,jg)*za(jn,jg)*za(jn,je)
     &    *zb(ju,je)*zb(jn,jc)*zb(je,jc) - 1.D0/2.D0*za(jd,jg)*za(jn,je
     &    )*za(je,jg)*zb(ju,je)*zb(je,jc)**2 )
      qcda(1,2,1) = qcda(1,2,1) + iza(jc,jg)*izb(jd,jg)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(jd,jg)*za(jn,jg)*zb(ju,jb)*zb(jd,je)
     &     + za(jn,jg)**2*zb(ju,jb)*zb(jn,je) )
      qcda(1,2,1) = qcda(1,2,1) + iza(jc,jg)*izb(jd,jg)*izb(jc,jg)*
     & sneu**(-1)*sbc**(-1)*fourrt2*m * (  - za(ju,jn)*za(jd,jg)*zb(ju,
     &    jb)*zb(ju,je)*zb(jd,jc) + za(jd,jg)*za(jn,je)*zb(ju,je)*zb(jd
     &    ,jc)*zb(je,jb) )
      qcda(1,2,1) = qcda(1,2,1) + iza(jc,jg)*izb(jd,jg)*izb(jc,jg)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * ( za(jd,jg)*za(jd,jn)*zb(ju,jb
     &    )*zb(jd,jc)*zb(jd,je) + za(jd,jn)*za(jn,jg)*zb(ju,jb)*zb(jd,
     &    jc)*zb(jn,je) )
      qcda(1,2,1) = qcda(1,2,1) + iza(jc,jg)*izb(jc,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jg)**2*za(jd,jn)*
     &    zb(ju,jb)*zb(jd,jc)*zb(jd,je) - za(jd,jg)*za(jd,jn)*za(jn,jg)
     &    *zb(ju,jb)*zb(jd,jc)*zb(jn,je) - za(jd,jg)*za(jd,jn)*za(jn,jg
     &    )*zb(ju,jb)*zb(jd,je)*zb(jn,jc) - za(jd,jg)*za(jd,jn)*za(je,
     &    jg)*zb(ju,jb)*zb(jd,je)*zb(je,jc) - za(jd,jn)*za(jn,jg)**2*
     &    zb(ju,jb)*zb(jn,jc)*zb(jn,je) - za(jd,jn)*za(jn,jg)*za(je,jg)
     &    *zb(ju,jb)*zb(jn,je)*zb(je,jc) )
      qcda(1,2,1) = qcda(1,2,1) + iza(jc,jg)*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * (  - 1.D0/2.D0*za(ju,jg)*za(jd
     &    ,jg)*za(jd,jn)*zb(ju,jb)*zb(ju,jc)*zb(jd,je) - 1.D0/2.D0*za(
     &    ju,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jb)*zb(ju,jc)*zb(jn,je) - 
     &    za(jd,jg)**2*za(jd,jn)*zb(ju,jb)*zb(jd,jc)*zb(jd,je) + 1.D0/2.
     &    D0*za(jd,jg)**2*za(jd,jn)*zb(ju,jc)*zb(jd,jb)*zb(jd,je) - za(
     &    jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jb)*zb(jd,jc)*zb(jn,je) - 
     &    za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jb)*zb(jd,je)*zb(jn,jc)
     &     + 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jc)*zb(jd,jb
     &    )*zb(jn,je) + 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,
     &    jc)*zb(jd,je)*zb(jn,jb) - za(jd,jg)*za(jd,jn)*za(je,jg)*zb(ju
     &    ,jb)*zb(jd,je)*zb(je,jc) + 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(
     &    je,jg)*zb(ju,jc)*zb(jd,je)*zb(je,jb) - za(jd,jn)*za(jn,jg)**2
     &    *zb(ju,jb)*zb(jn,jc)*zb(jn,je) + 1.D0/2.D0*za(jd,jn)*za(jn,jg
     &    )**2*zb(ju,jc)*zb(jn,jb)*zb(jn,je) - za(jd,jn)*za(jn,jg)*za(
     &    je,jg)*zb(ju,jb)*zb(jn,je)*zb(je,jc) )
      qcda(1,2,1) = qcda(1,2,1) + iza(jc,jg)*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * ( 1.D0/2.D0*za(jd,jn)*za(jn,jg
     &    )*za(je,jg)*zb(ju,jc)*zb(jn,je)*zb(je,jb) )
      qcda(1,2,1) = qcda(1,2,1) + iza(jc,jg)*izb(jc,jg)*sneu**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * ( 1.D0/2.D0*za(ju,jg)*za(ju,jn
     &    )*za(jd,jg)*zb(ju,jb)*zb(ju,jc)*zb(ju,je) + 1.D0/2.D0*za(ju,
     &    jg)*za(jd,jg)*za(jn,je)*zb(ju,jb)*zb(ju,je)*zb(je,jc) - za(ju
     &    ,jg)*za(jd,jg)*za(jn,je)*zb(ju,jc)*zb(ju,je)*zb(je,jb) + za(
     &    ju,jn)*za(jd,jg)**2*zb(ju,jb)*zb(ju,je)*zb(jd,jc) - 1.D0/2.D0
     &    *za(ju,jn)*za(jd,jg)**2*zb(ju,jc)*zb(ju,je)*zb(jd,jb) + za(ju
     &    ,jn)*za(jd,jg)*za(jn,jg)*zb(ju,jb)*zb(ju,je)*zb(jn,jc) - 1.D0/
     &    2.D0*za(ju,jn)*za(jd,jg)*za(jn,jg)*zb(ju,jc)*zb(ju,je)*zb(jn,
     &    jb) + za(ju,jn)*za(jd,jg)*za(je,jg)*zb(ju,jb)*zb(ju,je)*zb(je
     &    ,jc) - 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(je,jg)*zb(ju,jc)*zb(
     &    ju,je)*zb(je,jb) + 1.D0/2.D0*za(jd,jg)**2*za(jn,je)*zb(ju,je)
     &    *zb(jd,jb)*zb(je,jc) - za(jd,jg)**2*za(jn,je)*zb(ju,je)*zb(jd
     &    ,jc)*zb(je,jb) + 1.D0/2.D0*za(jd,jg)*za(jn,jg)*za(jn,je)*zb(
     &    ju,je)*zb(jn,jb)*zb(je,jc) - za(jd,jg)*za(jn,jg)*za(jn,je)*
     &    zb(ju,je)*zb(jn,jc)*zb(je,jb) )
      qcda(1,2,1) = qcda(1,2,1) + iza(jc,jg)*izb(jc,jg)*sneu**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * (  - 1.D0/2.D0*za(jd,jg)*za(jn
     &    ,je)*za(je,jg)*zb(ju,je)*zb(je,jb)*zb(je,jc) )

      qcda(1,1,2)= + izb(jd,jg)*izb(jb,jg)*sneu**(-1)*sbc**(-1)*fourrt2
     & *m * (  - za(ju,jn)*za(jc,jg)*zb(ju,jg)*zb(ju,je) + za(jc,jg)*
     &    za(jn,je)*zb(ju,je)*zb(je,jg) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jd,jg)*izb(jb,jg)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(jd,jc)*za(jn,jg)*zb(ju,jg)*zb(jd,je)
     &     + za(jc,jg)*za(jn,jg)*zb(ju,jg)*zb(je,jg) + za(jn,jg)*za(jn,
     &    jc)*zb(ju,jg)*zb(jn,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jd,jg)*izb(jb,jg)*izb(jc,jg)*
     & sneu**(-1)*sbc**(-1)*fourrt2*m * (  - za(ju,jn)*za(jd,jc)*zb(ju,
     &    jg)*zb(ju,je)*zb(jd,jc) + za(jd,jc)*za(jn,je)*zb(ju,je)*zb(jd
     &    ,jc)*zb(je,jg) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jd,jg)*izb(jb,jg)*izb(jc,jg)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * ( za(jd,jc)*za(jd,jn)*zb(ju,jg
     &    )*zb(jd,jc)*zb(jd,je) + za(jd,jn)*za(jc,jg)*zb(ju,jg)*zb(jd,
     &    jc)*zb(je,jg) + za(jd,jn)*za(jn,jc)*zb(ju,jg)*zb(jd,jc)*zb(jn
     &    ,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jd,jg)*izb(jc,jg)*sneu**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(ju,jn)*za(jb,jg)*zb(ju,jg)*zb(ju,
     &    je) + za(jb,jg)*za(jn,je)*zb(ju,je)*zb(je,jg) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jd,jg)*izb(jc,jg)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(jd,jb)*za(jn,jg)*zb(ju,jg)*zb(jd,je)
     &     + za(jb,jg)*za(jn,jg)*zb(ju,jg)*zb(je,jg) + za(jn,jg)*za(jn,
     &    jb)*zb(ju,jg)*zb(jn,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jd,jg)*izb(jc,jg)**2*sneu**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(ju,jn)*za(jd,jb)*zb(ju,jg)*zb(ju,
     &    je)*zb(jd,jc) + za(jd,jb)*za(jn,je)*zb(ju,je)*zb(jd,jc)*zb(je
     &    ,jg) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jd,jg)*izb(jc,jg)**2*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(jd,jb)*za(jd,jn)*zb(ju,jg)*zb(jd,jc)*
     &    zb(jd,je) + za(jd,jn)*za(jb,jg)*zb(ju,jg)*zb(jd,jc)*zb(je,jg)
     &     + za(jd,jn)*za(jn,jb)*zb(ju,jg)*zb(jd,jc)*zb(jn,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jb,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jg)*za(jd,jn)*za(jc,jg)*zb(ju,
     &    jg)*zb(jd,je) - za(jd,jn)*za(jc,jg)*za(jn,jg)*zb(ju,jg)*zb(jn
     &    ,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jb,jg)*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jg)*za(jd,jn)*za(jc,jg)*zb(ju,
     &    jg)*zb(jd,je) + 1.D0/2.D0*za(jd,jn)**2*za(jc,jg)*zb(ju,jd)*
     &    zb(jn,je) - 1.D0/2.D0*za(jd,jn)**2*za(jc,jg)*zb(ju,jn)*zb(jd,
     &    je) - 1.D0/2.D0*za(jd,jn)*za(jd,je)*za(jc,jg)*zb(ju,je)*zb(jd
     &    ,je) - za(jd,jn)*za(jc,jg)*za(jn,jg)*zb(ju,jg)*zb(jn,je) - 1.D
     &    0/2.D0*za(jd,jn)*za(jc,jg)*za(jn,je)*zb(ju,je)*zb(jn,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jb,jg)*sneu**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2*m * (  - 1.D0/2.D0*za(ju,jd)*za(jc,jg)*za(jn,
     &    je)*zb(ju,je)**2 + za(ju,jn)*za(jd,jg)*za(jc,jg)*zb(ju,jg)*
     &    zb(ju,je) + 1.D0/2.D0*za(ju,jn)*za(jd,jn)*za(jc,jg)*zb(ju,jn)
     &    *zb(ju,je) + 1.D0/2.D0*za(ju,jn)*za(jd,je)*za(jc,jg)*zb(ju,je
     &    )**2 - za(jd,jg)*za(jc,jg)*za(jn,je)*zb(ju,je)*zb(je,jg) + 1.D
     &    0/2.D0*za(jd,jn)*za(jc,jg)*za(jn,je)*zb(ju,je)*zb(jn,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jb,jg)*izb(jc,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jg)*za(jd,jc)*za(jd
     &    ,jn)*zb(ju,jg)*zb(jd,jc)*zb(jd,je) - za(jd,jg)*za(jd,jn)*za(
     &    jn,jc)*zb(ju,jg)*zb(jd,je)*zb(jn,jc) - za(jd,jg)*za(jd,jn)*
     &    za(je,jc)*zb(ju,jg)*zb(jd,je)*zb(je,jc) - za(jd,jc)*za(jd,jn)
     &    *za(jn,jg)*zb(ju,jg)*zb(jd,jc)*zb(jn,je) - za(jd,jn)*za(jn,jg
     &    )*za(jn,jc)*zb(ju,jg)*zb(jn,jc)*zb(jn,je) - za(jd,jn)*za(jn,
     &    jg)*za(je,jc)*zb(ju,jg)*zb(jn,je)*zb(je,jc) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jb,jg)*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * (  - za(ju,jg)*za(jd,jc)*za(jd
     &    ,jn)*zb(ju,jg)*zb(ju,jc)*zb(jd,je) - za(ju,jg)*za(jd,jn)*za(
     &    jn,jc)*zb(ju,jg)*zb(ju,jc)*zb(jn,je) + 1.D0/2.D0*za(ju,jc)*
     &    za(jd,jg)*za(jd,jn)*zb(ju,jg)*zb(ju,jc)*zb(jd,je) + 1.D0/2.D0
     &    *za(ju,jc)*za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(ju,jc)*zb(jn,je)
     &     - za(jd,jg)*za(jd,jc)*za(jd,jn)*zb(ju,jg)*zb(jd,jc)*zb(jd,je
     &    ) + 1.D0/2.D0*za(jd,jg)*za(jd,jc)*za(jd,jn)*zb(ju,jc)*zb(jd,
     &    jg)*zb(jd,je) - za(jd,jg)*za(jd,jn)*za(jn,jc)*zb(ju,jg)*zb(jd
     &    ,jc)*zb(jn,je) + 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(jn,jc)*zb(
     &    ju,jc)*zb(jd,je)*zb(jn,jg) + 1.D0/2.D0*za(jd,jg)*za(jd,jn)*
     &    za(je,jc)*zb(ju,jc)*zb(jd,je)*zb(je,jg) - za(jd,jc)*za(jd,jn)
     &    *za(jn,jg)*zb(ju,jg)*zb(jd,je)*zb(jn,jc) + 1.D0/2.D0*za(jd,jc
     &    )*za(jd,jn)*za(jn,jg)*zb(ju,jc)*zb(jd,jg)*zb(jn,je) - za(jd,
     &    jc)*za(jd,jn)*za(je,jg)*zb(ju,jg)*zb(jd,je)*zb(je,jc) - za(jd
     &    ,jn)*za(jn,jg)*za(jn,jc)*zb(ju,jg)*zb(jn,jc)*zb(jn,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jb,jg)*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * ( 1.D0/2.D0*za(jd,jn)*za(jn,jg
     &    )*za(jn,jc)*zb(ju,jc)*zb(jn,jg)*zb(jn,je) + 1.D0/2.D0*za(jd,
     &    jn)*za(jn,jg)*za(je,jc)*zb(ju,jc)*zb(jn,je)*zb(je,jg) - za(jd
     &    ,jn)*za(jn,jc)*za(je,jg)*zb(ju,jg)*zb(jn,je)*zb(je,jc) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jb,jg)*izb(jc,jg)*sned**(-1)*
     & sneud**(-1)*fourrt2*m * ( za(jd,jg)*za(jd,jn)*zb(ju,jg)*zb(jd,je
     &    ) + za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(jn,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jb,jg)*izb(jc,jg)*sneu**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * ( za(ju,jg)*za(ju,jn)*za(jd,jc
     &    )*zb(ju,jg)*zb(ju,jc)*zb(ju,je) - za(ju,jg)*za(jd,jc)*za(jn,
     &    je)*zb(ju,jc)*zb(ju,je)*zb(je,jg) - 1.D0/2.D0*za(ju,jc)*za(ju
     &    ,jn)*za(jd,jg)*zb(ju,jg)*zb(ju,jc)*zb(ju,je) + 1.D0/2.D0*za(
     &    ju,jc)*za(jd,jg)*za(jn,je)*zb(ju,jg)*zb(ju,je)*zb(je,jc) + 
     &    za(ju,jn)*za(jd,jg)*za(jd,jc)*zb(ju,jg)*zb(ju,je)*zb(jd,jc)
     &     - 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(jd,jc)*zb(ju,jc)*zb(ju,je
     &    )*zb(jd,jg) - 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(jn,jc)*zb(ju,
     &    jc)*zb(ju,je)*zb(jn,jg) - 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(je
     &    ,jc)*zb(ju,jc)*zb(ju,je)*zb(je,jg) + za(ju,jn)*za(jd,jc)*za(
     &    jn,jg)*zb(ju,jg)*zb(ju,je)*zb(jn,jc) + za(ju,jn)*za(jd,jc)*
     &    za(je,jg)*zb(ju,jg)*zb(ju,je)*zb(je,jc) + 1.D0/2.D0*za(jd,jg)
     &    *za(jd,jc)*za(jn,je)*zb(ju,je)*zb(jd,jg)*zb(je,jc) - za(jd,jg
     &    )*za(jd,jc)*za(jn,je)*zb(ju,je)*zb(jd,jc)*zb(je,jg) + 1.D0/2.D
     &    0*za(jd,jg)*za(jn,jc)*za(jn,je)*zb(ju,je)*zb(jn,jg)*zb(je,jc)
     &     )
      qcda(1,1,2) = qcda(1,1,2) + izb(jb,jg)*izb(jc,jg)*sneu**(-1)*
     & sneud**(-1)*sbc**(-1)*fourrt2*m * ( 1.D0/2.D0*za(jd,jg)*za(jn,je
     &    )*za(je,jc)*zb(ju,je)*zb(je,jg)*zb(je,jc) - za(jd,jc)*za(jn,
     &    jg)*za(jn,je)*zb(ju,je)*zb(jn,jc)*zb(je,jg) - za(jd,jc)*za(jn
     &    ,je)*za(je,jg)*zb(ju,je)*zb(je,jg)*zb(je,jc) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jb,jg)*izb(jc,jg)*sneu**(-1)*
     & sneud**(-1)*fourrt2*m * (  - za(ju,jn)*za(jd,jg)*zb(ju,jg)*zb(ju
     &    ,je) + za(jd,jg)*za(jn,je)*zb(ju,je)*zb(je,jg) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jc,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jg)*za(jd,jn)*za(jb,jg)*zb(ju,
     &    jg)*zb(jd,je) - za(jd,jn)*za(jb,jg)*za(jn,jg)*zb(ju,jg)*zb(jn
     &    ,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jc,jg)*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jg)*za(jd,jn)*za(jb,jg)*zb(ju,
     &    jg)*zb(jd,je) + 1.D0/2.D0*za(jd,jn)**2*za(jb,jg)*zb(ju,jd)*
     &    zb(jn,je) - 1.D0/2.D0*za(jd,jn)**2*za(jb,jg)*zb(ju,jn)*zb(jd,
     &    je) - 1.D0/2.D0*za(jd,jn)*za(jd,je)*za(jb,jg)*zb(ju,je)*zb(jd
     &    ,je) - za(jd,jn)*za(jb,jg)*za(jn,jg)*zb(ju,jg)*zb(jn,je) - 1.D
     &    0/2.D0*za(jd,jn)*za(jb,jg)*za(jn,je)*zb(ju,je)*zb(jn,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jc,jg)*sneu**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2*m * (  - 1.D0/2.D0*za(ju,jd)*za(jb,jg)*za(jn,
     &    je)*zb(ju,je)**2 + za(ju,jn)*za(jd,jg)*za(jb,jg)*zb(ju,jg)*
     &    zb(ju,je) + 1.D0/2.D0*za(ju,jn)*za(jd,jn)*za(jb,jg)*zb(ju,jn)
     &    *zb(ju,je) + 1.D0/2.D0*za(ju,jn)*za(jd,je)*za(jb,jg)*zb(ju,je
     &    )**2 - za(jd,jg)*za(jb,jg)*za(jn,je)*zb(ju,je)*zb(je,jg) + 1.D
     &    0/2.D0*za(jd,jn)*za(jb,jg)*za(jn,je)*zb(ju,je)*zb(jn,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jc,jg)**2*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jg)*za(jd,jb)*za(jd,jn)*zb(ju,
     &    jg)*zb(jd,jc)*zb(jd,je) - za(jd,jg)*za(jd,jn)*za(jn,jb)*zb(ju
     &    ,jg)*zb(jd,je)*zb(jn,jc) - za(jd,jg)*za(jd,jn)*za(je,jb)*zb(
     &    ju,jg)*zb(jd,je)*zb(je,jc) - za(jd,jb)*za(jd,jn)*za(jn,jg)*
     &    zb(ju,jg)*zb(jd,jc)*zb(jn,je) - za(jd,jn)*za(jn,jg)*za(jn,jb)
     &    *zb(ju,jg)*zb(jn,jc)*zb(jn,je) - za(jd,jn)*za(jn,jg)*za(je,jb
     &    )*zb(ju,jg)*zb(jn,je)*zb(je,jc) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jc,jg)**2*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(ju,jg)*za(jd,jb)*za(jd,jn)*zb(ju,
     &    jg)*zb(ju,jc)*zb(jd,je) - za(ju,jg)*za(jd,jn)*za(jn,jb)*zb(ju
     &    ,jg)*zb(ju,jc)*zb(jn,je) + 1.D0/2.D0*za(ju,jb)*za(jd,jg)*za(
     &    jd,jn)*zb(ju,jg)*zb(ju,jc)*zb(jd,je) + 1.D0/2.D0*za(ju,jb)*
     &    za(jd,jn)*za(jn,jg)*zb(ju,jg)*zb(ju,jc)*zb(jn,je) - za(jd,jg)
     &    *za(jd,jb)*za(jd,jn)*zb(ju,jg)*zb(jd,jc)*zb(jd,je) + 1.D0/2.D0
     &    *za(jd,jg)*za(jd,jb)*za(jd,jn)*zb(ju,jc)*zb(jd,jg)*zb(jd,je)
     &     - za(jd,jg)*za(jd,jn)*za(jn,jb)*zb(ju,jg)*zb(jd,jc)*zb(jn,je
     &    ) + 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(jn,jb)*zb(ju,jc)*zb(jd,
     &    je)*zb(jn,jg) + 1.D0/2.D0*za(jd,jg)*za(jd,jn)*za(je,jb)*zb(ju
     &    ,jc)*zb(jd,je)*zb(je,jg) - za(jd,jb)*za(jd,jn)*za(jn,jg)*zb(
     &    ju,jg)*zb(jd,je)*zb(jn,jc) + 1.D0/2.D0*za(jd,jb)*za(jd,jn)*
     &    za(jn,jg)*zb(ju,jc)*zb(jd,jg)*zb(jn,je) - za(jd,jb)*za(jd,jn)
     &    *za(je,jg)*zb(ju,jg)*zb(jd,je)*zb(je,jc) - za(jd,jn)*za(jn,jg
     &    )*za(jn,jb)*zb(ju,jg)*zb(jn,jc)*zb(jn,je) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jc,jg)**2*sned**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2*m * ( 1.D0/2.D0*za(jd,jn)*za(jn,jg)*za(jn,jb)*
     &    zb(ju,jc)*zb(jn,jg)*zb(jn,je) + 1.D0/2.D0*za(jd,jn)*za(jn,jg)
     &    *za(je,jb)*zb(ju,jc)*zb(jn,je)*zb(je,jg) - za(jd,jn)*za(jn,jb
     &    )*za(je,jg)*zb(ju,jg)*zb(jn,je)*zb(je,jc) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jc,jg)**2*sneu**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(ju,jg)*za(ju,jn)*za(jd,jb)*zb(ju,jg)*
     &    zb(ju,jc)*zb(ju,je) - za(ju,jg)*za(jd,jb)*za(jn,je)*zb(ju,jc)
     &    *zb(ju,je)*zb(je,jg) - 1.D0/2.D0*za(ju,jb)*za(ju,jn)*za(jd,jg
     &    )*zb(ju,jg)*zb(ju,jc)*zb(ju,je) + 1.D0/2.D0*za(ju,jb)*za(jd,
     &    jg)*za(jn,je)*zb(ju,jg)*zb(ju,je)*zb(je,jc) + za(ju,jn)*za(jd
     &    ,jg)*za(jd,jb)*zb(ju,jg)*zb(ju,je)*zb(jd,jc) - 1.D0/2.D0*za(
     &    ju,jn)*za(jd,jg)*za(jd,jb)*zb(ju,jc)*zb(ju,je)*zb(jd,jg) - 1.D
     &    0/2.D0*za(ju,jn)*za(jd,jg)*za(jn,jb)*zb(ju,jc)*zb(ju,je)*zb(
     &    jn,jg) - 1.D0/2.D0*za(ju,jn)*za(jd,jg)*za(je,jb)*zb(ju,jc)*
     &    zb(ju,je)*zb(je,jg) + za(ju,jn)*za(jd,jb)*za(jn,jg)*zb(ju,jg)
     &    *zb(ju,je)*zb(jn,jc) + za(ju,jn)*za(jd,jb)*za(je,jg)*zb(ju,jg
     &    )*zb(ju,je)*zb(je,jc) + 1.D0/2.D0*za(jd,jg)*za(jd,jb)*za(jn,
     &    je)*zb(ju,je)*zb(jd,jg)*zb(je,jc) - za(jd,jg)*za(jd,jb)*za(jn
     &    ,je)*zb(ju,je)*zb(jd,jc)*zb(je,jg) + 1.D0/2.D0*za(jd,jg)*za(
     &    jn,jb)*za(jn,je)*zb(ju,je)*zb(jn,jg)*zb(je,jc) )
      qcda(1,1,2) = qcda(1,1,2) + izb(jc,jg)**2*sneu**(-1)*sneud**(-1)*
     & sbc**(-1)*fourrt2*m * ( 1.D0/2.D0*za(jd,jg)*za(jn,je)*za(je,jb)*
     &    zb(ju,je)*zb(je,jg)*zb(je,jc) - za(jd,jb)*za(jn,jg)*za(jn,je)
     &    *zb(ju,je)*zb(jn,jc)*zb(je,jg) - za(jd,jb)*za(jn,je)*za(je,jg
     &    )*zb(ju,je)*zb(je,jg)*zb(je,jc) )

      return
      end
