      subroutine wbbgmassi(ju,jd,jg,jb,jc,jn,je,m,qedi)
      implicit none
      include 'constants.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      double precision m,sned,sneu,snedg,sneug,sbc
      double complex qedi(2,2,2)
      double complex iza(7,7),izb(7,7)
      integer i,j,ju,jd,jg,jb,jc,jn,je
C--   order of indices is gluon helicity hg,sb,sc
      sned=s(jn,je)+s(jn,jd)+s(je,jd)
      sneu=s(jn,je)+s(jn,ju)+s(je,ju)
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
      qedi(2,2,2)= + sneu**(-1)*sneug**(-1)*sbc**(-1)*fourrt2 * (  - 
     &    za(ju,jn)*za(jd,jc)*zb(ju,jg)*zb(ju,je)*zb(jb,jg) + za(jd,jc)
     &    *za(jn,je)*zb(ju,je)*zb(jb,jg)*zb(je,jg) )
      qedi(2,2,2) = qedi(2,2,2) + iza(ju,jg)*sned**(-1)*sbc**(-1)*
     & fourrt2 * (  - za(jd,jc)*za(jd,jn)*zb(jd,je)*zb(jb,jg) - za(jd,
     &    jn)*za(jn,jc)*zb(jb,jg)*zb(jn,je) )
      qedi(2,2,2) = qedi(2,2,2) + iza(ju,jg)*sneug**(-1)*sbc**(-1)*
     & fourrt2 * ( za(ju,jn)*za(jd,jc)*zb(ju,jb)*zb(je,jg) + za(jd,jc)*
     &    za(jn,jg)*zb(jb,jg)*zb(je,jg) - za(jd,jc)*za(jn,je)*zb(je,jg)
     &    *zb(je,jb) )
      qedi(2,2,2) = qedi(2,2,2) + iza(ju,jg)*iza(jd,jg)*sned**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(ju,jd)*za(jd,jc)*za(jd,jn)*zb(ju,jb)
     &    *zb(jd,je) - za(ju,jd)*za(jd,jn)*za(jn,jc)*zb(ju,jb)*zb(jn,je
     &    ) )
      qedi(2,2,2) = qedi(2,2,2) + iza(ju,jg)*iza(jd,jg)*sneug**(-1)*
     & sbc**(-1)*fourrt2 * ( za(ju,jd)*za(ju,jn)*za(jd,jc)*zb(ju,jb)*
     &    zb(ju,je) + za(ju,jd)*za(jd,jc)*za(jn,jg)*zb(ju,je)*zb(jb,jg)
     &     - za(ju,jd)*za(jd,jc)*za(jn,je)*zb(ju,je)*zb(je,jb) )
      qedi(2,2,2) = qedi(2,2,2) + iza(ju,jg)*iza(jd,jg)*iza(jb,jg)*izb(
     & jc,jg)*sned**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(ju,jd)*za(jd,
     &    jn)*za(jn,jg)*zb(ju,jg)*zb(jn,je) )
      qedi(2,2,2) = qedi(2,2,2) + iza(ju,jg)*iza(jb,jg)*izb(jc,jg)*
     & sned**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(ju,jd)*za(jd,jn)*zb(
     &    ju,jg)*zb(jd,je) )
      qedi(2,2,2) = qedi(2,2,2) + iza(ju,jg)*iza(jb,jg)*izb(jc,jg)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(ju,jd)*za(ju,jn)*zb(ju
     &    ,jg)*zb(ju,je) - za(ju,jd)*za(jn,je)*zb(ju,je)*zb(je,jg) + 
     &    za(ju,jn)*za(jd,jg)*zb(ju,jg)*zb(je,jg) - za(jd,jg)*za(jn,je)
     &    *zb(je,jg)**2 )
      qedi(2,2,2) = qedi(2,2,2) + iza(jd,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jc)*za(jd,jn)**2*zb(ju,jb)*zb(jd,jg)
     &    *zb(jn,je) + za(jd,jn)**2*za(jn,jc)*zb(ju,jb)*zb(jn,jg)*zb(jn
     &    ,je) + za(jd,jn)**2*za(je,jc)*zb(ju,jb)*zb(jn,je)*zb(je,jg) )
      qedi(2,2,2) = qedi(2,2,2) + iza(jd,jg)*sneu**(-1)*sneug**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(ju,jd)*za(ju,jn)*za(jd,jc)*zb(ju,jg)
     &    *zb(ju,jb)*zb(ju,je) + za(ju,jd)*za(jd,jc)*za(jn,je)*zb(ju,jb
     &    )*zb(ju,je)*zb(je,jg) + za(ju,jn)*za(jd,jc)*za(jd,jn)*zb(ju,
     &    jg)*zb(ju,je)*zb(jn,jb) + za(ju,jn)*za(jd,jc)*za(jd,je)*zb(ju
     &    ,jg)*zb(ju,je)*zb(je,jb) - za(jd,jc)*za(jd,jn)*za(jn,je)*zb(
     &    ju,je)*zb(jn,jb)*zb(je,jg) - za(jd,jc)*za(jd,je)*za(jn,je)*
     &    zb(ju,je)*zb(je,jg)*zb(je,jb) )
      qedi(2,2,2) = qedi(2,2,2) + iza(jd,jg)*iza(jb,jg)*izb(jc,jg)*
     & sned**(-1)*snedg**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jn)**2*
     &    za(jn,jg)*zb(ju,jg)*zb(jn,jg)*zb(jn,je) + za(jd,jn)**2*za(je,
     &    jg)*zb(ju,jg)*zb(jn,je)*zb(je,jg) )
      qedi(2,2,2) = qedi(2,2,2) + iza(jb,jg)*izb(jc,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jn)**2*zb(ju,jg)*
     &    zb(jd,jg)*zb(jn,je) )
      qedi(2,2,2) = qedi(2,2,2) + iza(jb,jg)*izb(jc,jg)*sneu**(-1)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(ju,jd)*za(ju,jn)*
     &    zb(ju,jg)**2*zb(ju,je) + za(ju,jd)*za(jn,je)*zb(ju,jg)*zb(ju,
     &    je)*zb(je,jg) + za(ju,jn)*za(jd,jn)*zb(ju,jg)*zb(ju,je)*zb(jn
     &    ,jg) + za(ju,jn)*za(jd,je)*zb(ju,jg)*zb(ju,je)*zb(je,jg) - 
     &    za(jd,jn)*za(jn,je)*zb(ju,je)*zb(jn,jg)*zb(je,jg) - za(jd,je)
     &    *za(jn,je)*zb(ju,je)*zb(je,jg)**2 )

      qedi(1,2,2)= + sned**(-1)*snedg**(-1)*sbc**(-1)*fourrt2 * ( za(jd
     &    ,jg)*za(jd,jn)*za(jc,jg)*zb(ju,jb)*zb(jd,je) + za(jd,jn)*za(
     &    jc,jg)*za(jn,jg)*zb(ju,jb)*zb(jn,je) )
      qedi(1,2,2) = qedi(1,2,2) + iza(jb,jg)*izb(ju,jg)*izb(jd,jg)*izb(
     & jc,jg)*sneu**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jg)*za(jn,je)
     &    *zb(ju,jd)*zb(ju,je)*zb(je,jg) )
      qedi(1,2,2) = qedi(1,2,2) + iza(jb,jg)*izb(ju,jg)*izb(jc,jg)*
     & sneu**(-1)*sneug**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jg)*za(
     &    jn,jg)*za(jn,je)*zb(ju,je)**2*zb(jn,jg) + za(jd,jg)*za(jn,je)
     &    *za(je,jg)*zb(ju,je)**2*zb(je,jg) )
      qedi(1,2,2) = qedi(1,2,2) + iza(jb,jg)*izb(jd,jg)*izb(jc,jg)*
     & sneu**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(ju,jn)*za(jd,jg)*zb(
     &    ju,jd)*zb(ju,je) )
      qedi(1,2,2) = qedi(1,2,2) + iza(jb,jg)*izb(jd,jg)*izb(jc,jg)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jg)*za(jd,jn)*zb(ju
     &    ,jd)*zb(jd,je) - za(jd,jg)*za(jn,jg)*zb(ju,jg)*zb(jd,je) + 
     &    za(jd,jn)*za(jn,jg)*zb(ju,jd)*zb(jn,je) - za(jn,jg)**2*zb(ju,
     &    jg)*zb(jn,je) )
      qedi(1,2,2) = qedi(1,2,2) + iza(jb,jg)*izb(jc,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(jd,jg)**2*za(jd,jn)
     &    *zb(ju,jd)*zb(jd,je) - za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jd
     &    )*zb(jn,je) - za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jn)*zb(jd,
     &    je) - za(jd,jg)*za(jd,jn)*za(je,jg)*zb(ju,je)*zb(jd,je) - za(
     &    jd,jn)*za(jn,jg)**2*zb(ju,jn)*zb(jn,je) - za(jd,jn)*za(jn,jg)
     &    *za(je,jg)*zb(ju,je)*zb(jn,je) )
      qedi(1,2,2) = qedi(1,2,2) + iza(jb,jg)*izb(jc,jg)*sneu**(-1)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(ju,jg)*za(jd,jg)*za(jn
     &    ,je)*zb(ju,je)**2 )
      qedi(1,2,2) = qedi(1,2,2) + izb(ju,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(jd,jg)*za(jd,jc)*za(jd,jn)*zb(ju,jd)
     &    *zb(ju,jb)*zb(jd,je) - za(jd,jg)*za(jd,jn)*za(jn,jc)*zb(ju,jb
     &    )*zb(ju,jn)*zb(jd,je) - za(jd,jg)*za(jd,jn)*za(je,jc)*zb(ju,
     &    jb)*zb(ju,je)*zb(jd,je) - za(jd,jc)*za(jd,jn)*za(jn,jg)*zb(ju
     &    ,jd)*zb(ju,jb)*zb(jn,je) - za(jd,jn)*za(jn,jg)*za(jn,jc)*zb(
     &    ju,jb)*zb(ju,jn)*zb(jn,je) - za(jd,jn)*za(jn,jg)*za(je,jc)*
     &    zb(ju,jb)*zb(ju,je)*zb(jn,je) )
      qedi(1,2,2) = qedi(1,2,2) + izb(ju,jg)*sneu**(-1)*sneug**(-1)*
     & sbc**(-1)*fourrt2 * ( za(ju,jg)*za(jd,jc)*za(jn,je)*zb(ju,jb)*
     &    zb(ju,je)**2 + za(jd,jc)*za(jn,jg)*za(jn,je)*zb(ju,je)**2*zb(
     &    jn,jb) + za(jd,jc)*za(jn,je)*za(je,jg)*zb(ju,je)**2*zb(je,jb)
     &     )
      qedi(1,2,2) = qedi(1,2,2) + izb(ju,jg)*izb(jd,jg)*sneu**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(ju,jn)*za(jd,jc)*zb(ju,jd)*zb(ju,jb)
     &    *zb(ju,je) + za(jd,jc)*za(jn,je)*zb(ju,jd)*zb(ju,je)*zb(je,jb
     &    ) )
      qedi(1,2,2) = qedi(1,2,2) + izb(ju,jg)*izb(jd,jg)*snedg**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jc)*za(jd,jn)*zb(ju,jd)*zb(ju,jb)*
     &    zb(jd,je) + za(jd,jn)*za(jc,jg)*zb(ju,jd)*zb(ju,jb)*zb(je,jg)
     &     + za(jd,jn)*za(jn,jc)*zb(ju,jd)*zb(ju,jb)*zb(jn,je) )
      qedi(1,2,2) = qedi(1,2,2) + izb(jd,jg)*sneu**(-1)*sbc**(-1)*
     & fourrt2 * ( za(ju,jn)*za(jc,jg)*zb(ju,jb)*zb(ju,je) - za(jc,jg)*
     &    za(jn,je)*zb(ju,je)*zb(je,jb) )
      qedi(1,2,2) = qedi(1,2,2) + izb(jd,jg)*snedg**(-1)*sbc**(-1)*
     & fourrt2 * (  - za(jd,jc)*za(jn,jg)*zb(ju,jb)*zb(jd,je) - za(jc,
     &    jg)*za(jn,jg)*zb(ju,jb)*zb(je,jg) - za(jn,jg)*za(jn,jc)*zb(ju
     &    ,jb)*zb(jn,je) )

      qedi(2,1,2)= + iza(ju,jg)*iza(jd,jg)*izb(jb,jg)*sned**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(jd,jc)*za(jd,jn)*zb(ju,jg)*
     &    zb(jd,je) + za(ju,jd)*za(jd,jn)*za(jn,jc)*zb(ju,jg)*zb(jn,je)
     &     )
      qedi(2,1,2) = qedi(2,1,2) + iza(ju,jg)*iza(jd,jg)*izb(jb,jg)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m * (  - za(ju,jd)*za(ju,jn)*za(jd
     &    ,jc)*zb(ju,jg)*zb(ju,je) + za(ju,jd)*za(jd,jc)*za(jn,je)*zb(
     &    ju,je)*zb(je,jg) )
      qedi(2,1,2) = qedi(2,1,2) + iza(ju,jg)*iza(jd,jg)*izb(jc,jg)*
     & sned**(-1)*sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(jd,jb)*za(jd,jn)
     &    *zb(ju,jg)*zb(jd,je) + za(ju,jd)*za(jd,jn)*za(jn,jb)*zb(ju,jg
     &    )*zb(jn,je) )
      qedi(2,1,2) = qedi(2,1,2) + iza(ju,jg)*iza(jd,jg)*izb(jc,jg)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m * (  - za(ju,jd)*za(ju,jn)*za(jd
     &    ,jb)*zb(ju,jg)*zb(ju,je) + za(ju,jd)*za(jd,jb)*za(jn,je)*zb(
     &    ju,je)*zb(je,jg) )
      qedi(2,1,2) = qedi(2,1,2) + iza(ju,jg)*izb(jb,jg)*sneug**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(ju,jn)*za(jd,jc)*zb(ju,jg)*zb(je,
     &    jg) + za(jd,jc)*za(jn,je)*zb(je,jg)**2 )
      qedi(2,1,2) = qedi(2,1,2) + iza(ju,jg)*izb(jc,jg)*sneug**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(ju,jn)*za(jd,jb)*zb(ju,jg)*zb(je,
     &    jg) + za(jd,jb)*za(jn,je)*zb(je,jg)**2 )
      qedi(2,1,2) = qedi(2,1,2) + iza(jd,jg)*izb(jb,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jc)*za(jd,jn)**2*
     &    zb(ju,jg)*zb(jd,jg)*zb(jn,je) - za(jd,jn)**2*za(jn,jc)*zb(ju,
     &    jg)*zb(jn,jg)*zb(jn,je) - za(jd,jn)**2*za(je,jc)*zb(ju,jg)*
     &    zb(jn,je)*zb(je,jg) )
      qedi(2,1,2) = qedi(2,1,2) + iza(jd,jg)*izb(jb,jg)*sneu**(-1)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(ju,jn)*za(jd,jc
     &    )*zb(ju,jg)**2*zb(ju,je) - za(ju,jd)*za(jd,jc)*za(jn,je)*zb(
     &    ju,jg)*zb(ju,je)*zb(je,jg) - za(ju,jn)*za(jd,jc)*za(jd,jn)*
     &    zb(ju,jg)*zb(ju,je)*zb(jn,jg) - za(ju,jn)*za(jd,jc)*za(jd,je)
     &    *zb(ju,jg)*zb(ju,je)*zb(je,jg) + za(jd,jc)*za(jd,jn)*za(jn,je
     &    )*zb(ju,je)*zb(jn,jg)*zb(je,jg) + za(jd,jc)*za(jd,je)*za(jn,
     &    je)*zb(ju,je)*zb(je,jg)**2 )
      qedi(2,1,2) = qedi(2,1,2) + iza(jd,jg)*izb(jc,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jb)*za(jd,jn)**2*
     &    zb(ju,jg)*zb(jd,jg)*zb(jn,je) - za(jd,jn)**2*za(jn,jb)*zb(ju,
     &    jg)*zb(jn,jg)*zb(jn,je) - za(jd,jn)**2*za(je,jb)*zb(ju,jg)*
     &    zb(jn,je)*zb(je,jg) )
      qedi(2,1,2) = qedi(2,1,2) + iza(jd,jg)*izb(jc,jg)*sneu**(-1)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(ju,jn)*za(jd,jb
     &    )*zb(ju,jg)**2*zb(ju,je) - za(ju,jd)*za(jd,jb)*za(jn,je)*zb(
     &    ju,jg)*zb(ju,je)*zb(je,jg) - za(ju,jn)*za(jd,jb)*za(jd,jn)*
     &    zb(ju,jg)*zb(ju,je)*zb(jn,jg) - za(ju,jn)*za(jd,jb)*za(jd,je)
     &    *zb(ju,jg)*zb(ju,je)*zb(je,jg) + za(jd,jb)*za(jd,jn)*za(jn,je
     &    )*zb(ju,je)*zb(jn,jg)*zb(je,jg) + za(jd,jb)*za(jd,je)*za(jn,
     &    je)*zb(ju,je)*zb(je,jg)**2 )

      qedi(2,2,1)= + iza(ju,jg)*iza(jd,jg)*iza(jb,jg)*sned**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(jd,jn)*za(jn,jg)*zb(ju,jc)*
     &    zb(jn,je) )
      qedi(2,2,1) = qedi(2,2,1) + iza(ju,jg)*iza(jd,jg)*iza(jc,jg)*
     & sned**(-1)*sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(jd,jn)*za(jn,jg)
     &    *zb(ju,jb)*zb(jn,je) )
      qedi(2,2,1) = qedi(2,2,1) + iza(ju,jg)*iza(jb,jg)*sned**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(jd,jn)*zb(ju,jc)*zb(jd,je)
     &     + za(jd,jg)*za(jd,jn)*zb(jd,je)*zb(jc,jg) + za(jd,jn)*za(jn,
     &    jg)*zb(jc,jg)*zb(jn,je) )
      qedi(2,2,1) = qedi(2,2,1) + iza(ju,jg)*iza(jb,jg)*sneug**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(ju,jd)*za(ju,jn)*zb(ju,jc)*zb(ju,
     &    je) - za(ju,jd)*za(jn,jg)*zb(ju,je)*zb(jc,jg) + za(ju,jd)*za(
     &    jn,je)*zb(ju,je)*zb(je,jc) - za(ju,jn)*za(jd,jg)*zb(ju,jc)*
     &    zb(je,jg) - za(jd,jg)*za(jn,jg)*zb(jc,jg)*zb(je,jg) + za(jd,
     &    jg)*za(jn,je)*zb(je,jg)*zb(je,jc) )
      qedi(2,2,1) = qedi(2,2,1) + iza(ju,jg)*iza(jc,jg)*sned**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(jd,jn)*zb(ju,jb)*zb(jd,je)
     &     + za(jd,jg)*za(jd,jn)*zb(jd,je)*zb(jb,jg) + za(jd,jn)*za(jn,
     &    jg)*zb(jb,jg)*zb(jn,je) )
      qedi(2,2,1) = qedi(2,2,1) + iza(ju,jg)*iza(jc,jg)*sneug**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(ju,jd)*za(ju,jn)*zb(ju,jb)*zb(ju,
     &    je) - za(ju,jd)*za(jn,jg)*zb(ju,je)*zb(jb,jg) + za(ju,jd)*za(
     &    jn,je)*zb(ju,je)*zb(je,jb) - za(ju,jn)*za(jd,jg)*zb(ju,jb)*
     &    zb(je,jg) - za(jd,jg)*za(jn,jg)*zb(jb,jg)*zb(je,jg) + za(jd,
     &    jg)*za(jn,je)*zb(je,jg)*zb(je,jb) )
      qedi(2,2,1) = qedi(2,2,1) + iza(jd,jg)*iza(jb,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jn)**2*za(jn,jg)*
     &    zb(ju,jc)*zb(jn,jg)*zb(jn,je) - za(jd,jn)**2*za(je,jg)*zb(ju,
     &    jc)*zb(jn,je)*zb(je,jg) )
      qedi(2,2,1) = qedi(2,2,1) + iza(jd,jg)*iza(jc,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jn)**2*za(jn,jg)*
     &    zb(ju,jb)*zb(jn,jg)*zb(jn,je) - za(jd,jn)**2*za(je,jg)*zb(ju,
     &    jb)*zb(jn,je)*zb(je,jg) )
      qedi(2,2,1) = qedi(2,2,1) + iza(jb,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jn)**2*zb(ju,jc)*zb(jd,jg)*zb(
     &    jn,je) )
      qedi(2,2,1) = qedi(2,2,1) + iza(jb,jg)*sneu**(-1)*sneug**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(ju,jn)*zb(ju,jg)*zb(ju,jc)*
     &    zb(ju,je) - za(ju,jd)*za(jn,je)*zb(ju,jc)*zb(ju,je)*zb(je,jg)
     &     + za(ju,jn)*za(jd,jg)*zb(ju,jg)*zb(ju,je)*zb(jc,jg) - za(ju,
     &    jn)*za(jd,jn)*zb(ju,jg)*zb(ju,je)*zb(jn,jc) - za(ju,jn)*za(jd
     &    ,je)*zb(ju,jg)*zb(ju,je)*zb(je,jc) - za(jd,jg)*za(jn,je)*zb(
     &    ju,je)*zb(jc,jg)*zb(je,jg) + za(jd,jn)*za(jn,je)*zb(ju,je)*
     &    zb(jn,jc)*zb(je,jg) + za(jd,je)*za(jn,je)*zb(ju,je)*zb(je,jg)
     &    *zb(je,jc) )
      qedi(2,2,1) = qedi(2,2,1) + iza(jc,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jn)**2*zb(ju,jb)*zb(jd,jg)*zb(
     &    jn,je) )
      qedi(2,2,1) = qedi(2,2,1) + iza(jc,jg)*sneu**(-1)*sneug**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(ju,jd)*za(ju,jn)*zb(ju,jg)*zb(ju,jb)*
     &    zb(ju,je) - za(ju,jd)*za(jn,je)*zb(ju,jb)*zb(ju,je)*zb(je,jg)
     &     + za(ju,jn)*za(jd,jg)*zb(ju,jg)*zb(ju,je)*zb(jb,jg) - za(ju,
     &    jn)*za(jd,jn)*zb(ju,jg)*zb(ju,je)*zb(jn,jb) - za(ju,jn)*za(jd
     &    ,je)*zb(ju,jg)*zb(ju,je)*zb(je,jb) - za(jd,jg)*za(jn,je)*zb(
     &    ju,je)*zb(jb,jg)*zb(je,jg) + za(jd,jn)*za(jn,je)*zb(ju,je)*
     &    zb(jn,jb)*zb(je,jg) + za(jd,je)*za(jn,je)*zb(ju,je)*zb(je,jg)
     &    *zb(je,jb) )

      qedi(1,1,1)= + sned**(-1)*snedg**(-1)*sbc**(-1)*fourrt2 * ( za(jd
     &    ,jg)*za(jd,jn)*za(jb,jg)*zb(ju,jc)*zb(jd,je) + za(jd,jn)*za(
     &    jb,jg)*za(jn,jg)*zb(ju,jc)*zb(jn,je) )
      qedi(1,1,1) = qedi(1,1,1) + iza(jc,jg)*izb(ju,jg)*izb(jd,jg)*izb(
     & jb,jg)*sneu**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jg)*za(jn,je)
     &    *zb(ju,jd)*zb(ju,je)*zb(je,jg) )
      qedi(1,1,1) = qedi(1,1,1) + iza(jc,jg)*izb(ju,jg)*izb(jb,jg)*
     & sneu**(-1)*sneug**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jg)*za(
     &    jn,jg)*za(jn,je)*zb(ju,je)**2*zb(jn,jg) + za(jd,jg)*za(jn,je)
     &    *za(je,jg)*zb(ju,je)**2*zb(je,jg) )
      qedi(1,1,1) = qedi(1,1,1) + iza(jc,jg)*izb(jd,jg)*izb(jb,jg)*
     & sneu**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(ju,jn)*za(jd,jg)*zb(
     &    ju,jd)*zb(ju,je) )
      qedi(1,1,1) = qedi(1,1,1) + iza(jc,jg)*izb(jd,jg)*izb(jb,jg)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jg)*za(jd,jn)*zb(ju
     &    ,jd)*zb(jd,je) - za(jd,jg)*za(jn,jg)*zb(ju,jg)*zb(jd,je) + 
     &    za(jd,jn)*za(jn,jg)*zb(ju,jd)*zb(jn,je) - za(jn,jg)**2*zb(ju,
     &    jg)*zb(jn,je) )
      qedi(1,1,1) = qedi(1,1,1) + iza(jc,jg)*izb(jb,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(jd,jg)**2*za(jd,jn)
     &    *zb(ju,jd)*zb(jd,je) - za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jd
     &    )*zb(jn,je) - za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jn)*zb(jd,
     &    je) - za(jd,jg)*za(jd,jn)*za(je,jg)*zb(ju,je)*zb(jd,je) - za(
     &    jd,jn)*za(jn,jg)**2*zb(ju,jn)*zb(jn,je) - za(jd,jn)*za(jn,jg)
     &    *za(je,jg)*zb(ju,je)*zb(jn,je) )
      qedi(1,1,1) = qedi(1,1,1) + iza(jc,jg)*izb(jb,jg)*sneu**(-1)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(ju,jg)*za(jd,jg)*za(jn
     &    ,je)*zb(ju,je)**2 )
      qedi(1,1,1) = qedi(1,1,1) + izb(ju,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(jd,jg)*za(jd,jb)*za(jd,jn)*zb(ju,jd)
     &    *zb(ju,jc)*zb(jd,je) - za(jd,jg)*za(jd,jn)*za(jn,jb)*zb(ju,jc
     &    )*zb(ju,jn)*zb(jd,je) - za(jd,jg)*za(jd,jn)*za(je,jb)*zb(ju,
     &    jc)*zb(ju,je)*zb(jd,je) - za(jd,jb)*za(jd,jn)*za(jn,jg)*zb(ju
     &    ,jd)*zb(ju,jc)*zb(jn,je) - za(jd,jn)*za(jn,jg)*za(jn,jb)*zb(
     &    ju,jc)*zb(ju,jn)*zb(jn,je) - za(jd,jn)*za(jn,jg)*za(je,jb)*
     &    zb(ju,jc)*zb(ju,je)*zb(jn,je) )
      qedi(1,1,1) = qedi(1,1,1) + izb(ju,jg)*sneu**(-1)*sneug**(-1)*
     & sbc**(-1)*fourrt2 * ( za(ju,jg)*za(jd,jb)*za(jn,je)*zb(ju,jc)*
     &    zb(ju,je)**2 + za(jd,jb)*za(jn,jg)*za(jn,je)*zb(ju,je)**2*zb(
     &    jn,jc) + za(jd,jb)*za(jn,je)*za(je,jg)*zb(ju,je)**2*zb(je,jc)
     &     )
      qedi(1,1,1) = qedi(1,1,1) + izb(ju,jg)*izb(jd,jg)*sneu**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(ju,jn)*za(jd,jb)*zb(ju,jd)*zb(ju,jc)
     &    *zb(ju,je) + za(jd,jb)*za(jn,je)*zb(ju,jd)*zb(ju,je)*zb(je,jc
     &    ) )
      qedi(1,1,1) = qedi(1,1,1) + izb(ju,jg)*izb(jd,jg)*snedg**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jb)*za(jd,jn)*zb(ju,jd)*zb(ju,jc)*
     &    zb(jd,je) + za(jd,jn)*za(jb,jg)*zb(ju,jd)*zb(ju,jc)*zb(je,jg)
     &     + za(jd,jn)*za(jn,jb)*zb(ju,jd)*zb(ju,jc)*zb(jn,je) )
      qedi(1,1,1) = qedi(1,1,1) + izb(jd,jg)*sneu**(-1)*sbc**(-1)*
     & fourrt2 * ( za(ju,jn)*za(jb,jg)*zb(ju,jc)*zb(ju,je) - za(jb,jg)*
     &    za(jn,je)*zb(ju,je)*zb(je,jc) )
      qedi(1,1,1) = qedi(1,1,1) + izb(jd,jg)*snedg**(-1)*sbc**(-1)*
     & fourrt2 * (  - za(jd,jb)*za(jn,jg)*zb(ju,jc)*zb(jd,je) - za(jb,
     &    jg)*za(jn,jg)*zb(ju,jc)*zb(je,jg) - za(jn,jg)*za(jn,jb)*zb(ju
     &    ,jc)*zb(jn,je) )

      qedi(2,1,1)= + sneu**(-1)*sneug**(-1)*sbc**(-1)*fourrt2 * (  - 
     &    za(ju,jn)*za(jd,jb)*zb(ju,jg)*zb(ju,je)*zb(jc,jg) + za(jd,jb)
     &    *za(jn,je)*zb(ju,je)*zb(jc,jg)*zb(je,jg) )
      qedi(2,1,1) = qedi(2,1,1) + iza(ju,jg)*sned**(-1)*sbc**(-1)*
     & fourrt2 * (  - za(jd,jb)*za(jd,jn)*zb(jd,je)*zb(jc,jg) - za(jd,
     &    jn)*za(jn,jb)*zb(jc,jg)*zb(jn,je) )
      qedi(2,1,1) = qedi(2,1,1) + iza(ju,jg)*sneug**(-1)*sbc**(-1)*
     & fourrt2 * ( za(ju,jn)*za(jd,jb)*zb(ju,jc)*zb(je,jg) + za(jd,jb)*
     &    za(jn,jg)*zb(jc,jg)*zb(je,jg) - za(jd,jb)*za(jn,je)*zb(je,jg)
     &    *zb(je,jc) )
      qedi(2,1,1) = qedi(2,1,1) + iza(ju,jg)*iza(jd,jg)*sned**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(ju,jd)*za(jd,jb)*za(jd,jn)*zb(ju,jc)
     &    *zb(jd,je) - za(ju,jd)*za(jd,jn)*za(jn,jb)*zb(ju,jc)*zb(jn,je
     &    ) )
      qedi(2,1,1) = qedi(2,1,1) + iza(ju,jg)*iza(jd,jg)*sneug**(-1)*
     & sbc**(-1)*fourrt2 * ( za(ju,jd)*za(ju,jn)*za(jd,jb)*zb(ju,jc)*
     &    zb(ju,je) + za(ju,jd)*za(jd,jb)*za(jn,jg)*zb(ju,je)*zb(jc,jg)
     &     - za(ju,jd)*za(jd,jb)*za(jn,je)*zb(ju,je)*zb(je,jc) )
      qedi(2,1,1) = qedi(2,1,1) + iza(ju,jg)*iza(jd,jg)*iza(jc,jg)*izb(
     & jb,jg)*sned**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(ju,jd)*za(jd,
     &    jn)*za(jn,jg)*zb(ju,jg)*zb(jn,je) )
      qedi(2,1,1) = qedi(2,1,1) + iza(ju,jg)*iza(jc,jg)*izb(jb,jg)*
     & sned**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(ju,jd)*za(jd,jn)*zb(
     &    ju,jg)*zb(jd,je) )
      qedi(2,1,1) = qedi(2,1,1) + iza(ju,jg)*iza(jc,jg)*izb(jb,jg)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(ju,jd)*za(ju,jn)*zb(ju
     &    ,jg)*zb(ju,je) - za(ju,jd)*za(jn,je)*zb(ju,je)*zb(je,jg) + 
     &    za(ju,jn)*za(jd,jg)*zb(ju,jg)*zb(je,jg) - za(jd,jg)*za(jn,je)
     &    *zb(je,jg)**2 )
      qedi(2,1,1) = qedi(2,1,1) + iza(jd,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2 * ( za(jd,jb)*za(jd,jn)**2*zb(ju,jc)*zb(jd,jg)
     &    *zb(jn,je) + za(jd,jn)**2*za(jn,jb)*zb(ju,jc)*zb(jn,jg)*zb(jn
     &    ,je) + za(jd,jn)**2*za(je,jb)*zb(ju,jc)*zb(jn,je)*zb(je,jg) )
      qedi(2,1,1) = qedi(2,1,1) + iza(jd,jg)*sneu**(-1)*sneug**(-1)*
     & sbc**(-1)*fourrt2 * (  - za(ju,jd)*za(ju,jn)*za(jd,jb)*zb(ju,jg)
     &    *zb(ju,jc)*zb(ju,je) + za(ju,jd)*za(jd,jb)*za(jn,je)*zb(ju,jc
     &    )*zb(ju,je)*zb(je,jg) + za(ju,jn)*za(jd,jb)*za(jd,jn)*zb(ju,
     &    jg)*zb(ju,je)*zb(jn,jc) + za(ju,jn)*za(jd,jb)*za(jd,je)*zb(ju
     &    ,jg)*zb(ju,je)*zb(je,jc) - za(jd,jb)*za(jd,jn)*za(jn,je)*zb(
     &    ju,je)*zb(jn,jc)*zb(je,jg) - za(jd,jb)*za(jd,je)*za(jn,je)*
     &    zb(ju,je)*zb(je,jg)*zb(je,jc) )
      qedi(2,1,1) = qedi(2,1,1) + iza(jd,jg)*iza(jc,jg)*izb(jb,jg)*
     & sned**(-1)*snedg**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jn)**2*
     &    za(jn,jg)*zb(ju,jg)*zb(jn,jg)*zb(jn,je) + za(jd,jn)**2*za(je,
     &    jg)*zb(ju,jg)*zb(jn,je)*zb(je,jg) )
      qedi(2,1,1) = qedi(2,1,1) + iza(jc,jg)*izb(jb,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m**2 * ( za(jd,jn)**2*zb(ju,jg)*
     &    zb(jd,jg)*zb(jn,je) )
      qedi(2,1,1) = qedi(2,1,1) + iza(jc,jg)*izb(jb,jg)*sneu**(-1)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m**2 * (  - za(ju,jd)*za(ju,jn)*
     &    zb(ju,jg)**2*zb(ju,je) + za(ju,jd)*za(jn,je)*zb(ju,jg)*zb(ju,
     &    je)*zb(je,jg) + za(ju,jn)*za(jd,jn)*zb(ju,jg)*zb(ju,je)*zb(jn
     &    ,jg) + za(ju,jn)*za(jd,je)*zb(ju,jg)*zb(ju,je)*zb(je,jg) - 
     &    za(jd,jn)*za(jn,je)*zb(ju,je)*zb(jn,jg)*zb(je,jg) - za(jd,je)
     &    *za(jn,je)*zb(ju,je)*zb(je,jg)**2 )

      qedi(1,2,1)= + iza(jb,jg)*izb(ju,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(jd,jg)**2*za(jd,jn)*zb(ju,jd)*zb(ju,
     &    jc)*zb(jd,je) + za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jd)*zb(ju
     &    ,jc)*zb(jn,je) + za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(ju,jc)*zb(
     &    ju,jn)*zb(jd,je) + za(jd,jg)*za(jd,jn)*za(je,jg)*zb(ju,jc)*
     &    zb(ju,je)*zb(jd,je) + za(jd,jn)*za(jn,jg)**2*zb(ju,jc)*zb(ju,
     &    jn)*zb(jn,je) + za(jd,jn)*za(jn,jg)*za(je,jg)*zb(ju,jc)*zb(ju
     &    ,je)*zb(jn,je) )
      qedi(1,2,1) = qedi(1,2,1) + iza(jb,jg)*izb(ju,jg)*sneu**(-1)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m * (  - za(ju,jg)*za(jd,jg)*za(jn
     &    ,je)*zb(ju,jc)*zb(ju,je)**2 - za(jd,jg)*za(jn,jg)*za(jn,je)*
     &    zb(ju,je)**2*zb(jn,jc) - za(jd,jg)*za(jn,je)*za(je,jg)*zb(ju,
     &    je)**2*zb(je,jc) )
      qedi(1,2,1) = qedi(1,2,1) + iza(jb,jg)*izb(ju,jg)*izb(jd,jg)*
     & sneu**(-1)*sbc**(-1)*fourrt2*m * ( za(ju,jn)*za(jd,jg)*zb(ju,jd)
     &    *zb(ju,jc)*zb(ju,je) - za(jd,jg)*za(jn,je)*zb(ju,jd)*zb(ju,je
     &    )*zb(je,jc) )
      qedi(1,2,1) = qedi(1,2,1) + iza(jb,jg)*izb(ju,jg)*izb(jd,jg)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jg)*za(jd,jn)*zb(ju
     &    ,jd)*zb(ju,jc)*zb(jd,je) - za(jd,jn)*za(jn,jg)*zb(ju,jd)*zb(
     &    ju,jc)*zb(jn,je) )
      qedi(1,2,1) = qedi(1,2,1) + iza(jb,jg)*izb(jd,jg)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(jd,jg)*za(jn,jg)*zb(ju,jc)*zb(jd,je)
     &     + za(jn,jg)**2*zb(ju,jc)*zb(jn,je) )
      qedi(1,2,1) = qedi(1,2,1) + iza(jc,jg)*izb(ju,jg)*sned**(-1)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * ( za(jd,jg)**2*za(jd,jn)*zb(ju
     &    ,jd)*zb(ju,jb)*zb(jd,je) + za(jd,jg)*za(jd,jn)*za(jn,jg)*zb(
     &    ju,jd)*zb(ju,jb)*zb(jn,je) + za(jd,jg)*za(jd,jn)*za(jn,jg)*
     &    zb(ju,jb)*zb(ju,jn)*zb(jd,je) + za(jd,jg)*za(jd,jn)*za(je,jg)
     &    *zb(ju,jb)*zb(ju,je)*zb(jd,je) + za(jd,jn)*za(jn,jg)**2*zb(ju
     &    ,jb)*zb(ju,jn)*zb(jn,je) + za(jd,jn)*za(jn,jg)*za(je,jg)*zb(
     &    ju,jb)*zb(ju,je)*zb(jn,je) )
      qedi(1,2,1) = qedi(1,2,1) + iza(jc,jg)*izb(ju,jg)*sneu**(-1)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m * (  - za(ju,jg)*za(jd,jg)*za(jn
     &    ,je)*zb(ju,jb)*zb(ju,je)**2 - za(jd,jg)*za(jn,jg)*za(jn,je)*
     &    zb(ju,je)**2*zb(jn,jb) - za(jd,jg)*za(jn,je)*za(je,jg)*zb(ju,
     &    je)**2*zb(je,jb) )
      qedi(1,2,1) = qedi(1,2,1) + iza(jc,jg)*izb(ju,jg)*izb(jd,jg)*
     & sneu**(-1)*sbc**(-1)*fourrt2*m * ( za(ju,jn)*za(jd,jg)*zb(ju,jd)
     &    *zb(ju,jb)*zb(ju,je) - za(jd,jg)*za(jn,je)*zb(ju,jd)*zb(ju,je
     &    )*zb(je,jb) )
      qedi(1,2,1) = qedi(1,2,1) + iza(jc,jg)*izb(ju,jg)*izb(jd,jg)*
     & snedg**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jg)*za(jd,jn)*zb(ju
     &    ,jd)*zb(ju,jb)*zb(jd,je) - za(jd,jn)*za(jn,jg)*zb(ju,jd)*zb(
     &    ju,jb)*zb(jn,je) )
      qedi(1,2,1) = qedi(1,2,1) + iza(jc,jg)*izb(jd,jg)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(jd,jg)*za(jn,jg)*zb(ju,jb)*zb(jd,je)
     &     + za(jn,jg)**2*zb(ju,jb)*zb(jn,je) )

      qedi(1,1,2)= + izb(ju,jg)*izb(jd,jg)*izb(jb,jg)*sneu**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jc)*za(jn,je)*zb(ju,jd)*zb(ju,
     &    je)*zb(je,jg) )
      qedi(1,1,2) = qedi(1,1,2) + izb(ju,jg)*izb(jd,jg)*izb(jc,jg)*
     & sneu**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jb)*za(jn,je)*zb(ju,
     &    jd)*zb(ju,je)*zb(je,jg) )
      qedi(1,1,2) = qedi(1,1,2) + izb(ju,jg)*izb(jb,jg)*sneu**(-1)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jc)*za(jn,jg)*za(jn
     &    ,je)*zb(ju,je)**2*zb(jn,jg) - za(jd,jc)*za(jn,je)*za(je,jg)*
     &    zb(ju,je)**2*zb(je,jg) )
      qedi(1,1,2) = qedi(1,1,2) + izb(ju,jg)*izb(jc,jg)*sneu**(-1)*
     & sneug**(-1)*sbc**(-1)*fourrt2*m * (  - za(jd,jb)*za(jn,jg)*za(jn
     &    ,je)*zb(ju,je)**2*zb(jn,jg) - za(jd,jb)*za(jn,je)*za(je,jg)*
     &    zb(ju,je)**2*zb(je,jg) )
      qedi(1,1,2) = qedi(1,1,2) + izb(jd,jg)*izb(jb,jg)*sneu**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(ju,jn)*za(jd,jc)*zb(ju,jd)*zb(ju,je)
     &     - za(ju,jn)*za(jc,jg)*zb(ju,jg)*zb(ju,je) + za(jc,jg)*za(jn,
     &    je)*zb(ju,je)*zb(je,jg) )
      qedi(1,1,2) = qedi(1,1,2) + izb(jd,jg)*izb(jb,jg)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jc)*za(jd,jn)*zb(ju,jd)*zb(jd,
     &    je) + za(jd,jc)*za(jn,jg)*zb(ju,jg)*zb(jd,je) - za(jd,jn)*za(
     &    jc,jg)*zb(ju,jd)*zb(je,jg) - za(jd,jn)*za(jn,jc)*zb(ju,jd)*
     &    zb(jn,je) + za(jc,jg)*za(jn,jg)*zb(ju,jg)*zb(je,jg) + za(jn,
     &    jg)*za(jn,jc)*zb(ju,jg)*zb(jn,je) )
      qedi(1,1,2) = qedi(1,1,2) + izb(jd,jg)*izb(jc,jg)*sneu**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(ju,jn)*za(jd,jb)*zb(ju,jd)*zb(ju,je)
     &     - za(ju,jn)*za(jb,jg)*zb(ju,jg)*zb(ju,je) + za(jb,jg)*za(jn,
     &    je)*zb(ju,je)*zb(je,jg) )
      qedi(1,1,2) = qedi(1,1,2) + izb(jd,jg)*izb(jc,jg)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(jd,jb)*za(jd,jn)*zb(ju,jd)*zb(jd,
     &    je) + za(jd,jb)*za(jn,jg)*zb(ju,jg)*zb(jd,je) - za(jd,jn)*za(
     &    jb,jg)*zb(ju,jd)*zb(je,jg) - za(jd,jn)*za(jn,jb)*zb(ju,jd)*
     &    zb(jn,je) + za(jb,jg)*za(jn,jg)*zb(ju,jg)*zb(je,jg) + za(jn,
     &    jg)*za(jn,jb)*zb(ju,jg)*zb(jn,je) )
      qedi(1,1,2) = qedi(1,1,2) + izb(jb,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(jd,jg)*za(jd,jc)*za(jd,jn)*zb(ju,jd)*
     &    zb(jd,je) - za(jd,jg)*za(jd,jn)*za(jc,jg)*zb(ju,jg)*zb(jd,je)
     &     + za(jd,jg)*za(jd,jn)*za(jn,jc)*zb(ju,jn)*zb(jd,je) + za(jd,
     &    jg)*za(jd,jn)*za(je,jc)*zb(ju,je)*zb(jd,je) + za(jd,jc)*za(jd
     &    ,jn)*za(jn,jg)*zb(ju,jd)*zb(jn,je) - za(jd,jn)*za(jc,jg)*za(
     &    jn,jg)*zb(ju,jg)*zb(jn,je) + za(jd,jn)*za(jn,jg)*za(jn,jc)*
     &    zb(ju,jn)*zb(jn,je) + za(jd,jn)*za(jn,jg)*za(je,jc)*zb(ju,je)
     &    *zb(jn,je) )
      qedi(1,1,2) = qedi(1,1,2) + izb(jb,jg)*sneu**(-1)*sneug**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(ju,jg)*za(jd,jc)*za(jn,je)*zb(ju,
     &    je)**2 )
      qedi(1,1,2) = qedi(1,1,2) + izb(jc,jg)*sned**(-1)*snedg**(-1)*
     & sbc**(-1)*fourrt2*m * ( za(jd,jg)*za(jd,jb)*za(jd,jn)*zb(ju,jd)*
     &    zb(jd,je) - za(jd,jg)*za(jd,jn)*za(jb,jg)*zb(ju,jg)*zb(jd,je)
     &     + za(jd,jg)*za(jd,jn)*za(jn,jb)*zb(ju,jn)*zb(jd,je) + za(jd,
     &    jg)*za(jd,jn)*za(je,jb)*zb(ju,je)*zb(jd,je) + za(jd,jb)*za(jd
     &    ,jn)*za(jn,jg)*zb(ju,jd)*zb(jn,je) - za(jd,jn)*za(jb,jg)*za(
     &    jn,jg)*zb(ju,jg)*zb(jn,je) + za(jd,jn)*za(jn,jg)*za(jn,jb)*
     &    zb(ju,jn)*zb(jn,je) + za(jd,jn)*za(jn,jg)*za(je,jb)*zb(ju,je)
     &    *zb(jn,je) )
      qedi(1,1,2) = qedi(1,1,2) + izb(jc,jg)*sneu**(-1)*sneug**(-1)*
     & sbc**(-1)*fourrt2*m * (  - za(ju,jg)*za(jd,jb)*za(jn,je)*zb(ju,
     &    je)**2 )

      return
      end
