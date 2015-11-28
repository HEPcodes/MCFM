!=========== virtual amplitude (LEADING COLOUR) for gamma gamma jet with the following helicity structure
!========== q^-(i1)+qb^+(i2)+g^+(i3)+gamma^+(i4)+gamma^-(i5)
!===== C. Williams March 2013
      double complex function virt_gmgmjt_GaMHV(i1,i2,i3,i4,i5,za,zb) 
      implicit none 
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f' 
      include 'scale.f'
      include 'epinv.f'
      integer i1,i2,i3,i4,i5
!======= box coefficients 
      double complex bxs13s45s25,bxs15s24s23,bxs23s45s15
      double complex bxs25s14s13,bxs45s13s23,box_sum
!======= bubble coefficients
      double complex c14,c25,c13,c45,bub_sum 
      double complex rat,comp_rat,rat_sum
!======== pole and finite pieces 
      double complex Vcc, Fcc
!======== basis functions
      double complex Lsm1,lnrat,l13,l23,l25
      double complex A5LO,amp_2gam1g
      double complex zab,zab2
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zab2(i1,i2,i3,i4)=zab(i1,i2,i4)+zab(i1,i3,i4)


      virt_gmgmjt_GaMHV=czip
      rat_sum=czip
      bub_sum=czip
      box_sum=czip
      A5LO=amp_2gam1g(i1,i2,i5,i4,i3,za,zb) 
      l13=lnrat(musq,-s(i1,i3))
      l23=lnrat(musq,-s(i2,i3))
      l25=lnrat(musq,-s(i2,i5))
!======== poles, 
      Vcc=(epinv**2+epinv*l13+0.5d0*l13**2)
     &   +(epinv**2+epinv*l23+0.5d0*l23**2)
     &   +3d0/2d0*(epinv+l25+2d0)
   
      Vcc=Vcc*A5LO

!======== box coefficients 
      
      bxs13s45s25=-((za(i1,i2)**3*za(i4,i5)**2)
     &     /(za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4)**3))

      bxs15s24s23=za(i1,i5)**2/(za(i1,i3)*za(i2,i4)*za(i3,i4))

      bxs23s45s15=-((za(i1,i2)*za(i1,i5)**2)
     &     /(za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4)))

      bxs25s14s13=-((za(i1,i3)**2*za(i4,i5)**2)
     &     /(za(i1,i4)*za(i2,i3)*za(i3,i4)**3))

      bxs45s13s23=-((za(i1,i2)*za(i1,i5)**2)
     &     /(za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4)))

!======= box summation 
      box_sum=
     &     +bxs13s45s25*Lsm1(-s(i4,i5),-s(i1,i3),-s(i2,i5),-s(i1,i3))
     &     +bxs15s24s23*Lsm1(-s(i2,i4),-s(i1,i5),-s(i2,i3),-s(i1,i5))
     &     +bxs23s45s15*Lsm1(-s(i4,i5),-s(i2,i3),-s(i1,i5),-s(i2,i3))
     &     +bxs25s14s13*Lsm1(-s(i1,i4),-s(i2,i5),-s(i1,i3),-s(i2,i5))
     &     +bxs45s13s23*Lsm1(-s(i1,i3),-s(i4,i5),-s(i2,i3),-s(i4,i5))
 
!=============== bubbles 
      c13=-((za(i1,i2)*za(i2,i5)**2*zb(i3,i2))/
     -     (za(i2,i3)*za(i2,i4)**2*zab2(i2,i1,i3,i2))) + 
     -  (za(i1,i3)*za(i4,i5)**2*zb(i4,i3))/
     -   (za(i2,i4)*za(i3,i4)**2*zab2(i4,i1,i3,i4)) + 
     -  (za(i1,i2)*za(i4,i5)**2*zb(i4,i3))/
     -   (za(i2,i4)**2*za(i3,i4)*zab2(i4,i1,i3,i4)) + 
     -  (za(i1,i3)*za(i4,i5)**2*zb(i4,i3)**2)/
     -   (2.*za(i2,i4)*za(i3,i4)*zab2(i4,i1,i3,i4)**2)

!      write(6,*) c13/A5LO*im
      c14=-(za(i1,i5)*za(i3,i5)*zb(i4,i3))/
     -   (2.*za(i2,i3)*za(i3,i4)*zab2(i3,i1,i4,i3)) + 
     -  (za(i1,i3)*za(i3,i5)*za(i4,i5)*zb(i4,i3))/
     -   (za(i2,i3)*za(i3,i4)**2*zab2(i3,i1,i4,i3)) + 
     -  (za(i1,i3)*za(i2,i5)*za(i3,i5)*zb(i3,i2)*zb(i4,i3))/
     -   (2.*za(i2,i3)*za(i3,i4)*zab2(i3,i1,i4,i3)**2)
      
!      write(6,*) c14/A5LO*im
      c45=(za(i1,i2)**2*za(i2,i5)*za(i4,i5)*zb(i4,i2))/
     -  (za(i1,i3)*za(i2,i3)*za(i2,i4)**2*zab2(i2,i4,i5,i2))

!      write(6,*) c45/A5LO*im
      c25=-c45-c13-c14
!      write(6,*) (c25+3d0/2d0*A5LO)*im/A5LO

      bub_sum=
     &     +c13*lnrat(musq,-s(i1,i3))
     &     +c14*lnrat(musq,-s(i1,i4))
     &     +c25*lnrat(musq,-s(i2,i5))
     &     +c45*lnrat(musq,-s(i4,i5))


!================ rational       
      rat=-(za(i1,i5)*za(i3,i5)*zb(i4,i3))/
     -   (4.*za(i1,i4)*za(i2,i3)*za(i3,i4)*zb(i4,i1)) - 
     -  (za(i1,i5)*za(i4,i5)*zb(i4,i3))/
     -   (4.*za(i2,i4)*za(i2,i5)*za(i3,i4)*zb(i5,i2)) - 
     -  (zb(i3,i2)*zb(i4,i3))/(4.*za(i2,i4)*zb(i5,i1)*zb(i5,i2)) - 
     -  (zb(i4,i2)*zb(i4,i3))/(4.*za(i2,i3)*zb(i5,i1)*zb(i5,i2)) - 
     -  (zb(i2,i1)**2*zb(i4,i3)**2)/
     -   (4.*za(i3,i4)*zb(i3,i1)*zb(i4,i1)*zb(i5,i1)*zb(i5,i2))

      comp_rat=-(za(i4,i5)**2*zb(i4,i3)**2*
     -     (za(i1,i3)*zb(i3,i1) + za(i2,i5)*zb(i5,i2)))/
     -  (4.*za(i2,i4)*za(i2,i5)*za(i3,i4)*zab2(i4,i1,i3,i4)*
     -    zb(i3,i1)*zb(i5,i2))

      comp_rat=comp_rat-(za(i1,i3)*za(i3,i5)*zb(i3,i2)*zb(i4,i3)*
     -     (za(i1,i4)*zb(i4,i1) + za(i2,i5)*zb(i5,i2)))/
     -  (4.*za(i1,i4)*za(i2,i3)*za(i3,i4)*zab2(i3,i1,i4,i3)*
     -    zb(i4,i1)*zb(i5,i2))

      rat_sum=rat+comp_rat

      virt_gmgmjt_GaMHV=Vcc+box_sum+rat_sum+bub_sum

!      write(6,*) 'Vcc (finite) ',Vcc
!      write(6,*) 'box sum ',box_sum
!      write(6,*) 'bub sum ',bub_sum
!      write(6,*) 'rat sum ',rat_sum
!      write(6,*) 'Check',im*(box_sum+Vcc+bub_sum+rat_sum)
      return 
      end
      
      
