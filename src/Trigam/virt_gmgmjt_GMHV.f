!=========== virtual amplitude (LEADING COLOUR) for gamma gamma jet with the following helicity structure
!========== q^-(i1)+qb^+(i2)+g^-(i3)+gamma^+(i4)+gamma^+(i5)
!===== C. Williams March 2013
      double complex function virt_gmgmjt_GMHV(i1,i2,i3,i4,i5,za,zb) 
      implicit none 
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f' 
      include 'scale.f'
      include 'epinv.f'
      integer i1,i2,i3,i4,i5
!======= box coefficients 
      double complex bxs14s25s23,bxs15s24s23,bxs24s15s13
      double complex bxs25s14s13,box_sum
!======= bubble coefficients
      double complex c14,c23,c15,bub_sum 
      double complex rat,comp_rat,rat_sum
!======== pole and finite pieces 
      double complex Vcc, Fcc
!======== basis functions
      double complex Lsm1,lnrat,l13,l23
      double complex A5LO,amp_2gam1g
      double complex zab,zab2
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zab2(i1,i2,i3,i4)=zab(i1,i2,i4)+zab(i1,i3,i4)


      virt_gmgmjt_GMHV=czip
      rat_sum=czip
      bub_sum=czip
      box_sum=czip
      A5LO=amp_2gam1g(i1,i2,i3,i4,i5,za,zb) 
      l13=lnrat(musq,-s(i1,i3))
      l23=lnrat(musq,-s(i2,i3))
!======== poles, 
      Vcc=(epinv**2+epinv*l13+0.5d0*l13**2)
     &   +(epinv**2+epinv*l23+0.5d0*l23**2)
     &   +3d0/2d0*(epinv+l23+2d0)
   
      Vcc=Vcc*A5LO
!======== box coefficients 
      
      bxs14s25s23=(za(i1,i3)**2/(za(i1,i4)*za(i2,i5)*za(i4,i5)))

!      write(6,*) 'bx s14 s25 s23 ',bxs14s25s23*im
!     &*(s(i2,i5)*s(i2,i3)/(2d0*A5LO))

      bxs15s24s23=-za(i1,i3)**2/(za(i1,i5)*za(i2,i4)*za(i4,i5))
!      write(6,*) 'bx s15 s24 s23 ',bxs15s24s23*im
!     &*(s(i2,i4)*s(i2,i3)/(2d0*A5LO))
      
      bxs24s15s13=-(za(i1,i3)**2/(za(i1,i5)*za(i2,i4)*za(i4,i5)))
!      write(6,*) 'bx s24 s15 s13 ',bxs24s15s13*im
!     &*(s(i1,i5)*s(i1,i3)/(2d0*A5LO))

      bxs25s14s13=za(i1,i3)**2/(za(i1,i4)*za(i2,i5)*za(i4,i5))
      
!      write(6,*) 'bx s25 s14 s13 ',bxs25s14s13*im
!     &*(s(i1,i4)*s(i1,i3)/(2d0*A5LO))

!======= box summation 
      box_sum=
     &     +bxs14s25s23*Lsm1(-s(i2,i5),-s(i1,i4),-s(i2,i3),-s(i1,i4))
     &     +bxs15s24s23*Lsm1(-s(i2,i4),-s(i1,i5),-s(i2,i3),-s(i1,i5))
     &     +bxs24s15s13*Lsm1(-s(i1,i5),-s(i2,i4),-s(i1,i3),-s(i2,i4))
     &     +bxs25s14s13*Lsm1(-s(i1,i4),-s(i2,i5),-s(i1,i3),-s(i2,i5))
 

!=============== bubbles 

      c14= -(-((za(i1,i3)*za(i3,i5)*zb(i5,i4))/
     &     (za(i2,i5)*za(i4,i5)*zab2(i5,i1,i4,i5))) - 
     &  (za(i1,i4)*za(i3,i5)**2*zb(i5,i4)**2)/
     &   (2.*za(i2,i5)*za(i4,i5)*zab2(i5,i1,i4,i5)**2))
      
      c15=-(-((za(i1,i3)*za(i3,i4)*zb(i5,i4))/
     &     (za(i2,i4)*za(i4,i5)*zab2(i4,i1,i5,i4))) + 
     &  (za(i1,i5)*za(i3,i4)**2*zb(i5,i4)**2)/
     &   (2.*za(i2,i4)*za(i4,i5)*zab2(i4,i1,i5,i4)**2))

      c23=-c15-c14
!      write(6,*) 'c14 ',im*c14/A5LO
!      write(6,*) 'c15 ',im*c15/A5LO

      bub_sum=
     &     +c14*lnrat(musq,-s(i1,i4))
     &     +c15*lnrat(musq,-s(i1,i5))
     &     +c23*lnrat(musq,-s(i2,i3))

      
      rat=  -(za(i1,i3)*(-(za(i2,i5)*za(i3,i4)) - 
     &        za(i2,i4)*za(i3,i5))*zb(i5,i4))/
     &   (4.*za(i2,i3)*za(i2,i4)*za(i2,i5)*za(i4,i5)*
     &     zb(i3,i2)) - 
     &  ((za(i2,i4)*zb(i4,i2) - za(i2,i5)*zb(i5,i2))*
     &     zb(i5,i4))/
     &   (4.*za(i2,i4)*za(i2,i5)*zb(i3,i1)*zb(i3,i2)) + 
     &  (zb(i2,i1)*(zb(i4,i2)*zb(i5,i1) + 
     &       zb(i4,i1)*zb(i5,i2))*zb(i5,i4))/
     &   (4.*za(i4,i5)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*
     &     zb(i5,i1))

      comp_rat=(za(i3,i4)**2*(za(i2,i3)*zb(i3,i2) + 
     -      za(i1,i5)*zb(i5,i1))*zb(i5,i4)**2)/
     -  (4.*za(i2,i3)*za(i2,i4)*za(i4,i5)*
     -    zab2(i4,i1,i5,i4)*zb(i3,i2)*zb(i5,i1))
      
      comp_rat=comp_rat -(za(i3,i5)**2*
     -     (za(i2,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i1))*
     -     zb(i5,i4)**2)/
     -  (4.*za(i2,i3)*za(i2,i5)*za(i4,i5)*
     -    zab2(i5,i1,i4,i5)*zb(i3,i2)*zb(i4,i1))

      rat_sum=rat+comp_rat
      virt_gmgmjt_GMHV=(box_sum+Vcc+bub_sum+rat_sum)

!      write(6,*) 'Vcc (finite) ',Vcc
!      write(6,*) 'box sum ',box_sum
!      write(6,*) 'bub sum ',bub_sum
!      write(6,*) 'rat sum ',rat_sum
!      write(6,*) 'Check',im*(box_sum+Vcc+bub_sum+rat_sum)
      return 
      end
      
      
