!===== nf loops for gamma gamma jet helicity structure : 
!===== q(i1)^-qb(i2)^+g(i3)^+gamma(i4)^+gamma(i5)^-

      double complex function virt_gmgmjt_nfGaMHV(i1,i2,i3,i4,i5
     &,za,zb)
      implicit none 
      include 'constants.f'
      include 'zprods_decl.f'
      include 'scale.f' 
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
!====== boxes 
      double complex box_sum,bxs12s34s45,Lsm1
!====== bubbles
      double complex bub_sum,c35,c45,c12,lnrat 
!====== rational 
      double complex rat,comp_rat,rat_sum 
      double complex zab,zab2
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zab2(i1,i2,i3,i4)=zab(i1,i2,i4)+zab(i1,i3,i4)
      
      
      virt_gmgmjt_nfGaMHV=czip
!======= boxes
      bxs12s34s45=-(za(i1,i4)**2*za(i3,i5)**2 
     &+ za(i1,i3)**2*za(i4,i5)**2)/(za(i1,i2)*za(i3,i4)**4)
      box_sum=
     & 2d0*bxs12s34s45*Lsm1(-s(i3,i5),-s(i1,i2),-s(i4,i5),-s(i1,i2))

!====== bubs
      c35= (-2*za(i1,i4)**2*za(i3,i5)**2*zb(i4,i3))/
     &   (za(i1,i2)*za(i3,i4)**3*zab2(i4,i3,i5,i4)) - 
     &  (4*za(i1,i3)*za(i1,i4)*za(i3,i5)*za(i4,i5)*zb(i4,i3))/
     &   (za(i1,i2)*za(i3,i4)**3*zab2(i4,i3,i5,i4)) + 
     &  (2*za(i1,i4)**2*za(i3,i5)**2*za(i4,i5)*zb(i4,i3)*zb(i5,i4))/
     &   (za(i1,i2)*za(i3,i4)**3*zab2(i4,i3,i5,i4)**2)

      c45= (-4*za(i1,i3)*za(i1,i4)*za(i3,i5)*za(i4,i5)*zb(i4,i3))/
     &   (za(i1,i2)*za(i3,i4)**3*zab2(i3,i4,i5,i3)) - 
     &  (2*za(i1,i3)**2*za(i4,i5)**2*zb(i4,i3))/
     &   (za(i1,i2)*za(i3,i4)**3*zab2(i3,i4,i5,i3)) + 
     &  (2*za(i1,i3)**2*za(i3,i5)*za(i4,i5)**2*zb(i4,i3)*zb(i5,i3))/
     &   (za(i1,i2)*za(i3,i4)**3*zab2(i3,i4,i5,i3)**2)

      c12=-c45-c35

      bub_sum=
     &    - c35*lnrat(musq,-s(i3,i5))
     &    - c45*lnrat(musq,-s(i4,i5))
     &    - c12*lnrat(musq,-s(i1,i2))
      
      rat= (-2*za(i1,i3)*za(i1,i4)*za(i3,i5)*za(i4,i5)*zb(i4,i3))/
     -   (za(i1,i2)**2*za(i3,i4)**3*zb(i2,i1)) - 
     -  (za(i3,i5)*za(i4,i5)*zb(i4,i3)*zb(i5,i2)**2)/
     -   (za(i3,i4)**3*zb(i2,i1)*zb(i5,i3)*zb(i5,i4))

      comp_rat= -((za(i1,i4)**2*za(i3,i5)*za(i4,i5)*zb(i4,i3)*
     -       (za(i1,i2)*zb(i2,i1) + za(i3,i5)*zb(i5,i3))*
     -       zb(i5,i4))/
     -     (za(i1,i2)**2*za(i3,i4)**3*zab2(i4,i1,i2,i4)*
     -       zb(i2,i1)*zb(i5,i3))) - 
     -  (za(i1,i3)**2*za(i3,i5)*za(i4,i5)*zb(i4,i3)*zb(i5,i3)*
     -     (za(i1,i2)*zb(i2,i1) + za(i4,i5)*zb(i5,i4)))/
     -   (za(i1,i2)**2*za(i3,i4)**3*zab2(i3,i1,i2,i3)*
     -     zb(i2,i1)*zb(i5,i4))

      rat_sum=comp_rat+rat

      virt_gmgmjt_nfGaMHV=bub_sum+box_sum+rat_sum

!      write(6,*) 'box coeff ', im*bxs12s34s45/(2d0*A5LO)
!     & *s(i3,i5)*s(i4,i5)
!      write(6,*) 'bub coeff 35 ', im*c35/A5LO
!      write(6,*) 'bub coeff 45 ', im*c45/A5LO
!      write(6,*) 'bub coeff 12 ', im*c12/A5LO
!      write(6,*) 'box sum ',box_sum 
!      write(6,*) 'bub sum ',bub_sum 
!      write(6,*) 'bub + box ',(box_sum+bub_sum)*im
!      write(6,*) 'rat sum',rat_sum*im 
!      write(6,*) 'Ans = ',virt_gmgmjt_nfGaMHV
     
      return 
      end 
