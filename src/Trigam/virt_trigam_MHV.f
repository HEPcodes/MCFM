!====== C. Williams March 2013
!====== Virtual Amplitude for the following tri-gamma process 
!===== q(i1)^- + qb(i2)^+ +gamma(i3)^+ + gamma(i4)^+ + gamma(i5)^- 

!===== version 1, list coefficients directly from mathematica with minimal changes 
!===== more specifically, only change = addition of Completion rational term 
!===== and factor of st/2 extracted from box functions to determine Lsm1 function  
      double complex function virt_trigam_MHV(i1,i2,i3,i4,i5,za,zb) 
      implicit none 
      include 'constants.f' 
      include 'sprods_com.f' 
      include 'zprods_decl.f' 
      include 'epinv.f' 
      include 'epinv2.f' 
      include 'scale.f'
      integer i1,i2,i3,i4,i5
      double complex Vcc,A5LO,Fcc
!====== box coefficients 
      double complex bxs12s35s45,bxs13s45s25,bxs14s35s25,bxs23s45s15,
     & bxs24s35s15,bxs45s13s12,bxs35s14s12,bxs35s12s24,bxs45s12s23
      double complex Box_sum,bub_sum,rat_sum
      double complex zab,zab2
      double complex comp_rat
!====== bubble coefficients
      double complex c13,c14,c45,c35,c25
!===== logs 
      double complex lnrat,l12,l25,Lsm1,qlLsm1
      double complex phase,bf,trigam
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zab2(i1,i2,i3,i4)=zab(i1,i2,i4)+zab(i1,i3,i4)

    
      virt_trigam_MHV=czip
      box_sum=czip
      bub_sum=czip
      rat_sum=czip
      comp_rat=czip
      l12=lnrat(musq,-s(i1,i2))
      l25=lnrat(musq,-s(i2,i5))

      A5LO=trigam(i1,i2,i3,i4,i5,za,zb)

!===== pole pieces
      Vcc=-(epinv**2+epinv*l12+0.5d0*l12**2)
     &     -(3d0/2d0*(epinv+l25+2d0))
      Vcc=Vcc*A5LO

    
!      write(6,*) 'epinv**2 pole'
!      write(6,*) -A5LO/phase*im 
!      write(6,*) 'epinv pole'
!      write(6,*) (-3d0/2d0-l12)*A5LO/phase*im     
!      write(6,*) 

!===============================================================
!==== BOXES
!===============================================================

      bxs12s35s45= (-(za(i1,i4)**3*za(i2,i3)*za(i3,i5)**2) + 
     &    za(i1,i3)**3*za(i2,i4)*za(i4,i5)**2)/
     &  (za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4)*za(i3,i4)**3)

!      write(6,*) -im*bxs12s35s45*s(i3,i5)*s(i4,i5)/2d0/A5LO
      
      bxs13s45s25=-((za(i1,i2)**2*za(i4,i5)**2)/
     &    (za(i1,i3)*za(i2,i4)**3*za(i3,i4)))

!      write(6,*) -im*bxs13s45s25*s(i2,i5)*s(i4,i5)/2d0/A5LO
      
      bxs14s35s25=(za(i1,i2)**2*za(i3,i5)**2)/
     &  (za(i1,i4)*za(i2,i3)**3*za(i3,i4))

!      write(6,*) -im*bxs14s35s25*s(i3,i5)*s(i2,i5)/2d0/A5LO
      
      bxs23s45s15=za(i1,i5)**2/(za(i1,i4)*za(i2,i3)*za(i3,i4))

!      write(6,*) -im*bxs23s45s15*s(i1,i5)*s(i4,i5)/2d0/A5LO
      
      bxs24s35s15=-(za(i1,i5)**2/(za(i1,i3)*za(i2,i4)*za(i3,i4)))

!      write(6,*) -im*bxs24s35s15*s(i3,i5)*s(i1,i5)/2d0/A5LO
      
      bxs45s13s12= -((za(i1,i2)**2*za(i3,i5)**2)/
     &    (za(i1,i3)*za(i2,i3)**2*za(i2,i4)*za(i3,i4)))

!      write(6,*) -im*bxs45s13s12*s(i3,i1)*s(i1,i2)/2d0/A5LO
      
      bxs35s12s24=-(za(i1,i5)**2/(za(i1,i3)*za(i2,i4)*za(i3,i4)))
      
!      write(6,*) -im*bxs35s12s24*s(i1,i2)*s(i2,i4)/2d0/A5LO
      
      bxs35s14s12= (za(i1,i2)**2*za(i4,i5)**2)/
     &  (za(i1,i4)*za(i2,i3)*za(i2,i4)**2*za(i3,i4))

!      write(6,*) -im*bxs35s14s12*s(i1,i2)*s(i1,i4)/2d0/A5LO
      
      bxs45s12s23=za(i1,i5)**2/(za(i1,i4)*za(i2,i3)*za(i3,i4))

!      write(6,*) -im*bxs45s12s23*s(i1,i2)*s(i2,i3)/2d0/A5LO
      
!===============================================================
!==== Box summation 
!===============================================================
      
      box_sum=box_sum
     & -bxs12s35s45*Lsm1(-s(i3,i5),-s(i1,i2),-s(i4,i5),-s(i1,i2))
     & -bxs13s45s25*Lsm1(-s(i4,i5),-s(i1,i3),-s(i2,i5),-s(i1,i3))
     & -bxs14s35s25*Lsm1(-s(i3,i5),-s(i1,i4),-s(i2,i5),-s(i1,i4))
     & -bxs23s45s15*Lsm1(-s(i4,i5),-s(i2,i3),-s(i1,i5),-s(i2,i3))
     & -bxs24s35s15*Lsm1(-s(i3,i5),-s(i2,i4),-s(i1,i5),-s(i2,i4))
     & -bxs35s14s12*Lsm1(-s(i1,i4),-s(i3,i5),-s(i1,i2),-s(i3,i5))
     & -bxs45s13s12*Lsm1(-s(i1,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))
     & -bxs45s12s23*Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))
     & -bxs35s12s24*Lsm1(-s(i1,i2),-s(i3,i5),-s(i2,i4),-s(i3,i5))
 
!===============================================================
!==== bubbles (UNCOMPLETED) v1 
!===============================================================
      c13=(za(i1,i2)*za(i2,i5)**2*zb(i3,i2))/
     &   (za(i2,i3)*za(i2,i4)**2*zab2(i2,i1,i3,i2)) - 
     &  (za(i1,i2)*za(i4,i5)**2*zb(i4,i3))/
     &   (za(i2,i4)**2*za(i3,i4)*zab2(i4,i1,i3,i4)) - 
     &  (za(i1,i3)*za(i4,i5)**2*zb(i4,i3)**2)/
     &   (2.*za(i2,i4)*za(i3,i4)*zab2(i4,i1,i3,i4)**2)
      
      c14=(za(i1,i2)*za(i2,i5)**2*zb(i4,i2))/
     &   (za(i2,i3)**2*za(i2,i4)*zab2(i2,i1,i4,i2)) - 
     &  (za(i1,i2)*za(i3,i5)**2*zb(i4,i3))/
     &   (za(i2,i3)**2*za(i3,i4)*zab2(i3,i1,i4,i3)) + 
     &  (za(i1,i4)*za(i3,i5)**2*zb(i4,i3)**2)/
     &   (2.*za(i2,i3)*za(i3,i4)*zab2(i3,i1,i4,i3)**2)

!==== note this c35 is not in correct format for completion 
      c35= -((za(i1,i2)**2*za(i2,i5)*za(i3,i5)*zb(i3,i2))/
     &    (za(i1,i4)*za(i2,i3)**2*za(i2,i4)*zab2(i2,i3,i5,i2))
     &    )
      c45= -((za(i1,i2)**2*za(i2,i5)*za(i4,i5)*zb(i4,i2))/
     &    (za(i1,i3)*za(i2,i3)*za(i2,i4)**2*zab2(i2,i4,i5,i2))
     &    )

!===== 3/2 tree is absorbed into poles pieces, so c25=-sum of rest
      c25=-c35-c45-c13-c14
!      write(6,*) 
!      write(6,*) -c14*im/A5LO
!      write(6,*) -c35*im/A5LO
!      write(6,*) -c13*im/A5LO
!      write(6,*) -c45*im/A5LO
!      write(6,*) -(c25+3d0/2d0*A5LO)*im/A5LO

!===============================================================
!==== bubbles summation 
!===============================================================
      
      bub_sum=
     &    -c13*lnrat(musq,-s(i1,i3))
     &    -c14*lnrat(musq,-s(i1,i4))
     &    -c45*lnrat(musq,-s(i4,i5))
     &    -c35*lnrat(musq,-s(i3,i5))
     &    -c25*lnrat(musq,-s(i2,i5))

!===============================================================
!==== rational 
!===============================================================

      rat_sum=(za(i1,i5)*(za(i2,i4)*za(i3,i5) + 
     &       za(i2,i3)*za(i4,i5))*zb(i4,i3))/
     &   (4d0*za(i2,i3)*za(i2,i4)*za(i2,i5)*za(i3,i4)*
     &     zb(i5,i2)) - ((-(za(i2,i3)*zb(i3,i2)) + 
     &       za(i2,i4)*zb(i4,i2))*zb(i4,i3))/
     &   (4d0*za(i2,i3)*za(i2,i4)*zb(i5,i1)*zb(i5,i2)) - 
     &  (zb(i2,i1)*(zb(i3,i2)*zb(i4,i1) + 
     &       zb(i3,i1)*zb(i4,i2))*zb(i4,i3))/
     &   (4d0*za(i3,i4)*zb(i3,i1)*zb(i4,i1)*zb(i5,i1)*
     &     zb(i5,i2))
      rat_sum=-rat_sum

!======= completion 
      comp_rat=-(za(i4,i5)**2*zb(i4,i3)**2*
     &    (za(i1,i3)*zb(i3,i1) + za(i2,i5)*zb(i5,i2)))/
     &  (4.*za(i2,i4)*za(i2,i5)*za(i3,i4)*zab2(i4,i1,i3,i4)*
     &    zb(i3,i1)*zb(i5,i2))
      
      comp_rat=comp_rat+(za(i3,i5)**2*zb(i4,i3)**2*
     &     (za(i1,i4)*zb(i4,i1) + za(i2,i5)*zb(i5,i2)))/
     &  (4.*za(i2,i3)*za(i2,i5)*za(i3,i4)*zab2(i3,i1,i4,i3)*
     &    zb(i4,i1)*zb(i5,i2))
      

!===============================================================
! TOTAL 
!===============================================================

      Fcc=box_sum+bub_sum+rat_sum+comp_rat
      virt_trigam_MHV=Vcc+Fcc
           
!      write(6,*) 'LO ',A5LO/phase
!      write(6,*) 'Pole pieces, epinv =',Vcc/phase,epinv
!      write(6,*) 'box_sum = ',box_sum/phase
!      write(6,*) 'bub_sum = ',bub_sum/phase
!      write(6,*) 'CR = ',comp_rat/phase
!      write(6,*) 'rat only',rat_sum/phase
!      write(6,*) 
!      write(6,*) 'rational sum = ',(rat_sum/phase+comp_rat/phase)*im
!      write(6,*) 'Total finite= ',virt_trigam_MHV/phase

 !     write(6,*) 'Check',im*(Vcc+box_sum+bub_sum+rat_sum+comp_rat)
 !    & /phase
 !     write(6,*) 'Check',(im*(Vcc+box_sum+bub_sum-0*rat_sum-0*comp_rat)
 !    & /phase-dcmplx(31.731336419133008d0,  24.557820072575367d0))/A5lo
     
  
      return 
      end 
      
