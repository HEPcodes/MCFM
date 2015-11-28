      subroutine A5NLO(j1,j2,j3,j4,j5,za,zb,A5LOm,A5NLOm)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1999.                                                      *
* Eqns (IV.3,IV.8,IV.5,IV.10) of hep-ph/9708239 of Bern,Dixon,Kosower  *
* Virtual terms are in units of                                        *
*  (as/4/pi)  (4 pi)^ep  Gamma(1+ep)*Gamma(1-ep)^2/Gamma(1-2*ep)       *
************************************************************************
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      integer j1,j2,j3,j4,j5
      double complex A5LOm,A5NLOm,A51,A52,xl15,xl25,xl12,xl34
      double complex fcc5,fsc5,Vcc,Vsc,Lnrat,L0,L1,Lsm1
      double complex Rst,Rsu,L0uq,L1uq,L0qu,L1qu,L0qs,L1qs,Ttemp
c      double complex Rtu

* As originally written these correspond to 
* 0 --> q_R(1)+qb_L(3)+g_R(2)+ebar_L(4)+e_R(5)
* with all RH couplings
* However we want it in our 
* standard form 
*       0--> qb_R(1)+q_L(2)++e_L(3)+ebar_R(4)+g_L(5)
* with all LH couplings 

* so we have made the changes 
*
*                    'q+g+qb-'
*                   (1 ---> 2)
*                   (2 ---> 5)
*                   (3 ---> 1)
*                   (4 ---> 4)
*                   (5 ---> 3)

*                    'q+qb-g+'
*                   (1 ---> 2)
*                   (2 ---> 1)
*                   (3 ---> 5)
*                   (4 ---> 4)
*                   (5 ---> 3)

*  and also exchanged za and zb.

C--- corresponds to (1V.1) times minus i
      A5LOm=-zb(j1,j4)**2/(zb(j2,j5)*zb(j5,j1)*zb(j4,j3))

      xl15=Lnrat(musq,-s(j2,j5))
      xl25=Lnrat(musq,-s(j1,j5))
      xl12=Lnrat(musq,-s(j2,j1))
      xl34=Lnrat(musq,-s(j3,j4))

      Rst=Lsm1(-s(j2,j1),-s(j3,j4),-s(j2,j5),-s(j3,j4))
c      Rtu=Lsm1(-s(j2,j5),-s(j3,j4),-s(j1,j5),-s(j3,j4))
      Rsu=Lsm1(-s(j2,j1),-s(j3,j4),-s(j1,j5),-s(j3,j4))
      L0uq=L0(-s(j5,j1),-s(j3,j4))
      L1uq=L1(-s(j5,j1),-s(j3,j4))
      L0qu=L0(-s(j3,j4),-s(j5,j1))
      L1qu=L1(-s(j3,j4),-s(j5,j1))
      L0qs=L0(-s(j3,j4),-s(j2,j1))
      L1qs=L1(-s(j3,j4),-s(j2,j1))
      Ttemp=zb(j1,j2)*za(j2,j3)*zb(j3,j4)


C--leading N
      Vcc=-four-two*(epinv+xl25)
     . -(epinv**2+xl15*epinv+half*xl15**2)
     . -(epinv**2+xl25*epinv+half*xl25**2)
      Vsc=one+half*(epinv+xl25)
      fcc5=-A5LOm*(Rst-two*Ttemp*L0uq/(zb(j1,j4)*s(j3,j4)))
      fsc5=Ttemp/(zb(j2,j5)*zb(j5,j1)*zb(j4,j3)*s(j3,j4))
     . *(zb(j1,j4)*L0uq+half*Ttemp*L1uq/s(j3,j4))
      A51=(Vcc+Vsc)*A5LOm+fcc5+fsc5

C--subleading N
C--reversed in sign to take account of different sign of A5LOm, Eqn (IV.6)
      Vcc=four+two*(epinv+xl34)
     . +(epinv**2+xl12*epinv+half*xl12**2)
      Vsc=-half*(one+epinv+xl34)

      fcc5=zb(j1,j4)/(zb(j1,j5)*zb(j5,j2)*zb(j4,j3))
     . *(-zb(j1,j4)*Rst
     .+(zb(j2,j1)*zb(j5,j4)-zb(j2,j4)*zb(j1,j5))*Rsu/zb(j5,j2)
     . +two*za(j2,j5)*zb(j2,j4)*zb(j5,j1)*L0uq/s(j3,j4))

      fsc5=zb(j2,j4)**2*zb(j1,j5)/(zb(j2,j5)**3*zb(j4,j3))*Rsu
     .-half*(zb(j4,j2)*za(j2,j5))**2*zb(j1,j5)*L1qu
     . /(zb(j2,j5)*zb(j4,j3)*s(j1,j5)**2)
     .+zb(j2,j4)**2*zb(j1,j5)*za(j5,j2)*L0qu
     ./(zb(j2,j5)**2*zb(j4,j3)*s(j1,j5))
     .+zb(j1,j2)*za(j2,j5)*zb(j4,j5)*za(j5,j3)*L1qs
     ./(zb(j2,j5)*s(j2,j1)**2)
     .-(zb(j1,j2)*za(j2,j5)*zb(j5,j4))*zb(j2,j4)
     . *L0qs/(zb(j2,j5)**2*zb(j4,j3)*s(j2,j1))
     . -half*za(j5,j3)/(za(j2,j1)*za(j1,j5)*zb(j2,j5)*za(j4,j3))
     . *(za(j2,j5)*za(j1,j3)+za(j1,j5)*za(j2,j3))
      A52=(Vcc+Vsc)*A5LOm+fcc5+fsc5

      A5NLOm=A51+A52/xnsq
  
      return
      end
