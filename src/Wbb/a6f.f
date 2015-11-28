      double complex function a6f(st,j1,j2,j3,j4,j5,j6,za,zb) 
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
      implicit none
c---Atreepm is the amplitude for
c---q-(-p4)+Q-(-p2)+l-(-p5) ---> q+(p1)+Q+(p3)+l+(p6)
c---All outgoing particles are right-handed
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      include 'epinv.f'
      include 'scale.f'
      integer j1,j2,j3,j4,j5,j6
      double complex atree,virt,Lnrat
      character*2 st 
      virt=epinv+Lnrat(musq,-s(j2,j3))+2d0
c---???continuation
      a6f=atree(st,j1,j2,j3,j4,j5,j6,za,zb)*virt 
      return
      end

 
