      double complex function a6loops(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'constants.f'
      include 'sprodx.f'
c---  DKS Eq. 3.15
      integer j1,j2,j3,j4,j5,j6
      double complex a6loopa
      a6loops=a6loopa(j1,j2,j3,j4,j5,j6,za,zb)
     .       +a6loopa(j1,j2,j6,j5,j4,j3,za,zb)
      return
      end
	
      double complex function a6loopa(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
c---  DKS Eq. 2.10 for alpha = a
c---  note that (-i) included in A(alpha),A(tree,alpha)
c---  so no factor of (+i) in front of F(alpha)
      include 'constants.f'
      include 'dprodx.f'
      include 'sprodx.f'
      integer j1,j2,j3,j4,j5,j6
      double complex tree,Vpole,a6treea,fa
      tree=a6treea(j1,j2,j3,j4,j5,j6,za,zb)
      a6loopa=tree*Vpole(s(1,2))+fa(j1,j2,j3,j4,j5,j6,za,zb)
      return 
      end

      double complex function a6loopb(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
c---  DKS Eq. 2.10 for alpha = b
c---  note that (-i) included in A(alpha),A(tree,alpha)
c---  so no factor of (+i) in front of F(alpha)
      include 'constants.f'
      include 'dprodx.f'
      include 'sprodx.f'
      integer j1,j2,j3,j4,j5,j6
      double complex tree,Vpole,a6treeb
	
      tree=a6treeb(j1,j2,j3,j4,j5,j6,za,zb)	
      a6loopb=tree*Vpole(s(j1,j2))

      return 
      end

      double complex function Vpole(sij)
      implicit none
c---  DKS Eq. 2.12
      include 'epinv.f'
      include 'scale.f'
      double precision sij
      double complex Lnrat,xl12
	
      xl12=Lnrat(-sij,musq)

      Vpole=-epinv**2+epinv*(-1.5d0+xl12)
     .   -0.5d0*xl12**2+1.5d0*xl12-3.5d0

      return
	end
	      

	
