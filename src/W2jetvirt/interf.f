      double precision function interf(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
C---hep-ph/9708239, Eqn 2.12
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprodx.f'
      integer j,j1,j2,j3,j4,j5,j6 
      double complex tree,loop,a6treeg,a61g,a63g
      double precision gw,fac
      character*9 st1(4),st2(4),st3(4)
      common/gw/gw     
      data st1/'q+g+g+qb-','q+g+g-qb-','q+g-g+qb-','q+g-g-qb-'/
      data st2/'q+g+g+qb-','q+g-g+qb-','q+g+g-qb-','q+g-g-qb-'/
      data st3/'q+qb-g+g+','q+qb-g+g-','q+qb-g-g+','q+qb-g-g-'/
      interf=0d0
      fac=8d0*gw**4*gsq**2
      do j=1,4
      tree=a6treeg(st1(j),j1,j2,j3,j4,j5,j6,za,zb)
      loop=+V*a61g(st1(j),j1,j2,j3,j4,j5,j6,za,zb)
     .       -a61g(st2(j),j1,j2,j3,j4,j5,j6,za,zb)
     .       +a63g(st3(j),j1,j2,j3,j4,j5,j6,za,zb)
      interf=interf+dble(Dconjg(tree)*loop)
      enddo
      interf=fac*interf
      return
      end
