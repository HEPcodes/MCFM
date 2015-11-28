      double complex function atrLLL(j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
c---atrLLL is the amplitude for
c---q+(-p4)+Q+(-p2)+l+(-p5) ---> q-(p1)+Q-(p3)+l-(p6)
c---All outgoing particles are left-handed, 
c---(obtained from atree(pm) all right-handed) 
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      include 'masses.f'
      integer j1,j2,j3,j4,j5,j6
      double complex atree
      double precision prop
      prop=s(5,6)/sqrt((s(5,6)-wmass**2)**2+(wmass*wwidth)**2)
c---Note interchange of za,zb to effect complex conjugation
c---to get from all right-handed to all left-handed

      atrLLL=atree('pm',j1,j2,j3,j4,j5,j6,zb,za)*prop
c      write(6,*) 'atrLLL',atrLLL
      return
      end

