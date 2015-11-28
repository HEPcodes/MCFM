      double complex function a61LRL(j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      include 'masses.f'
      integer j1,j2,j3,j4,j5,j6
      double complex a61
      double precision prop
      prop=s(5,6)/sqrt((s(5,6)-wmass**2)**2+(wmass*wwidth)**2)
c---Note interchange of za,zb to effect complex conjugation
      a61LRL=a61('pp',j1,j2,j3,j4,j5,j6,zb,za)*prop
      return
      end

