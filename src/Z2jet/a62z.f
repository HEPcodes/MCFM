      double complex function a62z(i1,i2,i3,j1,j2,j3,j4,j5,j6,za,zb)
c--- this function should be passed the (true) helicities of the particles,
c--- namely with i1 incoming and i2,i3 outgoing (1=left, 2=right)
c--- to compare with BDKW paper, first flip i1 to outgoing and then
c--- note that left corresponds to + and right to -
      implicit none
      include 'constants.f'
      include 'sprodx.f'
      include 'masses.f'
      integer i1,i2,i3,nhel
      integer j1,j2,j3,j4,j5,j6
      double complex a62

      nhel=(i1-1)+(i2-1)*2+(i3-1)*4

c--- note that reversing all helicities swaps za,zb
c--- and switching 3rd helicity swaps j5 and j6      
      if     (nhel .eq. 0) then
        a62z=+a62('pm',j1,j2,j3,j4,j5,j6,zb,za)
      elseif (nhel .eq. 1) then
        a62z=+a62('pp',j1,j2,j3,j4,j6,j5,za,zb)
      elseif (nhel .eq. 2) then
        a62z=+a62('pp',j1,j2,j3,j4,j5,j6,zb,za)
      elseif (nhel .eq. 3) then
        a62z=+a62('pm',j1,j2,j3,j4,j6,j5,za,zb)
      elseif (nhel .eq. 4) then
        a62z=+a62('pm',j1,j2,j3,j4,j6,j5,zb,za)      
      elseif (nhel .eq. 5) then
        a62z=+a62('pp',j1,j2,j3,j4,j5,j6,za,zb)
      elseif (nhel .eq. 6) then
        a62z=+a62('pp',j1,j2,j3,j4,j6,j5,zb,za)
      elseif (nhel .eq. 7) then
        a62z=+a62('pm',j1,j2,j3,j4,j5,j6,za,zb)
      else
        write(*,*) 'Unsupported helicity assignment'
        stop
      endif

      return
      end
     
