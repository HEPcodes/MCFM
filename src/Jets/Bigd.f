      double precision function Bigd()
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'qcdcouple.f'
      integer j1(12),j2(12),j3(12),j4(12),j5(12),j,k1(10),k2(10)
      double precision xnum,xdenom,temp
      data j1/1,1,1, 1,1,1, 1,1,1, 1,1,1/
      data j2/2,2,2, 2,2,2, 3,3,3, 3,4,4/
      data j3/3,3,4, 4,5,5, 2,2,4, 5,2,3/
      data j4/4,5,3, 5,3,4, 4,5,2, 2,3,2/
      data j5/5,4,5, 3,4,3, 5,4,5, 4,5,5/

      data k1/1,1,1, 1,2,2, 2,3,3, 4/
      data k2/2,3,4, 5,3,4, 5,4,5, 5/

      xdenom=1d0
      xnum=0d0
      do j=1,10
      xnum=xnum+s(k1(j),k2(j))**4
      xdenom=xdenom*s(k1(j),k2(j))
      enddo
      temp=xnum/xdenom

      xnum=0d0
      do j=1,12
          xnum=xnum+s(j1(j),j2(j))*s(j2(j),j3(j))
     .             *s(j3(j),j4(j))*s(j4(j),j5(j))*s(j5(j),j1(j))
      enddo
     
      Bigd=4d0*gsq**3*xn**3*V*xnum*temp
      return
      end

