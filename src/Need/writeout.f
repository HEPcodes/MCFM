      subroutine writeout(p)
      implicit none
      include 'constants.f'
      integer j,n
      double precision p(mxpart,4),dot,sum(4)
      write(6,*) 'In writeout'
      write(6,*) 'p1',p(1,1),p(1,2),p(1,3),p(1,4)
      write(6,*) 'p2',p(2,1),p(2,2),p(2,3),p(2,4)
      write(6,*) 'p3',p(3,1),p(3,2),p(3,3),p(3,4)
      write(6,*) 'p4',p(4,1),p(4,2),p(4,3),p(4,4)
      write(6,*) 'p5',p(5,1),p(5,2),p(5,3),p(5,4)
      write(6,*) 'p6',p(6,1),p(6,2),p(6,3),p(6,4)
      write(6,*) 'p7',p(7,1),p(7,2),p(7,3),p(7,4)
      write(6,*) 'p8',p(8,1),p(8,2),p(8,3),p(8,4)
      write(6,*) 'p9',p(9,1),p(9,2),p(9,3),p(9,4)
      write(6,*) 'p10',p(10,1),p(10,2),p(10,3),p(10,4)

      write(6,*) 's12',2d0*dot(p,1,2)
      write(6,*) 'sqrt(s34)',sqrt(2d0*dot(p,3,4))
      write(6,*) 's56',2d0*dot(p,5,6)
      write(6,*) 'sqrt(s345)', 
     .   sqrt(2d0*dot(p,3,4)+2d0*dot(p,3,5)+2d0*dot(p,4,5))
      write(6,*) 'sqrt(s3457)',
     .   sqrt(2d0*dot(p,3,4)+2d0*dot(p,3,5)+2d0*dot(p,3,7)
     .                                 +2d0*dot(p,4,5)+2d0*dot(p,4,7)
     .                                                +2d0*dot(p,5,7))

      do j=1,4
      sum(j)=p(1,j)+p(2,j)
      enddo
      do n=3,mxpart
      do j=1,4
      sum(j)=sum(j)+p(n,j)
      enddo
      enddo

      write(6,*) '     psum1',sum(1)
      write(6,*) '     psum2',sum(2)
      write(6,*) '     psum3',sum(3)
      write(6,*) '     psum4',sum(4)
c      do j=1,4
c      sum(j)=-p(1,j)-p(2,j)
c      enddo
c      do n=3,mxpart
c      do j=1,4
c      sum(j)=sum(j)+p(n,j)
c      enddo
c      enddo

c      write(6,*) '     msum1',sum(1)
c      write(6,*) '     msum2',sum(2)
c      write(6,*) '     msum3',sum(3)
c      write(6,*) '     msum4',sum(4)
      write(6,*) 'p1Dp1',dot(p,1,1)
      write(6,*) 'p2Dp2',dot(p,2,2)
      write(6,*) 'p3Dp3',dot(p,3,3)
      write(6,*) 'p4Dp4',dot(p,4,4)
      write(6,*) 'p5Dp5',dot(p,5,5)
      write(6,*) 'p6Dp6',dot(p,6,6)
      write(6,*) 'p7Dp7',dot(p,7,7)
      write(6,*) 'p8Dp8',dot(p,8,8)
      write(6,*) 'p9Dp9',dot(p,9,9)
      write(6,*) 'p10Dp10',dot(p,10,10)
      write(6,*)

      call flush(6)
      pause

      return
      end
