      subroutine gampsgh(p1,p2,p3,p4,t5,t6,gg_g,gg_h)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer p1,p2,p3,p4,t5,t6
      double complex 
     . gg_g(2,2,2,2,2),gg_h(2,2,2,2,2)
      double precision al5,al6,s125,s126,s34,s25,s16

      al5=mb**2/s(p2,t5)
      s125=(1d0+al5)*s(p1,p2)+s(p1,t5)+s(p2,t5)
      al6=mb**2/s(p1,t6)
      s126=(1d0+al6)*s(p1,p2)+s(p1,t6)+s(p2,t6)
      s34=s(p3,p4)
      s25=s(p2,t5)
      s16=s(p1,t6)

C   Notation is hz,h1,h2,h5,h6 
 
      include 'gg2m.f'
      include 'hh2m.f'
 
c      do hz=1,2      
c      do h1=1,2      
c      do h2=1,2      
c     do h5=1,2      
c      do h6=1,2      
c      write(6,*) 'gg_g(hz,h1,h2,h5,h6)+gg_h(hz,h2,h1,h5,h6)',
c     .            gg_g(hz,h1,h2,h5,h6)+gg_h(hz,h2,h1,h5,h6)     
c      enddo
c      enddo
c      enddo
c      enddo
c     enddo

c      pause
      return
      end
      
