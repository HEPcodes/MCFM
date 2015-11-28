      subroutine gampsabc(p1,p2,p3,p4,t5,t6,gg_a,gg_b,gg_c)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer p1,p2,p3,p4,t5,t6
      double complex 
     . gg_a(2,2,2,2,2),gg_b(2,2,2,2,2),gg_c(2,2,2,2,2)
      double precision al5,al6,s125,s126,s34,s12,s15,s25,s16,s26

      al5=mb**2/s(p2,t5)
      s125=(1d0+al5)*s(p1,p2)+s(p1,t5)+s(p2,t5)
      al6=mb**2/s(p1,t6)
      s126=(1d0+al6)*s(p1,p2)+s(p1,t6)+s(p2,t6)
      s12=s(p1,p2)
      s34=s(p3,p4)
      s15=s(p1,t5)+al5*s(p1,p2)
      s25=s(p2,t5)
      s16=s(p1,t6)
      s26=s(p2,t6)+al6*s(p1,p2)

C   Notation is hz,h1,h2,h5,h6 

      include 'aa2m.f' 
      include 'bb2m.f' 
      include 'cc2m.f' 
 
      return
      end
