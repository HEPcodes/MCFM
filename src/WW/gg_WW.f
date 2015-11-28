      subroutine gg_WW(p,msqgg)
      implicit none
c--- Author: J. M. Campbell, March 2011
c--- For now, work in the approximation of two massless isodoublets
c--- Box contributions are then complete
c--- Triangle (vector) pieces always vanish 
c--- Triangle (axial) pieces cancel for massless isodoublets

      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer h1,h2
      double precision p(mxpart,4),msqgg,fac
      double complex Avec(2,2)
      double complex a64v

c--- set up spinor products      
      call spinoru(6,p,za,zb)

c--- fill amplitudes
      Avec(2,2)=a64v('q+qb-g-g-',3,4,1,2,6,5,zb,za)*(-im)
      Avec(2,1)=a64v('q+qb-g-g+',3,4,1,2,6,5,zb,za)*(-im)
      Avec(1,2)=a64v('q+qb-g+g-',3,4,1,2,6,5,zb,za)*(-im)
      Avec(1,1)=a64v('q+qb-g+g+',3,4,1,2,6,5,zb,za)*(-im)

      msqgg=0d0
      do h1=1,2
      do h2=1,2
      msqgg=msqgg+cdabs(Avec(h1,h2))**2
      enddo
      enddo
      
c--- overall factor from diagrams
      fac=avegg*V*(2d0*gwsq*gsq/(16d0*pisq)*gwsq/2d0)**2
     & *s(3,4)**2/((s(3,4)-wmass**2)**2+(wwidth*wmass)**2)
     & *s(5,6)**2/((s(5,6)-wmass**2)**2+(wwidth*wmass)**2)
c--- two complete isodoublets
      fac=fac*(2d0)**2
      
      msqgg=msqgg*fac
      
      return
      end
      
      
