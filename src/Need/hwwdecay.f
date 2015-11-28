      subroutine hwwdecay(p,j1,j2,j3,j4,msq)
************************************************************************
*     Author: R.K. Ellis, September 2012                               *
*                                                                      *
*     matrix element squared for the process of                        *
*     Higgs decay  H --> WW --> ne(j1)+e+(j2)+e-(j3)+nu~(j4)           *
************************************************************************
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      integer j1,j2,j3,j4
      double precision p(mxpart,4),s1234,msq,dot,hdecay
      double precision s12,s13,s14,s23,s24,s34
      
      s12=2d0*dot(p,j1,j2)
      s13=2d0*dot(p,j1,j3)
      s14=2d0*dot(p,j1,j4)
      s23=2d0*dot(p,j2,j3)
      s24=2d0*dot(p,j2,j4)
      s34=2d0*dot(p,j3,j4)
      s1234=s12+s13+s14+s23+s24+s34

      hdecay=gwsq**3*wmass**2*s13*s24
      hdecay=hdecay
     .  /(((s12-wmass**2)**2+(wmass*wwidth)**2)
     .   *((s34-wmass**2)**2+(wmass*wwidth)**2))
      msq=hdecay
      
      return
      end
      
      
