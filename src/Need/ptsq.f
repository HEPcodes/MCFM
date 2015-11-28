      double precision function ptsq(j,p)
      implicit none
      include 'constants.f'
      integer j
      double precision p(mxpart,4),dot
      ptsq=two*dot(p,1,j)*dot(p,2,j)/dot(p,1,2)
      return
      end
