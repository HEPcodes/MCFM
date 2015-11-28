      double precision function r(p,i,j)
c----calculate the jets separation between p(i) and p(j)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),vs1,vs2,h1,h2,dely,delfi,pt1pt2,pipj
      integer i,j
      vs1=p(i,3)/dsqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)
      if (vs1 .ge. 1d0) then
        h1=0d0
      elseif (vs1 .le. -1d0) then
        h1=pi
      else
        h1=dacos(vs1)
      endif
      vs2=p(j,3)/dsqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      if (vs2 .ge. 1d0) then
        h2=0d0
      elseif (vs2 .le. -1d0) then
        h2=pi
      else
        h2=dacos(vs2)
      endif
      dely=-dlog(dabs(dtan(h1/two)/dtan(h2/two)))
      pt1pt2=dsqrt((p(i,1)**2+p(i,2)**2)*(p(j,1)**2+p(j,2)**2))
      pipj=(p(i,1)*p(j,1)+p(i,2)*p(j,2))/(pt1pt2)
      if (pipj .ge. 1d0) then
      delfi=0d0
      elseif (pipj .le. -1d0) then
      delfi=pi
      else
      delfi=dacos(pipj)
      endif
      r=dsqrt(dely**2+delfi**2)
      end
