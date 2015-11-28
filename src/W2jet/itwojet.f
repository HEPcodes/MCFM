      subroutine itwojet(p,nwz)      
      implicit none
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'     
      integer nwz
      double precision wmass,wwidth,zmass,zwidth
      double precision p(mxpart,4)
      double complex bwf
      common/vbmass/wmass,wwidth,zmass,zwidth
      common/bwf/bwf

c---calculate spinor and dot products

      call spinoru(6,p,za,zb)
    
C--- calculate the Breit-Wigner factor
      if (nwz .eq. 0) 
     & bwf=(s(3,4)-zmass**2+im*zmass*zwidth)/s(3,4)
c---initializes four quark matrix elements
      call initqqqq(za,zb)
c initializes two quark two gluon matrix elements
      call initqqgg(za,zb,nwz)
      return        
      end
