      double precision function vol(s,n)
c     volume of massless phase space including all (2*pi)'s
      implicit none
      include 'constants.f'
      double precision s
      integer n,fac
      vol=twopi**(4-3*n)*(pi/two)**(n-1)*s**(n-2)
      vol=vol/dfloat(fac(n-1)*fac(n-2))
c      if (n .eq. 6) vol=vol/(5.d0*4.d0*3.d0*2.d0)/(4.d0*3.d0*2.d0)
c      if (n .eq. 5) vol=vol/(4.d0*3.d0*2.d0)     /(3.d0*2.d0)
c      if (n .eq. 4) vol=vol/(3.d0*2.d0)          /(2.d0)
c      if (n .eq. 3) vol=vol/(2.d0)               /(1.d0)
      return
      end
      
      integer function fac(n)
      integer j,n
      j=1
      fac=1
 10   continue
      if (j .le. n) then
      fac=fac*j
      j=j+1
      else
      return
      endif
      go to 10
      end
