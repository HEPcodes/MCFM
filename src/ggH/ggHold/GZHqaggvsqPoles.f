      subroutine GZHqaggvsqPoles(p,sqres)
      implicit none 
      include 'constants.f' 
      include 'zprods_com.f'
      include 'sprods_com.f' 
      include 'scale.f' 
      include 'qcdcouple.f'
      include 'b0.f'
      double precision   p(4,5),sqres(-2:0)
      double precision   Hqagg,Hqagg_ab,Hqagg_ba,Hqagg_sym
      double precision   q(mxpart,4)
      integer i,j
      double precision   s12,s13,s14,s23,s24,s34
      double complex lnrat 


C --  reversal of GZ notation and restore dimension of the momenta 
      do i=1,4
         do j=1,4
            q(i,j) = p(j,i)*scale
         enddo
      enddo
      call spinoru(4,q,za,zb)

      call hqqgg(1,2,3,4,Hqagg,Hqagg_ab,Hqagg_ba,Hqagg_sym)

      s12=s(1,2)
      s13=s(1,3)
      s14=s(1,4)
      s23=s(2,3)
      s24=s(2,4)
      s34=s(3,4)
               
c--- JMC on 3/8/05 - removed factor of ason2pi which is restored
c--- in the wrapping routine
      sqres(-2) = (-3d0*xn+one/xn)*Hqagg
      sqres(-1) = (
     .     (-3d0*Cf+(11d0*xn/3d0-2d0/3d0*nf-1d0/xn*lnrat(-s12,musq)))
     .     *Hqagg
     .     +(xn*(lnrat(-s13,musq)+lnrat(-s34,musq)+lnrat(-s24,musq)))
     .       *Hqagg_ab 
     .     +(xn*(lnrat(-s14,musq)+lnrat(-s34,musq)+lnrat(-s23,musq)))
     .       *Hqagg_ba 
     .     +(xn*(-lnrat(-s12,musq)+lnrat(-s13,musq)+lnrat(-s23,musq)
     .       +lnrat(-s14,musq)+lnrat(-s24,musq)))*Hqagg_sym)
c--- debug: remove the b0 piece which is part of the renormalization
      sqres(-1) = (
     .     (-3d0*Cf+(-1d0/xn*lnrat(-s12,musq)))
     .     *Hqagg
     .     +(xn*(lnrat(-s13,musq)+lnrat(-s34,musq)+lnrat(-s24,musq)))
     .       *Hqagg_ab
     .     +(xn*(lnrat(-s14,musq)+lnrat(-s34,musq)+lnrat(-s23,musq)))
     .       *Hqagg_ba
     .     +(xn*(-lnrat(-s12,musq)+lnrat(-s13,musq)+lnrat(-s23,musq)
     .       +lnrat(-s14,musq)+lnrat(-s24,musq)))*Hqagg_sym)

c      write(6,*) 'Hqagg_ab ',Hqagg_ab
c      write(6,*) 'Hqagg_ba ',Hqagg_ba
c      write(6,*) 'Hqagg_sym',Hqagg_sym

      sqres(0)  = 0d0
c      write(*,*) 'GZHqaggvsqPoles: 1/e^2',sqres(-2)/ason2pi
c      write(*,*) 'GZHqaggvsqPoles: 1/e  ',sqres(-1)/ason2pi
      !write(*,*) 'Hqagg',Hqagg,Hqagg_ab,Hqagg_ba,Hqagg_sym

c--- Add on renormalization factor (JMC)
      sqres(-1)=sqres(-1)-2d0*b0*Hqagg

      return
      end
