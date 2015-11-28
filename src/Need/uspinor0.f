      subroutine uspinor0(q,i,f)
!-----subroutine for massless spinor
!     Weyl representation
      implicit none
      include 'swapxz.f'
      integer i
      double complex p(4),q(4),f(4),fc,czip
      double precision  px,py
      parameter(czip=(0d0,0d0))
C     Translate from MCFM notation (E=q(4)), 
      if (swapxz) then
C performing the swap (x<->z),(y->-y)
      p(1)=q(4)
      p(2)=q(3)
      p(3)=-q(2)
      p(4)=q(1)
      else
      p(1)=q(4)
      p(2)=q(1)
      p(3)=q(2)
      p(4)=q(3)
      endif
      px=+dreal(p(2))
      py=+dreal(p(3))
c      pz=+dreal(p(4))

      fc=sqrt(p(1)+p(4))

      if (abs(fc).gt.1D-8) then 

      if (i.eq.+1) then 
        f(1)=+fc
        f(2)=dcmplx(px,py)/fc
        f(3)=czip
        f(4)=czip
      elseif (i.eq.-1) then 
        f(1)=czip
        f(2)=czip
        f(3)=+dcmplx(px,-py)/fc
        f(4)=-fc
      endif 
      else
      if (i.eq.+1) then 
        f(1)=czip
        f(2)=sqrt(p(1))
        f(3)=czip
        f(4)=czip
      elseif (i.eq.-1) then 
        f(1)=czip
        f(2)=czip
        f(3)=sqrt(p(1))
        f(4)=czip
      endif 
      endif 
         
      return
      end

