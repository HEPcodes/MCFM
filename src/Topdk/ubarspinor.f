c      PROCEDURES FOR POLARIZATION VECTORS
C------Weyl representation
c      subroutine for ubar spinor
      subroutine ubarspinor(p,m,i,f)
      implicit none
      integer i
      double precision m
      double complex p(4)
      double complex f(4),fc,fc2,czip
      double precision  p0,px,py,pz
      parameter(czip=(0d0,0d0))

      if (m .eq. 0d0) then
        call ubar0spinor(p,i,f)
        return
      endif

      p0=dreal(p(1))
      px=dreal(p(2))
      py=dreal(p(3))
      pz=dreal(p(4))

      fc2=dcmplx(p0+pz,0d0)
      fc=sqrt(fc2)
      if (i.eq.1) then 
      f(1)=+dcmplx(m)/fc
      f(2)=czip
      f(3)=fc
      f(4)=dcmplx(px,-py)/fc

      elseif (i.eq.-1) then 
      f(1)=dcmplx(px,py)/fc
      f(2)=-fc
      f(3)=czip
      f(4)=-dcmplx(m)/fc
      endif 

      return
      end


