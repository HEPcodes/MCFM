       subroutine polarization(p,i,f)
       implicit none
       integer i,pol
       double complex p(4),f(4)
       double precision p0,px,py,pz,pv,ct,st,cphi,sphi

c--- for i=0 set the polarization vector equal to momentum vector,
c---   in order to perform gauge check:
       if (i .eq. 0) then
         do pol=1,4
         f(pol)=p(pol)
         enddo
	 return
       endif
       
       p0=dreal(p(4))
       px=dreal(p(1))
       py=dreal(p(2))
       pz=dreal(p(3))

       pv=dsqrt(dabs(p0**2))
       ct=pz/pv
       st=dsqrt(dabs(1d0-ct**2))

       if (st.lt.1d-8) then 
        cphi=1d0
        sphi=0d0
       else
        cphi= px/pv/st
        sphi= py/pv/st
       endif


c      the following if statement distinguishes between 
c      positive and negative energies
       if (p0.gt.0d0) then  
       pol=i
       else
       pol=-i
       endif


       f(4)=dcmplx(0d0,0d0)
       f(1)=dcmplx(ct*cphi/dsqrt(2d0),-pol*sphi/dsqrt(2d0))
       f(2)=dcmplx(ct*sphi/dsqrt(2d0), pol*cphi/dsqrt(2d0))
       f(3)=dcmplx(-st/dsqrt(2d0),0d0)

       return
       end

