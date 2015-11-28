      double precision function etarap(j,p)
      implicit none
C---returns the value of the pseudorapidity
      include 'constants.f'
      integer j
      double precision p(mxpart,4)
      etarap=dsqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      etarap=(etarap+p(j,3))/(etarap-p(j,3))
      if (etarap .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etarap=100d0
      else
      etarap=0.5d0*dlog(etarap)
      endif
      return
      end

      double precision function aetarap(j,p)
      implicit none
C---returns the absolute value of the pseudorapidity
      include 'constants.f'
      integer j
      double precision p(mxpart,4)
      aetarap=dsqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      aetarap=(aetarap+p(j,3))/(aetarap-p(j,3))
      if (aetarap .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      aetarap=100d0
      else
      aetarap=0.5d0*abs(dlog(aetarap))
      endif
      return
      end
 
      double precision function yrap(j,p)
      implicit none
C---returns the value of the rapidity
      include 'constants.f'
      integer j
      double precision p(mxpart,4)
      yrap=(p(j,4)+p(j,3))/(p(j,4)-p(j,3))
      if (yrap .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrap=100d0
      else
      yrap=0.5d0*dlog(yrap)
      endif
      return
      end

      double precision function ayrap(j,p)
      implicit none
C---returns the absolute value of the rapidity
      include 'constants.f'
      integer j
      double precision p(mxpart,4)
      ayrap=(p(j,4)+p(j,3))/(p(j,4)-p(j,3))
      if (ayrap .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      ayrap=100d0
      else
      ayrap=0.5d0*dabs(dlog(ayrap))
      endif
      return
      end
 
      double precision function pt(j,p)
      implicit none
      include 'constants.f'
      integer j
      double precision p(mxpart,4)
c--- This is the formula for pt
      pt=dsqrt(p(j,1)**2+p(j,2)**2)
c--- This is the formula for Et 
c      pt=dsqrt(p(j,1)**2+p(j,2)**2)
c     . *p(j,4)/dsqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      return
      end

      double precision function pttwo(j,k,p)
      implicit none
      include 'constants.f'
      integer j,k
      double precision p(mxpart,4)
      pttwo=dsqrt((p(j,1)+p(k,1))**2+(p(j,2)+p(k,2))**2)
      return
      end

c--- this is the rapidity of pair j,k
      double precision function yraptwo(j,k,p)
      implicit none
      include 'constants.f'
      integer j,k
      double precision p(mxpart,4)
      yraptwo=(p(j,4)+p(k,4)+p(j,3)+p(k,3))
     .       /(p(j,4)+p(k,4)-p(j,3)-p(k,3))
      if (yraptwo .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yraptwo=100d0
      else 
      yraptwo=0.5d0*dlog(yraptwo)
      endif
            
      return
      end

c--- this is the pseudo-rapidity
      double precision function etaraptwo(j,k,p)
      implicit none
      include 'constants.f'
      integer j,k
      double precision p(mxpart,4)
      
      etaraptwo=dsqrt((p(j,1)+p(k,1))**2+(p(j,2)+p(k,2))**2
     .               +(p(j,3)+p(k,3))**2)
      etaraptwo=(etaraptwo+p(j,3)+p(k,3))
     .         /(etaraptwo-p(j,3)-p(k,3))
      if (etaraptwo .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etaraptwo=100d0
      else 
      etaraptwo=0.5d0*dlog(etaraptwo)
      endif
      
      return
      end

