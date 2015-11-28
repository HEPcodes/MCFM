      subroutine gencol(x,xjac,xmin,emit,r)
      implicit none
c---Generate an x value and store it for later retrieval
      integer emit,lemit
      double precision x,xjac,xmin,xl,xljac,xlmin,r
      save xl,xljac,xlmin,lemit
      x=1-(1-xmin)*abs(1-2*r)
      xjac=2*(1-xmin)
      xl=x
      xljac=xjac
      xlmin=xmin
      lemit=emit
      return

      entry getcol(x,xjac,xmin,emit)
c---return the same values as last time
      x=xl
      xjac=xljac
      xmin=xlmin
      emit=lemit
      end
