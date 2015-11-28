************************************************************************
*     Author: J. M. Campbell                                           *
*     August, 1999  (updated April, 2001)                              *
*                                                                      *
*     Comments added October 15th 2001.                                *
*                                                                      *
*     Revised by R.K. Ellis, November 9th, 16th 2001.                  *
*                                                                      *
*     Routines which return various pieces of the integrated           *
*     subtraction terms, used in both _v and _z routines               *
************************************************************************

************************************************************************
*                                                                      *
*     The labelling of the routines is as follows:                     *
*     The collinear pair is assumed to be incoming,                    *
*     so a reversal has to be made for the final state cases           *
*                                                                      *
*              -------->------------>--------                          *
*                  j        /         i                                *
*                          /                                           *
*                         /                                            *
*                                                                      *
*                represented by {ii/if}_ij                             *
*                                                                      *
************************************************************************

***********************************************************************
*************************** INITIAL-INITIAL ***************************
***********************************************************************

***************************** Quark-Quark *****************************
      double precision function ii_qq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSBR
c-- Id,aqq=
c--  [delta(1-x)]*(epinv*(epinv-L)+1/2*L^2+3/2*epinv-[pi]^2/6)
c--  +(1-x)-(1+x)*(L+2*[ln(1-x)])-(1+x^2)*[ln(x)]/[1-x]
c--  +4*[ln(1-x)/(1-xp)]+2*L/[1-xp]

      
      if (vorz .eq. 1) then
        ii_qq=epinv*(epinv2-L)+0.5d0*L**2+1.5d0*epinv-pisqo6
     .   -epinv*1.5d0
        if (scheme .eq. 'msbr') then
           return
        elseif (scheme .eq. 'dred') then
           ii_qq=ii_qq-half
           return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        ii_qq=omx-(one+x)*(two*lomx+L-epinv)-(one+x**2)/omx*lx
        return
      endif
      
      ii_qq=two/omx*(two*lomx+L-epinv)
      
      return
      end

***************************** Quark-Gluon *****************************
      double precision function ii_qg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSBR
c-- Id,aqg=(1-2*x*(1-x))*(-[ln(x)]+L+2*[ln(1-x)])+2*x*(1-x)
      ii_qg=0d0
      if ((vorz .eq. 1) .or. (vorz .eq. 3)) return
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        ii_qg=(one-two*x*omx)*(two*lomx-lx+L-epinv)+two*x*omx
      endif
      return
      end
      
***************************** Gluon-Quark *****************************
      double precision function ii_gq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial quark-quark (--> gluon) antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSBR
c-- Id,agq=(1+(1-x)^2)/x*(-[ln(x)]+L+2*[ln(1-x)])+x

      
      ii_gq=0d0
      if ((vorz .eq. 1) .or. (vorz .eq. 3)) return
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        ii_gq=(one+omx**2)/x*(two*lomx-lx+L-epinv)+x
        return
      endif

      return
      end

***************************** Gluon-Gluon *****************************
      double precision function ii_gg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSBR
c-- Id,agg=(epinv*(epinv-L)+1/2*L^2+epinv*11/6-[pi]^2/6
c--  -nf/3/xn*epinv)*[delta(1-x)]
c--  -2*[ln(x)]/[1-x]
c--  +2*(-1+x*(1-x)+(1-x)/x)*(-[ln(x)]+L+2*[ln(1-x)])
c--  +(4*[ln(1-x)/(1-xp)]+2*L/[1-xp])
      
      if (vorz .eq. 1) then
        ii_gg=epinv*(epinv2-L)+half*L**2-pisqo6
     .       +(11d0-two*nf/xn)/6d0*epinv
     .       -(11d0-two*nf/xn)/6d0*epinv
        if (scheme .eq. 'msbr') then
        return
        elseif (scheme .eq. 'dred') then
        ii_gg=ii_gg-1d0/6d0
        return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      
      if (vorz .eq. 2) then
        lx=dlog(x)
        ii_gg=two*(omx/x+x*omx-one)*(two*lomx-lx+L-epinv)-two*lx/omx
        return
      endif
      
      ii_gg=two*(two*lomx+L-epinv)/omx
      
      return
      end

***********************************************************************
**************************** INITIAL-FINAL ****************************
***********************************************************************

***************************** Quark-Quark *****************************
      double precision function if_qq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,ltmx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- initial-final quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- MSBR
c-- Id,aqq=(epinv*(epinv-L)+1/2*L^2+3/2*epinv+[pi]^2/6)*[delta(1-x)]
c--  +(1-x-2/[1-x]*[ln(2-x)]
c--  -(1+x)*(L+[ln(1-x)])-(1+x^2)*[ln(x)]/[1-x]
c--  +4*[ln(1-x)/(1-xp)]+2*L/[1-xp]
      
      if (vorz .eq. 1) then
        if_qq=epinv*(epinv2-L)+half*L**2+1.5d0*epinv+pisqo6
     . -1.5d0*epinv
        if (scheme .eq. 'msbr') then
        return
        elseif (scheme .eq. 'dred') then
        if_qq=if_qq-half
        return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      
      if (vorz .eq. 2) then
        ltmx=dlog(two-x)
        lx=dlog(x)
        if_qq=omx-two*ltmx/omx-(one+x)*(lomx+L-epinv)-(1+x**2)/omx*lx
        return
      endif
      
      if_qq=two/omx*(two*lomx+L-epinv)
      
      return
      end

***************************** Quark-Gluon *****************************
      double precision function if_qg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- MSBR
c-- Id,aqg=(1-2*x*(1-x))*(L-[ln(x)]+[ln(1-x)])+2*x*(1-x)
      
      if_qg=0d0
      if ((vorz .eq. 1).or.(vorz .eq. 3)) return
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        if_qg=(one-two*x*omx)*(lomx-lx+L-epinv)+two*x*omx
      endif
      
      return
      end
      
***************************** Gluon-Quark *****************************
      double precision function if_gq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- MSBR
c-- Id,agq=(1+(1-x)^2)/x*(L-[ln(x)]+[ln(1-x)])+x
      
      if_gq=0d0
      if ((vorz .eq. 1).or.(vorz .eq. 3)) return
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        if_gq=(one+omx**2)/x*(lomx-lx+L-epinv)+x
      endif
      
      return
      end

***************************** Gluon-Gluon *****************************
      double precision function if_gg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,ltmx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- MSBR
c-- Id,agg=[delta(1-x)]*(
c--  epinv*(epinv-L)+1/2*L^2+11/6*epinv+[pi]^2/6-1/3*epinv*nf/xn)
c--  +2*(-1+(1-x)/x+x*(1-x))*(L-[ln(x)]+[ln(1-x)])
c--  -2*[ln(2-x)]/[1-x]-2*[ln(x)]/[1-x]
c--  +4*[ln(1-x)/(1-xp)]+2*L/[1-xp]
      
      if (vorz .eq. 1) then
        if_gg=epinv*(epinv2-L)+half*L**2+pisq/6d0
c     .       +(11d0-two*nf/xn)*epinv/6d0
        if (scheme .eq. 'msbr') then
        return
        elseif (scheme .eq. 'dred') then 
        if_gg=if_gg-1d0/6d0
        return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      
      if (vorz .eq. 2) then
        ltmx=dlog(two-x)
        lx=dlog(x)
        if_gg=two*((lomx-lx+L-epinv)*(omx/x+x*omx-one)-(ltmx+lx)/omx)
        return
      endif
      
       if_gg=two/omx*(two*lomx+L-epinv)
      
      return
      end

***********************************************************************
**************************** FINAL-INITIAL ****************************
***********************************************************************

***************************** Quark-Quark *****************************
      double precision function fi_qq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- final-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSbar
c-- Id,aqq=(epinv*(epinv-L)+1/2*L^2+3/2*(epinv-L)+7/2-[pi]^2/2)*[delta(1-x)]
c--  +2/[1-x]*[ln(2-x)]
c--  +(-2*[ln(1-x)/(1-xp)]-3/2/[1-xp])
      
      if (vorz .eq. 1) then
         fi_qq=epinv*(epinv2-L)+half*L**2+1.5d0*(epinv-L)
     .   +3.5d0-half*pisq
         if (scheme .eq. 'msbr') then
           return
         elseif (scheme .eq. 'dred') then
           fi_qq=fi_qq-half
         return
         endif
      endif
      
      omx=one-x
      
      if (vorz .eq. 2) then
        fi_qq=two*dlog(two-x)/omx
        return
      endif
      
      fi_qq=-(two*dlog(omx)+1.5d0)/omx
      
      return
      end

***************************** Gluon-Gluon *****************************
      double precision function fi_gg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSbar
c-- Id,aqg=(-2/3*(epinv-L)-10/9)*[delta(1-x)]
c--  +0
c--  +2/3/[1-xp]
c-- Id,agg=
c--  (2*epinv*(epinv-L)+L^2+(epinv-L)*11/3+67/9-[pi]^2)
c--  *[delta(1-x)]
c--  +4*[ln(2-x)]/[1-x]
c--  +2*(-2*[ln(1-x)/(1-xp)]-11/6/[1-xp])


      if (vorz .eq. 1) then
        fi_gg=two*epinv*(epinv2-L)+L**2+67d0/9d0-pisq+11d0*(epinv-L)/3d0
        fi_gg=fi_gg+dfloat(nf)/xn*(2d0/3d0*(-epinv+L)-10d0/9d0)
        if (scheme .eq. 'msbr') then
        return
        elseif (scheme .eq. 'dred') then
        fi_gg=fi_gg-dfloat(nf)/xn/3d0
        return
        endif
      endif
      
      omx=one-x
      
      if (vorz .eq. 2) then
        fi_gg=four*dlog(two-x)/omx
        return
      endif
      
      fi_gg=-(four*dlog(omx)+11d0/3d0-dfloat(nf)/xn*2d0/3d0)/omx
      return
      end




***************************** Quark-Gluon *****************************
c      double precision function fi_qg(x,L,vorz)
c      implicit none
c      integer vorz
c      double precision x,L,omx
c      include 'constants.f'
c      include 'epinv.f'
c      include 'epinv2.f'
c      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSbar
c--Id,aqg=(-2/3*(epinv-L)-10/9)*[delta(1-x)]
c-- +0
c-- +2/3/[1-xp]

      
c      if (vorz .eq. 1) then
c       fi_qg=2d0/3d0*(-epinv+L)-10d0/9d0
c       if (scheme .eq. 'msbr') then
c          return
c       elseif (scheme .eq. 'dred') then
c          fi_qg=fi_qg-1d0/3d0
c          return
c       endif
c      elseif (vorz .eq. 2) then
c        fi_qg=0d0
c      elseif (vorz .eq. 3) then
c        fi_qg=2d0/3d0/(one-x)
c      endif
c      return
c      end


***********************************************************************
***************************** FINAL-FINAL *****************************
***********************************************************************

***************************** Quark-Quark *****************************
      double precision function ff_qq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- final-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C --MSbar
c Id,aqq=epinv*(epinv-L)+1/2*L^2+3/2*(epinv-L)+5-[pi]^2/2
      
      ff_qq=0d0
      if (vorz .eq. 1) then
        ff_qq=epinv*(epinv2-L)+half*L**2+1.5d0*(epinv-L)+5d0-half*pisq
        if (scheme .eq. 'msbr') then
        return
        elseif (scheme .eq. 'dred') then
        ff_qq=ff_qq-half
        return
        endif
      endif
      return
      end

***************************** Quark-Gluon *****************************
c      double precision function ff_qg(x,L,vorz)
c      implicit none
c      integer vorz
c      double precision x,L
c      include 'constants.f'
c      include 'epinv.f'
c      include 'epinv2.f'
c      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
C --MSbar
c Id,aqg=-2/3*(epinv-L)-16/9
c      
c      ff_qg=0d0
c      if (vorz .eq. 1) then
c        ff_qg=-2d0/3d0*(epinv-L)-16d0/9d0
c        if (scheme .eq. 'msbr') then
c          return
c        elseif (scheme .eq. 'dred') then
c          ff_qg=ff_qg-1d0/3d0
c          return
c       endif
c      endif
c      return
c      end

***************************** Gluon-Gluon *****************************
      double precision function ff_gg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C --MSbar
c Id,aqg=-2/3*(epinv-L)-16/9
c Id,agg=2*epinv*(epinv-L)+L^2+11/3*(epinv-L)+100/9-[pi]^2
      
      ff_gg=0d0
      if (vorz .eq. 1) then
        ff_gg=two*epinv*(epinv2-L)+L**2+100d0/9d0-pisq
     .        +11d0/3d0*(epinv-L)
        ff_gg=ff_gg+dfloat(nf)/xn*(-2d0/3d0*(epinv-L)-16d0/9d0)
        if (scheme .eq. 'msbr') then
          return
        elseif (scheme .eq. 'dred') then
          ff_gg=ff_gg-dfloat(nf)/xn/3d0
          return
        endif
      endif
      return
      end

