*
* ADJUSTED TO REMOVE ALL 1/eps**2 PIECES IF EPINV2 IN EPINV2.F = 0d0
* ADJUSTED TO REMOVE ALL 1/eps**2 PIECES IF EPINV2 IN EPINV2.F = 0d0
* ADJUSTED TO REMOVE ALL 1/eps**2 PIECES IF EPINV2 IN EPINV2.F = 0d0
*
************************************************************************
*     Author: J. M. Campbell                                           *
*     August, 1999  (updated April, 2001)                              *
*     Routines which return various pieces of the integrated           *
*     subtraction terms, used in both _v and _z routines               *
************************************************************************

C**********************************************



C  An explanation of the notation qg gq gg would be helpful



***********************************************************************
*************************** INITIAL-INITIAL ***************************
***********************************************************************

      double precision function ii_qg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      if (vorz .eq. 1) then
c--- old version
c        ii_qg=epinv*(epinv-L)+0.5d0*L**2+1.5d0*(epinv-L)-pisqo6
c--- new version (2/22/00)
        ii_qg=epinv*(epinv2-L)+0.5d0*L**2+1.5d0*epinv-half-pisqo6
        return
      endif
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        ii_qg=omx-(one+x**2)/omx*lx-(one+x)*(L+two*lomx)
        return
      endif
      
      ii_qg=two/omx*(L+two*lomx)
      
      return
      end

      double precision function ii_gq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      if (vorz .eq. 1) then
        ii_gq=0d0
        return
      endif
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        ii_gq=two*x*omx+(one-two*x*omx)*(two*lomx-lx+L)
        return
      endif
      
      ii_gq=0d0
      
      return
      end
      
      double precision function ii_gg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      if (vorz .eq. 1) then
        ii_gg=epinv*(epinv2-L)+half*L**2-pisqo6-1d0/6d0
     .       +(11d0-two*nf/xn)/6d0*(epinv-L)

        return
      endif
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        ii_gg=two*(two*lomx-lx+L)*(omx/x+x*omx)-two*x*lx/omx
     &       -two*(two*lomx+L)
        return
      endif
      
      ii_gg=two*(L+two*lomx)/omx
      
      return
      end

***********************************************************************
**************************** INITIAL-FINAL ****************************
***********************************************************************

      double precision function if_qg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,ltmx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- initial-final quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      if (vorz .eq. 1) then
c--- old version
c        if_qg=epinv*(epinv-L)+half*L**2+1.5d0*(epinv-L)+pisqo6
c--- new version (2/22/00)
        if_qg=epinv*(epinv2-L)+half*L**2+1.5d0*epinv-half+pisqo6
        return
      endif
      
      omx=one-x
      lomx=dlog(omx)
      ltmx=dlog(two-x)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        if_qg=omx-two*ltmx/omx-(one+x)*(L-lx+lomx)-two/omx*lx
        return
      endif
      
      if_qg=two/omx*(L+two*lomx)
      
      return
      end

      double precision function if_gq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      if (vorz .eq. 1) then
        if_gq=0d0
        return
      endif
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        if_gq=(one-two*x*omx)*(L-lx+lomx)+two*x*omx
        return
      endif
      
      if_gq=0d0
      
      return
      end
      
      double precision function if_gg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,ltmx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      if (vorz .eq. 1) then
        if_gg=epinv*(epinv2-L)+half*L**2+pisq/6d0-1d0/6d0
     .       +(11d0-two*nf/xn)/6d0*(epinv-L)
        return
      endif
      
      omx=one-x
      lomx=dlog(omx)
      ltmx=dlog(two-x)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        if_gg=two*(L-lx+lomx)*(omx/x+x*omx-one)-two*ltmx/omx-two/omx*lx
        return
      endif
      
      if_gg=two/omx*(L+two*lomx)
      
      return
      end

***********************************************************************
**************************** FINAL-INITIAL ****************************
***********************************************************************

      double precision function fi_qg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lomx,ltmx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- final-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      if (vorz .eq. 1) then
       fi_qg=epinv*(epinv2-L)+half*L**2+1.5d0*(epinv-L)-half*pisq+3d0
         return
      endif
      
      omx=one-x
      lomx=dlog(omx)
      ltmx=dlog(two-x)
      
      if (vorz .eq. 2) then
        fi_qg=two*ltmx/omx
        return
      endif
      
      fi_qg=-two*lomx/omx-1.5d0/omx
      
      return
      end

      double precision function fi_gg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lomx,ltmx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      if (vorz .eq. 1) then
        fi_gg=two*(epinv*(epinv2-L)+half*L**2)-pisq+64d0/9d0
     .        +(11d0-two*nf/xn)/3d0*(epinv-L)
        return
      endif
      
      omx=one-x
      lomx=dlog(omx)
      ltmx=dlog(two-x)
      
      if (vorz .eq. 2) then
        fi_gg=four*ltmx/omx
        return
      endif
      
      fi_gg=-(four*lomx+(11d0-two*nf/xn)/3d0)/omx
      
      return
      end


      double precision function fi_gq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lomx,ltmx
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      if (vorz .eq. 1) then
        fi_gq=2d0/3d0*(-epinv+L)-10d0/9d0
        return
      endif
            
      if (vorz .eq. 2) then
        fi_gq=0d0
        return
      endif
      
      fi_gq=2d0/3d0/(one-x)
      
      return
      end


C---------------Really not checked yet below this point. RKE 9/1/99
      double precision function ff_qg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- final-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      ff_qg=0d0
      if (vorz .eq. 1) then
       ff_qg=epinv*(epinv2-L)+half*L**2+1.5d0*(epinv-L)+4.5d0-half*pisq
      endif
      return
      end

      double precision function ff_gq(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      ff_gq=0d0
      if (vorz .eq. 1) then
        ff_gq=-2d0/3d0*(epinv-L)-16d0/9d0
      endif
      return
      end

      double precision function ff_gg(x,L,vorz)
      implicit none
      integer vorz
      double precision x,L
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
      
      ff_gg=0d0
      if (vorz .eq. 1) then
        ff_gg=two*(epinv*(epinv2-L)+half*L**2)+97d0/9d0-pisq
     .        +(11d0-two*nf/xn)/3d0*(epinv-L)
      endif
      return
      end

