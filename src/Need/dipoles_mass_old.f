************************************************************************
*     Author: R.K. ELLIS                                               *
*     December 2002                                                    *
*                                                                      *
*     Routines which return various pieces of the integrated           *
*     massive subtraction terms, used in both _v and _z routines       *
*                                                                      *
*     These routines are for MASSIVE dipole subtractions               *
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
      double precision function ii_mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,mbar,Pqqreg,alfax
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--TH-V
c-- Id,aqq=
c--  [delta(1-x)]*(epinv*(epinv-L)+1/2*L^2+3/2*epinv-[pi]^2/6)
c--  +(1-x)-(1+x)*(L+2*[ln(1-x)])-(1+x^2)*[ln(x)]/[1-x]
c--  +4*[ln(1-x)/(1-xp)]+2*L/[1-xp]

      
      if (vorz .eq. 1) then
        ii_mqq=epinv*(epinv2-L)+0.5d0*L**2-pisqo6
        if (scheme .eq. 'tH-V') then
           return
        elseif (scheme .eq. 'dred') then
           ii_mqq=ii_mqq-half
           return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        Pqqreg=-one-x
        ii_mqq=omx+Pqqreg*(two*lomx+L-epinv)-(one+x**2)/omx*lx
        alfax=aii/omx
        if (alfax .lt. 1d0) ii_mqq=ii_mqq+(two/omx+Pqqreg)*dlog(alfax)
        return
      endif
      
      ii_mqq=two/omx*(two*lomx+L-epinv)
      
      return
      end

***************************** Quark-Gluon *****************************
      double precision function ii_mqg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pqgreg,alfax,mbar
      include 'constants.f'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--TH-V
c-- Id,aqg=(1-2*x*(1-x))*(-[ln(x)]+L+2*[ln(1-x)])+2*x*(1-x)
      ii_mqg=0d0
      if ((vorz .eq. 1) .or. (vorz .eq. 3)) return
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        Pqgreg=one-two*x*omx
        ii_mqg=Pqgreg*(two*lomx-lx+L-epinv)+two*x*omx
        alfax=aii/omx
        if (alfax .lt. 1d0) ii_mqg=ii_mqg+Pqgreg*dlog(alfax)
      endif
      return
      end
      
***************************** Gluon-Quark *****************************
      double precision function ii_mgq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pgqreg,alfax,mbar
      include 'constants.f'
      include 'epinv.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial quark-quark (--> gluon) antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--TH-V
c-- Id,agq=(1+(1-x)^2)/x*(-[ln(x)]+L+2*[ln(1-x)])+x

      
      ii_mgq=0d0
      if ((vorz .eq. 1) .or. (vorz .eq. 3)) return
      
      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)
      
      if (vorz .eq. 2) then
        Pgqreg=(one+omx**2)/x
        ii_mgq=Pgqreg*(two*lomx-lx+L-epinv)+x
        alfax=aii/omx
        if (alfax .lt. 1d0) ii_mgq=ii_mgq+Pgqreg*dlog(alfax)
        return
      endif

      return
      end

***************************** Gluon-Gluon *****************************
      double precision function ii_mgg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pggreg,alfax,mbar
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
c--- returns the integral of the subtraction term for an
c--- initial-initial gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--TH-V
c-- Id,agg=(epinv*(epinv-L)+1/2*L^2+epinv*11/6-[pi]^2/6
c--  -nf/3/xn*epinv)*[delta(1-x)]
c--  -2*[ln(x)]/[1-x]
c--  +2*(-1+x*(1-x)+(1-x)/x)*(-[ln(x)]+L+2*[ln(1-x)])
c--  +(4*[ln(1-x)/(1-xp)]+2*L/[1-xp])
      
      if (vorz .eq. 1) then
        ii_mgg=epinv*(epinv2-L)+half*L**2-pisqo6
        if (scheme .eq. 'tH-V') then
        return
        elseif (scheme .eq. 'dred') then
        ii_mgg=ii_mgg-1d0/6d0
        return
        endif
      endif
      
      omx=one-x
      lomx=dlog(omx)
      
      if (vorz .eq. 2) then
        Pggreg=omx/x+x*omx-one
        lx=dlog(x)
        ii_mgg=two*Pggreg*(two*lomx-lx+L-epinv)-two*lx/omx
        alfax=aii/omx
        if (alfax .lt. 1d0) ii_mgg=ii_mgg
     .   +two*(one/omx+Pggreg)*dlog(alfax)
        return
      endif
      
      ii_mgg=two*(two*lomx+L-epinv)/omx
      
      return
      end

***********************************************************************
**************************** INITIAL-FINAL ****************************
******************************Section 5.3******************************
***********************************************************************
***************************** Quark-Quark *****************************
      double precision function if_mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,omx,lx,lomx,Pqqreg,ddilog,zp
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'alfacut.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- initial-final quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     

c g zipifqq=colfac*(
c  del(omx)*((epinv+log(1+mbarsq))*(epinv-L)+1/2*L**2
c            +1/2*log(1+mbarsq)**2+2*dilog(1/(1+musq))-pisq/6)
c +Preg(q,q)/colfac*(-(epinv-L)+2*log(omx)-log(x)-log(x*mbarsq+omx))
c  +omx-2/omx*(log(x)+log(1+x*mbarsq+omx,1+mbarsq))
c +2/omxp*(-[epinv-L]+2*log(omx)-log(1+mbarsq))
c  )-ifqq;

      mbarsq=mbar**2
      if_mqq=0d0
      if (vorz .eq. 1) then
c         if_mqq=epinv*epinv2-epinv*L+0.5d0*L**2
c     .   +(epinv-L)*dlog(1d0+mbarsq)+0.5d0*dlog(1d0+mbarsq)**2 
c     .   +2d0*ddilog(1d0/(1d0+mbarsq))-pisq/6d0 
          if_mqq=(epinv+dlog(1d0+mbarsq))*(epinv-L)+0.5d0*L**2
     .     +0.5d0*dlog(1d0+mbarsq)**2
     .     +2d0*ddilog(1d0/(1d0+mbarsq))-pisqo6
          if (scheme .eq. 'tH-V') then
            return
          elseif (scheme .eq. 'dred') then
            if_mqq=if_mqq-half
            return
          endif
      endif
      omx=one-x
      lomx=dlog(omx)
      zp=omx/(omx+mbarsq)
      if (vorz .eq. 2) then
         Pqqreg=-(1d0+x)
         lx=dlog(x)
c         if_mqq=(2d0*lomx-lx-dlog(1d0-x+mbarsq)+L-epinv)*Pqqreg+Pqqpr
c     .   -two/omx*(dlog((2d0-x+mbarsq)/(1d0+mbarsq))+lx)
         if_mqq=Pqqreg*(-(epinv-L)+2d0*dlog(omx)-lx-dlog(x*mbarsq+omx))
     .   +omx-2d0/omx*(lx+dlog((1d0+x*mbarsq+omx)/(1d0+mbarsq)))
         if (aif .lt. zp)
     .   if_mqq=if_mqq-(two/omx*(dlog(zp*(omx+aif)/(aif*(omx+zp))))
     .   +Pqqreg*dlog(zp/aif))
         return
      elseif (vorz .eq. 3) then
c         if_mqq=two/omx*(two*lomx+L-epinv-dlog(1d0+mbarsq))
         if_mqq=2d0/omx*(-(epinv-L)+2d0*dlog(omx)-dlog(1d0+mbarsq))
      endif
      return
      end

***************************** Gluon-Gluon *****************************
      double precision function if_mgg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,mbar,mbarsq,Pggreg,ddilog,zp
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'alfacut.f'
      include 'scheme.f'

c g zipifgg=colfac*(
c  del(omx)*((epinv+log(1+mbarsq))*(epinv-L)+1/2*L**2
c            +1/2*log(1+mbarsq)**2+2*dilog(1/(1+musq))-pisq/6)
c +Preg(g,g)/colfac*(-(epinv-L)+2*log(omx)-log(x)-log(x*mbarsq+omx))
c  +2*mbarsq*log(x*mbarsq,x*mbarsq+omx)
c  -2/omx*(log(x)+log(1+x*mbarsq+omx,1+mbarsq))
c +2/omxp*(-[epinv-L]+2*log(omx)-log(1+mbarsq))
c  )-ifgg;

      mbarsq=mbar**2
      if_mgg=0d0
CDTS (5.88)
      if (vorz .eq. 1) then
c         if_mgg=epinv*epinv2-epinv*L+0.5d0*L**2
c     .   +2d0*ddilog(1d0/(1d0+mbarsq))-pisq/6d0
c     .   +(epinv-L)*(dlog(1d0+mbarsq))+0.5d0*(dlog(1d0+mbarsq))**2
         if_mgg=(epinv+dlog(1d0+mbarsq))*(epinv-L)+0.5d0*L**2
     .    +0.5d0*dlog(1d0+mbarsq)**2+2d0*ddilog(1d0/(1d0+mbarsq))-pisqo6
          if (scheme .eq. 'tH-V') then
            return
          elseif (scheme .eq. 'dred') then
            if_mgg=if_mgg-1d0/6d0
            return
          endif
         return
      endif
      omx=one-x
      lomx=dlog(omx)
      zp=omx/(omx+mbarsq)
      if (vorz .eq. 2) then
C regular contribution
         Pggreg=2d0*(omx/x-1d0+x*omx)
         lx=dlog(x)
c         if_mgg=(2d0*lomx-log(omx+mbarsq)+L-epinv-lx)*Pggreg
c     .   -two/omx*(log((two-x+mbarsq)/(1d0+mbarsq))+lx)
c     .   -2d0*mbarsq/x*log((omx+mbarsq)/mbarsq)
         if_mgg=Pggreg*(-(epinv-L)+2d0*dlog(omx)-lx-dlog(x*mbarsq+omx))
     .    +2d0*mbarsq*dlog(x*mbarsq/(x*mbarsq+omx))
     .    -2d0/omx*(lx+dlog((1d0+x*mbarsq+omx)/(1d0+mbarsq)))
c +Preg(g,g)/colfac*(-(epinv-L)+2*log(omx)-log(x)-log(x*mbarsq+omx))
c  +2*mbarsq*log(x*mbarsq,x*mbarsq+omx)
c  -2/omx*(log(x)+log(1+x*mbarsq+omx,1+mbarsq))
         if (aif .lt. zp)
     .   if_mgg=if_mgg-(two/omx*(dlog(zp*(omx+aif)/(aif*(omx+zp))))
     .   +Pggreg*dlog(zp/aif))
         return
      elseif (vorz .eq. 3) then
C plus contribution
c         if_mgg=two/omx*(two*lomx+L-epinv-log(1d0+mbarsq))
         if_mgg=2d0/omx*(-(epinv-L)+2d0*dlog(omx)-dlog(1d0+mbarsq))
         return 
      endif
      return
      end

***************************** Quark-Gluon *****************************
C--- Not necessary because for off-diagonal (no soft singularity)
C--- we always choose to use the initial spectator
c      double precision function if_mqg(x,L,mbar,vorz)
c      implicit none
c      integer vorz
c      double precision x,L,mbar,mbarsq,omx,lx,lomx,Pqgreg
c      include 'constants.f'
c      include 'epinv.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- TH-V
c-- Id,aqg=(1-2*x*(1-x))*(L-[ln(x)]+[ln(1-x)])+2*x*(1-x)
c      
c      if_mqg=0d0
c      if ((vorz .eq. 1).or.(vorz .eq. 3)) return
c      
c      mbarsq=mbar**2
c      omx=one-x
c      lomx=dlog(omx)
c      lx=dlog(x)
c      
c      if (vorz .eq. 2) then
c        Pqgreg=one-two*x*omx
c        if_mqg=Pqgreg*(lomx-lx+L-epinv+log(omx/(omx+mbarsq)))+two*x*omx
c      endif
c      
c      return
c      end
      
***************************** Gluon-Quark *****************************
C--- Not necessary because for off-diagonal (no soft singularity)
C--- we always choose to use the initial spectator
c      double precision function if_mgq(x,L,mbar,vorz)
c      implicit none
c      integer vorz
c      double precision x,L,mbar,mbarsq,omx,lx,lomx
c      include 'constants.f'
c      include 'epinv.f'
c--- returns the integral of the subtraction term for an
c--- initial-final gluon-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
c-- TH-V
c-- Id,agq=(1+(1-x)^2)/x*(L-[ln(x)]+[ln(1-x)])+x
c      
c      if_mgq=0d0
c      if ((vorz .eq. 1).or.(vorz .eq. 3)) return
      
c      omx=one-x
c      lomx=dlog(omx)
c      lx=dlog(x)
      
c      if (vorz .eq. 2) then
c        if_mgq=(one+omx**2)/x*(lomx-lx+L-epinv)+x
c      endif
      
c      return
c      end

***********************************************************************
**************************** FINAL-INITIAL ****************************
*****************************Section 5.2*******************************
***********************************************************************

***************************** Quark-Quark *****************************
      double precision function fi_mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,omx,ddilog,theta
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'alfacut.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- final-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSbar
c-- Id,aqq=(epinv*(epinv-L)+1/2*L^2+3/2*(epinv-L)+7/2-[pi]^2/2)*[delta(1-x)]
c--  +2/[1-x]*[ln(2-x)]
c--  +(-2*[ln(1-x)/(1-xp)]-3/2/[1-xp])
      
c g zipfiqq=colfac*(
c  del(omx)*((1+log(mbarsq,1+mbarsq))*(epinv-L)
c            +1/2*log(mbarsq)-1/2*log(mbarsq)**2
c            +1/2*log(1+mbarsq)+1/2*log(1+mbarsq)**2
c            -2*log(mbarsq)*log(1+mbarsq)
c            -4*dilog(-mbarsq)+mbarsq/2/(1+mbarsq)
c            +3/2-2/3*pisq)
c +2/omx*(log(1+x*mbarsq+omx,1+mbarsq))
c +omx/2/(x*mbarsq+omx)**2
c +2/omxp*(log(1+mbarsq,omx+x*mbarsq)-1)
c  )-fiqq;

c--- note that the modification for afi .ne. 1d0 is, at the moment,
c--- speculative (based on the massless case) and needs to be checked
      theta=0d0
      if (x .gt. 1d0-afi) theta=1d0      

      mbarsq=mbar**2
      if (vorz .eq. 1) then
CDST(5.59)
c      JaS=
c     .  +epinv*(epinv2-L)+half*L**2
c     .  -epinv*(epinv2-L)-half*L**2
c     .  +dlog(mbarsq)*(epinv-L)-0.5d0*dlog(mbarsq)**2
c     .  -0.5d0*(epinv-L-dlog(mbarsq))-pisq/6d0-2d0
c     . -(epinv-L)*(dlog(1d0+mbarsq))
c     . +1.5d0*(epinv-L)+3.5d0-half*pisq
CDST(5.60)
c      JaNS=0.5d0*dlog(1d0+mbarsq)**2-2d0*dlog(mbarsq)*dlog(1d0+mbarsq)
c     . -4d0*ddilog(-mbarsq)+0.5d0*dlog(1d0+mbarsq)
c     . +0.5d0*mbarsq/(1d0+mbarsq)
c      fi_mqq=JaS+JaNS
        fi_mqq=(1d0+dlog(mbarsq/(1d0+mbarsq)))*(epinv-L)
     .        +0.5d0*dlog(mbarsq)-0.5d0*dlog(mbarsq)**2
     .        +0.5d0*dlog(1d0+mbarsq)+0.5d0*dlog(1d0+mbarsq)**2
     .        -2d0*dlog(mbarsq)*dlog(1d0+mbarsq)
     .        -4d0*ddilog(-mbarsq)+mbarsq/2d0/(1d0+mbarsq)
     .        +1.5d0-2d0/3d0*pisq
         if (scheme .eq. 'tH-V') then
           return
         elseif (scheme .eq. 'dred') then
           fi_mqq=fi_mqq-half
           return
         endif
      endif
      
      omx=one-x
      
      if (vorz .eq. 2) then
c        fi_mqq=two*dlog(two-x)/omx
        fi_mqq=2d0/omx*(dlog((1d0+x*mbarsq+omx)/(1d0+mbarsq)))
     .   +omx/2d0/(x*mbarsq+omx)**2
        fi_mqq=fi_mqq*theta
        return
      endif
      
CDST(5.58)
      if (vorz .eq. 3) then
c      fi_mqq=-(two*dlog(omx)+1.5d0)/omx
        fi_mqq=theta*2d0/omx*(dlog((1d0+mbarsq)/(omx+x*mbarsq))-1d0)
      endif
      
      return
      end

***************************** Gluon-Gluon *****************************
      double precision function fi_mgg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,omx
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
        fi_mgg=two*epinv*(epinv2-L)+L**2+67d0/9d0-pisq
     .    +11d0*(epinv-L)/3d0
        if (scheme .eq. 'tH-V') then
        return
        elseif (scheme .eq. 'dred') then
        fi_mgg=fi_mgg-dfloat(nflav)/xn/3d0
        return
        endif
      endif
      
      omx=one-x
      
      if (vorz .eq. 2) then
        fi_mgg=four*dlog(two-x)/omx
        return
      endif
      
      fi_mgg=-(four*dlog(omx)+11d0/3d0)/omx
      return
      end




***************************** Quark-Gluon *****************************
      double precision function fi_mqg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,omx,rt,JaS,JaNS
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C--MSbar

      mbarsq=mbar**2
      if (1d0-4d0*mbarsq .lt. 0d0) then
      write(6,*) 'error in fi_mqg,(1d0-4d0*mbarsq .lt. 0d0)'
      stop
      else
      rt=dsqrt(1d0-4d0*mbarsq)
      endif
      if (vorz .eq. 1) then
CDTS 5.63
      JaS=-2d0/3d0*dlog(mbarsq)-10d0/9d0
CDTS 5.64
      JaNS=10d0/9d0*(1d0-rt)-8d0/9d0*mbarsq*rt
     . +4d0/3d0*dlog(0.5d0*(1d0+rt))
      fi_mqg=JaS+JaNS
       
      elseif (vorz .eq. 2) then
C regular
        fi_mqg=0d0
      elseif (vorz .eq. 3) then
C plus at x=x+
      omx=1d0-x
CDTS 5.62
      fi_mqg=2d0/3d0*(omx+2d0*mbarsq)/omx**2*sqrt(1d0-4d0*mbarsq/omx)
      endif
      return
      end


***********************************************************************
***************************** FINAL-FINAL *****************************
***********************************************************************

***************************** Quark-Quark *****************************
      double precision function ff_mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,Lro,ro,vtijk,arg,
     . Ieik,Icoll,ddilog

      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
C     mbarsq=mass**2/Qsq
C     L=log(Qsq/musq)
c--- returns the integral of the subtraction term for an
c--- final-initial quark-gluon antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)     
C --MSbar
c Id,aqq=epinv*(epinv-L)+1/2*L^2+3/2*(epinv-L)+5-[pi]^2/2

c g zipffqq=colfac*del(omx)*(
c  1/vtijk*((epinv-L)*Lro-1/2*Lro**2-2*Lro*log(1-4*mbarsq)
c           +4*dilog(-ro)-6*dilog(1-ro)+pisq/3)
c +epinv-L+3-2*log(1-2*mbar)+log(1-mbar)+1/2*log(mbarsq)
c  -2*mbar/(1-2*mbarsq)*(1-2*mbar)-mbar/(1-mbar)
c  -2*mbarsq/(1-2*mbarsq)*log(mbar/(1-mu))
c  )-ffqq;

      ff_mqq=0d0
      if (vorz .eq. 1) then
        mbarsq=mbar**2
        arg=1d0-4d0*mbarsq
        if (arg .lt. 0d0) then
            write(6,*) 'Threshold problem in ff_mqq'
            stop
        endif
        vtijk=dsqrt(arg)/(1d0-2d0*mbarsq)
        ro=dsqrt((1d0-vtijk)/(1d0+vtijk))
        Lro=dlog(ro)
        ff_mqq=1d0/vtijk*((epinv-L)*Lro-0.5d0*Lro**2-2d0*Lro*dlog(arg)
     .   +4d0*ddilog(-ro)-6d0*ddilog(1d0-ro)+pisq/3d0)
     .   +epinv-L+3d0
     .   -2d0*dlog(1d0-2d0*mbar)+dlog(1d0-mbar)+0.5d0*dlog(mbarsq)
     .   -2d0*mbar/(1d0-2d0*mbarsq)*(1d0-2d0*mbar)-mbar/(1d0-mbar)
     .   -2d0*mbarsq/(1d0-2d0*mbarsq)*dlog(mbar/(1d0-mbar))
c        Ieik=1d0/vtijk*(0.5d0*(epinv-L)*Lro-Lro*dlog(arg)
c     .  -dlog(roj)**2+pisq/6d0
c     .  +2d0*ddilog(-ro)-2d0*ddilog(1d0-ro)-ddilog(1d0-roj**2))
c        Icoll=(epinv-L)+0.5d0*Lmbarsq-2d0-2d0*dlog((1d0-mbar)**2-mbarsq)
c     . +dlog(1d0-mbar)-2d0*mbarsq/(1d0-2d0*mbarsq)*dlog(mbar/(1d0-mbar))
c     . +5d0-mbar/(1d0-mbar)-2d0*mbar*(1d0-2d0*mbar)/(1d0-2d0*mbarsq)
c        ff_mqq=2d0*Ieik+Icoll
        if (scheme .eq. 'tH-V') then
        return
        elseif (scheme .eq. 'dred') then
        ff_mqq=ff_mqq-half
        return
        endif
        return
      endif
      return
      end

***************************** Quark-Gluon *****************************
      double precision function ff_mqg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,ro,arg
C-----L=Log(Qsq/musq)
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
     
c--- returns the integral of the subtraction term for an
c--- final-initial gluon-quark antenna, either
c--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
C --MSbar
c Id,aqg=-2/3*(epinv-L)-16/9
c      
      ff_mqg=0d0
      if (vorz .eq. 1) then
        mbarsq=mbar**2
        arg=1d0-4d0*mbarsq
        if (arg .lt. 0d0) then
        write(6,*) 'Threshold problem in ff_mqg'
        stop
        else
        ro=dsqrt(arg)  
        ff_mqg=-2d0/3d0*(2d0*dlog(mbarsq)
     .   -2d0*dlog(0.5d0*(1d0+ro))+2d0/3d0*ro*(3d0+ro**2))
        return
        endif
      endif
      return
      end

***************************** Gluon-Gluon *****************************
      double precision function ff_mgg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,Ieik,Icoll
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
      
      ff_mgg=0d0
CDST Eq.5.32
      if (vorz .eq. 1) then
      Ieik=0.5d0*epinv*(epinv2-L)+0.5d0*L**2-pisq/4d0
      Icoll=11d0/6d0*(epinv-L)+50d0/9d0
CDST Eq.5.36 needs to be added too - need to fix this
      Icoll=Icoll+half*dfloat(nf)/xn*(-2d0/3d0*(epinv-L)-16d0/9d0)
      ff_mgg=2d0*(2d0*Ieik+Icoll)
        if (scheme .eq. 'tH-V') then
          return
        elseif (scheme .eq. 'dred') then
          ff_mgg=ff_mgg-dfloat(nflav)/xn/3d0
          return
        endif
      return
      endif
      return
      end

c
c * Now write these expressions in a neater form
c g zipffqq=colfac*del(omx)*(
c  1/vtijk*((epinv-L)*Lro-1/2*Lro**2-2*Lro*log(1-4*mbarsq)
c           +4*dilog(-ro)-6*dilog(1-ro)+pisq/3)
c +epinv-L+3-2*log(1-2*mbar)+log(1-mbar)+1/2*log(mbarsq)
c  -2*mbar/(1-2*mbarsq)*(1-2*mbar)-mbar/(1-mbar)
c  -2*mbarsq/(1-2*mbarsq)*log(mbar/(1-mu))
c  )-ffqq;

c g zipifqq=colfac*(
c  del(omx)*((epinv+log(1+mbarsq))*(epinv-L)+1/2*L**2
c            +1/2*log(1+mbarsq)**2+2*dilog(1/(1+musq))-pisq/6)
c +Preg(q,q)/colfac*(-(epinv-L)+2*log(omx)-log(x)-log(x*mbarsq+omx))
c  +omx-2/omx*(log(x)+log(1+x*mbarsq+omx,1+mbarsq))
c +2/omxp*(-[epinv-L]+2*log(omx)-log(1+mbarsq))
c  )-ifqq;

c g zipifgg=colfac*(
c  del(omx)*((epinv+log(1+mbarsq))*(epinv-L)+1/2*L**2
c            +1/2*log(1+mbarsq)**2+2*dilog(1/(1+musq))-pisq/6)
c +Preg(g,g)/colfac*(-(epinv-L)+2*log(omx)-log(x)-log(x*mbarsq+omx))
c  +2*mbarsq*log(x*mbarsq,x*mbarsq+omx)
c  -2/omx*(log(x)+log(1+x*mbarsq+omx,1+mbarsq))
c +2/omxp*(-[epinv-L]+2*log(omx)-log(1+mbarsq))
c  )-ifgg;

c g zipfiqq=colfac*(
c  del(omx)*((1+log(mbarsq,1+mbarsq))*(epinv-L)
c            +1/2*log(mbarsq)-1/2*log(mbarsq)**2
c            +1/2*log(1+mbarsq)+1/2*log(1+mbarsq)**2
c            -2*log(mbarsq)*log(1+mbarsq)
c            -4*dilog(-mbarsq)+mbarsq/2/(1+mbarsq)
c            +3/2-2/3*pisq)
c +2/omx*(log(1+x*mbarsq+omx,1+mbarsq))
c +omx/2/(x*mbarsq+omx)**2
c +2/omxp*(log(1+mbarsq,omx+x*mbarsq)-1)
c  )-fiqq;
