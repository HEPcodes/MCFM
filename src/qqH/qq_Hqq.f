      subroutine qq_Hqq(p,msq)
      implicit none 
c---Matrix element squared averaged over initial colors and spins
c
c     q(-p1)+q(-p2) -->  H(p3,p4)+q(p5)+q(p6) 
c                           |
c                           |
c                           |
c                           ---> b(p3)+bbar(p4)
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ckm.f'
      include 'sprods_com.f'
      integer j,k,m,n,pn(-nf:nf),x1(2),x2(2)
      double precision p(mxpart,4),fac,s34
      double precision msq(-nf:nf,-nf:nf),hdecay,
     . msqll,msqlr,msqzzin,msqwzin,msqwl,
     . msxll,msxlr,msxzzin,msxwzin,msxwl
      double precision msqx(fn:nf,fn:nf,fn:nf,fn:nf)
      common/msq_all/msqx
      logical includeall
      
      data pn/-1,-2,-1,-2,-1,0,1,2,1,2,1/
      data x1/1,2/,x2/2,1/

c--- This flag decides whether or not to include all types of diagrams:
c---  FALSE --> only diagrams that look like WBF
c---  TRUE --> extra "associated" diagrams e.g. qqb -> W --> W(->qqb)H
c--- The comparison with Madgraph verified by JMC on 2/4/03, as follows:
c---    FALSE --> exact agreement
c---    TRUE ---> some differences at 10^-4 level, presumably due to width
      includeall=.false.

      call dotem(6,p,s)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      do m=-nf,nf
      do n=-nf,nf
      msqx(j,k,m,n)=0d0
      enddo
      enddo
      enddo
      enddo

      s34=s(3,4)+2d0*mb**2
      hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s34-4d0*mb**2) 
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      fac=0.25d0*gwsq**3*hdecay
C Color cancels, 0.25d0 is spin average

C q-q and qbar-qbar
      call msqpieces(1,2,5,6,msqll,msqlr,msqzzin,msqwzin,msqwl)
      call msqpieces(1,2,6,5,msxll,msxlr,msxzzin,msxwzin,msxwl)

      do j=1,nf
      do k=1,nf
      do m=1,nf
      do n=1,nf

      if ((j.eq.m) .and. (k.eq.n)) then     
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*( 
     . +msqll*((L(j)*L(k))**2+(R(j)*R(k))**2)
     . +msqlr*((L(j)*R(k))**2+(R(j)*L(k))**2))
      endif
      if ((j.eq.n) .and. (k.eq.m)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*( 
     . +msxll*((L(j)*L(k))**2+(R(j)*R(k))**2)
     . +msxlr*((L(j)*R(k))**2+(R(j)*L(k))**2))
      endif
      if ((j.eq.k) .and. (j.eq.m) .and. (j.eq.n)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*( 
     . +msqzzin*((L(j)*L(k))**2+(R(j)*R(k))**2))
      endif
      if (  (pn(j)+pn(m) .eq. +3) .and. (pn(k)+pn(n) .eq. +3)
     ..and. (pn(j)+pn(k) .eq. +3) ) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*Vsq(j,-m)*Vsq(k,-n)*( 
     . +msxwzin*L(j)*L(k)+msxwl)
      endif 
      if (  (pn(j)+pn(n) .eq. +3) .and. (pn(k)+pn(m) .eq. +3)
     ..and. (pn(j)+pn(k) .eq. +3) ) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*Vsq(j,-n)*Vsq(k,-m)*( 
     . +msqwzin*L(j)*L(k)+msqwl)
      endif 
      if (j .eq. k) then
        msqx(j,k,m,n)=msqx(j,k,m,n)/2d0
      endif

      msqx(-j,-k,-m,-n)=msqx(j,k,m,n)

      enddo
      enddo
      enddo
      enddo

C q-qb
      call msqpieces(1,6,5,2,msqll,msqlr,msqzzin,msqwzin,msqwl)
      call msqpieces(1,6,2,5,msxll,msxlr,msxzzin,msxwzin,msxwl)

      do j=1,nf
      do k=-nf,-1
      do m=1,nf
      do n=-nf,-1
      
      if ((j.eq.m) .and. (k.eq.n)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*( 
     . +msqll*((L(j)*L(-k))**2+(R(j)*R(-k))**2)
     . +msqlr*((L(j)*R(-k))**2+(R(j)*L(-k))**2))
       if (includeall) then       
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*Vsq(j,k)*Vsq(m,n)*( 
     . +msqwzin*L(j)*L(-k)+msqwl)
         if (j .eq. -k) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*( 
     . +msxll*((L(j)*L(j))**2+(R(j)*R(j))**2)
     . +msxlr*((L(j)*R(j))**2+(R(j)*L(j))**2)
     . +msqzzin*((L(j)*L(j))**2+(R(j)*R(j))**2))
         endif
       endif
      endif
      
      if (   (pn(j)+pn(m) .eq. +3) .and. (pn(k)+pn(n) .eq. -3)
     . .and. (pn(j)+pn(k) .eq. pn(m)+pn(n)) ) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*Vsq(j,-m)*Vsq(k,-n)*msxwl
       if (includeall) then
         if ((j .eq. -k) .and. (m .eq. -n)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*( 
     . +msxll*((L(j)*L(m))**2+(R(j)*R(m))**2)
     . +msxlr*((L(j)*R(m))**2+(R(j)*L(m))**2)
     . +msxwzin*L(j)*L(m)*Vsq(j,-m)*Vsq(k,-n)
     . )
         endif
       endif
      endif
      
      msqx(-j,-k,-m,-n)=msqx(j,k,m,n)

      enddo
      enddo
      enddo
      enddo

C qb-q
      call msqpieces(6,2,1,5,msqll,msqlr,msqzzin,msqwzin,msqwl)
      call msqpieces(6,2,5,1,msxll,msxlr,msxzzin,msxwzin,msxwl)

      do j=-nf,1
      do k=1,nf
      do m=1,nf
      do n=-nf,-1
      
      if ((j.eq.n) .and. (k.eq.m)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*( 
     . +msqll*((L(-j)*L(k))**2+(R(-j)*R(k))**2)
     . +msqlr*((L(-j)*R(k))**2+(R(-j)*L(k))**2))
       if (includeall) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*Vsq(j,k)*Vsq(m,n)*( 
     . +msqwzin*L(-j)*L(k)+msqwl)
         if (j .eq. -k) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*( 
     . +msxll*((L(-j)*L(-j))**2+(R(-j)*R(-j))**2)
     . +msxlr*((L(-j)*R(-j))**2+(R(-j)*L(-j))**2)
     . +msqzzin*((L(-j)*L(-j))**2+(R(-j)*R(-j))**2))
         endif
       endif
      endif
      
      if (   (pn(j)+pn(n) .eq. -3) .and. (pn(k)+pn(m) .eq. +3)
     . .and. (pn(j)+pn(k) .eq. pn(m)+pn(n)) ) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*Vsq(j,-n)*Vsq(k,-m)*msxwl
       if (includeall) then
         if ((j .eq. -k) .and. (m .eq. -n)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+fac*( 
     . +msxll*((L(-j)*L(m))**2+(R(-j)*R(m))**2)
     . +msxlr*((L(-j)*R(m))**2+(R(-j)*L(m))**2)
     . +msxwzin*L(-j)*L(m)*Vsq(j,-n)*Vsq(k,-m)
     . )
         endif
       endif
      endif
      
      msqx(-j,-k,-m,-n)=msqx(j,k,m,n)

      enddo
      enddo
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf
      do m=-nf,nf
      do n=-nf,nf
c--- add in to total for msq(j,k)
c        if     ((j .gt. 0) .and. (k .lt. 0)) then
c          if (m.ge.n) msq(j,k) = msq(j,k)+msqx(j,k,m,n)
c        elseif ((j .lt. 0) .and. (k .gt. 0)) then
c          if (m.le.n) msq(j,k) = msq(j,k)+msqx(j,k,m,n)
c        elseif ((j .gt. 0) .and. (k .gt. 0)) then
c          if ((pn(j).eq.pn(n)) .and. (pn(k).eq.pn(m)))
c     .                msq(j,k) = msq(j,k)+msqx(j,k,m,n)
c        elseif ((j .lt. 0) .and. (k .lt. 0)) then
c          if ((pn(j).eq.pn(n)) .and. (pn(k).eq.pn(m)))
c     .                msq(j,k) = msq(j,k)+msqx(j,k,m,n)
c        endif 
        if (m.ne.n) then
          msq(j,k) = msq(j,k)+msqx(j,k,m,n)/2d0
        else
          msq(j,k) = msq(j,k)+msqx(j,k,m,n)
        endif 
c        if ((j.eq.m) .and. (k .eq. n)) then
c          msq(j,k) = msq(j,k)+msqx(j,k,m,n)
c        endif
      enddo
      enddo
      enddo
      enddo

      return
      end


      subroutine msqpieces(i1,i2,i5,i6,zll,zLR,zzLL,wzLL,wll)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      double precision zll,zLR,zzLL,wzll,wll,htheta
      double precision propw,propz,x
      integer i1,i2,i5,i6
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+sign(half,x)
      propw(i1,i2)=sign(one,(s(i1,i2)-wmass**2))*dsqrt(
     .((s(i1,i2)-wmass**2)**2+htheta(s(i1,i2))*(wmass*wwidth)**2)/wmass)
      propz(i1,i2)=sign(one,(s(i1,i2)-zmass**2))
     . *dsqrt(dsqrt(1d0-xw)/xw/2d0/zmass
     . *((s(i1,i2)-zmass**2)**2+htheta(s(i1,i2))*(zmass*zwidth)**2))
      zll=s(i1,i2)*s(i5,i6)/(propz(i1,i5)*propz(i2,i6))**2
      zlr=s(i1,i6)*s(i2,i5)/(propz(i1,i5)*propz(i2,i6))**2
      zzll=2d0/xn*s(i1,i2)*s(i5,i6)
     . /(propz(i1,i5)*propz(i2,i6)*propz(i1,i6)*propz(i2,i5))
      wzll=2d0/xn*s(i1,i2)*s(i5,i6)
     . /(propz(i1,i5)*propz(i2,i6)*propw(i1,i6)*propw(i2,i5))
      wll=s(i1,i2)*s(i5,i6)/(propw(i1,i6)*propw(i2,i5))**2
      return
      end
