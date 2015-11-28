      subroutine qqb_QQbdk_old(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 2008.                                                      *
*     calculate the element squared and subtraction terms              *
*     for the process                                                  *
*                                                                      *
*     q(-p1) +qbar(-p2)=bbar(p6)+e-(p7)+nubar(p8)+nu(p3)+e+(p4)+b(p5)  * 
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'msq_cs.f'
      include 'process.f'
      integer j,k,nu,cs
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),q(mxpart,4)
      double precision s,c1,c2,c6,c8,s1t,s2t,s12,mt2,fac,qqb
      double complex  sprod,amp(2),prop,loab(2,2),loba(2,2),loqed(2,2)
      s(j,k)=2d0
     & *(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      do cs=0,2
      msq_cs(cs,j,k)=0d0
      enddo
      enddo
      enddo

      s1t=s(1,3)+s(1,4)+s(1,5)
      s2t=s(2,3)+s(2,4)+s(2,5)
      s12=s(1,2)
      prop=dcmplx(s(3,4)-wmass**2,wmass*wwidth)
     .    *dcmplx(s(7,8)-wmass**2,wmass*wwidth)
     .    *dcmplx(zip,mt*twidth)**2
      fac=V*gwsq**4*gsq**2/abs(prop)**2*s(5,3)*s(6,8)
C--include factor for hadronic decays
      if ((case .eq. 'tt_bbh') .or. (case .eq. 'tt_hdk')) fac=2d0*xn*fac

      mt2=mt**2
      c1=mt2/(s(3,4)+s(4,5))
      c2=mt2/(s(6,7)+s(7,8))
      c6=mt2/s1t
      c8=mt2/s2t
      do nu=1,4
      q(1,nu)=p(1,nu)    
      q(2,nu)=p(2,nu)    
      q(3,nu)=+(p(3,nu)+(1d0-c1)*p(4,nu)+p(5,nu))
      q(4,nu)=p(4,nu)    
      q(5,nu)=-(p(6,nu)+(1d0-c2)*p(7,nu)+p(8,nu))
      q(6,nu)=p(3,nu)+p(4,nu)+p(5,nu)-c6*p(1,nu)    
      q(7,nu)=p(7,nu)    
      q(8,nu)=p(3,nu)+p(4,nu)+p(5,nu)-c8*p(2,nu)    
      enddo

      call spinoru(8,q,za,zb)
      sprod=zb(4,3)*za(5,7)
      amp(1)=(sprod*za(3,2)*zb(1,5)+mt2*zb(4,1)*za(2,7))/s12
      amp(2)=(sprod*za(3,1)*zb(2,5)+mt2*zb(4,2)*za(1,7))/s12
      qqb=fac*aveqq*(abs(amp(1))**2+abs(amp(2))**2)

      call ggttww1(s1t,s2t,s12,loab,loba)

c--- note that filling of colour structures here looks unnatural:
c---    1 <--> loba , 2 <--> loab
c---  but this does correspond to the filling in qqb_QQb.f
      do j=1,2
      do k=1,2
      loqed(j,k)=loab(j,k)+loba(j,k)
      msq_cs(1,0,0)=msq_cs(1,0,0)+fac*avegg*xn*abs(loba(j,k))**2
      msq_cs(2,0,0)=msq_cs(2,0,0)+fac*avegg*xn*abs(loab(j,k))**2
      msq_cs(0,0,0)=msq_cs(0,0,0)-fac*avegg/xn*abs(loqed(j,k))**2
      enddo
      enddo
      
C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if ((j .lt. 0) .or. (j .gt. 0)) then
          msq(j,-j)=qqb
C Division of quark into color structures is arbitrary
          msq_cs(1,j,-j)=qqb/3d0 
          msq_cs(2,j,-j)=qqb/3d0 
          msq_cs(0,j,-j)=qqb/3d0 
      elseif (j .eq. 0) then
          msq(0,0)=msq_cs(1,0,0)+msq_cs(2,0,0)+msq_cs(0,0,0)
      endif
      enddo

      return
      end
