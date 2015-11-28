      subroutine qqb_ttb(q,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     October, 1998.                                                   *
*     calculate the element squared and subtraction terms              *
*     for the process                                                  *
c----My notation                                                       *
C     q(-p1) +qbar(-p2)=bbar(p6)+e-(p7)+nubar(p8)+b(p5)+nu(p3)+e+(p4)  * 
*     This routine is nothing more than a wrapper for                  *
*     Kleiss and Stirling                                              *
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'process.f'
      integer j,k,nu
      double precision msq(-nf:nf,-nf:nf),q(mxpart,4),plab(4,10)
      double precision RMT,RGT,RMW,RGW,RMB,RMTLO,RMTUP
      double precision GW_KS,GS_KS,wtqqb,wtqbq,wtgg,fac
      COMMON/COUPS/GW_KS,GS_KS
      COMMON/PARS/RMT,RGT,RMW,RGW,RMB,RMTLO,RMTUP
      COMMON/MOM/plab
      logical first
      data first/.true./       
c---- Fill common blocks for Kleiss and Stirling   
      if (first) then
         first=.false.
         rmw=wmass
         rgw=wwidth
         rmt=mt
         rgt=twidth
         rmb=mb
         rmtlo=100d0
         rmtup=250d0
         GS_KS=sqrt(gsq)
         GW_KS=sqrt(gwsq/8d0)
         write(6,*) 
         write(6,*) 'rmw',rmw
         write(6,*) 'rmt',rmt
         write(6,*) 'rgw',rgw
         write(6,*) 'rgt',rgt
         write(6,*) 'GS_KS',GS_KS
         write(6,*) 'GW_KS',GW_KS
         write(6,*) 
      endif

c---step one change momentum notation into the notation of Kleiss and Stirling
c----My notation
C   q(-p1) +qbar(-p2)=bbar(p6)+e-(p7)+nubar(p8)+b(p5)+nu(p3)+e+(p4)  * 
c----KS notation
c---qbar(p2)+q(p1) = bbar(p3)+e-(p4)+nubar(p5)+b(p6)+nu(p7)+e+(p8)
c---        
      do nu=1,4
      plab(nu,1)=-q(2,nu)
      plab(nu,2)=-q(1,nu)
      plab(nu,3)=q(6,nu)
      plab(nu,4)=q(7,nu)
      plab(nu,5)=q(8,nu)
      plab(nu,6)=q(5,nu)
      plab(nu,7)=q(3,nu)
      plab(nu,8)=q(4,nu)
      enddo

C-----Matrix elements of Kleiss and Stirling include all averaging
C-----but are for leptonic decays
      fac=1d0
      if (case .eq. 'tt_bbh') then
      fac=2d0*xn      
      endif
      call ttbbww(1,2,wtqqb)
      wtqbq=wtqqb
c      call ttbbww(2,1,wtqbq)
      call ggttww(wtgg)
c      write(6,*) 'wtqqb',wtqqb
c      write(6,*) 'wtqbq',wtqbq
c      write(6,*) 'wtgg',wtgg

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if (j .lt. 0) then
          msq(j,-j)=fac*wtqbq
      elseif (j .eq. 0) then
          msq(j,j)=fac*wtgg
      elseif (j .gt. 0) then
          msq(j,-j)=fac*wtqqb
      endif
c      write(6,*) j,msq(j,j)
      enddo
c      pause
      return
      end
