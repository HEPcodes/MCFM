      subroutine qqb_QQbdk_v(p,msq)
      implicit none
c----Matrix element for tt~ production
C----averaged over initial colours and spins
C***********************************************************************
*     Author: R.K. Ellis                                               *
*     July, 2002.                                                      *
*     calculate the element squared                                    *
*     for the process                                                  *
*----My notation                                                       *
*     g(-p1) +g(-p2)=t(nu(p3)+e^+(p4)+b(p5))                           *
*                   +t~(b~(p6)+e^-(p7)+nu(p8))                         *
*     t and tbar are assumed to be on shell                            * 
************************************************************************

      include 'constants.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision wtgg,wtqqb,wtqbq
     

C---set all matrix elements equal to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      wtqbq=0d0
      wtgg=0d0
      call fundk(p,wtqqb,wtqbq,wtgg)

      do j=-nf,nf
      k=-j      
      if (j .lt.0) then
      msq(j,k)=aveqq*wtqbq
      elseif (j.gt.0) then
      msq(j,k)=aveqq*wtqqb
      elseif (j.eq.0) then
      msq(0,0)=avegg*wtgg
      endif
      enddo
      return
      end
 
