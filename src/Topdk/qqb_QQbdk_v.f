      subroutine qqb_QQbdk_v(q,msq)
      implicit none
************************************************************************
*     Author: J. M. Campbell                                           *
*     April 2010.                                                      *
*     Virtual matrix elements squared for the process                  *
*                                                                      *
*     q(-p1) +qbar(-p2)=bbar(p6)+e-(p7)+nubar(p8)+nu(p3)+e+(p4)+b(p5)  * 
************************************************************************
c--- Computed using the helicity amplitudes from:
c---
c--- \bibitem{Korner:2002hy}
c--- J.~G.~Korner and Z.~Merebashvili,
c--- %``One-loop corrections to four-point functions with two external massive
c--- %fermions and two external massless partons,''
c--- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c--- [arXiv:hep-ph/0207054].
c--- %%CITATION = PHRVA,D66,054023;%%
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'scale.f'
      include 'scheme.f'
      include 'epinv.f'
      include 'eplog.f'
      include 'spinorsw.f'
      include 'process.f'
      integer j,k,nu
      double precision msq(-nf:nf,-nf:nf),q(mxpart,4),fac
      double precision q1(4),q2(4),q3(4),q4(4),q5(4),q6(4),q7(4),q8(4)
      double precision resqqb0,resqbq0,resgg0,resqqb1,resqbq1,resgg1
c      double complex  zk1(4),zk2(4),zk3(4),zk4(4)

c--- set up common variable to indicate KM routines should be used
c--- to calculate spinor strings
      spinorsw='KM'      
      
      scheme='dred'

      do nu=1,4
      q1(nu)=q(1,nu)
      q2(nu)=q(2,nu)
      q3(nu)=q(3,nu)
      q4(nu)=q(4,nu)
      q5(nu)=q(5,nu)
      q6(nu)=q(6,nu)
      q7(nu)=q(7,nu)
      q8(nu)=q(8,nu)
      enddo

c--- set up variables that include the expansion of the KM overall
c--- factor in KM Eq. (2.3):
c---     epinv --> eplog ,  epinv**2 -> epsqlog
      eplog=epinv+log(musq/mt**2)
      epsqlog=epinv**2+epinv*log(musq/mt**2)
     &       +log(musq/mt**2)**2/2d0+pisqo6
     
      call KMampsqqqbQQbdk(q1,q2,q3,q4,q5,q6,q7,q8,
     & mt,wmass,resqqb0,resqqb1)
      call KMampsqqqbQQbdk(q2,q1,q3,q4,q5,q6,q7,q8,
     & mt,wmass,resqbq0,resqbq1)
      call KMampsqggQQbdk (q1,q2,q3,q4,q5,q6,q7,q8,
     & mt,wmass,resgg0,resgg1)

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      fac=ason2pi*gsq**2
C--include factor for hadronic decays
      if (case .eq. 'tt_bbh') fac=2d0*xn

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if (j .lt. 0) then
          msq(j,-j)=fac*aveqq*resqbq1
c          msq(j,-j)=gsq**2*aveqq*resqbq0
      elseif (j .eq. 0) then
          msq(j,j)=fac*avegg*resgg1
c          msq(j,j)=gsq**2*avegg*resgg0
      elseif (j .gt. 0) then
          msq(j,-j)=fac*aveqq*resqqb1
c          msq(j,-j)=gsq**2*aveqq*resqqb0
      endif
      enddo
      
      return
      end
