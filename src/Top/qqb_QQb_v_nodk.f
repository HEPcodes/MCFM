      subroutine qqb_QQb_v_nodk(p,msq)
      implicit none
************************************************************************
*     Author: J.M. Campbell                                            *
*     July, 2002.                                                      *
*     calculate the matrix element squared                             *
*     for the process                                                  *
*                                                                      *
*     g(-p1) + g(-p2) --> t(p3)+t~(p4)                                 *
*                                                                      *
*     t and tbar are assumed to be on shell                            * 
************************************************************************
      include 'constants.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'virtexp.f'
      include 'scheme.f'
      integer j,k,n2,n3
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),facgg,facqqb,
     . wtgg,wtqqb,wtqbq,t1_pass,t2_pass,msq0(-nf:nf,-nf:nf)
      double precision 
     . gnodkbasisae0,gnodkbasisae1,gnodkbasisae2,
     . gnodkbasisbe0,gnodkbasisbe1,gnodkbasisbe2,
     . gnodkbasisce0,gnodkbasisce1,gnodkbasisce2,
     . gnodkbasis1,gnodkbasis2,gnodkbasis3,gnodkbasis4,
     . gnodkbasis5,gnodkbasis6,gnodkbasis7,gnodkbasis8,
     . gnodkbasis9,gnodkbasis10,
     . gnodkbasise0,gnodkbasise1,
     . gcoeffae1,gcoeffae2,
     . gcoeffbe1,gcoeffbe2,
     . gcoeffce1,gcoeffce2,
     . gcoeff1,gcoeff2,gcoeff3,gcoeff4,gcoeff5,
     . gcoeff6,gcoeff7,gcoeff8,gcoeff9,gcoeff10,
     . gcoeffe1,gcoeffe2,
     . qnodkbasis0,qnodkbasis1,qnodkbasis2,qnodkbasis3,
     . qnodkbasis4,qnodkbasis5,qnodkbasis6,
     . qnodkbasise0,qnodkbasise1,
     . qcoeff0,qcoeff1,qcoeff2,qcoeff3,qcoeff4,
     . qcoeff5,qcoeff6,
     . qcoeffe1,qcoeffe2
      double precision 
     . gbasisae0,gbasisae1,gbasisae2,
     . gbasisbe0,gbasisbe1,gbasisbe2,
     . gbasisce0,gbasisce1,gbasisce2,
     . gbasis1,gbasis2,gbasis3,gbasis4,
     . gbasis5,gbasis6,gbasis7,gbasis8,
     . gbasis9,gbasis10,
     . gbasise0,gbasise1,
     . gcoefae1,gcoefae2,
     . gcoefbe1,gcoefbe2,
     . gcoefce1,gcoefce2,
     . gcoef1,gcoef2,gcoef3,gcoef4,gcoef5,
     . gcoef6,gcoef7,gcoef8,gcoef9,gcoef10,
     . gcoefe1,gcoefe2,
     . qbasis0,qbasis1,qbasis2,qbasis3,
     . qbasis4,qbasis5,qbasis6,
     . qbasise0,qbasise1,
     . qcoef0,qcoef1,qcoef2,qcoef3,qcoef4,
     . qcoef5,qcoef6,
     . qcoefe1,qcoefe2
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3

c--- choose the scheme here
      scheme='dred'

c--- set all matrix elements to zero
      do j=-nf,nf
        do k=-nf,nf
          msq(j,k)=0d0
        enddo
      enddo

c--- calculate momenta
      call dotem(4,p,s)

c--- calculate t1 and t2
      t1_pass=-s(2,4)/s(1,2)
      t2_pass=-s(2,3)/s(1,2)

c--- fill the special functions
      call virtexp_fill(t1_pass,t2_pass,mass2,s(1,2))

c--- calculate the gg piece
      call fill_gcoeff(gcoef1,gcoef2,gcoef3,gcoef4,gcoef5,
     . gcoef6,gcoef7,gcoef8,gcoef9,gcoef10)
      call fill_gcoeffe(gcoefae1,gcoefbe1,gcoefce1,
     . gcoefae2,gcoefbe2,gcoefce2)
      call fill_gnodkbasis(gbasis1,gbasis2,gbasis3,gbasis4,
     . gbasis5,gbasis6,gbasis7,gbasis8,gbasis9,gbasis10)
      call fill_gnodkbasise(gbasisae0,gbasisae1,gbasisae2,
     . gbasisbe0,gbasisbe1,gbasisbe2,gbasisce0,gbasisce1,gbasisce2)

c--- calculate the gg piece
      if (scheme .eq. 'tH-V') then
c--- lowest order has O(ep) and O(ep**2) pieces to interfere with poles
      wtgg=
     . +gbasisae0*(gcoefae2*EPINV**2+gcoefae1*EPINV)
c     . +gbasisae1*(gcoefae2*EPINV+gcoefae1)
c     . +gbasisae2*gcoefae2

     . +gbasisbe0*(gcoefbe2*EPINV**2+gcoefbe1*EPINV)
c     . +gbasisbe1*(gcoefbe2*EPINV+gcoefbe1)

     . +gbasisce0*(gcoefce2*EPINV**2+gcoefce1*EPINV)
c     . +gbasisce1*(gcoefce2*EPINV+gcoefce1)
      else
c--- scheme is 'dred' and lowest order is calculated in 4 dimensions 
      wtgg=
     . +gbasisae0*(gcoefae2*EPINV**2+gcoefae1*EPINV)
     . +gbasisbe0*(gcoefbe2*EPINV**2+gcoefbe1*EPINV)
     . +gbasisce0*(gcoefce2*EPINV**2+gcoefce1*EPINV)
      endif

c--- finite pieces, independent of scheme      
      wtgg=wtgg
     . +gbasis1*gcoef1
     . +gbasis2*gcoef2
     . +gbasis3*gcoef3
     . +gbasis4*gcoef4
     . +gbasis5*gcoef5
     . +gbasis6*gcoef6
     . +gbasis7*gcoef7
     . +gbasis8*gcoef8
     . +gbasis9*gcoef9
     . +gbasis10*gcoef10

c--- calculate the qqb piece
      call fill_qcoeff(qcoefe1,qcoefe2,qcoef0,qcoef1,qcoef2,
     . qcoef3,qcoef4,qcoef5,qcoef6)
      call fill_qnodkbasis(1,2,qbasise0,qbasise1,qbasis0,qbasis1,
     . qbasis2,qbasis3,qbasis4,qbasis5,qbasis6)

c--- calculate the qqb piece
      if (scheme .eq. 'tH-V') then
c--- lowest order has O(ep) and O(ep**2) pieces to interfere with poles
      wtqqb=
     . +qbasise0*(qcoefe2*EPINV**2+qcoefe1*EPINV)
c     . +qbasise1*(qcoefe2*EPINV+qcoefe1)
      else
c--- scheme is 'dred' and lowest order is calculated in 4 dimensions 
      wtqqb=
     . +qbasise0*(qcoefe2*EPINV**2+qcoefe1*EPINV)
      endif

c--- finite pieces, independent of scheme     
      wtqqb=wtqqb
     . +qbasis0*qcoef0
     . +qbasis1*qcoef1
     . +qbasis2*qcoef2
     . +qbasis3*qcoef3
     . +qbasis4*qcoef4
     . +qbasis5*qcoef5
     . +qbasis6*qcoef6

c--- calculate the qbq piece
      call virtexp_fill(t2_pass,t1_pass,mass2,s(1,2))

      call fill_qcoeff(qcoefe1,qcoefe2,qcoef0,qcoef1,qcoef2,
     . qcoef3,qcoef4,qcoef5,qcoef6)
      call fill_qnodkbasis(2,1,qbasise0,qbasise1,qbasis0,qbasis1,
     . qbasis2,qbasis3,qbasis4,qbasis5,qbasis6)

c--- calculate the qbq piece
      if (scheme .eq. 'tH-V') then
c--- lowest order has O(ep) and O(ep**2) pieces to interfere with poles
      wtqbq=
     . +qbasise0*(qcoefe2*EPINV**2+qcoefe1*EPINV)
c     . +qbasise1*(qcoefe2*EPINV+qcoefe1)
      else
c--- scheme is 'dred' and lowest order is calculated in 4 dimensions 
      wtqbq=
     . +qbasise0*(qcoefe2*EPINV**2+qcoefe1*EPINV)
      endif

c--- finite pieces, independent of scheme     
      wtqbq=wtqbq
     . +qbasis0*qcoef0
     . +qbasis1*qcoef1
     . +qbasis2*qcoef2
     . +qbasis3*qcoef3
     . +qbasis4*qcoef4
     . +qbasis5*qcoef5
     . +qbasis6*qcoef6

c--- this factor comes from nodksquare.frm and dksplit.frm, combined
      facgg=gsq**2*ason2pi*V/192d0/b/xnsq
c--- this factor comes from fortify.frm and hfortify.frm, combined
      facqqb=V/4d0*gsq**2*ason2pi/32d0

      wtgg=avegg*facgg*wtgg
      wtqqb=aveqq*facqqb*wtqqb
      wtqbq=aveqq*facqqb*wtqbq
 
c--- extra finite pieces in the DR scheme
      if (scheme .eq. 'dred') then
        call qqb_QQb(p,msq0)
        wtqqb=wtqqb+cf*ason2pi*msq0(+1,-1)
        wtqbq=wtqbq+cf*ason2pi*msq0(-1,+1)
        wtgg=wtgg+xn/3d0*ason2pi*msq0(0,0)
      endif

c-- fill the matrix elements
      do j=-nf,nf
        k=-j
        if     (j .lt. 0) then
          msq(j,k)=wtqbq
        elseif (j .gt. 0) then
          msq(j,k)=wtqqb
        elseif (j .eq. 0) then
          msq(j,k)=wtgg
        endif
      enddo

      return
      end

      
      
      
