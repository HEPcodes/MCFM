      subroutine qqb_QQb_v_dk(p,msq)
      implicit none
************************************************************************
*     Author: J.M. Campbell                                            *
*     July, 2002.                                                      *
*     calculate the matrix element squared                             *
*     for the process                                                  *
*                                                                      *
*     g(-p1) + g(-p2) --> t(nu(p3)+e^+(p4)+b(p5))                      *
*                        +t~(b~(p6)+e^-(p7)+nu(p8))                    *
*                                                                      *
*     t and tbar are assumed to be on shell                            * 
************************************************************************
      include 'constants.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'virtexp.f'
      include 'scheme.f'
      integer j,k,nu
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),facgg,facqqb,
     . facdk,
     . wtgg,wtqqb,wtqbq,t1_pass,t2_pass,ps(mxpart,4),
     . s34,s35,s68,s78,propw,propt
      double precision 
     . gdkbasisae0,gdkbasisae1,gdkbasisae2,
     . gdkbasisbe0,gdkbasisbe1,gdkbasisbe2,
     . gdkbasisce0,gdkbasisce1,gdkbasisce2,
     . gdkbasis1,gdkbasis2,gdkbasis3,gdkbasis4,
     . gdkbasis5,gdkbasis6,gdkbasis7,gdkbasis8,
     . gdkbasis9,gdkbasis10,
     . gdkbasise0,gdkbasise1,
     . gcoeffae1,gcoeffae2,
     . gcoeffbe1,gcoeffbe2,
     . gcoeffce1,gcoeffce2,
     . gcoeff1,gcoeff2,gcoeff3,gcoeff4,gcoeff5,
     . gcoeff6,gcoeff7,gcoeff8,gcoeff9,gcoeff10,
     . gcoeffe1,gcoeffe2,
     . qdkbasis0,qdkbasis1,qdkbasis2,qdkbasis3,
     . qdkbasis4,qdkbasis5,qdkbasis6,
     . qdkbasise0,qdkbasise1,
     . qcoeff0,qcoeff1,qcoeff2,qcoeff3,qcoeff4,
     . qcoeff5,qcoeff6,
     . qcoeffe1,qcoeffe2

c--- set all matrix elements to zero
      do j=-nf,nf
        do k=-nf,nf
          msq(j,k)=0d0
        enddo
      enddo

c--- calculate dot-products that will be used in propagators later
      s34=2d0*(p(3,4)*p(4,4)-p(3,3)*p(4,3)-p(3,2)*p(4,2)-p(3,1)*p(4,1))
      s78=2d0*(p(7,4)*p(8,4)-p(7,3)*p(8,3)-p(7,2)*p(8,2)-p(7,1)*p(8,1))
      s35=2d0*(p(3,4)*p(5,4)-p(3,3)*p(5,3)-p(3,2)*p(5,2)-p(3,1)*p(5,1))
      s68=2d0*(p(6,4)*p(8,4)-p(6,3)*p(8,3)-p(6,2)*p(8,2)-p(6,1)*p(8,1))

c--- transform momenta p --> ps
      do nu=1,4
c--- t momentum
        ps(3,nu)=p(3,nu)+p(4,nu)+p(5,nu)
c--- t-bar momentum
        ps(4,nu)=p(6,nu)+p(7,nu)+p(8,nu)
c--- positron
        ps(5,nu)=p(4,nu)
c--- electron
        ps(6,nu)=p(7,nu)
      enddo
           
c--- calculate dot-products with the new momenta, ps
      call dotem(6,ps,s)


c--- these are the propagator factors
      propw=((s34-wmass**2)**2+(wmass*wwidth)**2)
     .     *((s78-wmass**2)**2+(wmass*wwidth)**2)
      propt=(mt*twidth)**4
      facdk=gwsq**4*s35*s68/propw/propt
c--- include the factors for the decay here: amplitude has (gw/rt2) for
c--- each vertex, but two factors of 2 from the Fierz identity --> gw^8
c--- lepton correlations are (<35> [68])^2 = s35*s68      


c--- calculate t1 and t2
      t1_pass=-s(2,4)/s(1,2)
      t2_pass=-s(2,3)/s(1,2)

c--- fill the special functions
      call virtexp_fill(t1_pass,t2_pass,mt,s(1,2))

c--- calculate the gg piece
      if (scheme .eq. 'tH-V') then
c--- lowest order has O(ep) and O(ep**2) pieces to interfere with poles
      wtgg=
     . +gdkbasisae0*(gcoeffae2()*EPINV**2+gcoeffae1()*EPINV)
     . +gdkbasisae1*(gcoeffae2()*EPINV+gcoeffae1())
     . +gdkbasisae2*gcoeffae2()

     . +gdkbasisbe0*(gcoeffbe2()*EPINV**2+gcoeffbe1()*EPINV)
     . +gdkbasisbe1*(gcoeffbe2()*EPINV+gcoeffbe1())

     . +gdkbasisce0*(gcoeffce2()*EPINV**2+gcoeffce1()*EPINV)
     . +gdkbasisce1*(gcoeffce2()*EPINV+gcoeffce1())
      else
c--- scheme is 'dred' and lowest order is calculated in 4 dimensions 
      wtgg=
     . +gdkbasisae0*(gcoeffae2()*EPINV**2+gcoeffae1()*EPINV)
     . +gdkbasisbe0*(gcoeffbe2()*EPINV**2+gcoeffbe1()*EPINV)
     . +gdkbasisce0*(gcoeffce2()*EPINV**2+gcoeffce1()*EPINV)
      endif

c--- finite pieces, independent of scheme      
      wtgg=wtgg
     . +gdkbasis1()*gcoeff1()
     . +gdkbasis2()*gcoeff2()
     . +gdkbasis3()*gcoeff3()
     . +gdkbasis4()*gcoeff4()
     . +gdkbasis5()*gcoeff5()
     . +gdkbasis6()*gcoeff6()
     . +gdkbasis7()*gcoeff7()
     . +gdkbasis8()*gcoeff8()
     . +gdkbasis9()*gcoeff9()
     . +gdkbasis10()*gcoeff10()

c--- calculate the qqb piece
      if (scheme .eq. 'tH-V') then
c--- lowest order has O(ep) and O(ep**2) pieces to interfere with poles
      wtqqb=
     . +qdkbasise0(1,2)*(qcoeffe2()*EPINV**2+qcoeffe1()*EPINV)
     . +qdkbasise1(1,2)*(qcoeffe2()*EPINV+qcoeffe1())
      else
c--- scheme is 'dred' and lowest order is calculated in 4 dimensions 
      wtqqb=
     . +qdkbasise0(1,2)*(qcoeffe2()*EPINV**2+qcoeffe1()*EPINV)
      endif

c--- finite pieces, independent of scheme      
      wtqqb=wtqqb
     . +qdkbasis0(1,2)*qcoeff0()
     . +qdkbasis1(1,2)*qcoeff1()
     . +qdkbasis2(1,2)*qcoeff2()
     . +qdkbasis3(1,2)*qcoeff3()
     . +qdkbasis4(1,2)*qcoeff4()
     . +qdkbasis5(1,2)*qcoeff5()
     . +qdkbasis6(1,2)*qcoeff6()

c--- calculate the qbq piece
      call virtexp_fill(t2_pass,t1_pass,mt,s(1,2))
      if (scheme .eq. 'tH-V') then
c--- lowest order has O(ep) and O(ep**2) pieces to interfere with poles
      wtqbq=
     . +qdkbasise0(2,1)*(qcoeffe2()*EPINV**2+qcoeffe1()*EPINV)
     . +qdkbasise1(2,1)*(qcoeffe2()*EPINV+qcoeffe1())
      else
c--- scheme is 'dred' and lowest order is calculated in 4 dimensions 
      wtqbq=
     . +qdkbasise0(2,1)*(qcoeffe2()*EPINV**2+qcoeffe1()*EPINV)
      endif

c--- finite pieces, independent of scheme
      wtqbq=wtqbq      
     . +qdkbasis0(2,1)*qcoeff0()
     . +qdkbasis1(2,1)*qcoeff1()
     . +qdkbasis2(2,1)*qcoeff2()
     . +qdkbasis3(2,1)*qcoeff3()
     . +qdkbasis4(2,1)*qcoeff4()
     . +qdkbasis5(2,1)*qcoeff5()
     . +qdkbasis6(2,1)*qcoeff6()

c--- this factor comes from nodksquare.frm and dksplit.frm, combined
      facgg=gsq**2*ason2pi*V/192d0/b/xnsq
c--- this factor comes from fortify.frm and hfortify.frm, combined
      facqqb=V/4d0*gsq**2*ason2pi/32d0

      wtgg=facdk*facgg*wtgg
      wtqqb=facdk*facqqb*wtqqb
      wtqbq=facdk*facqqb*wtqbq
      
c-- fill the matrix elements
 
      do j=-nf,nf
        k=-j
        if     (j .lt. 0) then
          msq(j,k)=aveqq*wtqbq
        elseif (j .gt. 0) then
          msq(j,k)=aveqq*wtqqb
        elseif (j .eq. 0) then
          msq(j,k)=avegg*wtgg
        endif
      enddo

      write(6,*) 'qqb_QQb_v_dk: msq(0,0) = ',msq(0,0)
      pause
      
      return
      end

      
      
