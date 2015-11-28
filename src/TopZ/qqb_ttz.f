      subroutine qqb_ttz(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     December, 1999.                                                  *
*     calculate the element squared and subtraction terms              *
*     for the process                                                  *
c----My notation                                                       *
C     q(-p1) +qbar(-p2)=+nu(p3)+e+(p4)+b(p5)+bbar(p6)+e-(p7)+nubar(p8)
C     +Z(f(p9)+fb(p10))
C  
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'prods.f'
      
      integer j,k,nu
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),
     . pw1(4),pw2(4),p12(4),q(4),a(4),h(4),rr(4),bb(4),
     . w(4),z(4),ss(4),d1,d2,d3,d4,s12,sh,sw1,sw2,qDq,aDa,rDr,bDb,zDz,
     . wDw,ssDss,wtgg,qqbup,qqbdo,qbqup,qbqdo,densq,fac
      double complex T1(2,2),T2(2,2),T3(2,2),T4(2,2),
     . sumup(2,2),sumdo(2,2)


c--initialize to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do nu=1,4
      h(nu)=p(9,nu)+p(10,nu)
      p12(nu)=p(1,nu)+p(2,nu)
      pw1(nu)=p(3,nu)+p(4,nu)
      pw2(nu)=p(7,nu)+p(8,nu)
      q(nu)=+p(3,nu)+p(4,nu)+p(5,nu)
      rr(nu)=q(nu)+h(nu)
      a(nu)=-p(6,nu)-p(7,nu)-p(8,nu)
      bb(nu)=a(nu)-h(nu)
      ss(nu)=p12(nu)+h(nu)
      w(nu)=+p(2,nu)+h(nu)
      z(nu)=-p(1,nu)-h(nu)
      enddo      


      s12=(p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2)
      ssDss=(ss(4)**2-ss(1)**2-ss(2)**2-ss(3)**2)
      sh=(h(4)**2-h(1)**2-h(2)**2-h(3)**2)
      sw1=(pw1(4)**2-pw1(1)**2-pw1(2)**2-pw1(3)**2)
      sw2=(pw2(4)**2-pw2(1)**2-pw2(2)**2-pw2(3)**2)
      qDq=(q(4)**2-q(1)**2-q(2)**2-q(3)**2)
      aDa=(a(4)**2-a(1)**2-a(2)**2-a(3)**2)
      wDw=(w(4)**2-w(1)**2-w(2)**2-w(3)**2)
      zDz=(z(4)**2-z(1)**2-z(2)**2-z(3)**2)
      rDr=(rr(4)**2-rr(1)**2-rr(2)**2-rr(3)**2)
      bDb=(bb(4)**2-bb(1)**2-bb(2)**2-bb(3)**2)


C---Fill spinor products
      call spinoru(10,p,za,zb)

      densq=      ((sw1-wmass**2)**2+wmass**2*wwidth**2)
      densq=densq*((sw2-wmass**2)**2+wmass**2*wwidth**2)
      densq=densq*((sh-zmass**2)**2+zmass**2*zwidth**2)
      densq=densq*((qDq-mt**2)**2+mt**2*twidth**2)
      densq=densq*((aDa-mt**2)**2+mt**2*twidth**2)
      densq=densq*(1d0+(mt*twidth/(rDr-mt**2))**2)
      densq=densq*(1d0+(mt*twidth/(bDb-mt**2))**2)


      fac=V/4d0/densq*(gwsq/2d0)**4*gsq**2*esq**2

      d1=1d0/wDw/ssDss
      d2=1d0/zDz/ssDss
      d3=1d0/s12/(bDb-mt**2)
      d4=1d0/s12/(rDr-mt**2)

      call diag1(1,2,3,4,5,6,7,8,9,10,T1)
      call diag2(1,2,3,4,5,6,7,8,9,10,T2)
      call diag3(1,2,3,4,5,6,7,8,9,10,T3)
      call diag4(1,2,3,4,5,6,7,8,9,10,T4)


      sumup(1,1)=l1*(L(2)
     . *(d1*T1(1,1)+d2*T2(1,1))+L(2)*(d3*T3(1,1)+d4*T4(1,1)))
      sumup(1,2)=l1*(R(2)
     . *(d1*t1(1,2)+d2*T2(1,2))+L(2)*(d3*T3(1,2)+d4*T4(1,2)))
      sumup(2,1)=r1*(L(2)
     . *(d1*t1(2,1)+d2*T2(2,1))+L(2)*(d3*T3(2,1)+d4*T4(2,1)))
      sumup(2,2)=r1*(R(2)
     . *(d1*t1(2,2)+d2*T2(2,2))+L(2)*(d3*T3(2,2)+d4*T4(2,2)))

      sumdo(1,1)=l1
     . *(L(1)*(d1*t1(1,1)+d2*T2(1,1))+L(2)*(d3*T3(1,1)+d4*T4(1,1)))
      sumdo(1,2)=l1
     . *(R(1)*(d1*t1(1,2)+d2*T2(1,2))+L(2)*(d3*T3(1,2)+d4*T4(1,2)))
      sumdo(2,1)=r1
     . *(L(1)*(d1*t1(2,1)+d2*T2(2,1))+L(2)*(d3*T3(2,1)+d4*T4(2,1)))
      sumdo(2,2)=r1
     . *(R(1)*(d1*t1(2,2)+d2*T2(2,2))+L(2)*(d3*T3(2,2)+d4*T4(2,2)))
      qqbup=abs(sumup(1,1))**2+abs(sumup(2,1))**2
     .     +abs(sumup(1,2))**2+abs(sumup(2,2))**2
      qqbdo=abs(sumdo(1,1))**2+abs(sumdo(2,1))**2
     .     +abs(sumdo(1,2))**2+abs(sumdo(2,2))**2
c      write(6,*) 'sumup(1,1)',sumup(1,1)
c      write(6,*) 'sumup(2,1)',sumup(2,1)
c      write(6,*) 'sumup(1,2)',sumup(1,2)
c      write(6,*) 'sumup(2,2)',sumup(2,2)
c      write(6,*) 'sumdo(1,1)',sumdo(1,1)
c      write(6,*) 'sumdo(2,1)',sumdo(2,1)
c      write(6,*) 'sumdo(1,2)',sumdo(1,2)
c      write(6,*) 'sumdo(2,2)',sumdo(2,2)

      d2=1d0/wDw/ssDss
      d1=1d0/zDz/ssDss
      d3=1d0/s12/(bDb-mt**2)
      d4=1d0/s12/(rDr-mt**2)

      call diag1(2,1,3,4,5,6,7,8,9,10,T1)
      call diag2(2,1,3,4,5,6,7,8,9,10,T2)
      call diag3(2,1,3,4,5,6,7,8,9,10,T3)
      call diag4(2,1,3,4,5,6,7,8,9,10,T4)


      sumup(1,1)=l1*(L(2)
     . *(d1*T1(1,1)+d2*T2(1,1))+L(2)*(d3*T3(1,1)+d4*T4(1,1)))
      sumup(1,2)=l1*(R(2)
     . *(d1*t1(1,2)+d2*T2(1,2))+L(2)*(d3*T3(1,2)+d4*T4(1,2)))
      sumup(2,1)=r1*(L(2)
     . *(d1*t1(2,1)+d2*T2(2,1))+L(2)*(d3*T3(2,1)+d4*T4(2,1)))
      sumup(2,2)=r1*(R(2)
     . *(d1*t1(2,2)+d2*T2(2,2))+L(2)*(d3*T3(2,2)+d4*T4(2,2)))

      sumdo(1,1)=l1
     . *(L(1)*(d1*t1(1,1)+d2*T2(1,1))+L(2)*(d3*T3(1,1)+d4*T4(1,1)))
      sumdo(1,2)=l1
     . *(R(1)*(d1*t1(1,2)+d2*T2(1,2))+L(2)*(d3*T3(1,2)+d4*T4(1,2)))
      sumdo(2,1)=r1
     . *(L(1)*(d1*t1(2,1)+d2*T2(2,1))+L(2)*(d3*T3(2,1)+d4*T4(2,1)))
      sumdo(2,2)=r1
     . *(R(1)*(d1*t1(2,2)+d2*T2(2,2))+L(2)*(d3*T3(2,2)+d4*T4(2,2)))
      qbqup=abs(sumup(1,1))**2+abs(sumup(2,1))**2
     .     +abs(sumup(1,2))**2+abs(sumup(2,2))**2
      qbqdo=abs(sumdo(1,1))**2+abs(sumdo(2,1))**2
     .     +abs(sumdo(1,2))**2+abs(sumdo(2,2))**2

c      write(6,*)
c      write(6,*) 'sumup(1,1)',sumup(1,1)
c      write(6,*) 'sumup(2,1)',sumup(2,1)
c      write(6,*) 'sumup(1,2)',sumup(1,2)
c      write(6,*) 'sumup(2,2)',sumup(2,2)
c      write(6,*) 'sumdo(1,1)',sumdo(1,1)
c      write(6,*) 'sumdo(2,1)',sumdo(2,1)
c      write(6,*) 'sumdo(1,2)',sumdo(1,2)
c      write(6,*) 'sumdo(2,2)',sumdo(2,2)
c      pause

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      wtgg=0d0

      do j=-nf,nf
      k=-j
          if ((j .eq. 0) .and. (k .eq. 0)) then
           msq(j,k)=avegg*fac*wtgg
          elseif ((j .eq. 1) .or.(j .eq. 3)  .or.(j .eq. 5)) then
           msq(j,k)=aveqq*fac*qqbdo
          elseif ((j .eq. 2) .or.(j .eq. 4)) then
           msq(j,k)=aveqq*fac*qqbup
          elseif ((j .eq. -1) .or.(j .eq. -3)  .or.(j .eq. -5)) then
           msq(j,k)=aveqq*fac*qbqdo
          elseif ((j .eq. -2) .or.(j .eq. -4)) then
           msq(j,k)=aveqq*fac*qbqup
          endif
      enddo

      return
      end
