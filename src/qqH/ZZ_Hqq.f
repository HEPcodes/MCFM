      subroutine ZZ_Hqq(p,msq)
      implicit none 
c--- Weak Bosion Fusion by Z-Z exchange only
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
     . ud_ud_LL,udb_udb_LL,ud_ud_LR,udb_udb_LR
      
      data pn/-1,-2,-1,-2,-1,0,1,2,1,2,1/
      data x1/1,2/,x2/2,1/

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      
      call dotem(6,p,s)

      s34=s(3,4)+2d0*mb**2
      hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s34-4d0*mb**2) 
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      fac=0.25d0*gwsq**3*hdecay
C Color cancels, 0.25d0 is spin average

C q-q and qbar-qbar
c--- u(1)+d(2) -> u(5)+d(6)
c--- ub(1)+db(2) -> ub(5)+db(6)
      call msqpieces_zz(1,2,5,6,ud_ud_LL,ud_ud_LR)

C q-qbar and qbar-q
c--- u(1)+db(2) -> u(5)+db(6)
c--- ub(1)+d(2) -> ub(5)+d(6)
      call msqpieces_zz(1,6,5,2,udb_udb_LL,udb_udb_LR)

      do j=-nf,nf
      do k=-nf,nf
        if     ((j .gt. 0) .and. (k .lt. 0)) then
          msq(j,k)=fac*(
     .     +udb_udb_LL*((L(+j)*L(-k))**2+(R(+j)*R(-k))**2)
     .     +udb_udb_LR*((L(+j)*R(-k))**2+(R(+j)*L(-k))**2))
        elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=fac*(
     .     +udb_udb_LL*((L(-j)*L(k))**2+(R(-j)*R(k))**2)
     .     +udb_udb_LR*((L(-j)*R(k))**2+(R(-j)*L(k))**2))
        elseif ((j .gt. 0) .and. (k .gt. 0)) then
          msq(j,k)=fac*(
     .     +ud_ud_LL*((L(+j)*L(+k))**2+(R(+j)*R(+k))**2)
     .     +ud_ud_LR*((L(+j)*R(+k))**2+(R(+j)*L(+k))**2))
        elseif ((j .lt. 0) .and. (k .lt. 0)) then
          msq(j,k)=fac*(
     .     +ud_ud_LL*((L(-j)*L(-k))**2+(R(-j)*R(-k))**2)
     .     +ud_ud_LR*((L(-j)*R(-k))**2+(R(-j)*L(-k))**2))
        endif 
      enddo
      enddo

      return
      end


      subroutine msqpieces_zz(i1,i2,i5,i6,zll,zlr)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      double precision zll,zLR,htheta
      double precision propw,propz,x
      integer i1,i2,i5,i6
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+sign(half,x)

      propz(i1,i2)=sign(one,(s(i1,i2)-zmass**2))
     . *dsqrt(dsqrt(1d0-xw)/xw/2d0/zmass
     . *((s(i1,i2)-zmass**2)**2+htheta(s(i1,i2))*(zmass*zwidth)**2))

      zll=s(i1,i2)*s(i5,i6)/(propz(i1,i5)*propz(i2,i6))**2
      zlr=s(i1,i6)*s(i2,i5)/(propz(i1,i5)*propz(i2,i6))**2
      return
      end
