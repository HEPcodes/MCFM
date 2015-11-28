      double complex function amp_qqgg(i1,h1,i2,h2,i3,h3,i4,lh)
C--Results taken from Nagy and Trocsanyi-hep-ph/9806317
C--Appendix A
C  lh is the helicity of lepton pair
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'prods.f'
      integer i1,i2,i3,i4,i5,i6,j,k,h1,h2,h3,lh,g1,g2,g3,lg
      double complex xa(mxpart,mxpart),xb(mxpart,mxpart),t2a
      double precision t123,t234
      
      t2a(i1,i2,i3,i4)=xa(i1,i2)*xb(i2,i4)+xa(i1,i3)*xb(i3,i4)
      
C---A(hq,h2,h3,lh)      
C---h=1 LH
C---h=2 RH

      i5=5
      i6=6
      g1=h1
      g2=h2
      g3=h3
      lg=lh

      if (g1*lg.eq.2)then
         i5=6
         i6=5
         lg=3-lg  
      endif 

      if (g1.eq.2) then
          do j=1,mxpart
          do k=1,mxpart
             xa(j,k)=za(j,k)
             xb(j,k)=zb(j,k)
          enddo
          enddo
      elseif (g1.eq.1) then
          g1=3-g1
          g2=3-g2
          g3=3-g3
          lg=3-lg
          do j=1,mxpart
          do k=1,mxpart
             xa(j,k)=zb(k,j)
             xb(j,k)=za(k,j)
          enddo
          enddo
      endif


C--A37
c      A(2,2,2,2) =
      if ((g1.eq.2).and.(g2.eq.2).and.(g3.eq.2).and.(lg.eq.2)) then
      amp_qqgg=
     .-xa(i4,i5)**2*xb(i5,i6)/(xa(i1,i2)*xa(i2,i3)*xa(i3,i4)*s(i5,i6))
 
C--A38
c      A(2,1,1,2) =
      elseif ((g1.eq.2).and.(g2.eq.1).and.(g3.eq.1).and.(lg.eq.2)) then
      amp_qqgg=
     . -xb(i1,i6)**2*xa(i5,i6)/(xb(i1,i2)*xb(i2,i3)*xb(i3,i4)*s(i5,i6))


C--A39
c      A(2,2,1,2)=
      elseif ((g1.eq.2).and.(g2.eq.2).and.(g3.eq.1).and.(lg.eq.2)) then
      t123=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t234=s(i2,i3)+s(i3,i4)+s(i4,i2)
      amp_qqgg=
     . -xa(i3,i1)*xb(i1,i2)*xa(i4,i5)*t2a(i3,i1,i2,i6)
     . /(xa(i1,i2)*s(i2,i3)*t123*s(i5,i6)) 
     . +xa(i3,i4)*xb(i4,i2)*xb(i1,i6)*t2a(i5,i3,i4,i2)
     . /(xb(i3,i4)*s(i2,i3)*t234*s(i5,i6))
     . +t2a(i5,i3,i4,i2)*t2a(i3,i1,i2,i6)
     . /(xa(i1,i2)*xb(i3,i4)*s(i2,i3)*s(i5,i6)) 
C---A40
c      A(2,1,2,2) =
      elseif ((g1.eq.2).and.(g2.eq.1).and.(g3.eq.2).and.(lg.eq.2)) then
      t123=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t234=s(i2,i3)+s(i3,i4)+s(i4,i2)
      amp_qqgg=
     . +xb(i1,i3)**2*xa(i4,i5)*t2a(i2,i1,i3,i6)
     . /(xb(i1,i2)*s(i2,i3)*t123*s(i5,i6))
     . -xa(i2,i4)**2*xb(i1,i6)*t2a(i5,i2,i4,i3)
     . /(xa(i3,i4)*s(i2,i3)*t234*s(i5,i6))
     . -xb(i1,i3)*xa(i2,i4)*xb(i1,i6)*xa(i4,i5)
     . /(xb(i1,i2)*xa(i3,i4)*s(i2,i3)*s(i5,i6))
      endif

      return 
      end
