      subroutine chkC5_47(k,l,i1,i2,DetGr,f,Gtwiddle,Gtt,
     . Shat3,Shat4,S00ii,Shat4zz,N0)
      implicit none
      include 'constants.f' 
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,n,m,i1,i2,np
      parameter(np=2)
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np),f(np)
      double complex S00ii(z2max,-2:0),Shat4zz(np,z1max,-2:0),
     . Shat3(np,z2max,-2:0),Shat4(np,z3max,-2:0),bit,pole,diff

      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat3(m,z2(i1,i2),ep)
     . +2*(delta(n,i1)*Shat4zz(m,i2,ep)+delta(n,i2)*Shat4zz(m,i1,ep))
     . -f(n)*f(m)*Cv(cii(z2(i1,i2))+N0,ep)
     . -2*(f(n)*delta(m,i1)+f(m)*delta(n,i1))*Cv(czzi(i2)+N0,ep)
     . -2*(f(n)*delta(m,i2)+f(m)*delta(n,i2))*Cv(czzi(i1)+N0,ep)
     . -4*(delta(n,i1)*delta(m,i2)+delta(n,i2)*delta(m,i1))
     . *Cv(cc0000+N0,ep))
      enddo
      enddo
      pole=czip
      if (ep .gt. -2) pole=-4*Gtwiddle(k,l)*Cv(czzii(z2(i1,i2))+N0,ep-1)
      diff=
     . 14*Gtwiddle(k,l)*Cv(czzii(z2(i1,i2))+N0,ep)
     . +pole
     . +DetGr*Cv(ciiii(z4(k,l,i1,i2))+N0,ep)
     . -Gtwiddle(k,l)*S00ii(z2(i1,i2),ep)
     . -Gtwiddle(1,l)*Shat4(1,z3(k,i1,i2),ep)
     . -Gtwiddle(2,l)*Shat4(2,z3(k,i1,i2),ep)
     . +Gtwiddle(k,l)
     . *(Shat4(1,z3(1,i1,i2),ep)+Shat4(2,z3(2,i1,i2),ep))
     . +bit
      if (abs(diff) .gt. weenumber) 
     . write(6,*) 'chkC5_47',k,l,i1,i2,diff
      enddo


      return
      end
  



