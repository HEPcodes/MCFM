      subroutine chkC5_44(k,l,i1,DetGr,f,Gtwiddle,Gtt,
     . Shat2,Shat3,Shat3zz,S00i,N0)
      implicit none
      include 'constants.f'  
      include 'Cnames.f'  
      include 'Cv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      include 'weenumber.f'  
      integer ep,N0,k,l,n,m,i1,np
      parameter(np=2)
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np),f(np)
      double complex S00i(np,-2:0),Shat3zz(np,-2:0),
     . Shat2(np,np,-2:0),Shat3(np,z2max,-2:0),pole,bit,diff
 
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat2(m,i1,ep)
     . +2d0*delta(n,i1)*Shat3zz(m,ep)-f(n)*f(m)*Cv(ci(i1)+N0,ep)
     . -2d0*(f(n)*delta(m,i1)+f(m)*delta(n,i1))*Cv(cc00+N0,ep))
      enddo
      enddo
      pole=0d0
      if (ep.ge.-1) pole=-4d0*Gtwiddle(k,l)*Cv(czzi(i1)+N0,ep-1)

      diff=
     . 10*Gtwiddle(k,l)*Cv(czzi(i1)+N0,ep)
     . +pole
     . +DetGr*Cv(ciii(z3(k,l,i1))+N0,ep)
     . -Gtwiddle(k,l)*S00i(i1,ep)
     . -Gtwiddle(1,l)*Shat3(1,z2(k,i1),ep)
     . -Gtwiddle(2,l)*Shat3(2,z2(k,i1),ep)
     . +Gtwiddle(k,l)*(Shat3(1,z2(1,i1),ep)
     .                +Shat3(2,z2(2,i1),ep))
     . +bit

      if (abs(diff).gt.weenumber) 
     . write(6,*) 'chkC5_44:',k,l,i1,diff
      enddo

      return
      end
  



