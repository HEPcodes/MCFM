      subroutine chkC_0000ii(k,l,i1,i2,DetGr,f,Gtwiddle,Gtt, 
     . Shat5zz,Shat6zzzz,S0000ii,Shat6zz,N0) 
      implicit none 
      include 'constants.f'   
      include 'Cnames.f'   
      include 'Cv.f'   
      include 'Carraydef.f'   
      include 'Carrays.f'   
      include 'weenumber.f'   
      integer ep,N0,k,l,n,m,i1,i2,np 
      parameter(np=2) 
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np), 
     . f(np) 
      double complex Shat5zz(np,z2max,-2:0),Shat6zzzz(np,z1max,-2:0), 
     . S0000ii(z2max,-2:0),Shat6zz(np,z3max,-2:0),bit,pole,diff
        
      do ep=-2,0
      bit=czip 
      do n=1,np 
      do m=1,np 
      bit=bit 
     . +Gtt(k,n,l,m)*(f(n)*Shat5zz(m,z2(i1,i2),ep) 
     . +2*(delta(n,i1)*Shat6zzzz(m,i2,ep)
     .    +delta(n,i2)*Shat6zzzz(m,i1,ep)) 
     . -f(n)*f(m)*Cv(czzii(z2(i1,i2))+N0,ep) 
     . -2*(f(n)*delta(m,i1)+f(m)*delta(n,i1))*Cv(czzzzi(i2)+N0,ep)
     . -2*(f(n)*delta(m,i2)+f(m)*delta(n,i2))*Cv(czzzzi(i1)+N0,ep)
     . -4*(delta(n,i1)*delta(m,i2)+delta(m,i1)*delta(n,i2))
     . *Cv(cc000000+N0,ep))
      enddo 
      enddo 

      pole=czip 
      if (ep .gt. -2)  
     . pole=-4*Gtwiddle(k,l)*Cv(czzzzii(z2(i1,i2))+N0,ep-1) 
      diff= 
     . 18d0*Gtwiddle(k,l)*Cv(czzzzii(z2(i1,i2))+N0,ep) 
     . +pole 
     . +DetGr*Cv(czziiii(z4(k,l,i1,i2))+N0,ep) 
     . -Gtwiddle(k,l)*S0000ii(z2(i1,i2),ep) 
     . -Gtwiddle(1,l)*Shat6zz(1,z3(k,i1,i2),ep) 
     . -Gtwiddle(2,l)*Shat6zz(2,z3(k,i1,i2),ep) 
     . +Gtwiddle(k,l) 
     . *(Shat6zz(1,z3(1,i1,i2),ep) 
     .  +Shat6zz(2,z3(2,i1,i2),ep)) 
     . +bit 
     
      if ((abs(diff) .gt. weenumber) .and. (Gsing .eqv. .false.))  
     . write(6,*) 'chkC_0000ii:ep,k,l,i1,i2,diff',ep,k,l,i1,i2,diff 
      enddo 
      
      return
      end
