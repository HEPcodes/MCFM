      subroutine chkC5_43(j,i1,DetGr,Xtwiddle0,Gtwiddle,Shat2,N0)
      implicit none
      include 'constants.f'  
      include 'Cnames.f'  
      include 'Cv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      include 'weenumber.f'  
      integer ep,N0,j,i1,n,np
      parameter(np=2)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat2(np,np,-2:0),Shat2s(np,np,-2:0),diff,bit
       
      do ep=-2,0
      bit=czip
      do n=1,np
      Shat2s(n,i1,ep)=Shat2(n,i1,ep)-2d0*delta(n,i1)*Cv(N0+cc00,ep)
      bit=bit+Gtwiddle(j,n)*Shat2s(n,i1,ep)
      enddo
      diff=bit
     . +Xtwiddle0(j)*Cv(ci(i1)+N0,ep)-DetGr*Cv(cii(z2(j,i1))+N0,ep)

      if (abs(diff) .gt. weenumber) write(6,*) 
     . 'chk5_43:np,j,i1,diff',np,j,i1,diff 
      enddo
      return
      end
