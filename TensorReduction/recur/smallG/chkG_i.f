      subroutine chkG_i(j,i1,DetGr,Xtwiddle0,Gtwiddle,Shat2,N0)
      implicit none
C--   DD 5.43
      include 'constants.f'  
      include 'Dnames.f'  
      include 'Dv.f'  
      include 'Darraydef.f'  
      include 'Darrays.f'  
      include 'weenumber.f'  
      integer ep,N0,j,i1,n,np
      parameter(np=3)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat2(np,np,-2:0),Shat2s(np,np,-2:0),diff,bit
       
      do ep=-2,0
      bit=czip
      do n=1,np
      Shat2s(n,i1,ep)=Shat2(n,i1,ep)-2d0*delta(n,i1)*Dv(N0+dd00,ep)
      bit=bit+Gtwiddle(j,n)*Shat2s(n,i1,ep)
      enddo
      diff=bit
     . +Xtwiddle0(j)*Dv(di(i1)+N0,ep)-DetGr*Dv(dii(z2(j,i1))+N0,ep)
      if ((abs(diff) .gt. weenumber)  .and. (Gsing .eqv. .false.))
     . write(6,*) 'chkG_i:np,j,i1,diff',np,j,i1,diff 
      enddo
      return
      end
