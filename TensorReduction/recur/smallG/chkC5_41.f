      subroutine chkC5_41(j,DetGr,Xtwiddle0,Gtwiddle,Shat1,N0)
      implicit none
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f'
      integer ep,N0,j,np
      parameter(np=2)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat1(np,-2:0),diff
       
      do ep=-2,0
      diff=Xtwiddle0(j)*Cv(cc0+N0,ep)
     . +Gtwiddle(j,1)*Shat1(1,ep)
     . +Gtwiddle(j,2)*Shat1(2,ep)
     . -DetGr*Cv(ci(j)+N0,ep)
      if (abs(diff) .gt. weenumber) write(6,*) 'chkC5_41',j,diff
      enddo
      return
      end
