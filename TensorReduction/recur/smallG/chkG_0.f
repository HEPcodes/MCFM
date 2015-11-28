      subroutine chkG_0(j,DetGr,Xtwiddle0,Gtwiddle,Shat1,N0)
      implicit none
C-----DD Eq. 5.41
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f'
      integer ep,N0,j,np
      parameter(np=3)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat1(np,-2:0),diff

      do ep=-2,0
      diff=Xtwiddle0(j)*Dv(dd0+N0,ep)
     . +Gtwiddle(j,1)*Shat1(1,ep)
     . +Gtwiddle(j,2)*Shat1(2,ep)
     . +Gtwiddle(j,3)*Shat1(3,ep)
     . -DetGr*Dv(di(j)+N0,ep)

      if ((abs(diff) .gt. weenumber) .and. (Gsing .eqv. .false.)) 
     . write(6,*) 'chkG_0',j,diff
      enddo
      return
      end
