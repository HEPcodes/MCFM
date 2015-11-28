      subroutine tinyb(s,t,u,msqx)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      double precision s,t,u,ssq,tsq,usq,msqx(0:2),fac
      ssq=s**2
      tsq=t**2
      usq=u**2
      fac=gsq**2

      msqx(0)=4d0*fac*((ssq+usq)/tsq+(ssq+tsq)/usq-(xn+1d0/xn)*ssq/u/t)
      msqx(1)=2d0*fac*(xn**2*(ssq+tsq)/usq-(ssq+usq)/tsq
     . -2d0*(ssq+tsq)/usq+2d0/xn*ssq/u/t)

      msqx(2)=2d0*fac*(xn**2*(ssq+usq)/tsq-(ssq+tsq)/usq
     . -2d0*(ssq+usq)/tsq+2d0/xn*ssq/t/u)
      return
      end

