      subroutine tinyc(s,t,u,msqx)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      double precision s,t,u,ssq,tsq,usq,msqx(0:2),fac
      ssq=s**2
      tsq=t**2
      usq=u**2
      fac=gsq**2

      msqx(0)=-2d0*fac*V/xn*(t/u+u/t)
      msqx(1)=2d0*fac*V*xn*(t/u-2d0*tsq/ssq)
      msqx(2)=2d0*fac*V*xn*(u/t-2d0*usq/ssq)

      return
      end
