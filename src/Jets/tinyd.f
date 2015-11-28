      subroutine tinyd(s,t,u,msqx)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      double precision s,t,u,ssq,tsq,usq,msqx(0:2),fac
      ssq=s**2
      tsq=t**2
      usq=u**2
      fac=gsq**2
      msqx(0)=16d0*fac*V*xn**2*(2d0-t*u/ssq-ssq/(t*u))
      msqx(1)=16d0*fac*V*xn**2*(2d0-s*u/tsq-tsq/(s*u))
      msqx(2)=16d0*fac*V*xn**2*(2d0-s*t/usq-usq/(t*s))

      return
      end
