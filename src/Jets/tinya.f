      subroutine tinya(s,t,u,msqx)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      double precision s,t,u,ssq,tsq,usq,msqx(0:2),fac
      ssq=s**2
      tsq=t**2
      usq=u**2
      fac=gsq**2
      msqx(0)=2d0*fac*V*(ssq+usq)/tsq
      msqx(1)=0d0
      msqx(2)=0d0
      return
      end

