C  Result for the easy box is
C  2*Rgamma/(s*t-m1^2*m3^2)*
C     (1/ep^2*(1-ep*lnrat(s/m1sq)-ep*lnrat(t/m3sq)
C     +0.5*log^2(-s)+0.5*log^2(-t)-0.5*log^2(-m1^2)+0.5*log^2(-m3^2)
C     +Lsm1_2me(s,t,m1sq,m3sq))
C 
C        m1sq
C        \\          /
C         \\--------/
C          |   a1   |
C          |a2    a4|
C          |   a3   |
c          /--------\\
c         /          \\m3^2

C in units of r_Gamma (1/msq)^ep 
C        \sim 1/Gamma(1-ep)*(1/msq)^ep 
      double complex function BDKBoxe(s,t,m1sq,m3sq,msq)
      implicit none
      include 'constants.f'
      include 'epinv.f'
      double precision s,t,m1sq,m3sq,msq
      double complex lnrat,wlogs,wlogt,Lsm1_2me
      wlogs=lnrat(-s,-m1sq)
      wlogt=lnrat(-t,-m3sq)
      BDKBoxe=2d0*(epinv**2-(wlogs+wlogt)*epinv
     . +0.5d0*wlogs**2+0.5d0*wlogt**2
     . +wlogs*lnrat(m1sq,msq)+wlogt*lnrat(m3sq,msq)
     . +Lsm1_2me(s,t,m1sq,m3sq))
      return
      end     
