      double complex function I41(s,t,m0sq,m1sq)
c     Def : 
c     IntX  =  Mu^2*eps Integral[(1/DenominatorX) d^Dq1/((2*Pi))^D]
c     D = 4 - 2*eps
c     funX = IntX /((4*Pi*Mu^2/m^2)^eps*(I/16*Pi^2)/Gamma[1 - eps])
c     FIntX = Finite part


c     Three mass box I41
c     All momenta outgoing
c     Denominator4: q1^2*(q1-d)^2*[(q1-d-c)^2-m^2]*(q1-d-c-b)^2
c     Overall factor: 1/(s*(t-m1sq))
C     s=(b+c)^2,t=(d+c)^2
c     Divergences: epinv^2/2 + epinv*Log(m1sq*(-m0sq)/(-s*(m1sq-t)) )
C     FP-Denominator
C     a2*a4*(-s)+a1*a3*(-t)+a1*a4*(-m0sq)+a3*(1-a2-a4)*m1sq-i*ep)
C     Hence for s,t,m0sq<0,m1sq>0 it should be real

c                 d     c(m1sq)
c                 \ _____// 
c                  |     || 
c                  |     ||  <--s 
c                  |_____||
c                ///      \\
c         (m0sq)(u+q)   b(m1sq)
c
c                     ^
c                     |
c                     |
c                     t

C
c This integral is derived from Eq.(34) of hep/0211390 
c with a4/a5=s15/s23=m0sq/s   (Eq.(17))
c with m^2*M^2*a4^2
c  =-m^4*s51/s23/sbar34*s45/sbar12*(-s51/s23/sbar34*sbar12/s45)
c  =(m^2*s51/s23/sbar34)^2=m1sq*(-m0sq)/(-s)/(-2d.c)   Eqs(17) and (20)
c Overall factor = 1/rt(a2*a3*a4*a5)=1/s23/sbar34=1/s/(t-m1sq)
c Hence answer is
c 1/2/ep^2+1/ep*Log(m1sq*(-m0sq)/(-s)/(m1sq-t))
c +Log(m1sq*(-m0sq)/(-s)/(m1sq-t))^2-log(-m0sq/-s)-2*Li2(1-m0sq/s)+pisq/6

      implicit none
      include 'constants.f'
      double precision m0sq,t,s,m1sq,x1,y1,r1,omr1,DDILOG
      double complex Lnrat,dilog1,wlog1
      wlog1=Lnrat(-m0sq,-s)
      r1=m0sq/s
      omr1=1d0-m0sq/s
      if (omr1 .gt. one) then 
         dilog1=dcmplx(pisqo6-DDILOG(r1))-wlog1*dcmplx(log(omr1))
      else
         dilog1=dcmplx(DDILOG(omr1))
      endif
      I41 =(wlog1+Lnrat(m1sq,m1sq-t))**2-wlog1**2-2d0*dilog1+pisq/4d0
      return
      end


      double complex function I42(s,t,m0sq,m1sq,m2sq)
c    Three mass box FInt8
c    Denominator4: q1^2*[(q1-u-b)^2-m^2]*(q1-u-b-c)^2*(q1-u-b-c-d)^2
c    Overall factor: 1/(u*(v2-m2sq))
c    Divergences: epinv^2/2 + epinv*Log( Q2*m2sq/(u*(m2sq-v2)) )
C    FP-Denominator
C    +a1*a3*(-s)+a2*a4*(-t)+a4*a1*(-m0sq)+a2*a1*(-m1sq)+a2*(1-a3)*m2sq-i*ep
C    Hence for s,t,m0sq,m1sq<0,m2sq>0 it should be real


c                 d(4:0)     c(3:m2sq)
c                  \------//
c                  |      || 
c                  |      ||   <--s
c                  /------\\
c                 q(5:m0sq)  (u+b)(12:m1sq)     
c                  
c                     ^
c                     |
c                      t
c
c This integral is derived from Eq.(46) of hep/0211390 
c lambda=q2*s23/s45/s51
C lambda*M^2=q2*m^2/sbar12/sbar34=m0sq*m2sq/(m1sq-m2sq)/(t-m2sq)
C m*lambda*M*a4
C              =m*(q2*s23/s45/s51)*rt(-m^2*s45*s51/sbar12/s23/sbar34)
C              *rt(-s51*sbar12/(s23*sbar34*s45))
C              =m^2*q2/s45/sbar34
C              =m2sq*(-m0sq)/(-s)/(m2sq-t)
c a4/a5=s15/s23   (Eq.(17))
c lambda*a4/a5=q2/s45=m0sq/s
C so that function is
C   0.5/ep^2 +1/ep*log((m2sq*(-m0sq))/(-s)*(t-m2sq)))
C               +log^2((m2sq*(-m0sq))/((-s)*(t-m2sq)))
C   -1/2*log(m2sq*(-m0sq)/(-m1sq+m2sq)/(m2sq-t))
C   -log^2(-m0sq/-s)-dilog(1-m0sq*(-m2sq)/(-m1sq+m2sq)/(m2sq-t))
C   -2*dilog(1-m0sq/s)

      implicit none
      include 'constants.f'
      double precision s,t,m0sq,m1sq,m2sq,DDILOG
      double precision r1,r2,x2,y1,y2,omr1,omr2
      double complex Lnrat,dilog1,dilog2,wlogs,wlogt,wlogm,wlog1,wlog2
      wlogs=lnrat(-m0sq,-s)
      wlogt=lnrat(m2sq,m2sq-t)
      wlogm=lnrat(-m0sq,-m1sq+m2sq)
      wlog1=wlogs+wlogt
      wlog2=wlogm+wlogt

      r1=m0sq/s
      omr1=1d0-r1
      if (omr1 .gt. one) then 
         dilog1=dcmplx(pisqo6-DDILOG(r1))-wlogs*dcmplx(log(omr1))
      else
         dilog1=dcmplx(DDILOG(omr1))
      endif
      r2=-m0sq*m2sq/((m2sq-m1sq)*(m2sq-t))
      omr2=1d0-r2
      if (omr2 .gt. one) then 
         dilog2=dcmplx(pisqo6-DDILOG(r2))-wlog2*dcmplx(log(omr2))
      else
         dilog2=dcmplx(DDILOG(omr2))
      endif
      I42=wlog1**2-0.5d0*wlog2**2-wlogs**2-2d0*dilog1-dilog2+pisq/12d0
      return
      end




