      double precision function virtc(s,t,u)
      implicit none
      include 'constants.f'
      include 'epinv.f'
      include 'scale.f'
C     Virtual corrections taken from Ellis and Sexton
C  %\cite{Ellis:1985er}
C  \bibitem{Ellis:1985er}
C   R.~K.~Ellis and J.~C.~Sexton,
C  %``QCD Radiative Corrections To Parton Parton Scattering,''
C  Nucl.\ Phys.\ B {\bf 269}, 445 (1986).
C  %%CITATION = NUPHA,B269,445;%%
C     In units of 
C     as/2/pi*(4 \pi)^epsilon*Gamma(1+ep)*Gamma(1-ep)^2/Gamma(1-2*ep) 
C      \approx as/2/pi*(4 \pi)^epsilon /Gamma(1-ep)
      double precision s,t,u,Ws,Wt,Wu,Ws2,Wt2,Wu2,
     . Wmu,theta,C4,FC,HCS,HCT,HCU,x,Q2,Tf
c===statement functions
      FC(S,T,U,WS,WT,WU,WS2,WT2)=
     . xn*V*(    
     . +WS2*(-(xn+1d0/xn)*(T**2+U**2)/S**2
     . +(5d0/2d0+2d0*(T**2+U**2)/T/U)/xn-V/xn**3*S**2/T/U)
     . +WT2*(xn*(S/T-4d0*U/S-1d0)-1d0/xn/U/S*(S**2+T**2)
     . +1d0/xn**3/U/T*(S**2+T**2))
     . +WT*WU*2d0*(T**2+U**2)/T/U/xn
     . +WT*WS*(xn*(4d0*(2d0*T**2+S**2+2*S*T)*T+2d0*S**2*(S+T))/T/S**2
     . +2d0/xn/U/S*(S**2+T**2)-2d0/xn**3/U/T*(S**2+T**2))
     . +PISQ/2d0*(xn*(T**2+U**2)*(1d0/T/U-4d0/S**2)
     . +1d0/xn+1d0/xn**3*(3d0*(T**2+U**2)/U/T+4d0))
     . +WT*(V/xn*(3d0*S*U+5d0*U*T)/T/S-4d0/xn*T/S+V/xn**3*((T+3d0*S)/T))
     . +WS*(V/2/xn*(12d0*U*T/S**2-1d0)-2d0/xn**3)
     . +(xn+1d0/xn)*(0.5d0-(T**2+U**2)/S**2))
      theta(x)=half+sign(half,x)
c===statement functions


      Tf=TR*dfloat(nf)
      q2=musq
      Ws=log(abs(s/q2))
      Wt=log(abs(t/q2))
      Wu=log(abs(u/q2))
      Wmu=log(abs(musq/q2))

      Ws2=Ws**2-theta(s)*pisq
      Wt2=Wt**2-theta(t)*pisq
      Wu2=Wu**2-theta(u)*pisq
      hcs=2d0*(t**2+u**2)/t/u
      hct=2d0*(t**2+u**2)/s**2*t/u
      hcu=2d0*(t**2+u**2)/s**2*u/t
      c4=2d0*V/xn*(V/(u*t)-2d0*xn**2/s**2)*(u**2+t**2)


      virtc=
     . +(CF*(-2d0*EPINV**2-3d0*EPINV)
     . +xn*(-2d0*EPINV**2-11d0/3d0*EPINV)       
     . +TF*4d0/3d0*EPINV)*C4
     . +V*WS*EPINV*(-HCT-HCU+HCS+HCS/xn**2+xn**2*(HCT+HCU))
     . +2d0*V*WT*EPINV*(xn**2*HCU-HCS)
     . +2d0*V*WU*EPINV*(xn**2*HCT-HCS)
      virtc=virtc
     . +(xn*11d0/3d0*WMU-4d0/3d0*TF*WMU-7d0*CF)*C4
     . +FC(S,T,U,WS,WT,WU,WS2,WT2)
     .  +FC(S,U,T,WS,WU,WT,WS2,WU2) 
      return
      end


