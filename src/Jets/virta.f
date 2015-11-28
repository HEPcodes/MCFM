      double precision function virta(s,t,u)
      implicit none
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
      include 'constants.f'
      include 'epinv.f'
      include 'scale.f'
      double precision a4,theta,x,s,t,u,Ws,Wt,Wu,Ws2,WT2,Wu2,Wmu,Q2,TF
C===statement functions
      a4(s,t,u) =2D0*V*(s**2+u**2)/t**2      
      theta(x)=half+sign(half,x)
C===statement functions

      Tf=TR*dfloat(nf)
      Q2=musq
      Ws=log(abs(s/q2))
      Wt=log(abs(t/q2))
      Wu=log(abs(u/q2))
      Wmu=log(abs(musq/q2))

      Ws2=Ws**2-theta(s)*pisq
      Wt2=Wt**2-theta(t)*pisq
      Wu2=Wu**2-theta(u)*pisq

      virta=
     . +(CF*(-4d0*EPINV**2-EPINV*(6d0+8d0*WS-8d0*WU-4d0*WT)) 
     . +xn*EPINV*(4d0*WS-2d0*WU-2d0*WT))*A4(S,T,U)
      virta=virta
     . +(CF*(-16d0-2d0*WT2+WT*(6d0+8d0*WS-8d0*WU))
     . +xn*(85d0/9d0+PISQ+2d0*WT*(WT+WU-2d0*WS) 
     . +11d0/3d0*(WMU-WT))
     . +TF*(4d0/3d0*(WT-WMU)-20d0/9d0))*A4(S,T,U)
     . -4d0*V*CF*(S**2-U**2)/T**2
     . *(2d0*PISQ+2d0*WT2+WS2+WU2-2d0*WS*WT-2d0*WT*WU)
     . +xn*V*((S**2-U**2)/T**2
     . *(3d0*PISQ+3d0*WT2+2d0*WS2+WU2-4d0*WS*WT-2d0*WT*WU))
     . +4d0*V*(CF*((WU-WS)-(U-S)/T*(2d0*WT-WS-WU))
     . +xn*(-0.5d0*S/T*(WT-WU)+U/T*(WT-WS)))


 
      return
      end
