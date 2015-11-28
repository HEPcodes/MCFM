      double precision function virtd(s,t,u)
      implicit none
      include 'constants.f'
      include 'epinv.f'
      include 'scale.f'
      double precision s,t,u,Ws,Wt,Wu,Wmu,Ws2,Wt2,Wu2,FD,D4,HD,
     . theta,x,Q2,Tf
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
c===statement functions
      D4(S,T,U)=16d0*V*xn**2*(3d0-U*T/S**2-S*T/U**2-U*S/T**2)
      HD(S,T,U)=2d0*V*xn**3
     . *(12d0-8d0*U*T/S**2+4d0*(T**4+U**4)/(U*T)**2)
      FD(S,T,U,WS,WT,WU,WS2)=  
     . xn*(2d0*(U**2+T**2)/U/T*WS2
     . +(4d0*S*(U**3+T**3)/U**2/T**2-6d0)*WT*WU      
     . +(22d0/3d0*(U**2+T**2)/S**2+8d0/3d0
     . -14d0/3d0*(T**2+U**2)/U/T)*WS-1d0-PISQ)
     . +TF*((10d0/3d0*(U**2+T**2)/U/T+16d0/3d0*U*T/S**2-2d0)*WS   
     . -2d0*(U**2+T**2)/U/T*WU*WT-(S**2+U*T)/U/T*WS2+2d0-PISQ)   

      theta(x)=half+sign(half,x)
c===statement functions

      Tf=TR*dfloat(nf)
      Q2=musq
      WS=log(abs(s/Q2))
      WT=log(abs(t/Q2))
      WU=log(abs(u/Q2))
      WMU=log(abs(musq/Q2))
      WS2=WS**2-theta(s)*pisq
      WT2=WT**2-theta(t)*pisq
      WU2=WU**2-theta(u)*pisq


      virtd=D4(S,T,U)
     . *(-4d0*xn*EPINV**2-22d0*xn/3d0*EPINV+8d0*TF/3d0*EPINV)
     . +2d0*EPINV*(+WS*HD(S,T,U)+WT*HD(T,S,U)+WU*HD(U,T,S))
     . +D4(S,T,U)*(-67d0/9d0*xn+20d0/9d0*TF+xn*PISQ
     . +11d0*xn/3d0*WMU-4d0*TF/3d0*WMU)
     . +4d0*V*xn**2*(FD(S,T,U,WS,WT,WU,WS2)
     .              +FD(T,S,U,WT,WS,WU,WT2)
     .              +FD(U,T,S,WU,WT,WS,WU2))
      return
      end

