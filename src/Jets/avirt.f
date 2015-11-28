      double precision function avirt(s,t,u)
      implicit none
      include 'constants.f'
      include 'epinv.f'
      double precision s,t,u,Ws,Wt,Wu,Wmu,A4,U2PS2
      
c A4(S,T,U)*U2PS2**-1=2*V/T**2         
c A4(S,T,U) =2*V*(S**2+U**2-EP*T**2)/T**2        
      U2PS2=u**2+s**2      
      avirt=(CF*(-4d0*EPINV**2-EPINV*(6d0+8d0*WS-8d0*WU-4d0*WT)) 
     . +xn*EPINV*(4d0*WS-2d0*WU-2d0*WT)      
     . +CF*(-16d0-2d0*WT**2+WT*(6d0+8d0*WS-8d0*WU)  
     . -2d0*(S**2-U**2)/U2PS2*(2d0*pisq+(WT-WS)**2+(WT-WU)**2)       
     . +2d0*(S+U)/U2PS2*((S+U)*(WU-WS)+(U-S)*(2*WT-WS-WU)))         
     . +xn*(85d0/9d0+pisq+2d0*WT*(WT+WU-2d0*WS)      
     . +(S**2-U**2)/2d0/U2PS2*(3d0*pisq+2d0*(WT-WS)**2+(WT-WU)**2)     
     . -S*T/U2PS2*(WT-WU)+2d0*U*T/U2PS2*(WT-WS)+11d0/3d0*(WMU-WT))      
     . +TR*(4d0/3d0*(WT-WMU)-20d0/9d0))*A4(S,T,U)     
      return
      end
