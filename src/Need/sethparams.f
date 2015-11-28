      subroutine sethparams(br,wwbr,zzbr,tautaubr)
c--- set up the necessary parameters for a Standard Model Higgs boson
c---   hwidth : either the NLO value from Spira,
c---                or the LO value valid for low Higgs masses only
c---   br,wwbr,zzbr,tautaubr : the LO calculated values
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      double precision br,gamgambr,wwbr,zzbr,tautaubr,x_w,x_z
      logical spira
      common/spira/spira

c--- calculate the value of hwidth, depending on "spira"
c--- note that the branching ratios "br","wwbr","zzbr" will be over-written   
      if (spira) then
        call higgsp(br,gamgambr,wwbr,zzbr)
      else
        call higgsw(br)
      endif 

c--- branching ratio H->bb : note that the mass is mostly retained in the
c--- coupling only, since mb=0 usually (the mass in the phase space)
      br=xn*gwsq/32d0/pi*mbsq*hmass/wmass**2
     . *(1d0-4d0*mb**2/hmass**2)**1.5d0/hwidth
c--- branching ratio H->tau^- tau^+ : this is obtained by a rescaling of
c--- the H->bb BR since the tau mass is not included in the phase space
      tautaubr=br*(mtausq/mbsq/xn)
      
      x_w=4d0*wmass**2/hmass**2
      x_z=4d0*zmass**2/hmass**2
c--- branching ratio H-> WW
      if (x_w .lt. 1d0) then
        wwbr=gwsq/64d0/pi*hmass**3/wmass**2
     .   *dsqrt(1d0-x_w)*(1d0-x_w+0.75d0*x_w**2)/hwidth
      else
        wwbr=0d0
      endif
c--- branching ratio H-> ZZ
      if (x_z .lt. 1d0) then
        zzbr=gwsq/128d0/pi*hmass**3/wmass**2
     .   *dsqrt(1d0-x_z)*(1d0-x_z+0.75d0*x_z**2)/hwidth
      else
        zzbr=0d0
      endif
                
      write(6,99) hmass,hwidth,br,tautaubr,wwbr,zzbr
      
      return

 99   format(/,
     .       ' ****************** Higgs parameters ****************'/, 
     .       ' *                                                  *'/, 
     .       ' *   mass(H) = ',f7.2,'      width(H) = ',e12.5,' *'/,
     .       ' *                                                  *'/, 
     .       ' *              Br( H -> b bbar)  = ',f8.4,'        *'/,
     .       ' *              Br( H -> tau tau) = ',f8.4,'        *'/,
     .       ' *              Br( H -> W W)     = ',f8.4,'        *'/,
     .       ' *              Br( H -> Z Z)     = ',f8.4,'        *'/,
     .       ' ****************************************************')

      end
      
