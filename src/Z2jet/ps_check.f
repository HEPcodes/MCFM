      subroutine ps_check(p,n)
c--- Subroutine containing phase space points used for checking
c--- loop amplitudes for Z+2 jets against Madloop
      implicit none
      include 'constants.f'
      double precision p(mxpart,4)
      integer n
      
      if     (n .eq. 1) then
c pu = (79.343740010234328, 0, 0, 79.343740010234328)
c pg = (330.04970916921303, 0, 0, -330.04970916921303)
c pe- = (44.502332440251230, -9.0128674723142481, 20.893959782906148, -38.244846151337029)
c pe+ = (173.55124786955787, -83.712846866171589, -85.405837095621621, -125.76978133334636)
c pu = (142.62173792637984, 51.789721902948564, 71.086870231086095, -112.27395831226823)
c pg = (48.718130943258451, 40.935992435537273, -6.5749929183706302, 25.582616637972905)
        p(1,4) = -79.343740010234328d0
        p(1,1) = 0d0
        p(1,2) = 0d0
        p(1,3) = -79.343740010234328d0
        p(2,4) = -330.04970916921303d0
        p(2,1) = 0d0
        p(2,2) = 0d0
        p(2,3) = 330.04970916921303d0
        p(3,4) = 44.502332440251230d0
        p(3,1) = -9.0128674723142481d0
        p(3,2) = 20.893959782906148d0
        p(3,3) = -38.244846151337029d0
        p(4,4) = 173.55124786955787d0
        p(4,1) = -83.712846866171589d0
        p(4,2) = -85.405837095621621d0
        p(4,3) = -125.76978133334636d0
        p(5,4) = 142.62173792637984d0
        p(5,1) = 51.789721902948564d0
        p(5,2) = 71.086870231086095d0
        p(5,3) = -112.27395831226823d0
        p(6,4) = 48.718130943258451d0
        p(6,1) = 40.935992435537273d0
        p(6,2) = -6.5749929183706302d0
        p(6,3) = 25.582616637972905d0
      elseif (n .eq. 2) then
c 0.100000000000000E+04    0.000000000000000E+00    0.000000000000000E+00    0.100000000000000E+04
c 0.100000000000000E+04    0.000000000000000E+00    0.000000000000000E+00   -0.100000000000000E+04
c 0.235495416274342E+03   -0.121411843692723E+03    0.142462091724653E+03    0.142904890464937E+03
c 0.701908346155938E+03   -0.635763334052215E+02    0.167878857386347E+03    0.678564709866912E+03
c 0.698664457147580E+03    0.368018607991527E+03   -0.103055595878474E+03   -0.584870816515437E+03
c 0.363931780422140E+03   -0.183030430893582E+03   -0.207285353232526E+03   -0.236598783816412E+03
        p(1,4) = -0.100000000000000d+04
        p(1,1) = 0d0
        p(1,2) = 0d0
        p(1,3) = -0.100000000000000d+04
        p(2,4) = -0.100000000000000d+04
        p(2,1) = 0d0
        p(2,2) = 0d0
        p(2,3) = 0.100000000000000d+04
        p(3,4) = 0.235495416274342d+03
        p(3,1) = -0.121411843692723d+03
        p(3,2) = 0.142462091724653d+03
        p(3,3) = 0.142904890464937d+03
        p(4,4) = 0.701908346155938d+03
        p(4,1) = -0.635763334052215d+02
        p(4,2) = 0.167878857386347d+03
        p(4,3) = 0.678564709866912d+03
        p(5,4) = 0.698664457147580d+03
        p(5,1) = 0.368018607991527d+03
        p(5,2) = -0.103055595878474d+03
        p(5,3) = -0.584870816515437d+03
        p(6,4) = 0.363931780422140d+03
        p(6,1) = -0.183030430893582d+03
        p(6,2) = -0.207285353232526d+03
        p(6,3) = -0.236598783816412d+03
      elseif (n .eq. 3) then
c 0.100000000000000E+04    0.000000000000000E+00    0.000000000000000E+00    0.100000000000000E+04
c 0.100000000000000E+04    0.000000000000000E+00    0.000000000000000E+00   -0.100000000000000E+04
c 0.350815654600186E+03    0.274947585251735E+03   -0.177057379208442E+03   -0.126988713453734E+03
c 0.602328810143063E+02    0.156351883341439E+02   -0.472677797496510E+02   -0.339012955294915E+02
c 0.910714300630089E+03   -0.757455408085941E+03   -0.146796400460784E+03    0.483851897738626E+03
c 0.678237163755404E+03    0.466872634500051E+03    0.371121559418886E+03   -0.322961888755395E+03
        p(1,4) = -0.100000000000000d+04
        p(1,1) = 0d0
        p(1,2) = 0d0
        p(1,3) = -0.100000000000000d+04
        p(2,4) = -0.100000000000000d+04
        p(2,1) = 0d0
        p(2,2) = 0d0
        p(2,3) = 0.100000000000000d+04
        p(3,4) = 0.350815654600186d+03
        p(3,1) = 0.274947585251735d+03
        p(3,2) = -0.177057379208442d+03
        p(3,3) = -0.126988713453734d+03
        p(4,4) = 0.602328810143063d+02
        p(4,1) = 0.156351883341439d+02
        p(4,2) = -0.472677797496510d+02
        p(4,3) = -0.339012955294915d+02
        p(5,4) = 0.910714300630089d+03
        p(5,1) = -0.757455408085941d+03
        p(5,2) = -0.146796400460784d+03
        p(5,3) = 0.483851897738626d+03
        p(6,4) = 0.678237163755404d+03
        p(6,1) = 0.466872634500051d+03
        p(6,2) = 0.371121559418886d+03
        p(6,3) = -0.322961888755395d+03
      elseif (n .eq. 4) then
c 0.100000000000000E+04    0.000000000000000E+00    0.000000000000000E+00    0.100000000000000E+04
c 0.100000000000000E+04    0.000000000000000E+00    0.000000000000000E+00   -0.100000000000000E+04
c 0.682544080941019E+03    0.296005212070423E+03    0.428723040630645E+03   -0.440957924622268E+03
c 0.309882149139339E+03    0.766567197312570E+02    0.209302020640822E+03   -0.215275074801816E+03
c 0.204924765556479E+03   -0.140944376091923E+03   -0.867294692272388E+02    0.120858766972024E+03
c 0.802649004363150E+03   -0.231717555709762E+03   -0.551295592044237E+03    0.535374232452070E+03
        p(1,4) =-0.100000000000000d+04 
        p(1,1) =0d0		       
        p(1,2) =0d0		       
        p(1,3) =-0.100000000000000d+04 
        p(2,4) =-0.100000000000000d+04 
        p(2,1) =0d0		       
        p(2,2) =0d0		       
        p(2,3) =0.100000000000000d+04  
        p(3,4) = 0.682544080941019d+03
        p(3,1) = 0.296005212070423d+03
        p(3,2) = 0.428723040630645d+03
        p(3,3) = -0.440957924622268d+03
        p(4,4) = 0.309882149139339d+03
        p(4,1) = 0.766567197312570d+02
        p(4,2) = 0.209302020640822d+03
        p(4,3) = -0.215275074801816d+03
        p(5,4) = 0.204924765556479d+03
        p(5,1) = -0.140944376091923d+03
        p(5,2) = -0.867294692272388d+02
        p(5,3) = 0.120858766972024d+03
        p(6,4) = 0.802649004363150d+03
        p(6,1) = -0.231717555709762d+03
        p(6,2) = -0.551295592044237d+03
        p(6,3) = 0.535374232452070d+03
      else
        write(6,*) 'No such phase space point in ps_check.f'
	stop
      endif
      
      call writeout(p)
      write(6,*) 'Madloop check: PS point # ',n

      return
      end
      
      
