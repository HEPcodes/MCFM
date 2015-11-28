      double complex function a6(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
*     st is a string to choose between pp,pm or sl
*     implementation of Eq. (3.3) of BDKW hep-ph/9610370
*     character string st can take the value pp,pm or sl
************************************************************************
      include 'constants.f'
      include 'zprods_decl.f'
      character*2 st
      integer j1,j2,j3,j4,j5,j6
      double complex atree,vv,ff

      a6=atree(st,j1,j2,j3,j4,j5,j6,za,zb)*vv(st,j1,j2,j3,j4,j5,j6)
      a6=a6+ff(st,j1,j2,j3,j4,j5,j6,za,zb)

      end

      double complex function ff(st,j1,j2,j3,j4,j5,j6,za,zb) 
c-----switchyard function to direct to pm,pp or st
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      character*2 st
      integer j1,j2,j3,j4,j5,j6
      double complex fpp,fpm,fsl
      if     (st .eq. 'pp') then
           ff=fpp(j1,j2,j3,j4,j5,j6,za,zb) 
      elseif (st .eq. 'pm') then
           ff=fpm(j1,j2,j3,j4,j5,j6,za,zb) 
      elseif (st .eq. 'sl') then
           ff=fsl(j1,j2,j3,j4,j5,j6,za,zb) 
      endif
      end



