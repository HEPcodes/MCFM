      subroutine chooser(nproc)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'vegas_common.f'
      include 'zerowidth.f'
      include 'bbproc.f'
      include 'lc.f'
      include 'nwz.f'
      character*4 part
      common/part/part
      double precision rtsmin
      common/rtsmin/rtsmin
      double precision bbbr,gamgambr,wwbr,zzbr,x,y,lambda
      double precision br,BrnRat,brwen,brzee,brtau,brtop
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      integer nproc,nnproc,mproc,j,n2,n3,nqcdjets,nqcdstart,isub
      logical spira
      character*66 pname
      character*2 plabel(mxpart)
      character*6 case
      common/process/case
      common/nproc/nnproc
      common/BrnRat/BrnRat
      common/spira/spira
      common/nqcdjets/nqcdjets,nqcdstart
      common/plabel/plabel
      common/isub/isub
      do j=1,mxpart
      plabel(j)=''
      enddo

c set-up twidth
      x=wmass/mt
      y=mb/mt
      lambda=sqrt((1d0-(x+y)**2)*(1d0-(x-y)**2))
      twidth=gf*mt**3/(8d0*pi*sqrt(2d0))
     . *((1d0-x**2)*(1+2d0*x**2)-y**2*(2d0-x**2-y**2))*lambda

      open(unit=21,file='process.DAT',status='old',err=43)
      call checkversion(21,'process.DAT')
      
      nnproc=nproc
c      write(6,*) 'Chooser:process chosen by nproc=',nproc
      do j=1,300
      read(21,*,err=44) mproc,pname
      if (nproc .lt. 0) then 
      write(6,*) mproc,pname 
      endif

      if (mproc .eq. nproc) go to 42
      enddo

      goto 44

c 42   write(6,*) mproc
 42   write(6,*)
      write(6,*) '**************************************************'//
     . '********************'
      write(6,*) '* ',pname,' *'
      write(6,*) '**************************************************'//
     . '********************'
      write(6,*)

      close(unit=21)

      plabel(1)='pp'
      plabel(2)='pp'

      BrnRat=1d0
      call coupling

      nqcdjets=0
      isub=0
      bbproc=.false.

      if (nproc/10 .eq. 0) then
             case='W_only'
             mass3=wmass
             width3=wwidth
             n3=1
             ndim=4
             BrnRat=1d0
             nqcdjets=0
c---W^+
             if (nproc .lt. 5) then
              plabel(3)='nl'
              plabel(4)='ea'
              plabel(5)='pp'
              nwz=1
c---W^-
             elseif (nproc .ge. 5) then
              plabel(3)='el'
              plabel(4)='na'
              plabel(5)='pp'
              nwz=-1
             endif
c---total cross-section
             if ((nproc .eq. 0) .or. (nproc .eq. 5)) then
             call branch(brwen,brzee,brtau,brtop)
             BrnRat=brwen
             endif

      elseif ((nproc .ge. 10) .and. (nproc .le. 16)) then
             nqcdjets=1
             BrnRat=1d0
             case='W_1jet'
             ndim=7
             mb=0
             n2=0
             n3=1
             mass3=wmass
             width3=wwidth
             plabel(5)='pp'
             plabel(6)='pp'
             BrnRat=1d0
             if (nproc .eq. 10) then
             plabel(3)='nl'
             plabel(4)='ea'
             nwz=1
             call branch(brwen,brzee,brtau,brtop)
             BrnRat=brwen
             elseif (nproc .eq. 11) then
             nwz=1
             plabel(3)='nl'
             plabel(4)='ea'
             elseif (nproc .eq. 16) then
             plabel(3)='el'
             plabel(4)='na'
             nwz=-1
             endif
      elseif ((nproc .eq. 18) .or. (nproc .eq. 19)) then
             nqcdjets=0
             BrnRat=1d0
             case='Wgamma'
             ndim=7
             mb=0
             n2=0
             n3=1
             mass3=wmass
             width3=wwidth
             if (nproc .eq. 18) then
             nwz=1
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='ga'
             plabel(6)='pp'
             elseif (nproc .eq. 19) then
             nwz=-1
             plabel(3)='el'
             plabel(4)='na'
             plabel(5)='ga'
             plabel(6)='pp'
             write(6,*) 'chooser',case
             endif

      elseif (nproc/10 .eq. 2) then
             nqcdjets=2
             if (nproc .eq. 20) then
             case='Wbbmas'
             write(6,*) 'mb=',mb
             bbproc=.true.
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='bq'
             plabel(6)='ba'
             plabel(7)='pp'
             nwz=1
             ndim=10
             n2=0
             n3=1
             mass3=wmass
             width3=wwidth
             endif
             if (nproc .eq. 21) then
             case='Wbbbar'
             bbproc=.true.
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='bq'
             plabel(6)='ba'
             plabel(7)='pp'
             nwz=1
             ndim=10
             mb=0
             n2=0
             n3=1
             mass3=wmass
             width3=wwidth
             endif
             if (nproc .eq. 22) then
             case='W_2jet'
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='pp'
             plabel(6)='pp'
             plabel(7)='pp'
             nwz=1
             ndim=10
             mb=0
             n2=0
             n3=1
             mass3=wmass
             width3=wwidth
             endif
             if (nproc .eq. 23) then
             nqcdjets=3
             case='W_3jet'
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='pp'
             plabel(6)='pp'
             plabel(7)='pp'
             nwz=1
             ndim=13
             mb=0
             n2=0
             n3=1
             mass3=wmass
             width3=wwidth
             endif

             if (nproc .eq. 26) then
             case='Wbbbar'
             bbproc=.true.
             nwz=-1
             plabel(3)='el'
             plabel(4)='na'
             plabel(5)='bq'
             plabel(6)='ba'
             plabel(7)='pp'
             ndim=10
             mb=0
             n2=0
             n3=1
             mass3=wmass
             width3=wwidth
             endif
             if (nproc .eq. 27) then
             case='W_2jet'
             plabel(3)='el'
             plabel(4)='na'
             plabel(5)='pp'
             plabel(6)='pp'
             plabel(7)='pp'
             nwz=-1
             ndim=10
             mb=0
             n2=0
             n3=1
             mass3=wmass
             width3=wwidth
             endif

      elseif (nproc/10 .eq. 3) then
             nqcdjets=0
             case='Z_only'
             nwz=0
             mass3=zmass
             width3=zwidth
             n3=1
             ndim=4
             BrnRat=1d0
             if (nproc .eq. 30) then
             plabel(3)='el'
             plabel(4)='ea'
             plabel(5)='pp'
             q1=-1d0
             l1=le
             r1=re
             call branch(brwen,brzee,brtau,brtop)
             BrnRat=brzee
             elseif (nproc .eq. 31) then
             plabel(3)='el'
             plabel(4)='ea'
             plabel(5)='pp'
             q1=-1d0
             l1=le
             r1=re
             elseif (nproc .eq. 32) then
             plabel(3)='nl'
             plabel(4)='na'
             plabel(5)='pp'
             q1=0d0
             l1=ln*sqrt(3d0)
             r1=rn*sqrt(3d0)
             elseif (nproc .eq. 33) then
             plabel(3)='qb'
             plabel(4)='ab'
             plabel(5)='pp'
             q1=Q(5)*sqrt(xn)
             l1=l(5)*sqrt(xn)
             r1=r(5)*sqrt(xn)
             endif 

      elseif (nproc/10 .eq. 4) then
             nqcdjets=1
             case='Z_1jet'
             nwz=0
             ndim=7
             mb=0
             n2=0
             n3=1
             mass3=zmass
             width3=zwidth
             if (nproc .eq. 40) then
             plabel(3)='el'
             plabel(4)='ea'
             plabel(5)='pp'
             plabel(6)='pp'
             q1=-1d0
             l1=le
             r1=re
             call branch(brwen,brzee,brtau,brtop)
             BrnRat=brzee
             elseif (nproc .eq. 41) then
             nqcdjets=1
             plabel(3)='el'
             plabel(4)='ea'
             plabel(5)='pp'
             plabel(6)='pp'
             q1=-1d0
             l1=le
             r1=re
             elseif (nproc .eq. 42) then
             nqcdjets=1
             plabel(3)='nl'
             plabel(4)='na'
             plabel(5)='pp'
             plabel(6)='pp'
             q1=0d0
             l1=ln*sqrt(3d0)
             r1=rn*sqrt(3d0)
             elseif (nproc .eq. 43) then
             plabel(3)='qb'
             plabel(4)='ab'
             plabel(5)='pp'
             plabel(6)='pp'
             q1=Q(5)*sqrt(xn)
             l1=l(5)*sqrt(xn)
             r1=r(5)*sqrt(xn)
             elseif (nproc .eq. 44) then
             case='Z_2jet'
             ndim=10
             nqcdjets=2
             plabel(3)='el'
             plabel(4)='ea'
             plabel(5)='pp'
             plabel(6)='pp'
             plabel(7)='pp'
             q1=-1d0
             l1=le
             r1=re
             elseif (nproc .eq. 45) then
             nqcdjets=3
             case='Z_3jet'
             ndim=13
             nqcdjets=3
             plabel(3)='el'
             plabel(4)='ea'
             plabel(5)='pp'
             plabel(6)='pp'
             plabel(7)='pp'
             q1=-1d0
             l1=le
             r1=re
             elseif (nproc .eq. 48) then
             nqcdjets=0
             BrnRat=1d0
             case='Zgamma'
             q1=-1d0
             l1=le
             r1=re
             ndim=7
             mb=0
             n2=0
             n3=1
             mass3=zmass
             width3=zwidth
             nwz=0
             plabel(3)='el'
             plabel(4)='ea'
             plabel(5)='ga'
             plabel(6)='pp'
             write(6,*) 'chooser',case
             endif
      elseif (nproc/10 .eq. 5) then
             case='Zbbbar'
             colourchoice=0
             ndim=10
             n2=0
             n3=1
             mass3=zmass
             width3=zwidth
             if     (nproc .eq. 50) then
             case='Zbbmas'
               write(6,*) 'mb=',mb
               bbproc=.true.
               nqcdjets=2
               plabel(3)='el'
               plabel(4)='ea'
               plabel(5)='bq'
               plabel(6)='ba'
               plabel(7)='pp'
               q1=-1d0
               l1=le
               r1=re
             elseif (nproc .eq. 51) then
               bbproc=.true.
               mb=0d0
               nqcdjets=2
               plabel(3)='el'
               plabel(4)='ea'
               plabel(5)='bq'
               plabel(6)='ba'
               plabel(7)='pp'
               q1=-1d0
               l1=le
               r1=re
             elseif (nproc .eq. 52) then
               bbproc=.true.
               mb=0d0
               nqcdjets=2
               plabel(3)='nl'
               plabel(4)='na'
               plabel(5)='bq'
               plabel(6)='ba'
               plabel(7)='pp'
               q1=0d0
               l1=ln*sqrt(3d0)
               r1=rn*sqrt(3d0)
             elseif (nproc .eq. 53) then
               bbproc=.true.
               mb=0d0
               nqcdjets=2
               plabel(3)='qb'
               plabel(4)='ab'
               plabel(5)='bq'
               plabel(6)='ba'
               plabel(7)='pp'
               q1=Q(5)*sqrt(xn)
               l1=l(5)*sqrt(xn)
               r1=r(5)*sqrt(xn)
             else
             write(6,*) 'z+2jet'
             write(6,*) 'not implemented yet'
             pause
             endif

c---case WW
      elseif (nproc/10 .eq. 6) then
             call readcoup
             nqcdjets=0  
             case='WWqqbr'
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='el'
             plabel(6)='na'
             plabel(7)='pp'
             nwz=1
             ndim=10
             mb=0
             n2=1
             n3=1
             mass2=wmass
             width2=wwidth
             mass3=wmass
             width3=wwidth
             l1=1d0
             if    (nproc .eq. 60) then
               call branch(brwen,brzee,brtau,brtop)
               BrnRat=brwen**2
             elseif (nproc .eq. 62) then
c--- note: scattering diagrams are NOT included, only couplings change
               plabel(3)='nl'
               plabel(4)='ea'
               plabel(5)='qq'
               plabel(6)='aa'
               plabel(7)='pp'
               l1=dsqrt(xn*2d0)
             elseif (nproc .eq. 63) then
c--- note: scattering diagrams are NOT included, only couplings change
               plabel(3)='qq'
               plabel(4)='aa'
               plabel(5)='el'
               plabel(6)='na'
               plabel(7)='pp'
               l1=dsqrt(xn*2d0)
             endif
c---case WZ
      elseif (nproc/10 .eq. 7) then
             call readcoup
             nqcdjets=0
             case='WZbbar'
             ndim=10
             mb=0
             n2=1
             n3=1
             mass2=zmass
             width2=zwidth
             mass3=wmass
             width3=wwidth
             if (nproc .lt. 75) then
c---case W^+Z
             nwz=+1
             if (nproc .eq. 70) then
             q1=-1d0
             l1=le
             r1=re
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='ml'
             plabel(6)='ma'
             plabel(7)='pp'
             call branch(brwen,brzee,brtau,brtop)
             BrnRat=brwen*brzee
             elseif (nproc .eq. 71) then             
               plabel(3)='nl'
               plabel(4)='ea'
               plabel(5)='ml'
               plabel(6)='ma'
               plabel(7)='pp'
               q1=-1d0
               l1=le
               r1=re
             elseif (nproc .eq. 72) then
               plabel(3)='nl'
               plabel(4)='ea'
               plabel(5)='nl'
               plabel(6)='na'
               plabel(7)='pp'
               q1=0d0
               l1=ln
               r1=rn
            elseif (nproc .eq. 73) then
               bbproc=.true.
               plabel(3)='nl'
               plabel(4)='ea'
               plabel(5)='qb'
               plabel(6)='ab'
               plabel(7)='pp'
               q1=Q(5)*sqrt(xn)
               l1=l(5)*sqrt(xn)
               r1=r(5)*sqrt(xn)
             endif 
             elseif (nproc .ge. 75) then
c---case W^-Z
             nwz=-1
             if (nproc .eq. 75) then
               plabel(3)='el'
               plabel(4)='na'
               plabel(5)='ml'
               plabel(6)='ma'
               plabel(7)='pp'
               q1=-1d0
               l1=le
               r1=re
             call branch(brwen,brzee,brtau,brtop)
             BrnRat=brwen*brzee
             elseif (nproc .eq. 76) then
               plabel(3)='el'
               plabel(4)='na'
               plabel(5)='ml'
               plabel(6)='ma'
               plabel(7)='pp'
               q1=-1d0
               l1=le
               r1=re
             elseif (nproc .eq. 77) then
               plabel(3)='el'
               plabel(4)='na'
               plabel(5)='nl'
               plabel(6)='na'
               plabel(7)='pp'
               q1=0d0
               l1=ln
               r1=rn
             elseif (nproc .eq. 78) then
               bbproc=.true.
               plabel(3)='el'
               plabel(4)='na'
               plabel(5)='qb'
               plabel(6)='ab'
               plabel(7)='pp'
               q1=Q(5)*sqrt(xn)
               l1=l(5)*sqrt(xn)
               r1=r(5)*sqrt(xn)
             endif 

             endif
            
c---case ZZ
      elseif (nproc/10 .eq. 8) then
             call readcoup
             nqcdjets=0
             case='ZZlept'
             nwz=1
             ndim=10
             mb=0
             n2=1
             n3=1
             mass2=zmass
             width2=zwidth
             mass3=zmass
             width3=zwidth
             q1=-1d0
             l1=le
             r1=re
             if (nproc .eq. 80 .or. nproc .eq. 85) then
               plabel(3)='el'
               plabel(4)='ea'
               plabel(5)='ml'
               plabel(6)='ma'
               plabel(7)='pp'
               q2=-1d0
               l2=le
               r2=re
               call branch(brwen,brzee,brtau,brtop)
               BrnRat=2d0*brzee**2
             elseif (nproc .eq. 81 .or. nproc .eq. 86) then
               plabel(3)='el'
               plabel(4)='ea'
               plabel(5)='ml'
               plabel(6)='ma'
               plabel(7)='pp'
               q2=-1d0
               l2=le
               r2=re
             elseif (nproc .eq. 82 .or. nproc .eq. 87) then
               plabel(3)='el'
               plabel(4)='ea'
               plabel(5)='nl'
               plabel(6)='na'
               plabel(7)='pp'
               q2=0d0
               l2=ln*sqrt(3d0)
               r2=rn*sqrt(3d0)
             elseif (nproc .eq. 83 .or. nproc .eq. 88) then
               plabel(3)='el'
               plabel(4)='ea'
               plabel(5)='qb'
               plabel(6)='ab'
               plabel(7)='pp'
               q2=Q(5)*sqrt(xn)
               l2=l(5)*sqrt(xn)
               r2=r(5)*sqrt(xn)
             elseif ((nproc .eq. 84) .or. (nproc .eq. 89))  then
               bbproc=.true.
               plabel(3)='nl'
               plabel(4)='na'
               plabel(5)='qb'
               plabel(6)='ab'
               plabel(7)='pp'
               q2=0d0
               l2=ln*sqrt(3d0)
               r2=rn*sqrt(3d0)
               q1=Q(5)*sqrt(xn)
               l1=l(5)*sqrt(xn)
               r1=r(5)*sqrt(xn)
             endif 

             if (nproc .ge. 85) then
               q1=0d0
             endif
c---case WH
      elseif (nproc/10 .eq. 9) then
             nqcdjets=2
             nqcdstart=5
             case='WHbbar'
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='bq'
             plabel(6)='ba'
             plabel(7)='pp'

             if (spira) then
                 call higgsp(br,gamgambr,wwbr,zzbr)
             else
                 call higgsw(br)
             endif 
c--reset higgs width so that whole gives correct BR
             hwidth=xn*sqrt(2d0)/8d0/pi*gf*hmass*mbsq
     .         *(1d0-4d0*mbsq/hmass**2)**1.5d0/br

             ndim=10
             mb=0
             n2=1
             n3=1
             mass2=hmass
             width2=hwidth
             mass3=wmass
             width3=wwidth

             if (nproc .eq. 90) then
               nwz=1
               call branch(brwen,brzee,brtau,brtop)
               BrnRat=brwen*br
             elseif (nproc .eq. 91) then
               bbproc=.true.
               nqcdjets=2
               nwz=1
             elseif (nproc .eq. 96) then
               bbproc=.true.
               nqcdjets=2
               nwz=-1
             endif
             write(6,99) hmass,hwidth,br
c---case ZH
      elseif (nproc/10 .eq. 10) then
             case='ZHbbar'
             nqcdjets=2
             nqcdstart=5
             plabel(3)='el'
             plabel(4)='ea'
             plabel(5)='bq'
             plabel(6)='ba'
             plabel(7)='pp'
             if (spira) then
                 call higgsp(br,gamgambr,wwbr,zzbr)
C--reset higgs width so that whole gives coreect BR
             else
                 call higgsw(br)
             endif 
             hwidth=xn*sqrt(2d0)/8d0/pi*gf*hmass*mbsq
     .         *(1d0-4d0*mbsq/hmass**2)**1.5d0/br
             write(6,99) hmass,hwidth,br

             ndim=10
             nwz=0
             mb=0
             n2=1
             n3=1
             mass2=hmass
             width2=hwidth
             mass3=zmass
             width3=zwidth

             if (nproc .eq. 100) then
               plabel(3)='el'
               plabel(4)='ea'
               plabel(5)='bq'
               plabel(6)='ba'
               plabel(7)='pp'
               nqcdjets=2
               q1=-1d0
               l1=le
               r1=re
               call branch(brwen,brzee,brtau,brtop)
               BrnRat=brzee*br
             elseif (nproc .eq. 101) then
               bbproc=.true.
               plabel(3)='el'
               plabel(4)='ea'
               plabel(5)='bq'
               plabel(6)='ba'
               plabel(7)='pp'
               nqcdjets=2
               q1=-1d0
               l1=le
               r1=re
             elseif (nproc .eq. 102) then
               bbproc=.true.
               plabel(3)='nl'
               plabel(4)='na'
               plabel(5)='bq'
               plabel(6)='ba'
               plabel(7)='pp'
               nqcdjets=2
               q1=0d0
               l1=ln*sqrt(3d0)
               r1=rn*sqrt(3d0)
             elseif (nproc .eq. 103) then
C--- NEED TO LOOK AT THIS TO CLUSTER 4 JETS STARTING FROM 3
               nqcdjets=4
               plabel(3)='bq'
               plabel(4)='ba'
               plabel(5)='bq'
               plabel(6)='ba'
               plabel(7)='pp'
               q1=Q(5)*sqrt(xn)
               l1=l(5)*sqrt(xn)
               r1=r(5)*sqrt(xn)
             endif 

      elseif (nproc/10 .eq. 11) then
             if (spira) then
                 call higgsp(br,gamgambr,wwbr,zzbr)
             else
                 call higgsw(br)
             endif 
             write(6,99) hmass,hwidth,br

             if ((nproc .eq. 110) .or. (nproc .eq. 111)) then
c--reset higgs width so that whole gives correct BR
             hwidth=xn*sqrt(2d0)/8d0/pi*gf*hmass*mbsq
     .         *(1d0-4d0*mbsq/hmass**2)**1.5d0/br

             case='Hbbbar'
             plabel(3)='bq'
             plabel(4)='ba'
             nqcdjets=2
             nqcdstart=3
             mass3=hmass
             width3=hwidth
             n3=1
             ndim=4
             endif
             if (nproc .eq. 110) BrnRat=br

             if (nproc .ge. 112) then
             case='HWW_4l'             
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='el'
             plabel(6)='na'
             plabel(7)='pp'
             nqcdjets=0
             nqcdstart=7
             write(6,*) 'nproc=',nproc
             write(6,*) 'case=',case
             ndim=10
             n2=1
             n3=1
             mass2=wmass
             width2=wwidth
             mass3=wmass
             width3=wwidth
             endif
             if (nproc .eq. 112) then
             call branch(brwen,brzee,brtau,brtop)
             BrnRat=brwen**2
             endif
      elseif (nproc/10 .eq. 12) then
             case='HZZ_4l'
             plabel(3)='el'
             plabel(4)='ea'
             plabel(5)='ml'
             plabel(6)='ma'
             plabel(7)='pp'
             nqcdjets=0
             nqcdstart=7
             write(6,*) 'nproc=',nproc
             write(6,*) 'case=',case
             ndim=10
             n2=1
             n3=1
             mass2=zmass
             width2=zwidth
             mass3=zmass
             width3=zwidth
             if (spira) then
                 call higgsp(bbbr,gamgambr,wwbr,br)
             else
                 call higgsw(br)
             endif 
             write(6,99) hmass,hwidth,br

             if (nproc .eq. 121) then
               plabel(3)='el'
               plabel(4)='ea'
               plabel(5)='ml'
               plabel(6)='ma'
               plabel(7)='pp'
               l1=le
               r1=re
               l2=le
               r2=re
             elseif (nproc .eq. 122) then
               plabel(3)='nl'
               plabel(4)='na'
               plabel(5)='ml'
               plabel(6)='ma'
               plabel(7)='pp'
               l1=le
               r1=re
               l2=ln*sqrt(3d0)
               r2=rn*sqrt(3d0)
             elseif (nproc .eq. 123) then
               nqcdjets=2
               plabel(3)='ml'
               plabel(4)='ma'
               plabel(5)='bq'
               plabel(6)='ba'
               plabel(7)='pp'
               l1=le
               r1=re
               l2=l(5)*sqrt(3d0)
               r2=r(5)*sqrt(3d0)
             endif
      elseif (nproc/10 .eq. 14) then
             case='H_1jet'
             ndim=7
             nqcdstart=5
             plabel(3)='ml'
             plabel(4)='ma'
             plabel(5)='bq'
             nqcdjets=1

             if (spira) then
                 call higgsp(bbbr,gamgambr,wwbr,br)
             else
                 call higgsw(br)
             endif
C--reset higgs width so that whole gives correct BR, (if we have set
C mb=0, keeps mbsq in the coupling but ignores it in the phase space) 
             hwidth=xn*sqrt(2d0)/8d0/pi*gf*hmass*mbsq
     .         *(1d0-4d0*mb**2/hmass**2)**1.5d0/br
             write(6,99) hmass,hwidth,br

             n2=0
             n3=1
             mass3=hmass
             width3=hwidth
             BrnRat=br
             
             if     (nproc .eq. 141) then
             isub=1
             plabel(6)='pp'
             elseif (nproc .eq. 142) then
             isub=2
             plabel(5)='ba'
             plabel(6)='bq'
             elseif (nproc .eq. 143) then
             isub=3
             plabel(5)='ba'
             plabel(6)='qj'
             elseif (nproc .eq. 144) then
             isub=4
             plabel(5)='bq'
             plabel(6)='ba'
             elseif (nproc .eq. 145) then
             isub=2
             plabel(5)='ba'
             plabel(6)='bq'
             nqcdjets=2
             endif
             
      elseif (nproc .eq. 146) then
             case='H_1jet'
             ndim=7
             nqcdstart=3
             plabel(3)='bq'
             plabel(4)='ba'
             plabel(5)='bq'
             plabel(6)='pp'
             nqcdjets=3

             if (spira) then
                 call higgsp(bbbr,gamgambr,wwbr,br)
             else
                 call higgsw(br)
             endif
C--reset higgs width so that whole gives correct BR, (if we have set
C mb=0, keeps mbsq in the coupling but ignores it in the phase space) 
             hwidth=xn*sqrt(2d0)/8d0/pi*gf*hmass*mbsq
     .         *(1d0-4d0*mb**2/hmass**2)**1.5d0/br
             write(6,99) hmass,hwidth,br

             n2=0
             n3=1
             mass3=hmass
             width3=hwidth
             BrnRat=1d0

      elseif ((nproc .ge. 150) .and. (nproc .le. 152)) then
             nwz=1
             ndim=16
             n2=1
             n3=1
             mass2=mt
             width2=twidth
             mass3=mt
             width3=twidth
             if (nproc .eq. 150) then
                plabel(3)='nl'
                plabel(4)='ea'
                plabel(5)='bq'
                plabel(6)='ba'
                plabel(7)='el'
                plabel(8)='na'
                case='tt_bbl'
                write(6,*) 'nproc=',nproc
                write(6,*) 'case=',case
                call branch(brwen,brzee,brtau,brtop)
                BrnRat=(brwen*brtop)**2
             elseif (nproc .eq. 151) then
                bbproc=.true.
                plabel(3)='nl'
                plabel(4)='ea'
                plabel(5)='ea'
                plabel(6)='ea'
                plabel(7)='el'
                plabel(8)='na'
                plabel(9)='el'
c----debug
                nqcdjets=0
                nqcdstart=5
                case='tt_bbl'
                write(6,*) 'nproc=',nproc
                write(6,*) 'case=',case
             elseif (nproc .eq. 152) then
                bbproc=.true.
                plabel(3)='nl'
                plabel(4)='ea'
                plabel(5)='bq'
                plabel(6)='ba'
                plabel(7)='pj'
                plabel(8)='pj'
                nqcdjets=4
                nqcdstart=5
c--need to extend--not quite correct
                case='tt_bbh'
                write(6,*) 'nproc=',nproc
                write(6,*) 'case=',case
             endif 
      elseif ((nproc .ge. 153) .and. (nproc .le. 154)) then
             nwz=1
             ndim=14
             n2=1
             n3=1
             mass2=mt
             width2=twidth
             mass3=mt
             width3=twidth
             if (nproc .eq. 153) then
                bbproc=.true.
                plabel(3)='nl'
                plabel(4)='ea'
                plabel(5)='ea'
                plabel(6)='ea'
                plabel(7)='el'
                plabel(8)='na'
                plabel(9)='el'
c----debug
                nqcdjets=0
                nqcdstart=5
                case='ttbdkl'
                write(6,*) 'nproc=',nproc
                write(6,*) 'case=',case
             elseif (nproc .eq. 154) then
                bbproc=.true.
                plabel(3)='nl'
                plabel(4)='ea'
                plabel(5)='bq'
                plabel(6)='ba'
                plabel(7)='pj'
                plabel(8)='pj'
                nqcdjets=4
                nqcdstart=5
c--need to extend--not quite correct
                case='ttbdkh'
                write(6,*) 'nproc=',nproc
                write(6,*) 'case=',case
             endif
      elseif (nproc .eq. 156) then
             nwz=1
             ndim=19
             mb=0
c            nqcdjets=3
c            nqcdstart=5
             case='qq_ttg'
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='bq'
             plabel(6)='ba'
             plabel(7)='el'
             plabel(8)='na'
             plabel(9)='pj'
      elseif (nproc .eq. 157) then
c--- for now, pretend the heavy quarks are electrons, to avoid
c--- any jet cuts taking place
             case='tt_tot'
             nqcdjets=0
             ndim=4
             mass2=mt
             plabel(3)='ea'
             plabel(4)='ea'
             plabel(5)='pp'
      elseif (nproc .eq. 158) then
c--- for now, pretend the heavy quarks are electrons, to avoid
c--- any jet cuts taking place
             case='bb_tot'
             nqcdjets=0
             ndim=4
             mass2=mb
             plabel(3)='ea'
             plabel(4)='ea'
             plabel(5)='pp'
      elseif (nproc .eq. 159) then
c--- for now, pretend the heavy quarks are electrons, to avoid
c--- any jet cuts taking place
             case='cc_tot'
             nqcdjets=0
             ndim=4
             mass2=mc
             plabel(3)='ea'
             plabel(4)='ea'
             plabel(5)='pp'
 
      elseif (nproc/10 .eq. 16) then
             bbproc=.true.
             nqcdjets=3
             nqcdstart=5
             case='qg_tbb'
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='bq'
             plabel(6)='ba'
             plabel(7)='pj'
             write(6,*) 'nproc=',nproc
             write(6,*) 'case=',case
             nwz=1
             ndim=13
             mb=0
             mass2=mt
             width2=twidth
             mass3=mt
             width3=twidth
      elseif (nproc/10 .eq. 17) then
             bbproc=.true.
             nqcdjets=2
             nqcdstart=5
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='bq'
             plabel(6)='ba'
             case='t_bbar'
             rtsmin=100d0
             write(6,*) 'Setting rtsmin to',rtsmin
             nwz=1
             ndim=10
             mb=0
             n2=0
             n3=1
             mass3=mt
             width3=twidth
      elseif (nproc/10 .eq. 18) then
             case='tautau'
C----not correct yet!
             plabel(3)='nl'
             plabel(4)='ea'
             plabel(5)='na'
             plabel(6)='na'
             plabel(7)='na'
             plabel(8)='na'
             write(6,*) 'nproc=',nproc
             write(6,*) 'case=',case
             nqcdjets=0
             nwz=1
             ndim=16
             n2=1
             n3=1
             mass2=mtau
             width2=tauwidth
             mass3=mtau
             width3=tauwidth
             if (nproc .eq. 180) then
             call branch(brwen,brzee,brtau,brtop)
             BrnRat=brtau**2
             endif
      elseif (nproc .eq. 190) then
             case='tottth'
             n2=0
             n3=0
             write(6,*) 'nproc=',nproc
             write(6,*) 'case=',case
                nwz=1
                ndim=7
      elseif ((nproc .eq. 191) .or. (nproc .eq. 192)) then
             case='qq_tth'
                plabel(3)='nl'
                plabel(4)='ea'
                plabel(5)='bq'
                plabel(6)='ba'
                plabel(7)='el'
                plabel(8)='na'
                plabel(9)='bq'
                plabel(10)='ba'

             write(6,*) 'nproc=',nproc
             write(6,*) 'case=',case
                nwz=1
                ndim=22
                n2=1
                n3=1
             mass2=mt
             width2=twidth
             mass3=mt
             width3=twidth

             if (spira) then
                 call higgsp(br,gamgambr,wwbr,zzbr)
             else
                 call higgsw(br)
             endif 
C--reset higgs width so that whole gives correct BR, (if we have set
C mb=0, keeps mbsq in the coupling but ignores it in the phase space) 
             hwidth=xn*sqrt(2d0)/8d0/pi*gf*hmass*mbsq
     .         *(1d0-4d0*mb**2/hmass**2)**1.5d0/br
             write(6,99) hmass,hwidth,br
             write(6,*) 'chooser:mbsq',mbsq
             if (nproc .eq. 191) then 
             call branch(brwen,brzee,brtau,brtop)
             BrnRat=brwen**2*br
             endif


      elseif (nproc .eq. 193) then
             case='tottth'
             n2=0
             n3=0
             write(6,*) 'nproc=',nproc
             write(6,*) 'case=',case
                nwz=1
                ndim=7
      elseif ((nproc .eq. 194)
     .   .or. (nproc .eq. 195)
     .   .or. (nproc .eq. 196)) then
             case='qq_ttz'
                plabel(3)='nl'
                plabel(4)='ea'
                plabel(5)='bq'
                plabel(6)='ba'
                plabel(7)='el'
                plabel(8)='na'

             write(6,*) 'nproc=',nproc
             write(6,*) 'case=',case
                nwz=1
                ndim=22
                n2=1
                n3=1
             mass2=mt
             width2=twidth
             mass3=mt
             width3=twidth

             if (nproc .eq. 194) then 
                plabel(9)='el'
                plabel(10)='ea'
                q1=-1d0
                l1=le
                r1=re
                call branch(brwen,brzee,brtau,brtop)
                BrnRat=brwen**2*brzee
             elseif (nproc .eq. 195) then 
                plabel(9)='el'
                plabel(10)='ea'
                q1=-1d0
                l1=le
                r1=re
             elseif (nproc .eq. 196) then 
                plabel(9)='bq'
                plabel(10)='ba'
                q1=Q(5)*sqrt(xn)
                l1=l(5)*sqrt(xn)
                r1=r(5)*sqrt(xn)
             endif
      elseif (nproc .eq. 198) then
             case='qqttbb'
             ndim=10
      elseif (nproc .eq. 199) then
             case='qqttgg'
             ndim=10
      elseif (nproc/10 .eq. 20) then
             case='httjet'
             nqcdjets=1
             nqcdstart=5
             plabel(3)='tl'
             plabel(4)='ta'
             plabel(5)='pj'
             ndim=7
             mb=0
             n2=0 
             n3=1
             write(6,*) 'nproc=',nproc
             write(6,*) 'case=',case
             if (spira) then
                 call higgsp(br,gamgambr,wwbr,zzbr)
             else
                 call higgsw(br)
             endif
C--reset higgs width so that whole gives correct BR, (if we have set
C mb=0, keeps mbsq in the coupling but ignores it in the phase space) 
             hwidth=xn*sqrt(2d0)/8d0/pi*gf*hmass*mbsq
     .         *(1d0-4d0*mbsq/hmass**2)**1.5d0/br
             write(6,99) hmass,hwidth,br
             write(6,*) 'chooser:mbsq',mbsq
             mass3=hmass
             width3=hwidth
             if     (nproc .eq. 200) then
               BrnRat=br
             elseif (nproc .eq. 201) then
c---first work out width into a tau pair (EHSV, 1.3)
               BrnRat=gwsq/4d0/pi*mtausq*hmass/8d0/wmass**2*
     .                 (1d0-4d0*mtausq/hmass**2)**1.5d0
               Brnrat=BrnRat/hwidth
               write(6,*) 'Branching ratio of higgs -> tau tau: ',Brnrat
c---now compensate for b-bbar branching ratio
               Brnrat=br/Brnrat
             elseif (nproc .eq. 202) then
             ndim=19
             case='hlljet'
             plabel(3)='nl'
             plabel(4)='na'
             plabel(5)='el'
             plabel(6)='na'
             plabel(7)='nl'
             plabel(8)='ea'
             plabel(9)='pj'
             nqcdjets=1
             nqcdstart=9
             Brnrat=1d0
             elseif (nproc .eq. 205) then
               case='attjet'
               BrnRat=br
             elseif (nproc .eq. 206) then
               case='attjet'
c---first work out width into a tau pair (EHSV, 1.3)
               BrnRat=gwsq/4d0/pi*mtausq*hmass/8d0/wmass**2*
     .                 (1d0-4d0*mtausq/hmass**2)**1.5d0
               Brnrat=BrnRat/hwidth
               write(6,*) 'Branching ratio of higgs -> tau tau: ',Brnrat
c---now compensate for b-bbar branching ratio
               Brnrat=br/Brnrat
             endif
      elseif (nproc .eq. 240) then
             ndim=7 
             case='threeb'
             plabel(3)='bq'
             plabel(4)='bq'
             plabel(5)='ba'
      elseif (nproc/10 .ge. 30) then
             write(6,*) 'Setting part to lord and zerowidth to false'
             zerowidth=.false.
             part='lord'
          if (nproc .eq. 302) then
             case='vlchk2'
             nwz=1
             ndim=4
             n3=1
             mass3=wmass
             width3=wwidth
          elseif (nproc .eq. 303) then
             case='vlchk3'
             nwz=1
             ndim=7
             n3=1
             mass3=wmass
             width3=wwidth
          elseif (nproc .eq. 304) then
             case='vlchk4'
             nwz=1
             ndim=10
             n2=1
             n3=1
             mass2=hmass
             width2=hwidth
             mass3=wmass
             width3=wwidth
          elseif (nproc .eq. 305) then
             case='vlchk5'
             nwz=1
             ndim=13
             n2=1
             n3=1
             mass2=hmass
             width2=hwidth
             mass3=wmass
             width3=wwidth
          elseif (nproc .eq. 306) then
             case='vlchk6'
             nwz=1
             ndim=16
             n2=1
             n3=1
             mass2=mt
             width2=twidth
             mass3=mt
             width3=twidth
          elseif (nproc .eq. 308) then
             case='vlchk8'
             nwz=1
             ndim=22
             n2=1
             n3=1
             mass2=mt
             width2=twidth
             mass3=mt
             width3=twidth
          elseif (nproc .eq. 309) then
             case='vlchkm'
             write(6,*) 'mb=',mb
             nwz=1
             ndim=10
             n2=1
             n3=1
             mass2=hmass
             width2=hwidth
             mass3=wmass
             width3=wwidth
          elseif (nproc .eq. 310) then
             case='vlchm3'
             write(6,*) 'mt=',mt
             nwz=1
             ndim=7
             n2=0
             n3=0
             mass2=mt
             width2=twidth
             mass3=mt
             width3=twidth
          endif
      else 
      write(6,*) 'chooser:Unimplemented case'
      write(6,*) 'nproc=',nproc
      stop
      endif

      call ckmfill(nwz)

      return

 43   write(6,*) 'problems opening process.DAT'
      stop

 44   write(6,*) 'Unimplemented process number'
      stop
      
 99   format(/,
     .       ' ****************** Higgs parameters ****************'/, 
     .       ' *                                                  *'/, 
     .       ' *   mass(H) = ',f7.2,'      width(H) = ',e12.5,' *'/,
     .       ' *              Br( H -> b bbar) = ',f8.4,'         *'/,
     .       ' ****************************************************')
     
      end
