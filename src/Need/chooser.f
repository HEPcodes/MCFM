      subroutine chooser
c---- Note added 4/21/03
c---- plabel set to 'ig' (for 'ignore') means that this
c---- particle should not be subject to any cuts, so that the
c---- total cross-section comes out correctly when the BR is removed
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'vegas_common.f'
      include 'zerowidth.f'
      include 'removebr.f'
      include 'bbproc.f'
      include 'lc.f'
      include 'nwz.f'
      include 'process.f'
      include 'flags.f'
      include 'heavyflav.f'
      include 'nflav.f'
      character*4 part
      common/part/part
      double precision wwbr,zzbr,tautaubr
      double precision br,BrnRat,brwen,brzee,brtau,brtop
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      integer nproc,mproc,j,n2,n3,nqcdjets,nqcdstart,isub,notag
      character*83 pname
      character*2 plabel(mxpart)
      character*72 string
      character*30 runstring
      common/runstring/runstring
      common/nproc/nproc
      common/BrnRat/BrnRat
      common/nqcdjets/nqcdjets,nqcdstart
      common/plabel/plabel
      common/isub/isub
      common/notag/notag

      do j=1,mxpart
      plabel(j)=''
      enddo

      string='process.DAT' 
      open(unit=21,file=string,status='old',err=43)
      call checkversion(21,string)
      
      write(6,*) 'Chooser:process chosen by nproc=',nproc

      do j=1,400
      read(21,*,err=44) mproc,pname
      
      if (nproc .lt. 0) then 
      write(6,*) mproc,pname 
      endif

      if (mproc .eq. nproc) go to 42
      if (pname .eq. 'EOF') go to 44
      enddo
      goto 44

 42   write(6,*)
      write(6,*) '*************************** f(p1)+f(p2) --> *****'//
     . '********************'
      write(6,*) '* ',pname(19:83),' *'
      write(6,*) '*************************************************'//
     . '********************'
      write(6,*)

      close(unit=21)

      plabel(1)='pp'
      plabel(2)='pp'

c--- the default behaviour is to remove no branching ratio
      BrnRat=1d0

      call coupling

      notag=0
      nqcdjets=0
      isub=0
      bbproc=.false.

c-----------------------------------------------------------------------

      if (nproc/10 .eq. 0) then
        case='W_only'
        mass3=wmass
        width3=wwidth
        n3=1
        ndim=4
        nqcdjets=0
c---W^+
        if     (nproc .eq. 1) then
C-- 1  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))'
C--    '  f(p1)+f(p2) --> W^+ (for total Xsect)' (removebr=.true.)
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='pp'
          nwz=1
c---W^-
        elseif (nproc .eq. 6) then
c-- 6  '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))'
c--    '  f(p1)+f(p2) --> W^- (for total Xsect)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='pp'
          nwz=-1
        else
          call nprocinvalid()
        endif

c--- total cross-section
        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 11) .or. (nproc .eq. 16)) then
        case='W_1jet'
        nqcdjets=1
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        plabel(5)='pp'
        plabel(6)='pp'

        if     (nproc .eq. 11) then
C-- 11 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+f(p5)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) + f(p5)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 16) then
c-- 16 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+f(p5)'
c--    '  f(p1)+f(p2) --> W^- (no BR) + f(p5)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 12) .or. (nproc .eq. 17)) then
        case='Wgamma'
        nqcdjets=0
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        plabel(5)='ga'
        plabel(6)='pp'

        if     (nproc .eq. 12) then
c-- 12 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+gamma(p5)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) + gamma(p5)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 17) then
c-- 17 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~+(p4))+gamma(p5)'
c--    '  f(p1)+f(p2) --> W^- (no BR) + gamma(p5)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
      
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 13) .or. (nproc .eq. 18)) then
        case='W_cjet'
        nqcdjets=1
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        mass2=mc
        nflav=3
        plabel(5)='bq'
        plabel(6)='pp'

        if     (nproc .eq. 13) then
c-- 13 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+cbar(p5)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) + cbar(p5)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 18) then
c-- 18 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+c(p5)'
c--    '  f(p1)+f(p2) --> W^- (no BR) + c(p5)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 14) .or. (nproc .eq. 19)) then
        case='Wcjet0'
        nqcdjets=1
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        mass2=0d0
        nflav=3
        plabel(5)='bq'
        plabel(6)='pp'

        if     (nproc .eq. 14) then
c-- 13 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+cbar(p5) [massless]'
c--    '  f(p1)+f(p2) --> W^+ (no BR) + cbar(p5) [massless]' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 19) then
c-- 18 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+c(p5) [massless]'
c--    '  f(p1)+f(p2) --> W^- (no BR) + c(p5) [massless]' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 20) .or. (nproc .eq. 25)) then
        case='Wbbmas'
        write(6,*) 'mb=',mb
        nqcdjets=2
        flav=5
        bbproc=.true.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
 
        if     (nproc .eq. 20) then
c-- 20 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+b~(p6) [massive]'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +b(p5)+b~(p6) [massive]' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 25) then
c-- 25 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5)+b~(p6) [massive]'
c--    '  f(p1)+f(p2) --> W^- (no BR) +b(p5)+b~(p6) [massive]' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
 
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 21) .or. (nproc .eq. 26)) then
        case='Wbbbar'
        write(6,*) 'mb=',mb
        nqcdjets=2
        bbproc=.true.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        mb=0
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 21) then
c-- 21 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+b~(p6)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +b(p5)+b~(p6)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 26) then
c-- 26 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5)+b~(p6)'
c--    '  f(p1)+f(p2) --> W^- (no BR) +b(p5)+b~(p6)' (removebr=.true.)
          nwz=-1 
          plabel(3)='el'
          plabel(4)='na'
        endif
             
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 22) .or. (nproc .eq. 27)) then
        case='W_2jet'
        nqcdjets=2
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 22) then
c-- 22 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+f(p5)+f(p6)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +f(p5)+f(p6)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 27) then
c-- 27 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + f(p5)+f(p6)'
c--    '  f(p1)+f(p2) --> W^- (no BR) +f(p5)+f(p6)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 23) .or. (nproc .eq. 28)) then
        case='W_3jet'
        nqcdjets=3
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=13
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 23) then
c-- 23 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+f(p5)+f(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +f(p5)+f(p6)+f(p7)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 28) then
c-- 28 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + f(p5)+f(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> W^- (no BR) +f(p5)+f(p6)+f(p7)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 24) .or. (nproc .eq. 29)) then
        case='Wbbjet'
        write(6,*) 'mb=',mb
        nqcdjets=3
        bbproc=.true.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        mb=0d0
        ndim=13
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 24) then
c-- 24 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+b~(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +b(p5)+b~(p6)+f(p7)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 29) then
c-- 29 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5)+b~(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> W^- (no BR) +b(p5)+b~(p6)+f(p7)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
        endif
        
c-----------------------------------------------------------------------

      elseif (nproc/10 .eq. 3) then
        case='Z_only'
        call checkminzmass(1)
        nqcdjets=0
        nwz=0
        mass3=zmass
        width3=zwidth
        n3=1
        ndim=4
        plabel(5)='pp'

        if     (nproc .eq. 31) then
c-- 31 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))'
c--    '  f(p1)+f(p2) --> Z^0 (for total Xsect)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='ea'
          q1=-1d0
          l1=le
          r1=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            call branch(brwen,brzee,brtau,brtop)
            BrnRat=brzee
          endif
        elseif (nproc .eq. 32) then
c-- 32 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4)))'
          plabel(3)='nl'
          plabel(4)='na'
          q1=0d0
          l1=ln*dsqrt(3d0)
          r1=rn*dsqrt(3d0)
        elseif (nproc .eq. 33) then
c-- 33 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4))'
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
        else
          call nprocinvalid()
        endif 

c-----------------------------------------------------------------------
      
      elseif ((nproc .ge. 41) .and. (nproc .le. 43)) then
        case='Z_1jet'
        call checkminzmass(1)
        nqcdjets=1
        nwz=0
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        plabel(5)='pp'
        plabel(6)='pp'

        if     (nproc .eq. 41) then
c-- 41 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+f(p5)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +f(p5)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='ea'
          q1=-1d0
          l1=le
          r1=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            call branch(brwen,brzee,brtau,brtop)
            BrnRat=brzee
          endif
        elseif (nproc .eq. 42) then
c-- 42 '  f(p1)+f(p2) --> Z_0(-->3*(nu(p3)+nu~(p4)))-(sum over 3 nu)+f(p5)'
            plabel(3)='nl'
            plabel(4)='na'
            q1=0d0
            l1=ln*dsqrt(3d0)
            r1=rn*dsqrt(3d0)
        elseif (nproc .eq. 43) then
c-- 43 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4))+f(p5)'
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
        endif
      
c-----------------------------------------------------------------------
      
      elseif (nproc .eq. 44) then
c-- 44 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+f(p5)+f(p6)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +f(p5)+f(p6)' (removebr=.true.)
        case='Z_2jet'
        call checkminzmass(1)
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
        nwz=0   
        mass3=zmass
        width3=zwidth
       
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brzee
        endif

      elseif (nproc .eq. 45) then
c-- 45 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+f(p5)+f(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +f(p5)+f(p6)+f(p7)' (removebr=.true.)
        case='Z_3jet'
        call checkminzmass(1)
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
        nwz=0   
        mass3=zmass
        width3=zwidth
       
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brzee
        endif

c-----------------------------------------------------------------------

        elseif ((nproc .eq. 48) .or. (nproc .eq. 49)) then
          case='Zgamma'
          call checkminzmass(1)
          nqcdjets=0
          ndim=7
          n2=0
          n3=1
          mass3=zmass
          width3=zwidth
          nwz=0
          plabel(5)='ga'
          plabel(6)='pp'
          
          if     (nproc .eq. 48) then
c-- 48 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+gamma(p5)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +gamma(p5)' (removebr=.true.)
            plabel(3)='el'
            plabel(4)='ea'
            q1=-1d0
            l1=le
            r1=re
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              call branch(brwen,brzee,brtau,brtop)
              BrnRat=brzee
            endif
          elseif (nproc .eq. 49) then
c-- 49 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4)))-(sum over 3 nu)+gamma(p5)'
            plabel(3)='nl'
            plabel(4)='na'
            q1=0d0
            l1=ln*dsqrt(3d0)
            r1=rn*dsqrt(3d0)
          endif

c-----------------------------------------------------------------------
          
      elseif (nproc .eq. 50) then
c-- 50 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b~(p5)+b(p6) (massive)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +b~(p5)+b(p6) (massive)' (removebr=.true.)
        case='Zbbmas'
        call checkminzmass(1)
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
        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brzee
        endif

      elseif ((nproc .ge. 51) .and. (nproc .le. 53)) then
        case='Zbbbar'
        call checkminzmass(1)
        bbproc=.true.
        mb=0d0
        nqcdjets=2
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        q1=-1d0
        l1=le
        r1=re
        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

        if     (nproc .eq. 51) then
c-- 51 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+b~(p6)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +b(p5)+b~(p6)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='ea'
          q1=-1d0
          l1=le
          r1=re
c--- total cross-section             
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            call branch(brwen,brzee,brtau,brtop)
            BrnRat=brzee
          endif
        elseif (nproc .eq. 52) then
c-- 52 '  f(p1)+f(p2) --> Z_0(-->3*(nu(p3)+nu~(p4)))+b(p5)+b~(p6)'
          plabel(3)='nl'
          plabel(4)='na'
          q1=0d0
          l1=ln*dsqrt(3d0)
          r1=rn*dsqrt(3d0)
        elseif (nproc .eq. 53) then
c-- 53 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4))+b(p5)+b~(p6)'
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
        endif
        
      elseif (nproc .eq. 56) then
c-- 56 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+b~(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +b(p5)+b~(p6)+f(p7)' (removebr=.true.)
        case='Zbbjet'
        ndim=13
        bbproc=.true.
        mb=0d0
        nqcdjets=3
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        q1=-1d0
        l1=le
        r1=re
        mass3=zmass
        width3=zwidth

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brzee
        endif

c-----------------------------------------------------------------------
          
      elseif (nproc/10 .eq. 6) then
        case='WWqqbr'
        call readcoup
        nqcdjets=0  
        plabel(7)='pp'
        nwz=1
        ndim=10
        mb=0d0
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

        if    ((nproc .eq. 61) .or. (nproc .eq. 64)) then
c--  61 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6))'
c--     '  f(p1)+f(p2) --> W^+ + W^- (for total Xsect)' (removebr=.true.)
          if (nproc .eq. 64) then
c--  64 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6)) [no pol]'
            case='WWnpol'
          endif
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='el'
          plabel(6)='na'
          l1=1d0
c--- total cross-section             
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brtau,brtop)
            BrnRat=brwen**2
          endif
        elseif (nproc .eq. 62) then
c--  62 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->q(p5)+q~(p6))'
c--- note: scattering diagrams are NOT included, only couplings change
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='qq'
          plabel(6)='aa'
          plabel(7)='pp'
          l1=dsqrt(xn*2d0)
        elseif (nproc .eq. 63) then
c--  63 '  f(p1)+f(p2) --> W^+(--> q(p3)+ q~(p4)) +W^-(-->e^-(p5)+nu~(p6))'
c--- note: scattering diagrams are NOT included, only couplings change
          plabel(3)='qq'
          plabel(4)='aa'
          plabel(5)='el'
          plabel(6)='na'
          plabel(7)='pp'
          l1=dsqrt(xn*2d0)
        endif

c-----------------------------------------------------------------------

      elseif (nproc/10 .eq. 7) then
        case='WZbbar'
        call checkminzmass(2)
        call readcoup
        nqcdjets=0
        plabel(7)='pp'
        ndim=10
        mb=0
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=wmass
        width3=wwidth

        if (nproc .lt. 75) then
c-- W^+Z
          nwz=+1

          if     (nproc .eq. 71) then             
c--  71 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->e^-(p5)+e^+(p6))'
c--     '  f(p1)+f(p2) --> W^+ (for total Xsect) + Z^0 ' (removebr=.true.)
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='ml'
            plabel(6)='ma'
            q1=-1d0
            l1=le
            r1=re
c--- total cross-section             
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              plabel(5)='ig'
              plabel(6)='ig'
              call branch(brwen,brzee,brtau,brtop)
              BrnRat=brwen*brzee
            endif
          elseif (nproc .eq. 72) then
c--  72 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->nu_e(p5)+nu~_e(p6))'
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='nl'
            plabel(6)='na'
            plabel(7)='pp'
            q1=0d0
            l1=ln
            r1=rn
          elseif (nproc .eq. 73) then
c--  73 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->b(p5)+b~(p6))'
            bbproc=.true.
            nqcdjets=2
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='bq'
            plabel(6)='ba'
            plabel(7)='pp'
            q1=Q(5)*dsqrt(xn)
            l1=l(5)*dsqrt(xn)
            r1=r(5)*dsqrt(xn)
          else
            call nprocinvalid()
          endif 

        elseif (nproc .ge. 75) then
c-- W^-Z
          nwz=-1

          if     (nproc .eq. 76) then
c--  76 '  f(p1)+f(p2) --> W^-(-->mu^-(p3)+nu~(p4))+Z^0(-->e^-(p5)+e^+(p6))'
c--     '  f(p1)+f(p2) --> W^- + Z^0 (for total Xsect)' (removebr=.true.)
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='ml'
            plabel(6)='ma'
            q1=-1d0
            l1=le
            r1=re
c--- total cross-section             
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              plabel(5)='ig'
              plabel(6)='ig'
              call branch(brwen,brzee,brtau,brtop)
              BrnRat=brwen*brzee
            endif
          elseif (nproc .eq. 77) then
c--  77 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->nu(p5)+nu~(p6))'
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='nl'
            plabel(6)='na'
            q1=0d0
            l1=ln
            r1=rn
          elseif (nproc .eq. 78) then
c--  78 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->b(p5)+b~(p6))'
            bbproc=.true.
            nqcdjets=2
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='bq'
            plabel(6)='ba'
            q1=Q(5)*dsqrt(xn)
            l1=l(5)*dsqrt(xn)
            r1=r(5)*dsqrt(xn)
          else
            call nprocinvalid()
          endif 

        endif
            
c-----------------------------------------------------------------------

      elseif (nproc/10 .eq. 8) then
        case='ZZlept'
        call checkminzmass(1)
        call checkminzmass(2)
        call readcoup
        plabel(7)='pp'
        nqcdjets=0
        nwz=0
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
        
        if (nproc .eq. 81 .or. nproc .eq. 86) then
c--  81 '  f(p1)+f(p2) --> Z^0(-->mu^-(p3)+mu^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))'
c--     '  f(p1)+f(p2) --> Z^0 + Z^0 (for total Xsect)' (removebr=.true.)
c--  86 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->mu^-(p3)+mu^+(p4)) (NO GAMMA*)'
c--     '  f(p1)+f(p2) --> Z^0 + Z^0 (for total Xsect) (NO GAMMA*)' (removebr=.true.)
          plabel(3)='el' 
          plabel(4)='ea'
          plabel(5)='ml'
          plabel(6)='ma'
          q2=-1d0
          l2=le
          r2=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brtau,brtop)
            BrnRat=2d0*brzee**2
          endif
        elseif (nproc .eq. 82 .or. nproc .eq. 87) then
c--  82 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + Z^0(-->e^-(p5)+e^+(p6))'
c--  87 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->3*(nu(p3)+nu~(p4))) (NO GAMMA*)'
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='nl'
          plabel(6)='na'
          q2=0d0
          l2=ln*dsqrt(3d0)
          r2=rn*dsqrt(3d0)
        elseif (nproc .eq. 83 .or. nproc .eq. 88) then
c--  83 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6)) + Z^0(-->b(p3)+b~(p4))'
c--  88 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->b(p3)+b~(p4)) (NO GAMMA*)'
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='qb'
          plabel(6)='ab'
          q2=Q(5)*dsqrt(xn)
          l2=l(5)*dsqrt(xn)
          r2=r(5)*dsqrt(xn)
        elseif ((nproc .eq. 84) .or. (nproc .eq. 89))  then
c--  84 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + Z^0(-->b(p5)+bb(p6))'
c--  89 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + Z^0(-->b(p5)+bb(p6))(NO GAMMA*)'
          bbproc=.true.
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='na'
          plabel(5)='bq'
          plabel(6)='ba'
          q2=0d0
          l2=ln*dsqrt(3d0)
          r2=rn*dsqrt(3d0)
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
        else
          call nprocinvalid()
        endif 

c-- remove gamma^* if necessary
        if (nproc .ge. 85) then
          q1=0d0
          q2=0d0
        endif
        
c-----------------------------------------------------------------------

      elseif (nproc/10 .eq. 9) then
        case='WHbbar'
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr)
        nqcdjets=2
        bbproc=.true.
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        
        ndim=10
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 91) then
c--  91 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) H(-->b(p5)+b~(p6)) '
c--     '  f(p1)+f(p2) --> W+ + H (for total Xsect)' (removebr=.true.)
          nwz=1
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            call branch(brwen,brzee,brtau,brtop)
            BrnRat=brwen*br
            bbproc=.false.
            nqcdjets=0
          endif
        elseif (nproc .eq. 96) then
c--  96 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+ H(-->b(p5)+b~(p6))' 
c--     '  f(p1)+f(p2) --> W- + H (for total Xsect)' (removebr=.true.)
          nwz=-1
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            call branch(brwen,brzee,brtau,brtop)
            BrnRat=brwen*br
            bbproc=.false.
            nqcdjets=0
          endif
        else
            call nprocinvalid()
        endif

c-----------------------------------------------------------------------

      elseif (nproc/10 .eq. 10) then
        case='ZHbbar'
        call checkminzmass(1)
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr)
        nqcdjets=2
        bbproc=.true.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'

        ndim=10
        nwz=0
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=zmass
        width3=zwidth

        if (nproc .eq. 101) then
c--  101 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + H(-->b(p5)+b~(p6))'
c--      '  f(p1)+f(p2) --> H + Z0 (for total Xsect)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='ea'
          q1=-1d0
          l1=le
          r1=re
          if (removebr) then
            call branch(brwen,brzee,brtau,brtop)
            BrnRat=brzee*br
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'    
            bbproc=.false.
            nqcdjets=0           
          endif
        elseif (nproc .eq. 102) then
c--  102 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + H(-->b(p5)+b~(p6))'
          plabel(3)='nl'
          plabel(4)='na'
          q1=0d0
          l1=ln*dsqrt(3d0)
          r1=rn*dsqrt(3d0)
        elseif (nproc .eq. 103) then
c--  103 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + H(-->b(p5)+b~(p6))'     
          nqcdjets=4
          plabel(3)='bq'
          plabel(4)='ba'
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
        else
          call nprocinvalid()
        endif 

c-----------------------------------------------------------------------

       elseif ((nproc .eq. 111) .or. (nproc .eq. 112)) then
        case='ggfus0'
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr)
        plabel(5)='pp'
        ndim=4
      
        n2=0
        n3=1
        mass3=hmass
        width3=hwidth

        if     (nproc .eq. 111) then
c--  111 '  f(p1)+f(p2) --> H(-->b(p3)+bbar(p4))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)       
          plabel(3)='bq'
          plabel(4)='ba'
          nqcdjets=2
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=0
            BrnRat=br
          endif
          
        elseif (nproc .eq. 112) then
c--  112 '  f(p1)+f(p2) --> H(-->tau^-(p3)+tau^+(p4))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)       
          plabel(3)='tl'
          plabel(4)='ta'
          nqcdjets=0
          Brnrat=br/tautaubr
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            BrnRat=br
          endif
        endif

        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          nqcdjets=0
          BrnRat=br
        endif

      elseif (nproc .eq. 113) then
c--  113 '  f(p1)+f(p2) --> H (--> W^+(nu(p3)+e^+(p4)) + W^-(e^-(p5)+nu~(p6)))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)
        case='HWW_4l'             
        call sethparams(br,wwbr,zzbr,tautaubr)
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='pp'
        nqcdjets=0
        ndim=10
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

c--- print warning if we're below threshold
        if (hmass .lt. 2d0*wmass) then
          write(6,*)
          write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
          write(6,*) 'may not yield sensible results - check the number'
          write(6,*) 'of integration points and the value of zerowidth'
        endif
        
        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen**2*wwbr
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'               
        endif
        
      elseif ((nproc .ge. 114) .and. (nproc .le. 116)) then
        case='HZZ_4l'
        call sethparams(br,wwbr,zzbr,tautaubr)
        plabel(7)='pp'
        nqcdjets=0
        ndim=10
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth

c--- print warning if we're below threshold
        if (hmass .lt. 2d0*zmass) then
          write(6,*)
          write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
          write(6,*) 'may not yield sensible results - check the number'
          write(6,*) 'of integration points and the value of zerowidth'
        endif
        
        if     (nproc .eq. 114) then
c--  114 '  f(p1)+f(p2) --> H(-->Z^0(mu^-(p3)+mu^+(p4)) + Z^0(e^-(p5)+e^+(p6))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='ml'
          plabel(6)='ma'
          l1=le
          r1=re
          l2=le
          r2=re
          if (removebr) then
            call branch(brwen,brzee,brtau,brtop)
            BrnRat=brzee**2*zzbr
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'             
          endif
        elseif (nproc .eq. 115) then
c--  115 '  f(p1)+f(p2) --> H(-->Z^0(3*(nu(p3)+nu~(p4)))+ Z^0(e^-(p5)+e^+(p6))'
          plabel(3)='nl'
          plabel(4)='na'
          plabel(5)='ml'
          plabel(6)='ma'
          l1=le
          r1=re
          l2=ln*dsqrt(3d0)
          r2=rn*dsqrt(3d0)      
        elseif (nproc .eq. 116) then
c--  116 '  f(p1)+f(p2) --> H(-->Z^0(mu^-(p3)+mu^+(p4)) + Z^0(b(p5)+b~(p6))'
          nqcdjets=2
          plabel(3)='ml'
          plabel(4)='ma'
          plabel(5)='bq'
          plabel(6)='ba'
          l1=le
          r1=re
          l2=l(5)*dsqrt(xn)
          r2=r(5)*dsqrt(xn)
        else
          call nprocinvalid()
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .ge. 141) .and. (nproc .le. 143)) then
        case='H_1jet'
        call sethparams(br,wwbr,zzbr,tautaubr)
        call setmb_msbar
        ndim=7
        plabel(3)='bq'
        plabel(4)='ba'
        plabel(5)='bq'
        nqcdjets=3

        n2=0
        n3=1
        mass3=hmass
        width3=hwidth
        
        if ( (nproc .eq. 142) .and.  
     .       ((part .eq. 'virt') .or. (part .eq. 'tota')) ) then
          write(6,*) 'This process number is not suitable for the'
          write(6,*) 'NLO calculation. Please run processes'
          write(6,*) '141 (virtual+real) and 142 (real) separately.'
          stop
        endif
        if ( (nproc .eq. 143) .and. (part .ne. 'real') ) then
          write(6,*) 'This process number is not suitable for such a'
          write(6,*) 'calculation. Please run process 143 (real) only.'
          stop
        endif
             
        if     (nproc .eq. 141) then
c--  141 '  f(p1)+f(p2) --> H (no BR) + b(p5) [+g(p6)]'
          isub=1
          plabel(6)='pp'
        elseif (nproc .eq. 142) then
c--  142 '(p1)+f(p2) --> H (no BR) + b~(p5) [+b(p6)]'
          isub=2
          plabel(5)='ba'
          plabel(6)='bq'
        elseif (nproc .eq. 143) then
c--  143 '  f(p1)+f(p2) --> H (no BR) + b(p5) + b~(p6) [both observed]'
          isub=2
          plabel(5)='ba'
          plabel(6)='bq'
          nqcdjets=4
        endif
        
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          BrnRat=br
          nqcdjets=nqcdjets-2
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 151) .or. (nproc .eq. 152)) then
        nwz=1
        ndim=16
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth
        bbproc=.true.
        
        if (nproc .eq. 151) then
c--  151 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+e^-(p7)+nu~(p8)'
c--      '  f(p1)+f(p2) --> t t~ (with BR for total Xsect)' (removebr=.true.)
          case='tt_bbl'
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          if (removebr) then
            call branch(brwen,brzee,brtau,brtop)
            BrnRat=(brwen*brtop)**2
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            plabel(7)='ig'
            plabel(8)='ig'
            nqcdjets=0
            bbproc=.false.
          endif
        elseif (nproc .eq. 152) then
c--  152 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+q(p7)+q~(p8)'
          case='tt_bbh'
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='pp'
          plabel(8)='pp'
          nqcdjets=4
          if (runstring(1:4) .eq. 'stop') then
c--- For single top study, we only want to see 2 of the jets
            notag=2   
          endif             
        endif 

c-----------------------------------------------------------------------

      elseif (nproc .eq. 156) then
c--  156 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->nu~(p7)+e^-(p8)+b~(p6))+g(p9)'
c--      '  f(p1)+f(p2)-->t(p345)+t~(p678)+g(p9)' (removebr=.true.)
        case='qq_ttg'
        zerowidth=.true.
        Write(6,*) 'Setting zerowidth to true, obligatory for nproc=156' 
        nwz=1
        ndim=19
        mb=0
        nqcdjets=3
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='el'
        plabel(8)='na'
        plabel(9)='pp'

c--- total cross-section             
        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=(brwen*brtop)**2
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'               
          plabel(7)='ig'
          plabel(8)='ig'
          nqcdjets=1
        endif

c-----------------------------------------------------------------------

      elseif (nproc .eq. 157) then
c--  157 '  f(p1)+f(p2) --> t t~ (for total Xsect)'
        case='tt_tot'
        nqcdjets=0
        ndim=4
        mass2=mt
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
      elseif (nproc .eq. 158) then
c--  158 '  f(p1)+f(p2) --> b b~ (for total Xsect)'
        case='bb_tot'
        nqcdjets=0
        ndim=4
        mass2=mb
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
      elseif (nproc .eq. 159) then
c--  159 '  f(p1)+f(p2) --> c c~ (for total Xsect)'
        case='cc_tot'
        nqcdjets=0
        ndim=4
        mass2=mc
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
 
c-----------------------------------------------------------------------

      elseif (nproc .eq. 161) then
c--  161 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+q(p6) [t-channel]'
c--      '  f(p1)+f(p2) --> t(no BR) + q(p6)' (removebr=.true.)
        case='bq_tpq'
        isub=1
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='qj'
        plabel(7)='pp'
        nwz=+1

c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=0
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
        endif

      elseif (nproc .eq. 162) then
c--  162 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+q(p6)(decay)'
        case='ttdkay'
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='qj'
        plabel(7)='pp'
        nwz=+1

        if (part .eq. 'lord') then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '161 (lord) or process 162 (virt+real).'
          stop
        endif
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=0
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
        
      elseif (nproc .eq. 166) then
c--  166 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+q(p6) [t-channel]''
c--      '  f(p1)+f(p2) --> t~(no BR) + q(p6)' (removebr=.true.)
        case='bq_tpq'
        isub=1
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ba'
        plabel(6)='qj'
        plabel(7)='pp'
        nwz=-1

c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=0
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
        endif
        
      elseif (nproc .eq. 167) then
c--  167 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+q(p6)(decay)'
        case='ttdkay'
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ba'
        plabel(6)='qj'
        plabel(7)='pp'
        nwz=-1

        if (part .eq. 'lord') then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '166 (lord) or process 167 (virt+real).'
          stop
        endif
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=0
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
             
c-----------------------------------------------------------------------

      elseif (nproc .eq. 171) then
c--  171 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6)) [s-channel]'
c--      '  f(p1)+f(p2) --> t(no BR) + b~(p6)' (removebr=.true.)
        case='t_bbar'
        isub=2
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        nwz=1
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=0
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
        
        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
        endif
        
      elseif (nproc .eq. 172) then
c--  172 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))(decay)'
        case='tdecay'
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        nwz=1
        
        if (part .eq. 'lord') then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '171 (lord) or process 172 (virt+real).'
          stop
        endif
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=0
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

      elseif (nproc .eq. 176) then
c--  176 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+b(p6)) [s-channel]'
c--      '  f(p1)+f(p2) --> t~(no BR) + b(p6)' (removebr=.true.)
        case='t_bbar'
        isub=2
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ba'
        plabel(6)='bq'
        plabel(7)='pp'
        nwz=-1
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=0
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
        
         if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          nqcdjets=0           
        endif
             
      elseif (nproc .eq. 177) then
c--  177 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+b(p6))(decay)'
        case='tdecay'
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ba'
        plabel(6)='bq'
        plabel(7)='pp'
        nwz=-1
        
        if (part .eq. 'lord') then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '176 (lord) or process 177 (virt+real).'
          stop
        endif
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=0
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

c-----------------------------------------------------------------------

      elseif (nproc .eq. 180) then
c--  180 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(p5)'
        case='W_tndk'
        nqcdjets=1
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='bq'
        mass2=mt
        nflav=5
        nwz=-1

        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        
        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
             
      elseif (nproc .eq. 181) then
c--  181 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(nu(p5)+e^+(p6)+b(p7))'
        case='W_twdk'
        nqcdjets=1
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ne'
        plabel(6)='ea'
        plabel(7)='bq'
        nflav=5
        nwz=-1

        ndim=13
        mb=0
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=mt
        width3=twidth
             
        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen*brtop
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          nqcdjets=0           
        endif
             
      elseif (nproc .eq. 185) then
c--  185 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+tbar(p5)'
        case='W_tndk'
        nqcdjets=1
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='ba'
        mass2=mt
        nflav=5
        nwz=+1
        
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
             
        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
      elseif (nproc .eq. 186) then
c--  186 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+t~(e^-(p5)+nu~(p6)+bbar(p7)'
        case='W_twdk'
        nqcdjets=1
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='ba'
        nflav=5
        nwz=+1
        
        ndim=13
        mb=0
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=mt
        width3=twidth

        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen*brtop
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          nqcdjets=0           
        endif
             
c-----------------------------------------------------------------------

      elseif (nproc .eq. 190) then
c--  190 '  f(p1)+f(p2)-->t(p3)+t~(p4)+H(p5)'
        case='tottth'
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        nwz=1
        n2=0
        n3=0
        ndim=7
        
      elseif (nproc .eq. 191) then
c--  191 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))
c--         +t~(-->nu~(p7)+e^-(p8)+b~(p6))+H(p9+p10)'
c--      '  f(p1)+f(p2)-->t(p3+p4+p5)+t~(p6+p7+p8)+H(p9+p10)' (removebr=.true.)
        case='qq_tth'
        call sethparams(br,wwbr,zzbr,tautaubr)
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='el'
        plabel(8)='na'
        plabel(9)='bq'
        plabel(10)='ba'

        nwz=1
        ndim=22
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth

        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=(brtop*brwen)**2*br
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          plabel(8)='ig'
          plabel(9)='ig'
          plabel(10)='ig'
        endif

c-----------------------------------------------------------------------

        elseif ((nproc .eq. 196) .or. (nproc .eq. 197)) then
        case='qq_ttz'
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='el'
        plabel(8)='na'
        nwz=1

        ndim=22
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth

        if     (nproc .eq. 196) then 
c--  196 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->nu~(p7)+e^-(p8)+b~(p6))+Z(e(p9),e~(p10))'
          plabel(9)='el'
          plabel(10)='ea'
          q1=-1d0
          l1=le
          r1=re
          if (removebr) then
            call branch(brwen,brzee,brtau,brtop)
            BrnRat=(brtop*brwen)**2*brzee
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            plabel(7)='ig'
            plabel(8)='ig'
            plabel(9)='ig'
            plabel(10)='ig'
          endif
        elseif (nproc .eq. 197) then 
c--  197 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->nu~(p7)+e^-(p8)+b~(p6))+Z(b(p9),b~(p10))'
          plabel(9)='bq'
          plabel(10)='ba'
          q1=Q(5)*dsqrt(xn)
          l1=l(5)*dsqrt(xn)
          r1=r(5)*dsqrt(xn)
        endif

c-----------------------------------------------------------------------

      elseif (nproc/10 .eq. 20) then
        case='httjet'
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr)
        nqcdjets=1
        plabel(3)='tl'
        plabel(4)='ta'
        plabel(5)='pp'
 
        ndim=7
        n2=0
        n3=1
        mass3=hmass
        width3=hwidth
        
        mass3=hmass
        width3=hwidth
        
        if     (nproc .eq. 201) then
c--  201 '  f(p1)+f(p2)--> H(-->b(p3)+b~(p4)) + f(p5) [full mt dep.]'
c--      '  f(p1)+f(p2)--> H(p3+p4) + f(p5) (for total Xsect)' (removebr=.true.)
          plabel(3)='bq'
          plabel(4)='ba'  
          nqcdjets=3
          if (removebr) then        
            BrnRat=br
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
          endif

        elseif (nproc .eq. 202) then
c--  202 '  f(p1)+f(p2)--> H (-> tau(p3) tau~(p4)) + f(p5) [full mt dep.]'
          plabel(3)='tl'
          plabel(4)='ta'
          Brnrat=br/tautaubr

        elseif ((nproc .eq. 203) .or. (nproc .eq. 204)) then
          case='ggfus1'
          nqcdjets=1
          mb=0
          plabel(5)='pp'
          plabel(6)='pp'
          ndim=7
      
          n2=0
          n3=1

          if     (nproc .eq. 203) then
c--  203 '  f(p1)+f(p2) -->H(-->b(p3)+b~(p4)) + f(p5)'
c--      '  f(p1)+f(p2)--> H(p3+p4) + f(p5) (for total Xsect)' (removebr=.true.)
            plabel(3)='bq'
            plabel(4)='ba'  
            nqcdjets=3
            if (removebr) then        
              BrnRat=br
              plabel(3)='ig'
              plabel(4)='ig'
              nqcdjets=1
            endif
          elseif (nproc .eq. 204) then
c--  204 '  f(p1)+f(p2) -->H(-->tau^-(p3)+tau^+(p4)) + f(p5)'
            plabel(3)='tl'
            plabel(4)='ta'
            Brnrat=br/tautaubr
          endif
        
        elseif (nproc .eq. 206) then
c--  206 '  f(p1)+f(p2)--> A(-->b(p3)+b~(p4)) + f(p5) [full mt dep.]'
c--      '  f(p1)+f(p2)--> A(p3+p4) + f(p5) (for total Xsect)' (removebr=.true.)
          case='attjet'
          plabel(3)='bq'
          plabel(4)='ba'  
          nqcdjets=3
          if (removebr) then        
            BrnRat=br
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
          endif

        elseif (nproc .eq. 207) then
c--  207 '  f(p1)+f(p2)--> A (--> tau(p3) tau~(p4)) + f(p5) [full mt dep.]'
         case='attjet'
          plabel(3)='tl'
          plabel(4)='ta'
          Brnrat=br/tautaubr
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 211) .or. (nproc .eq. 212)) then
        case='qq_Hqq'
        mb=0d0
        call sethparams(br,wwbr,zzbr,tautaubr)
        nwz=2
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=10
        n2=0
        n3=1

        mass3=hmass
        width3=hwidth

        if     (nproc .eq. 211) then
c--  211 '  f(p1)+f(p2)--> H(-->b(p3)+b~(p4))+f(p5)+f(p6) [WBF]'
c--      '  f(p1)+f(p2)--> H(p3+p4)+f(p5)+f(p6) [WBF]' (removebr=.true.)
          plabel(3)='bq'
          plabel(4)='ba'
          nqcdjets=4
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=2
            BrnRat=br
          endif
        elseif (nproc .eq. 212) then
c--  212 '  f(p1)+f(p2)--> H(-->tau-(p3)+tau+(p4))+f(p5)+f(p6) [WBF]'
c--      '  f(p1)+f(p2)--> H(p3+p4)+f(p5)+f(p6) [WBF]' (removebr=.true.)
          plabel(3)='tl'
          plabel(4)='ta'
          nqcdjets=2
          Brnrat=br/tautaubr
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=0
            BrnRat=br
          endif
        endif
          
      elseif ((nproc .eq. 216) .or. (nproc .eq. 217)) then
        case='qqHqqg'
        mb=0d0
        call sethparams(br,wwbr,zzbr,tautaubr)
        nwz=2
        plabel(3)='bq'
        plabel(4)='ba'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=13
        n2=0
        n3=1

        mass3=hmass
        width3=hwidth

        if     (nproc .eq. 216) then
c-- 216 '  f(p1)+f(p2)--> H(-->b(p3)+b~(p4))+f(p5)+f(p6)+f(p7) [WBF+jet]'
c--     '  f(p1)+f(p2)--> H(p3+p4)+f(p5)+f(p6)+f(p7) [WBF+jet]' (removebr=.true.)
          plabel(3)='bq'
          plabel(4)='ba'
          nqcdjets=5
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=3
            BrnRat=br
          endif
        elseif (nproc .eq. 217) then
c-- 217 '  f(p1)+f(p2)--> H(-->tau-(p3)+tau+(p4))+f(p5)+f(p6)+f(p7) [WBF+jet]'
          plabel(3)='tl'
          plabel(4)='ta'
          nqcdjets=3
          Brnrat=br/tautaubr
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
            BrnRat=br
          endif
        endif

c-----------------------------------------------------------------------

      elseif (nproc .eq. 221) then
        case='tautau'
c--  221 '  f(p1)+f(p2)--> tau^-(-->e^-(p3)+nu~_e(p4)+nu_tau(p5))+tau^+(-->nu~_tau(p6)+nu_e(p7)+e^+(p8))'
c--      '  f(p1)+f(p2)--> tau tau~ [for total Xsect]' (removebr=.true.)
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='nl'
        plabel(6)='na'
        plabel(7)='nl'
        plabel(8)='ea'
        nqcdjets=0
        nwz=1
        ndim=16
        n2=1
        n3=1
        mass2=mtau
        width2=tauwidth
        mass3=mtau
        width3=tauwidth

        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brtau**2
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          plabel(8)='ig'
        endif

c-----------------------------------------------------------------------

      elseif ((nproc .eq. 261) .or. (nproc .eq. 266)) then
c--  261 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)'
c--  266 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)[+b~(p6)]'
        case='gQ__ZQ'
        nqcdjets=1
        flav=5
        nwz=0
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        isub=1+(nproc-261)/5
        if (nproc .eq. 261) then
          plabel(6)='pp'
        else
          plabel(6)='ba'
        endif
        q1=-1d0
        l1=le
        r1=re

        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
      elseif ((nproc .eq. 262) .or. (nproc .eq. 267)) then
c--  262 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c(p5)'
c--  267 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c(p5)[+c~(p6)]'
        case='gQ__ZQ'
        nqcdjets=1
        flav=4
        nwz=0
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        isub=1+(nproc-262)/5
        if (nproc .eq. 262) then
          plabel(6)='pp'
        else
          plabel(6)='ba'
        endif
        q1=-1d0
        l1=le
        r1=re
        
        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
      elseif (nproc .eq. 263) then
c--  263 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b~(p5)+b(p6) (1 b-tag)'
        case='Zbbmas'
        nqcdjets=2
        notag=1

        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

        write(6,*) 'mb=',mb
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        q1=-1d0
        l1=le
        r1=re
        
        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
      elseif (nproc .eq. 264) then
c--  264 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c~(p5)+c(p6) (1 c-tag)'
        case='Zccmas'
        nqcdjets=2
        notag=1
        
        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        
        mb=mc
        write(6,*) 'mc=',mb
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        q1=-1d0
        l1=le
        r1=re

        if (removebr) then
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif      
           
c-----------------------------------------------------------------------
          
      elseif ((nproc .eq. 271) .or. (nproc .eq. 272)) then
        case='ggfus2'
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr)
        plabel(3)='tl'
        plabel(4)='ta'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=10
      
        n2=0
        n3=1

        mass3=hmass
        width3=hwidth
        
        if     (nproc .eq. 271) then
c-- 271 '  f(p1)+f(p2) --> H(b(p3)+b~(p4))+f(p5)+f(p6)[in heavy top limit]'
c--     '  f(p1)+f(p2) --> H(no BR)+f(p5)+f(p6)[in heavy top limit]' (removebr=.true.)
          plabel(3)='bq'
          plabel(4)='ba'
          nqcdjets=4
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=2
            BrnRat=br
          endif
          
        elseif (nproc .eq. 272) then
c-- 272 '  f(p1)+f(p2) --> H(tau-(p3)+tau+(p4))+f(p5)+f(p6)[in heavy top limit]'
c--     '  f(p1)+f(p2) --> H(no BR)+f(p5)+f(p6)[in heavy top limit]' (removebr=.true.)
          plabel(3)='tl'
          plabel(4)='ta'
          nqcdjets=2
          Brnrat=br/tautaubr
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            BrnRat=br
          endif
        endif
                
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 311) .or. (nproc .eq. 316)) then
        case='W_bjet'
        nqcdjets=2
        flav=5
        isub=1
        
        nflav=5
        mb=0d0
        plabel(5)='bq'
        plabel(6)='pp'
        plabel(7)='pp'
        
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 311) then
c--  311 '  f(p1)+b(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+f(p6)'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 316) then
c--  316 '  f(p1)+b(p2) --> W^-(-->e^-(p3)+nu~(p4))+b(p5)+f(p6)'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
        
        if (removebr) then
c--      '  f(p1)+b(p2) --> W(no BR)+b(p5)+f(p6)' (removebr=.true.)
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 321) .or. (nproc .eq. 326)) then
        case='W_bjet'
        nqcdjets=2
        flav=4
        isub=1
        
        nflav=4
        mb=0d0
        plabel(5)='bq'
        plabel(6)='pp'
        plabel(7)='pp'
        
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        if     (nproc .eq. 321) then
c--  321 '  f(p1)+b(p2) --> W^+(-->nu(p3)+e^+(p4))+c(p5)+f(p6)'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 326) then
c--  326 '  f(p1)+b(p2) --> W^-(-->e^-(p3)+nu~(p4))+c(p5)+f(p6)'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
        
        if (removebr) then
c--      '  f(p1)+b(p2) --> W(no BR)+c(p5)+f(p6)' (removebr=.true.)
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc .eq. 331) .or. (nproc .eq. 336)) then
        case='Wcjetg'
        nqcdjets=2
        nflav=3
        
        plabel(5)='bq'
        plabel(6)='pp'

        ndim=10
        mb=0
        n2=0
        n3=1
        mass2=0d0
        mass3=wmass
        width3=wwidth
        
        if     (nproc .eq. 331) then
c--  331 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+c(p5)+f(p6) [c-s interaction]'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc .eq. 336) then
c--  336 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+c(p5)+f(p6) [c-s interaction]'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
        
        if (removebr) then
c--      '  f(p1)+f(p2) --> W(no BR)+c(p5)+f(p6) [c-s interaction]' (removebr=.true.)
          call branch(brwen,brzee,brtau,brtop)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif (nproc/10 .ge. 90) then
        write(6,*) 'Setting part to lord and zerowidth to false'
        zerowidth=.false.
        part='lord'
        if (nproc .eq. 902) then
          case='vlchk2'
          nwz=1
          ndim=4
          n3=1
          mass3=wmass
          width3=wwidth
        elseif (nproc .eq. 903) then
          case='vlchk3'
          nwz=1
          ndim=7
          n3=1
          mass3=wmass
          width3=wwidth
        elseif (nproc .eq. 904) then
          case='vlchk4'
          nwz=1
          ndim=10
          n2=1
          n3=1
          mass2=zmass
          width2=zwidth
          mass3=wmass
          width3=wwidth
        elseif (nproc .eq. 905) then
          case='vlchk5'
          nwz=1
          ndim=13
          n2=1
          n3=1
          mass2=hmass
          width2=hwidth
          mass3=wmass
          width3=wwidth
        elseif (nproc .eq. 906) then
          case='vlchk6'
          nwz=1
          ndim=16
          n2=1
          n3=1
          mass2=mt
          width2=twidth
          mass3=mt
          width3=twidth
        elseif (nproc .eq. 908) then
          case='vlchk8'
          nwz=1
          ndim=22
          n2=1
          n3=1
          mass2=mt
          width2=twidth
          mass3=mt
          width3=twidth
        elseif (nproc .eq. 909) then
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
        elseif (nproc .eq. 910) then
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
        call nprocinvalid()
      endif

c--- report on the removed BR, if necessary
      if (removebr) then
        write(6,*)'****************************************************'
        write(6,98) BrnRat
        write(6,*)'****************************************************'
      endif

c--- check that calculation can be performed
      call checkorder()

c--- fill up CKM matrix
      call ckmfill(nwz)

c--- set flags to true unless we're doing W+2 jet or Z+2 jet
      if ((case .ne. 'W_2jet') .and. (case .ne. 'Z_2jet')) then
        Qflag=.true.
        Gflag=.true.
      endif
      
      return

 43   write(6,*) 'problems opening process.DAT'
      stop

 44   write(6,*) 'Unimplemented process number, nproc = ',nproc, 
     . ' mcfm halted'
      stop
 
 98   format(' *              Brn.Rat. removed = ',  f8.4, '         *')
     
      end

      subroutine nprocinvalid()
      implicit none
      integer nproc
      common/nproc/nproc

      write(6,*) 'chooser: Unimplemented case'
      write(6,*) 'nproc=',nproc      
      stop
      
      return 
      end
      
      subroutine checkminzmass(i)
c--- Checks that the minimum invariant mass specified in the options
c--- file is not zero for boson 34 (i=1) or boson 56 (i=2)
      implicit none
      include 'limits.f'
      integer i
      
      if ((i .eq. 1) .and. (wsqmin .eq. 0d0)) then
        write(6,*)
        write(6,*) 'Please set m34min not equal to zero to'
        write(6,*) 'prevent the virtual photon from becoming real.'
        stop
      endif

      if ((i .eq. 2) .and. (bbsqmin .eq. 0d0)) then
        write(6,*)
        write(6,*) 'Please set m56min not equal to zero to'
        write(6,*) 'prevent the virtual photon from becoming real.'
        stop
      endif
      
      return
      end
      
      
