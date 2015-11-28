      subroutine compare(p)
      implicit none
      
      include 'constants.f'
      include 'zerowidth.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
 
    
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),
     . msqa(-nf:nf,-nf:nf)
      integer i,j
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3),P7(0:3)                     
      real*8 aemmz,ans,Vfac,ee
      real*8 SUUB_VEEPDUB,SUUB_EMVEUDB,SUUB_VEEPSCB,SUDB_VEEPSSB,
     .       SUDB_VEEPUUB,SUDB_VEEPDDB,SUSB_VEEPUCB,SCDB_VEEPSDB,
     .       SUBU_VEEPSCB,SUBU_VEEPDUB,SUBU_EMVEUDB,
     .       SUU_VEEPDU,SDU_EMVEUU,SUD_EMVEUU,SSU_EMVECU,SUS_EMVEUC,
     .       SUBUB_EMVEDBUB,SSBUB_VEEPCBUB,SUBSB_VEEPUBCB,
     .       SDBUB_VEEPUBUB,SUBDB_VEEPUBUB,SUDB_VEEPGGG,SDBU_VEEPGGG,
     .       SUDB_VEEPGG,SDBU_VEEPGG,SUG_VEEPDG,
     .       SDBG_VEEPUBG,SGU_VEEPDG,SGDB_VEEPUBG,SGG_VEEPUBD,
     . sub_emepub,suu_emepuu,suc_emepuc,
     . subcb_emepubcb,subub_emepubub,sucb_emepucb,subc_emepubc,
     . suub_emepuub,subu_emepubu,
     . subu_emepccb,subu_emepbbb,suub_emepccb,suub_emepbbb,
     . subu_emepcbc,suub_emepggg,sddb_emepggg,subu_emepggg,sdbd_emepggg,
     . sgg_emepsbsg,sgg_emepcbcg,
     . sgd_emepdgg,sgu_emepugg,sgdb_emepdbgg,sgub_emepubgg,
     . sdg_emepdgg,sug_emepugg,sdbg_emepdbgg,subg_emepubgg,
     . sddb_emepbbbg,sdbd_emepbbbg,sgdb_emepbbbdb,sdbg_emepbbbdb,
     . sgd_emepbbbd,sdg_emepbbbd,sdb_emepbdg,sdbbb_emepbbdbg,
     . suub_emepbbbg,subu_emepbbbg,sdd_emepddg,sdbdb_emepdbdbg
      integer nwz,n1,n2
      common/nwz/nwz

c--- implement the momentum exchange      
      do i=1,4
        if (i.lt.4) then
          j=i
        else
          j=0
        endif 
        p1(j)=-p(1,i)
        p2(j)=-p(2,i)
        p3(j)=p(3,i)
        p4(j)=p(4,i)
        p5(j)=p(5,i)
        p6(j)=p(6,i)
        p7(j)=p(7,i)
      enddo

c--- make call to MCFM routines
      aemmz=1d0/128d0      
      ee=sqrt(fourpi*aemmz)
      write(6,*) 'ee',ee
      gwsq=fourpi*aemmz/xw
      gw=sqrt(gwsq)
      gsq=1d0
      call ckmfill(nwz)
      
      write(*,*) 'zmass, wmass',zmass,wmass
      write(*,*) 'new xw,gw',xw,gw

      write(*,*)

      call qqb_z2jet_g(p,msq)
      
      call initialize

      ans=sddb_emepbbbg(p1,p2,p3,p4,p5,p6,p7)
      write(*,*) 'madgraph (1,-1)',ans,ans/msq(1,-1)
      ans=suub_emepbbbg(p1,p2,p3,p4,p5,p6,p7)
      write(*,*) 'madgraph (2,-2)',ans,ans/msq(2,-2)
      ans=sdbd_emepbbbg(p1,p2,p3,p4,p5,p6,p7)
      write(*,*) 'madgraph (-1,1)',ans,ans/msq(-1,1)
      ans=subu_emepbbbg(p1,p2,p3,p4,p5,p6,p7)
      write(*,*) 'madgraph (-2,2)',ans,ans/msq(-2,2)
      ans=sgdb_emepbbbdb(p1,p2,p3,p4,p5,p6,p7)
      write(*,*) 'madgraph (0,-1)',ans,ans/msq(0,-1)
      ans=sdbg_emepbbbdb(p1,p2,p3,p4,p5,p6,p7)
      write(*,*) 'madgraph (-1,0)',ans,ans/msq(-1,0)
      ans=sgd_emepbbbd(p1,p2,p3,p4,p5,p6,p7)
      write(*,*) 'madgraph (0,1)',ans,ans/msq(0,1)
      ans=sdg_emepbbbd(p1,p2,p3,p4,p5,p6,p7)
      write(*,*) 'madgraph (1,0)',ans,ans/msq(1,0)
      ans=sdbbb_emepbbdbg(p1,p2,p3,p4,p5,p6,p7)
      write(*,*) 'madgraph (-1,-5)',ans,ans/msq(-1,-5)
      ans=sdb_emepbdg(p1,p2,p3,p4,p5,p6,p7)
      write(*,*) 'madgraph (1,5)',ans,ans/msq(1,5)
      ans=sdd_emepddg(p1,p2,p3,p4,p5,p6,p7)
      write(*,*) 'madgraph (1,1)',ans,ans/msq(1,1)
      ans=sdbdb_emepdbdbg(p1,p2,p3,p4,p5,p6,p7)
      write(*,*) 'madgraph (-1,-1)',ans,ans/msq(-1,-1)

c      call qqb_zbb_g(p,msqa)
c      
c      write(*,*) '        MCFM new         zbb_g            ratio'
c      write(*,99) '(d db)',msq(1,-1),msqa(1,-1),msq(1,-1)/msqa(1,-1)
c      write(*,99) '(db d)',msq(-1,1),msqa(-1,1),msq(-1,1)/msqa(-1,1)
c      write(*,99) '(g db)',msq(0,-1),msqa(0,-1),msq(0,-1)/msqa(0,-1)
c      write(*,99) '(db g)',msq(-1,0),msqa(-1,0),msq(-1,0)/msqa(-1,0)
c      write(*,99) '(g d)',msq(0,1),msqa(0,1),msq(0,1)/msqa(0,1)
c      write(*,99) '(d g)',msq(1,0),msqa(1,0),msq(1,0)/msqa(1,0)
c      write(*,99) '(d b)',msq(1,5),msqa(1,5),msq(1,5)/msqa(1,5)
c      write(*,99)'(db bb)',msq(-1,-5),msqa(-1,-5),msq(-1,-5)/msqa(-1,-5)
c
c      pause

   99 format(a7,2e17.9,f13.9)

c------------checked
c      ans=2d0*sgg_emepcbcg(p1,p2,p3,p4,p5,p6,p7)
c     .   +3d0*sgg_emepsbsg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'MCFM new (gg)',msq(0,0)      
c      write(*,*) 'ans/msq(0,0)',ans/msq(0,0)
c      ans=suub_emepggg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'ans/msq(2,-2)',ans/msq(2,-2)
c      ans=sddb_emepggg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'ans/msq(1,-1)',ans/msq(1,-1)
c
c      ans=subu_emepggg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'ans/msq(2,-2)',ans/msq(-2,2)
c      ans=sdbd_emepggg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'ans/msq(1,-1)',ans/msq(-1,1)

c      ans=sgd_emepdgg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'ans/msq(0,1)',ans/msq(0,1)
c      ans=sgu_emepugg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'ans/msq(0,2)',ans/msq(0,2)
c      ans=sgdb_emepdbgg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'ans/msq(0,-1)',ans/msq(0,-1)
c      ans=sgub_emepubgg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'ans/msq(0,-2)',ans/msq(0,-2)

c      ans=sdg_emepdgg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'ans/msq(1,0)',ans/msq(1,0)
c      ans=sug_emepugg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'ans/msq(2,0)',ans/msq(2,0)
c      ans=sdbg_emepdbgg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'ans/msq(-1,0)',ans/msq(-1,0)
c      ans=subg_emepubgg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'ans/msq(-2,0)',ans/msq(-2,0)
c      ans=subub_emepubub(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'MCFM new (ub ub)',msq(-2,-2)      
c      write(*,*) 'Madgraph ub ub -> e-, e+, ub, ub',ans
c      ans=subcb_emepubcb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'MCFM new (ub cb)',msq(-2,-4)      
c      write(*,*) 'Madgraph ub cb -> e-, e+, ub, cb',ans
c      ans=sucb_emepucb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'MCFM new (u cb)',msq(2,-4)      
c      write(*,*) 'Madgraph u cb -> e-, e+, u, cb',ans
c      ans=subc_emepubc(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'MCFM new (ub  c)',msq(-2,4)      
c      write(*,*) 'Madgraph ub c -> e-, e+, ub, c',ans
c      ans=suc_emepuc(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'MCFM new (2  4)',msq(2,4)      
c      write(*,*) 'Madgraph u c -> e-, e+, u, c',ans
c      ans=suub_emepuub(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'MCFM new (2  -2)',msq(2,-2)      
c      write(*,*) 'Madgraph u ub -> e-, e+, u, ub',ans
c------checked

c      ans=subu_emepbbb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ub u -> e-, e+, bb, b',ans
c      ans=subu_emepcbc(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ub u -> e-, e+, cb, c',ans

c      ans=suub_emepbbb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph u ub -> e-, e+, b, bb',ans
c      ans=suub_emepccb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph u ub -> e-, e+, c, cb',ans



c      ans=sdbu_veepgg(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph db u -> ve, e+, g, g',ans
c      ans=sgg_veepubd(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph  g g -> ve, e+, ub, d',ans
c      ans=sug_veepdg(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph  u g -> ve, e+, d, g',ans
c      ans=sdbg_veepubg(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph db g -> ve, e+, ub, g',ans
c      ans=sgu_veepdg(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph  g u -> ve, e+, d, g',ans
c      ans=sgdb_veepubg(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph g db -> ve, e+, ub, g',ans
      
      return
c      call qqb_w2jet(p,msq)
c      write(*,*) 'MCFM old (u ub)',msq(2,-2)
c      write(*,*) 'MCFM old (u db)',msq(2,-1)
c      write(*,*) 'MCFM old (u sb)',msq(2,-3)
c      write(*,*) 'MCFM old (ub u)',msq(-2,2)
c      write(*,*) 'MCFM old (ub d)',msq(-2,1)
c      write(*,*) 'MCFM old (ub s)',msq(-2,3)
c      write(*,*) 'MCFM old (u u)',msq(2,2)
c      write(*,*) 'MCFM old (s u)',msq(3,2)
c      write(*,*) 'MCFM old (u s)',msq(2,3)
c      write(*,*) 'MCFM old (d u)',msq(1,2)
c      write(*,*) 'MCFM old (u d)',msq(2,1)
c      write(*,*) 'MCFM old (ub ub)',msq(-2,-2)
c      write(*,*) 'MCFM old (ub sb)',msq(-2,-3)
c      write(*,*) 'MCFM old (sb ub)',msq(-3,-2)
c      write(*,*) 'MCFM old (ub db)',msq(-2,-1)
c      write(*,*) 'MCFM old (db ub)',msq(-1,-2)
c      write(*,*)
c      call qqb_w2jet_2000(p,msq)
c      write(*,*) 'MCFM new (u ub)',msq(2,-2)
c      write(*,*) 'MCFM new (u db)',msq(2,-1)
c      write(*,*) 'MCFM new (u sb)',msq(2,-3)
c      write(*,*) 'MCFM new (ub u)',msq(-2,2)
c      write(*,*) 'MCFM new (ub d)',msq(-2,1)
c      write(*,*) 'MCFM new (ub s)',msq(-2,3)
c      write(*,*) 'MCFM new (u u)',msq(2,2)
c      write(*,*) 'MCFM new (s u)',msq(3,2)
c      write(*,*) 'MCFM new (u s)',msq(2,3)
c      write(*,*) 'MCFM new (d u)',msq(1,2)
c      write(*,*) 'MCFM new (u d)',msq(2,1)
c      write(*,*) 'MCFM new (ub ub)',msq(-2,-2)
c      write(*,*) 'MCFM new (ub sb)',msq(-2,-3)
c      write(*,*) 'MCFM new (sb ub)',msq(-3,-2)
c      write(*,*) 'MCFM new (ub db)',msq(-2,-1)
c      write(*,*) 'MCFM new (db ub)',msq(-1,-2)
       call qqb_w2jet_g(p,msq)
       write(*,*) 'MCFM new (u db)',msq(2,-1)/Vsq(2,-1)
       write(*,*) 'MCFM new (db u)',msq(-1,2)/Vsq(-1,2)
      
      call initialize
c      ans=suub_veepdub(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ve, e+, d, ub',ans
c      ans=suub_emveudb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph e-, ve, u, db',ans
c      ans=suub_veepscb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ve, e+, s, cb',ans
c      ans=sudb_veepssb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ve, e+, s, sb',ans
c      ans=sudb_veepuub(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ve, e+, u, ub',ans
c      ans=sudb_veepddb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ve, e+, d, db',ans
c      ans=susb_veepucb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ve, e+, u, cb',ans
c      ans=scdb_veepsdb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ve, e+, s, db',ans
c      ans=subu_veepscb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ub u -> ve, e+, s, cb',ans
c      ans=subu_veepdub(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ub u -> ve, e+, d, ub',ans
c      ans=subu_emveudb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ub u -> e-, ve, u, db',ans
c      ans=suu_veepdu(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph u u -> ve, e+, d, u',ans
c      ans=sdu_emveuu(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph d u -> e-, ve, u, u',ans
c      ans=sud_emveuu(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph u d -> e-, ve, u, u',ans
c      ans=ssu_emvecu(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph s u -> e-, ve, c, u',ans
c      ans=sus_emveuc(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph u s -> e-, ve, u, c',ans
c      ans=subub_emvedbub(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ub ub -> e-, ve, db, ub',ans
c      ans=ssbub_veepcbub(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph sb ub -> ve, e+, cb, ub',ans
c      ans=subsb_veepubcb(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ub sb -> ve, e+, ub, cb',ans
c      ans=sdbub_veepubub(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph db ub -> ve, e+, ub, ub',ans
c      ans=subdb_veepubub(p1,p2,p3,p4,p5,p6)
c      write(*,*) 'Madgraph ub db -> ve, e+, ub, ub',ans
c      ans=sudb_veepggg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'Madgraph u db -> ve, e+, g, g, g',ans
c      ans=sdbu_veepggg(p1,p2,p3,p4,p5,p6,p7)
c      write(*,*) 'Madgraph db u -> ve, e+, g, g, g',ans
c      pause
      
      return
      end
