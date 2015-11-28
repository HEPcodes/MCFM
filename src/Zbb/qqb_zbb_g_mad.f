      subroutine qqb_zbb_g_mad(p,msq)
************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2001.                                                      *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  Z^0 + bbar(p5)+b(p6)+p(p7)
c                           |
c                            --> e^-(p3)+e^+(p4)
c                           
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'prods.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'hardscale.f'
      include 'flags.f'
      integer i,j,k,nu,hq,Qh,hg,lh,n1,n2,nquark,j1,k2
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf)
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3),P7(0:3)
      double precision suu_emepuug,suc_emepucg,
     . sdd_emepddg,sds_emepdsg,
     . sud_emepudg,sdu_emepdug,
     . subub_emepububg,subcb_emepubcbg,
     . sdbdb_emepdbdbg,sdbsb_emepdbsbg,
     . subdb_emepubdbg,sdbub_emepdbubg,
     . sddb_emepddbg,suub_emepuubg,
     . sdsb_emepdsbg,sucb_emepucbg,
     . sudb_emepudbg,sdub_emepdubg,
     . sdbd_emepdbdg,subu_emepubug,
     . sdbs_emepdbsg,subc_emepubcg,
     . subd_emepubdg,sdbu_emepdbug,
     . suub_emepddbg,suub_emepccbg,
     . sddb_emepuubg,sddb_emepssbg,
     . sdbd_emepubug,sdbd_emepsbsg,
     . subu_emepdbdg,subu_emepcbcg,
     . sug_emepuubu,sug_emepddbu,sug_emepccbu,
     . sdg_emepddbd,sdg_emepuubd,sdg_emepssbd,
     . sgu_emepuubu,sgu_emepddbu,sgu_emepccbu,
     . sgd_emepddbd,sgd_emepuubd,sgd_emepssbd,
     . subg_emepuubub,subg_emepddbub,subg_emepccbub,
     . sdbg_emepddbdb,sdbg_emepuubdb,sdbg_emepssbdb,
     . sgub_emepuubub,sgub_emepddbub,sgub_emepccbub,
     . sgdb_emepddbdb,sgdb_emepuubdb,sgdb_emepssbdb,
     . sddb_emepggg,sdbd_emepggg,suub_emepggg,subu_emepggg,
     . sdg_emepdgg,sdbg_emepdbgg,sgd_emepdgg,sgdb_emepdbgg,
     . sug_emepugg,subg_emepubgg,sgu_emepugg,sgub_emepubgg,
     . sgg_emepdbdg,sgg_emepubug
      integer jj(-nf:nf),kk(-nf:nf)
      data jj/-1,-2,-1,-2,-1,0,1,2,1,2,1/
      data kk/-1,-2,-1,-2,-1,0,1,2,1,2,1/

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

      call initialize

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf
       
      if ((j .eq. 0) .and. (k .eq. 0)) then
c-gg
              msq(j,k)=
     .           +sgg_emepdbdg(p1,p2,p3,p4,p5,p6,p7)
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
c-qqb
           if (k.eq.-j) then
              if ((jj(j).eq.1) .and. (kk(k).eq.-1)) then
              msq(j,k)=msq(j,k)
     .           +sddb_emepssbg(p1,p2,p3,p4,p6,p5,p7)
              elseif ((jj(j).eq.2) .and. (kk(k).eq.-2)) then
              msq(j,k)=msq(j,k)
     .           +suub_emepddbg(p1,p2,p3,p4,p6,p5,p7)
              endif
          endif
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
c-qbq
          if (j .eq.-k) then
              if ((jj(j).eq.-1) .and. (kk(k).eq.1)) then
              msq(j,k)=msq(j,k)
     .           +sdbd_emepsbsg(p1,p2,p3,p4,p5,p6,p7)
              elseif ((jj(j).eq.-2) .and. (kk(k).eq.2)) then
              msq(j,k)=msq(j,k)
     .           +subu_emepdbdg(p1,p2,p3,p4,p5,p6,p7)
              endif
          endif
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
C-qg
            if( jj(j) .eq. 1) then
              msq(j,k)=msq(j,k)
     .        +sdg_emepssbd(p1,p2,p3,p4,p6,p5,p7)
            elseif(jj(j) .eq. 2) then
              msq(j,k)=msq(j,k)
     .        +sug_emepddbu(p1,p2,p3,p4,p6,p5,p7)
             endif
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
c-qbg
            if( jj(j) .eq. -1) then
              msq(j,k)=msq(j,k)
     .        +sdbg_emepssbdb(p1,p2,p3,p4,p6,p5,p7)
            elseif(jj(j) .eq. -2) then
              msq(j,k)=msq(j,k)
     .        +subg_emepddbub(p1,p2,p3,p4,p6,p5,p7)
             endif
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
c-gq
            if( kk(k) .eq. 1) then
              msq(j,k)=msq(j,k)
     .        +sgd_emepssbd(p1,p2,p3,p4,p6,p5,p7)
            elseif(kk(k) .eq. 2) then
              msq(j,k)=msq(j,k)
     .        +sgu_emepddbu(p1,p2,p3,p4,p6,p5,p7)
             endif
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
c-gqb
            if( kk(k) .eq. -1) then
              msq(j,k)=msq(j,k)
     .        +sgdb_emepssbdb(p1,p2,p3,p4,p6,p5,p7)
            elseif(kk(k) .eq. -2) then
              msq(j,k)=msq(j,k)
     .        +sgub_emepddbub(p1,p2,p3,p4,p6,p5,p7)
             endif
      endif
      enddo
      enddo




      return
      end
