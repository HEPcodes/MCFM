      subroutine qqb_2jetx(p,msq,msqx)
      implicit none
c--- matrix element squared and averaged over initial colours and spins
c     f(-p1) + f(-p2) --> f(p5) + f(p6)

      include 'constants.f'
      include 'sprods_com.f'
      include 'msq_cs.f'
      double precision msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
c      double precision mqq(0:2,fn:nf,fn:nf)
c      double precision msqx_cs(0:2,-nf:nf,-nf:nf)
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),
     .  qqjk_jkx(0:2),aajk_jkx(0:2),qqjk_kjx(0:2),aajk_kjx(0:2),
     .  qajk_jkx(0:2),aqjk_jkx(0:2),qajj_kkx(0:2),aqjj_kkx(0:2),
     .  qajk_kjx(0:2),aqjk_kjx(0:2),aqkk_kkx(0:2),
     .  qqjj_jjx(0:2),aajj_jjx(0:2),aqjj_jjx(0:2),qajj_jjx(0:2),
     .  gq_gqx(0:2),ga_gax(0:2),qg_qgx(0:2),ag_agx(0:2),
     .  gq_qgx(0:2),ga_agx(0:2),qg_gqx(0:2),ag_gax(0:2),
     .  aq_ggx(0:2),gg_ggx(0:2),gg_aqx(0:2),qa_ggx(0:2),
     .  ss,tt,uu

      integer j,k,l,m,jf,cs

      call dotem(3,p,s)

      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      call tinya(ss,tt,uu,qqjk_jkx)
      call tinya(ss,tt,uu,aajk_jkx)
      call tinya(tt,ss,uu,qajj_kkx)
      call tinya(ss,tt,uu,qajk_jkx)
      call tinya(uu,ss,tt,aqjj_kkx)
      call tinya(uu,tt,ss,aqjk_jkx)

      call tinya(ss,uu,tt,qqjk_kjx)
      call tinya(ss,uu,tt,aajk_kjx)
      call tinya(ss,uu,tt,qajk_kjx)
      call tinya(tt,uu,ss,aqjk_kjx)

      call tinyb(ss,tt,uu,qqjj_jjx)
      call tinyb(ss,tt,uu,aajj_jjx)
      call tinyb(uu,tt,ss,aqjj_jjx)
      call tinyb(uu,tt,ss,qajj_jjx)

      call tinyb(tt,uu,ss,aqkk_kkx)

      call tinyc(ss,uu,tt,aq_ggx)

      call tinyc(tt,ss,uu,qg_qgx)
      call tinyc(tt,uu,ss,ag_agx)
      call tinyc(tt,uu,ss,gq_gqx)
      call tinyc(tt,ss,uu,ga_gax)

      call tinyc(uu,ss,tt,qg_gqx)
      call tinyc(uu,tt,ss,ag_gax)
      call tinyc(uu,tt,ss,gq_qgx)
      call tinyc(uu,ss,tt,ga_agx)

      call tinyc(ss,tt,uu,qa_ggx)
      call tinyc(ss,uu,tt,gg_aqx)

      call tinyd(ss,tt,uu,gg_ggx)

c      do j=-nf,nf
c      do k=-nf,nf
      do j=-2,2
      do k=-2,2
c--set msqx=0 to initalize
      do cs=0,2
c          do l=-nf,nf
c          do m=-nf,nf
          do l=-2,2
          do m=-2,2
              msqx(cs,j,k,l,m)=0d0
          enddo
          enddo
C--qq      
      if ((j .gt. 0) .and. (k .gt. 0)) then
          if (j .eq. k) then
              msqx(cs,j,j,j,j)=aveqq*qqjj_jjx(cs)*half
          else
              msqx(cs,j,k,j,k)=aveqq*qqjk_jkx(cs)
c--- note that this spoils adding up for total msq, below
              msqx(cs,j,k,k,j)=aveqq*qqjk_kjx(cs)
          endif

C--qa      
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
            do jf=0,nf
            if (jf .eq. 0) then
              msqx(cs,j,k,0,0)=aveqq*qa_ggx(cs)*half
            elseif (jf .eq. j) then 
              msqx(cs,j,k,jf,-jf)=aveqq*qajj_jjx(cs)
            elseif (jf.ne.j) then 
              msqx(cs,j,k,jf,-jf)=aveqq*qajj_kkx(cs)
            endif  
            enddo
          else
              msqx(cs,j,k,j,k)=aveqq*qajk_jkx(cs)
c--- note that this spoils adding up for total msq, below
              msqx(cs,j,k,k,j)=aveqq*qajk_kjx(cs)
          endif

C--aa      
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (j .eq. k) then
              msqx(cs,j,j,j,j)=aveqq*aajj_jjx(cs)*half
          else
              msqx(cs,j,k,j,k)=aveqq*aajk_jkx(cs)
c--- note that this spoils adding up for total msq, below
              msqx(cs,j,k,k,j)=aveqq*aajk_kjx(cs)
          endif

C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) then
            do jf=0,nf
            if (jf .eq. 0) then
            msqx(cs,j,k,0,0)=aveqq*aq_ggx(cs)*half
            elseif (jf .eq. -j) then 
            msqx(cs,j,k,-jf,jf)=aveqq*aqjj_jjx(cs)
c--- note that this spoils adding up for total msq, below
            msqx(cs,j,k,jf,-jf)=aveqq*aqkk_kkx(cs)
            elseif (jf.ne.-j) then 
            msqx(cs,j,k,-jf,jf)=aveqq*aqjj_kkx(cs)
c--- note that this spoils adding up for total msq, below
            msqx(cs,j,k,jf,-jf)=aveqq*aqjj_kkx(cs)
            endif  
            enddo
          else
            msqx(cs,j,k,j,k)=aveqq*aqjk_jkx(cs)
c--- note that this spoils adding up for total msq, below
              msqx(cs,j,k,k,j)=aveqq*aqjk_kjx(cs)
          endif

C--qg_qg      
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msqx(cs,j,k,j,k)=-aveqg*qg_qgx(cs)
c--- note that this spoils adding up for total msq, below
            msqx(cs,j,k,k,j)=-aveqg*qg_gqx(cs)
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msqx(cs,j,k,j,k)=-aveqg*ag_agx(cs)
c--- note that this spoils adding up for total msq, below
            msqx(cs,j,k,k,j)=-aveqg*ag_gax(cs)
C--gq_gq      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msqx(cs,j,k,j,k)=-aveqg*gq_gqx(cs)
c--- note that this spoils adding up for total msq, below
            msqx(cs,j,k,k,j)=-aveqg*gq_qgx(cs)
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msqx(cs,j,k,j,k)=-aveqg*ga_gax(cs)
c--- note that this spoils adding up for total msq, below
            msqx(cs,j,k,k,j)=-aveqg*ga_agx(cs)
C--gg      
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
            do jf=0,nf
            if (jf .eq. 0) then
            msqx(cs,j,k,0,0)=avegg*gg_ggx(cs)*half
            else
            msqx(cs,j,k,jf,-jf)=avegg*gg_aqx(cs)
            endif
            enddo
      endif

      enddo
      enddo
      enddo

c--- now fill up msq_cs, which will be used in the virtual routine
      do j=-nf,nf
      do k=-nf,nf
      do cs=0,2
        msq_cs(cs,j,k)=0d0
        msq_cs(cs,0,0)=avegg*(gg_ggx(cs)*half+dfloat(nf)*gg_aqx(cs))
      enddo
      enddo
      enddo
      
c      do j=-nf,nf
c      do k=-nf,nf
c      msq(j,k)=0d0
c      do cs=0,2
c      do l=-nf,nf
c      do m=-nf,nf
c      msq(j,k)=msq(j,k)+msqx(cs,j,k,l,m)
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo

      return
      end



     
