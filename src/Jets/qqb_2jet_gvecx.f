      subroutine qqb_2jet_gvecx(p,n,in,msqvx)
      implicit none
c----Matrix element for 2jet production
C----averaged over initial colours and spins
c    line contracted with the vector n(mu)

      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'

C ip is the label of the emitter
      integer i,j,k,l,in,cs,m(-nf:nf)
      double precision p(mxpart,4)
      double precision fac,prop,n(4)
      double precision msqvx(0:2,-1:1,-1:1,-1:1,-1:1)
      double precision 
     .  aq_gg(0:2),qa_gg(0:2),gg_gg(0:2),gg_qa(0:2),
     .  gq_gq(0:2),gq_qg(0:2),qg_qg(0:2),qg_gq(0:2),
     .  ga_ga(0:2),ga_ag(0:2),ag_ag(0:2),ag_ga(0:2),
     .  ss,tt,uu,nDn,nDp4 
      data m/-1,-1,-1,-1,-1,0,1,1,1,1,1/

      do cs=0,2
      do i=-1,1
      do j=-1,1
      do k=-1,1
      do l=-1,1
        msqvx(cs,i,j,k,l)=0d0
      enddo
      enddo
      enddo
      enddo
      enddo

C--in is the label of the parton dotted with n

      call dotem(3,p,s)
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      nDn=n(4)**2-n(1)**2-n(2)**2-n(3)**2
      fac=-0.5d0*nDn
      nDp4=n(4)*p(in,4)-n(3)*p(in,3)-n(2)*p(in,2)-n(1)*p(in,1)

c--- appropriate scale is approx 1d-3*energy(incoming)
c--- so of order(1) for the Tevatron
      if (abs(nDp4).gt.1d-3*abs(p(1,4))) then 
         write(*,*) 'Error for in=',in
         write(*,*) 'cutoff',1d-3*abs(p(in,4))
         write(6,*) 'nDp4',nDp4
         call flush(6)
         stop
      endif


      call tinyc(tt,ss,uu,qg_qg)
      
      call tinyc(uu,ss,tt,qg_gq)
      call tinyc(tt,uu,ss,ag_ag)
      call tinyc(uu,tt,ss,ag_ga)

      call tinyc(tt,uu,ss,gq_gq)
      call tinyc(uu,tt,ss,gq_qg)
      call tinyc(tt,ss,uu,ga_ga)
      call tinyc(uu,ss,tt,ga_ag)

      call tinyc(uu,ss,tt,qg_gq)
      call tinyc(uu,tt,ss,ag_ga)
      call tinyc(uu,tt,ss,gq_qg)
      call tinyc(uu,ss,tt,ga_ag)

      call tinyc(ss,tt,uu,aq_gg)
      call tinyc(ss,tt,uu,qa_gg)
      call tinyc(ss,uu,tt,gg_qa)

      call tinyd(ss,tt,uu,gg_gg)

      do cs=0,2
      qg_qg(cs)=-aveqg*fac*qg_qg(cs)
      ag_ag(cs)=-aveqg*fac*ag_ag(cs)
      gq_gq(cs)=-aveqg*fac*gq_gq(cs)
      ga_ga(cs)=-aveqg*fac*ga_ga(cs)
      aq_gg(cs)=+aveqq*fac*aq_gg(cs)*half
      qa_gg(cs)=+aveqq*fac*qa_gg(cs)*half
      gg_qa(cs)=+avegg*fac*gg_qa(cs)
      gg_gg(cs)=+avegg*fac*gg_gg(cs)*half

      gq_qg(cs)=-aveqg*fac*gq_qg(cs)
      ga_ag(cs)=-aveqg*fac*ga_ag(cs)
      qg_gq(cs)=-aveqg*fac*qg_gq(cs)
      ag_ga(cs)=-aveqg*fac*ag_ga(cs)

c      do j=-nf,nf
c      do k=-nf,nf
      do j=-1,1
      do k=-1,1

c--- contracted with initial state momentum 1
      if    (in.eq.1) then
         if (j.eq.0) then
              if (k.eq.0) then
                msqvx(cs,m(j),m(k),0,0)=gg_gg(cs)
                msqvx(cs,m(j),m(k),1,-1)=gg_qa(cs)
              endif
              if (k.gt.0) then
                msqvx(cs,m(j),m(k),m(j),m(k))=gq_gq(cs)
                msqvx(cs,m(j),m(k),m(k),m(j))=gq_qg(cs)
              endif
              if (k.lt.0) then
                msqvx(cs,m(j),m(k),m(j),m(k))=ga_ga(cs)
                msqvx(cs,m(j),m(k),m(k),m(j))=ga_ag(cs)
              endif
         endif

c--- contracted with initial state momentum 2
      elseif (in.eq.2) then
         if (k.eq.0) then
              if (j.eq.0) then
                msqvx(cs,m(j),m(k),0,0)=gg_gg(cs)
                msqvx(cs,m(j),m(k),1,-1)=gg_qa(cs)
              endif
              if (j.gt.0) then
                msqvx(cs,m(j),m(k),m(j),m(k))=qg_qg(cs)
                msqvx(cs,m(j),m(k),m(k),m(j))=qg_gq(cs)
              endif
              if (j.lt.0) then
                msqvx(cs,m(j),m(k),m(j),m(k))=ag_ag(cs)
                msqvx(cs,m(j),m(k),m(k),m(j))=ag_ga(cs)
              endif
         endif

c--- contracted with final state momentum 3
      elseif (in.eq.3) then
C--qa      
         if    ((j .gt. 0) .and. (k .eq. -j)) then
              msqvx(cs,m(j),m(k),0,0)=qa_gg(cs)
C--aq      
         elseif ((j .lt. 0) .and. (k .eq. -j)) then
              msqvx(cs,m(j),m(k),0,0)=aq_gg(cs)
C--gq_gq      
         elseif ((j .eq. 0) .and. (k .gt. 0)) then
              msqvx(cs,m(j),m(k),m(j),m(k))=gq_gq(cs)
C--qg_gq      
         elseif ((k .eq. 0) .and. (j .gt. 0)) then
              msqvx(cs,m(j),m(k),m(k),m(j))=qg_gq(cs)
C--ga_ga      
         elseif ((j .eq. 0) .and. (k .lt. 0)) then
              msqvx(cs,m(j),m(k),m(j),m(k))=ga_ga(cs)
C--ag_ga      
         elseif ((k .eq. 0) .and. (j .lt. 0)) then
              msqvx(cs,m(j),m(k),m(k),m(j))=ag_ga(cs)
C--gg_gg     
         elseif ((j .eq. 0) .and. (k .eq. 0)) then
              msqvx(cs,m(j),m(k),m(j),m(k))=gg_gg(cs)
         endif

c--- contracted with final state momentum 4
      elseif (in.eq.4) then
C--qa      
         if    ((j .gt. 0) .and. (k .eq. -j)) then
              msqvx(cs,m(j),m(k),0,0)=qa_gg(cs)
C--aq      
         elseif ((j .lt. 0) .and. (k .eq. -j)) then
              msqvx(cs,m(j),m(k),0,0)=aq_gg(cs)
C--qg_qg      
         elseif ((j .gt. 0) .and. (k .eq. 0)) then
              msqvx(cs,m(j),m(k),m(j),m(k))=qg_qg(cs)
C--gq_qg      
         elseif ((k .gt. 0) .and. (j .eq. 0)) then
              msqvx(cs,m(j),m(k),m(k),m(j))=gq_qg(cs)
C--ag_ag      
         elseif ((j .lt. 0) .and. (k .eq. 0)) then
              msqvx(cs,m(j),m(k),m(j),m(k))=ag_ag(cs)
C--ga_ag      
         elseif ((k .lt. 0) .and. (j .eq. 0)) then
              msqvx(cs,m(j),m(k),m(k),m(j))=ga_ag(cs)
C--gg_gg      
         elseif ((j .eq. 0) .and. (k .eq. 0)) then
              msqvx(cs,m(j),m(k),m(j),m(k))=gg_gg(cs)
         endif

      endif

      enddo
      enddo

      enddo

      return
      end
