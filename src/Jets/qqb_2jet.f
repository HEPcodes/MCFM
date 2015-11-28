      subroutine qqb_2jet(p,msq)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     .  qqij_ij,aaij_ij,
     .  qaij_ij,aqij_ij,aqii_jj,qaii_jj,
     .  qqii_ii,aaii_ii,aqii_ii,qaii_ii,
     .  aq_gg,gq_gq,ga_ga,qg_qg,ag_ag,gg_gg,gg_qa,qa_gg,ss,tt,uu,
     .  smalla,smallb,smallc,smalld 
      call dotem(4,p,s)
      fac=gsq**2
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      qqij_ij=fac*aveqq*smalla(ss,tt,uu)
      aaij_ij=fac*aveqq*smalla(ss,tt,uu)
      qaii_jj=fac*aveqq*smalla(tt,ss,uu)
      qaij_ij=fac*aveqq*smalla(ss,tt,uu)
      aqii_jj=fac*aveqq*smalla(tt,ss,uu)
      aqij_ij=fac*aveqq*smalla(uu,tt,ss)

      qqii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aaii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aqii_ii=fac*aveqq*smallb(uu,tt,ss)
      qaii_ii=fac*aveqq*smallb(uu,tt,ss)

      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half

      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      gq_gq=-fac*aveqg*smallc(tt,uu,ss)
      ga_ga=-fac*aveqg*smallc(tt,ss,uu)

      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)

      gg_gg=fac*avegg*smalld(ss,tt,uu)*half

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
C--qq      
      if ((j .gt. 0) .and. (k .gt. 0)) then
          if (j .eq. k) then
            msq(j,k)=qqii_ii
          else
            msq(j,k)=qqij_ij
          endif

C--qa      
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=qaii_ii+dfloat(nf-1)*qaii_jj+qa_gg
          else
            msq(j,k)=qaij_ij
          endif

C--aa      
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (j .eq. k) then
            msq(j,k)=aaii_ii
          else
            msq(j,k)=aaij_ij
          endif

C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=aqii_ii+dfloat(nf-1)*aqii_jj+aq_gg
          else
            msq(j,k)=aqij_ij
          endif

C--qg_qg      
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=qg_qg
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=ag_ag
C--gq_gq      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=gq_gq
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=ga_ga
C--gg      
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=gg_gg+dfloat(nf)*gg_qa
      endif

      enddo
      enddo


      return
      end



