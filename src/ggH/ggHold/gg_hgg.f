      subroutine gg_hgg(p,msq)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c
c     g(-p1)+g(-p2) -->  H(p3)+g(p_iglue1=5)+g(p_iglue2=6) 

      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'msq_struc.f'
      integer j,k,iglue1,iglue2
      double precision p(mxpart,4),Asq,fac
      double precision Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625,
     .                       Hgggg_1652,Hgggg_1562,Hgggg_1526
      double precision Hqagg,Haqgg,Hgqqg,Hgaag,Hqgqg,Hagag,Hggqa
      double precision Hggqa_ab,Hggqa_ba,Hggqa_sym
      double precision Hqgqg_ab,Hqgqg_ba,Hqgqg_sym
      double precision Hgqqg_ab,Hgqqg_ba,Hgqqg_sym
      double precision Hagag_ab,Hagag_ba,Hagag_sym
      double precision Hgaag_ab,Hgaag_ba,Hgaag_sym
      double precision Hqagg_ab,Hqagg_ba,Hqagg_sym
      double precision Haqgg_ab,Haqgg_ba,Haqgg_sym
      double precision Hqqqq_a,Hqqqq_b,Hqqqq_i
      double precision Hqaqa_a,Hqaqa_b,Hqaqa_i
      double precision Haqaq_a,Haqaq_b,Haqaq_i
      double precision Hqaaq_a,Hqaaq_b,Hqaaq_i
      double precision 
     . Hqrqr,Hqqqq,
     . Habab,Haaaa,
     . Hqarb,Hqaqa,Hqbqb,
     . Haqbr,Haqaq,Hbqbq,
     . Hqaaq
      double precision msq(-nf:nf,-nf:nf),hdecay,s34

      parameter(iglue1=5,iglue2=6)


C---fill spinor products up to maximum number
      call spinoru(iglue2,p,za,zb)  

C   Deal with Higgs decay to b-bbar
      s34=s(3,4)+2d0*mb**2
      hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s34-4d0*mb**2) 
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      Asq=(as/(3d0*pi))**2/vevsq
      fac=gsq**2*Asq*hdecay
C--four gluon terms
      call h4g(1,2,iglue1,iglue2,Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625)
      
C--two quark two gluon terms
      call hqqgg(1,2,iglue1,iglue2,Hqagg,Hqagg_ab,Hqagg_ba,Hqagg_sym)
      call hqqgg(2,1,iglue1,iglue2,Haqgg,Haqgg_ab,Haqgg_ba,Haqgg_sym)
C====symmetric in first two arguments, but not the ab, ba terms
c      Haqgg=Hqagg

      call hqqgg(1,iglue1,2,iglue2,Hqgqg,Hqgqg_ab,Hqgqg_ba,Hqgqg_sym)
      call hqqgg(iglue1,1,2,iglue2,Hagag,Hagag_ab,Hagag_ba,Hagag_sym)
C====symmetric in first two arguments
c      Hagag=Hqgqg

      call hqqgg(2,iglue1,1,iglue2,Hgqqg,Hgqqg_ab,Hgqqg_ba,Hgqqg_sym)
      call hqqgg(iglue1,2,1,iglue2,Hgaag,Hgaag_ab,Hgaag_ba,Hgaag_sym)
C====symmetric in first two arguments
c      Hgaag=Hgqqg

      call hqqgg(iglue2,iglue1,1,2,Hggqa,Hggqa_ab,Hggqa_ba,Hggqa_sym)
      
C---four quark terms
      call H4qn(1,2,iglue1,iglue2,Hqrqr)
      call H4qi(1,2,iglue1,iglue2,Hqqqq,Hqqqq_a,Hqqqq_b,Hqqqq_i)
C---four anti-quark terms
c      call H4qn(iglue1,iglue2,1,2,Habab)
c      call H4qi(iglue1,iglue2,1,2,Haaaa)
      Habab=Hqrqr
      Haaaa=Hqqqq

C-qqb
      call H4qn(1,iglue2,2,iglue1,Hqarb)
      call H4qi(1,iglue2,iglue1,2,Hqaqa,Hqaqa_a,Hqaqa_b,Hqaqa_i)
      call H4qn(1,iglue2,iglue1,2,Hqbqb)
c      write(6,*) 'Hqaqa',Hqaqa_a,Hqaqa_b,Hqaqa_i
c      write(6,*) 'Hqbqb',Hqbqb
C-qbq
      Haqbr=Hqarb
      
      Haqaq=Hqaqa
      Haqaq_a=Hqaqa_a
      Haqaq_b=Hqaqa_b
      Haqaq_i=Hqaqa_i
      Hbqbq=Hqbqb

c--- subtraction matrix elements use qa->aq; calculate this and
c--- artificially store it in msq_struc(iqr,0,0), which is not
c--- used for anything else
      call H4qi(1,iglue1,iglue2,2,Hqaaq,Hqaaq_a,Hqaaq_b,Hqaaq_i)
      msq_struc(iqr,0,0)=aveqq*fac*Hqaaq
      
      do j=fn,nf
      do k=fn,nf
      msq(j,k)=0d0

      if ((j.gt.0).and.(k.gt.0)) then 
        if (j.eq.k) then
          msq(j,k)=0.5d0*aveqq*fac*Hqqqq
          msq_struc(iqq_a,j,k)=0.5d0*aveqq*fac*Hqqqq_a
          msq_struc(iqq_b,j,k)=0.5d0*aveqq*fac*Hqqqq_b
          msq_struc(iqq_i,j,k)=0.5d0*aveqq*fac*Hqqqq_i
        else
          msq(j,k)=aveqq*fac*Hqrqr
          msq_struc(iqq_a,j,k)=msq(j,k)
        endif
      endif
      
      if ((j.lt.0).and.(k.lt.0)) then 
        if (j.eq.k) then
          msq(j,k)=0.5d0*aveqq*fac*Haaaa
        else
          msq(j,k)=aveqq*fac*Habab
          msq_struc(iqq_a,j,k)=msq(j,k)
        endif
      endif

      if ((j.gt.0).and.(k.lt.0)) then
        if (j.eq.-k) then
          msq(j,k)=aveqq*fac*(0.5d0*Hqagg+Hqaqa+dfloat(nf-1)*Hqarb)
          msq_struc(iqr,j,k)=aveqq*fac*dfloat(nf-1)*Hqarb
          msq_struc(iqq_a,j,k)=aveqq*fac*Hqaqa_a
          msq_struc(iqq_b,j,k)=aveqq*fac*Hqaqa_b
          msq_struc(iqq_i,j,k)=aveqq*fac*Hqaqa_i
          msq_struc(igg_ab,j,k)=aveqq*fac*0.5d0*Hqagg_ab
          msq_struc(igg_ba,j,k)=aveqq*fac*0.5d0*Hqagg_ba
          msq_struc(igg_sym,j,k)=aveqq*fac*0.5d0*Hqagg_sym
        else
          msq(j,k)=aveqq*fac*Hqbqb
          msq_struc(iqq_a,j,k)=msq(j,k)
        endif
      endif

      if ((j.lt.0).and.(k.gt.0)) then
        if (j.eq.-k) then
          msq(j,k)=aveqq*fac*(0.5d0*Haqgg+Haqaq+dfloat(nf-1)*Haqbr)
          msq_struc(iqr,j,k)=aveqq*fac*dfloat(nf-1)*Haqbr
          msq_struc(iqq_a,j,k)=aveqq*fac*Haqaq_a
          msq_struc(iqq_b,j,k)=aveqq*fac*Haqaq_b
          msq_struc(iqq_i,j,k)=aveqq*fac*Haqaq_i
          msq_struc(igg_ab,j,k)=aveqq*fac*0.5d0*Haqgg_ab
          msq_struc(igg_ba,j,k)=aveqq*fac*0.5d0*Haqgg_ba
          msq_struc(igg_sym,j,k)=aveqq*fac*0.5d0*Haqgg_sym
        else
          msq(j,k)=aveqq*fac*Hbqbq
          msq_struc(iqq_a,j,k)=msq(j,k)
        endif
      endif

      if ((j.gt.0).and.(k.eq.0)) then
        msq(j,0)=aveqg*fac*Hqgqg
        msq_struc(igg_ab,j,0)=aveqg*fac*Hqgqg_ab
        msq_struc(igg_ba,j,0)=aveqg*fac*Hqgqg_ba
        msq_struc(igg_sym,j,0)=aveqg*fac*Hqgqg_sym
      endif
      
      if ((j.lt.0).and.(k.eq.0)) then
        msq(j,0)=aveqg*fac*Hagag
        msq_struc(igg_ab,j,0)=aveqg*fac*Hagag_ab
        msq_struc(igg_ba,j,0)=aveqg*fac*Hagag_ba
        msq_struc(igg_sym,j,0)=aveqg*fac*Hagag_sym
      endif

      if ((j.eq.0).and.(k.gt.0)) then
        msq(0,k)=aveqg*fac*Hgqqg
        msq_struc(igg_ab,0,k)=aveqg*fac*Hgqqg_ab
        msq_struc(igg_ba,0,k)=aveqg*fac*Hgqqg_ba
        msq_struc(igg_sym,0,k)=aveqg*fac*Hgqqg_sym
      endif

      if ((j.eq.0).and.(k.lt.0)) then
        msq(0,k)=aveqg*fac*Hgaag
        msq_struc(igg_ab,0,k)=aveqg*fac*Hgaag_ab
        msq_struc(igg_ba,0,k)=aveqg*fac*Hgaag_ba
        msq_struc(igg_sym,0,k)=aveqg*fac*Hgaag_sym
      endif

      if ((j.eq.0).and.(k.eq.0)) then
        msq(0,0)=avegg*fac*(0.5d0*Hgggg+dfloat(nf)*Hggqa)
        msq_struc(igg_ab,0,0)=avegg*fac*dfloat(nf)*Hggqa_ab
        msq_struc(igg_ba,0,0)=avegg*fac*dfloat(nf)*Hggqa_ba
        msq_struc(igg_sym,0,0)=avegg*fac*dfloat(nf)*Hggqa_sym
        msq_struc(igggg_a,0,0)=avegg*fac*0.5d0*Hgggg_1256
        msq_struc(igggg_b,0,0)=avegg*fac*0.5d0*Hgggg_1625
        msq_struc(igggg_c,0,0)=avegg*fac*0.5d0*Hgggg_1265
      endif
      
      enddo
      enddo

      return
      end

 
