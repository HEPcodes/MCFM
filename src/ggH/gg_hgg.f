      subroutine gg_hgg(p,msq)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c
c     g(-p1)+g(-p2) -->  H(p3)+g(p_iglue1=5)+g(p_iglue2=6) 

      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j,k,h1,h2,h3,h4,iglue1,iglue2
      double precision p(mxpart,4),Asq,fac
      double precision Hgggg,Hqagg,Haqgg,Hgqgq,Hgaga,Hqgqg,Hagag,Hggqa
      double precision 
     . Hqrqr,Hqqqq,
     . Habab,Haaaa,
     . Hqarb,Hqaqa,Hqbqb,
     . Haqbr,Haqaq,Hbqbq
      double precision msq(-nf:nf,-nf:nf),hdecay,s34
      double complex amp(6,2,2,2,2)
      parameter(iglue1=5,iglue2=6)

C---fill spinor products upto maximum number
      call spinoru(iglue2,p,za,zb)  

C   Deal withb Higgs decay to b-bbar
      s34=s(3,4)+2d0*mb**2
      hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s34-4d0*mb**2) 
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      Asq=(as/(3d0*pi))**2/vevsq

C--four gluon terms
      call h4g(1,2,iglue1,iglue2,Hgggg)

C--two quark two gluon terms
      call hqqgg(1,2,iglue1,iglue2,Hqagg)
      call hqqgg(2,1,iglue1,iglue2,Haqgg)
      call hqqgg(1,iglue1,2,iglue2,Hqgqg)
      call hqqgg(iglue1,1,2,iglue2,Hagag)
      call hqqgg(2,iglue1,1,iglue2,Hgqgq)
      call hqqgg(iglue1,2,1,iglue2,Hgaga)
      call hqqgg(iglue2,iglue1,1,2,Hggqa)

C---four quark terms
      call H4q(1,2,iglue1,iglue2,Hqrqr,Hqqqq)
C---four anti-quark terms
      call H4q(iglue1,iglue2,1,2,Habab,Haaaa)

C-qqb
      call H4q(1,iglue2,2,iglue1,Hqarb,Hqaqa)
      call H4q(1,iglue2,iglue1,2,Hqbqb,Hqaqa)

C-qbq
      call H4q(2,iglue2,1,iglue1,Haqbr,Haqaq)
      call H4q(2,iglue2,iglue1,1,Hbqbq,Haqaq)

      fac=gsq**2*Asq*hdecay

      do j=fn,nf
      do k=fn,nf
      msq(j,k)=0d0
      if ((j.gt.0).and.(k.gt.0)) then 
      if (j.eq.k) then
      msq(j,k)=0.5d0*aveqq*fac*Hqqqq
      else
      msq(j,k)=aveqq*fac*Hqrqr
      endif
      endif
      if ((j.lt.0).and.(k.lt.0)) then 
      if (j.eq.k) then
      msq(j,k)=0.5d0*aveqq*fac*Haaaa
      else
      msq(j,k)=aveqq*fac*Habab
      endif
      endif

      if ((j.gt.0).and.(k.lt.0)) then
      if (j.eq.-k) then
      msq(j,k)=aveqq*fac*(0.5d0*Hqagg+Hqaqa+(nf-1)*Hqarb)
      else
      msq(j,k)=aveqq*fac*Hqbqb
      endif
      endif

      if ((j.lt.0).and.(k.gt.0)) then
      if (j.eq.-k) then
      msq(j,k)=aveqq*fac*(0.5d0*Haqgg+Haqaq+dfloat(nf-1)*Haqbr)
      else
      msq(j,k)=aveqq*fac*Hbqbq
      endif
      endif

      if ((j.gt.0).and.(k.eq.0)) msq(j,0)=aveqg*fac*Hqgqg
      if ((j.lt.0).and.(k.eq.0)) msq(j,0)=aveqg*fac*Hagag

      if ((j.eq.0).and.(k.gt.0)) msq(0,k)=aveqg*fac*Hgqgq
      if ((j.eq.0).and.(k.lt.0)) msq(0,k)=aveqg*fac*Hgaga

      if ((j.eq.0).and.(k.eq.0)) then
      msq(0,0)=avegg*fac*(0.5d0*Hgggg+dfloat(nf)*Hggqa)
      endif
      
      enddo
      enddo

      return
      end

 
