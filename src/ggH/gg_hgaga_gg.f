      subroutine gg_hgaga_gg(p,msq)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c
c     g(-p1)+g(-p2) -->  H(p3)+g(p_iglue1=4)+g(p_iglue2=5) 

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
      parameter(iglue1=4,iglue2=5)

C---fill spinor products upto maximum number
      call spinoru(iglue2,p,za,zb)  

C   Deal withb Higgs decay to b-bbar
c      s34=s(3,4)+2d0*mb**2
c      hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s34-4d0*mb**2) 
c      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      hdecay=1d0
      Asq=(as/(3d0*pi))**2/vevsq


ccc----zero amplitudes
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do j=1,6
      amp(j,h1,h2,h3,h4)=czip
      enddo
      enddo
      enddo
      enddo
      enddo



C--four gluon terms
      call makeall(1,2,iglue1,iglue2,amp)
      Hgggg=0d0
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      write(6,*) h1,h2,h3,h4,
     . +amp(1,h1,h2,h3,h4)
     . +amp(2,h1,h2,h3,h4)
     . +amp(3,h1,h2,h3,h4)
      write(6,*) h1,h2,h3,h4,
     . +amp(4,h1,h2,h3,h4)
     . +amp(5,h1,h2,h3,h4)
     . +amp(6,h1,h2,h3,h4)
      do j=1,6
      Hgggg=Hgggg+cdabs(amp(j,h1,h2,h3,h4))**2
      enddo
      enddo
      enddo
      enddo
      enddo
      Hgggg=xn**2*V*Hgggg/4d0



C--two quark two gluon terms
      call hqqgg(1,2,iglue1,iglue2,Hqagg)
      call hqqgg(2,1,iglue1,iglue2,Haqgg)
      call hqqgg(1,iglue1,2,iglue2,Hqgqg)
      call hqqgg(iglue1,1,2,iglue2,Hagag)
      call hqqgg(2,iglue2,1,iglue1,Hgqgq)
      call hqqgg(iglue2,2,1,iglue1,Hgaga)
      call hqqgg(iglue2,iglue1,1,2,Hggqa)

C---four quark terms
      call H4q(1,2,iglue1,iglue2,Hqrqr,Hqqqq)
C---four anti-quark terms
      call H4q(iglue1,iglue2,1,2,Habab,Haaaa)
c      write(6,*) 'Hqrqr',Hqrqr
c      write(6,*) 'Habab',Habab
c      write(6,*) 'Hqqqq',Hqqqq
c      write(6,*) 'Haaaa',Haaaa

C-qqb
      call H4q(1,iglue2,2,iglue1,Hqarb,Hqaqa)
      call H4q(1,iglue2,iglue1,2,Hqbqb,Hqaqa)
c      write(6,*) 'Hqaqa',Hqaqa


C-qbq
      call H4q(2,iglue1,1,iglue2,Haqbr,Haqaq)
      call H4q(2,iglue1,iglue2,1,Hbqbq,Haqaq)
c      write(6,*) 'Haqaq',Haqaq

c      write(6,*) 'Hqarb',Hqarb
c      write(6,*) 'Haqbr',Haqbr

c      write(6,*) 'Hqbqb',Hqbqb
c      write(6,*) 'Hbqbq',Hbqbq


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

      if ((j.eq.0).and.(k.eq.0)) msq(0,0)=avegg*fac*(0.5d0*Hgggg+Hggqa)
      enddo
      enddo

      return
      end

 
      subroutine makeall(p1,p2,p3,p4,amp)
      implicit none
      include 'constants.f'
      include 'zprods_com.f'
      integer j,p1,p2,p3,p4
      double complex amp(6,2,2,2,2),
     .  amppp(3),apmpp(3),appmp(3),apppm(3),
     .  apppp(3),
     .  ammpp(3),ampmp(3),amppm(3),apmmp(3),apmpm(3),appmm(3)

      call makepppp(p1,p2,p3,p4,za,apppp)
      call makemppp(p1,p2,p3,p4,za,zb,amppp,apmpp,appmp,apppm)
      call makemmpp(p1,p2,p3,p4,za,zb,
     . ammpp,ampmp,amppm,apmmp,apmpm,appmm)

      do j=1,3
      amp(j,  2,2,2,2)=apppp(j)
      amp(j+3,2,2,2,2)=apppp(j)
      amp(j,  1,2,2,2)=amppp(j)
      amp(j+3,1,2,2,2)=amppp(j)
      amp(j,  2,1,2,2)=apmpp(j)
      amp(j+3,2,1,2,2)=apmpp(j)
      amp(j,  2,2,1,2)=appmp(j)
      amp(j+3,2,2,1,2)=appmp(j)
      amp(j,  2,2,2,1)=apppm(j)
      amp(j+3,2,2,2,1)=apppm(j)

      amp(j,  1,1,2,2)=ammpp(j)
      amp(j+3,1,1,2,2)=ammpp(j)
      amp(j,  1,2,1,2)=ampmp(j)
      amp(j+3,1,2,1,2)=ampmp(j)
      amp(j,  1,2,2,1)=amppm(j)
      amp(j+3,1,2,2,1)=amppm(j)
      amp(j,  2,1,1,2)=apmmp(j)
      amp(j+3,2,1,1,2)=apmmp(j)
      amp(j,  2,1,2,1)=apmpm(j)
      amp(j+3,2,1,2,1)=apmpm(j)
      amp(j,  2,2,1,1)=appmm(j)
      amp(j+3,2,2,1,1)=appmm(j)
      enddo


      call makepppp(p1,p2,p3,p4,zb,apppp)
      call makemppp(p1,p2,p3,p4,zb,za,amppp,apmpp,appmp,apppm)
      call makemmpp(p1,p2,p3,p4,zb,za,
     . ammpp,ampmp,amppm,apmmp,apmpm,appmm)

      do j=1,3
      amp(j,  1,1,1,1)=apppp(j)
      amp(j+3,1,1,1,1)=apppp(j)
      amp(j,  2,1,1,1)=amppp(j)
      amp(j+3,2,1,1,1)=amppp(j)
      amp(j,  1,2,1,1)=apmpp(j)
      amp(j+3,1,2,1,1)=apmpp(j)
      amp(j,  1,1,2,1)=appmp(j)
      amp(j+3,1,1,2,1)=appmp(j)
      amp(j,  1,1,1,2)=apppm(j)
      amp(j+3,1,1,1,2)=apppm(j)

      amp(j,  1,1,2,2)=ammpp(j)
      amp(j+3,1,1,2,2)=ammpp(j)
      amp(j,  1,2,1,2)=ampmp(j)
      amp(j+3,1,2,1,2)=ampmp(j)
      amp(j,  1,2,2,1)=amppm(j)
      amp(j+3,1,2,2,1)=amppm(j)
      amp(j,  2,1,1,2)=apmmp(j)
      amp(j+3,2,1,1,2)=apmmp(j)
      amp(j,  2,1,2,1)=apmpm(j)
      amp(j+3,2,1,2,1)=apmpm(j)
      amp(j,  2,2,1,1)=appmm(j)
      amp(j+3,2,2,1,1)=appmm(j)
      enddo

      return
      end

      subroutine makepppp(p1,p2,p3,p4,za,apppp)
      implicit none
C     Taken from Kauffman hep-ph/9903330
C     and (older formula)
C     %\cite{Kauffman:1996ix}
C     \bibitem{Kauffman:1996ix}
C     R.~P.~Kauffman, S.~V.~Desai and D.~Risal,
C     %``Production of a Higgs boson plus two jets in hadronic collisions,''
C     Phys.\ Rev.\ D {\bf 55}, 4005 (1997)
C     [Erratum-ibid.\ D {\bf 58}, 119901 (1998)]
C     [arXiv:hep-ph/9610541].
C     %%CITATION = HEP-PH 9610541;%%
      include 'constants.f'
      include 'masses.f'
      include 'zprods_decl.f'
      integer j,p1,p2,p3,p4,i1(4),i2(4),i3(4),i4(4)
      double complex apppp(3)
      do j=1,3
      i1(j)=p1
            if (j.eq.1) then 
            i2(j)=p2
            i3(j)=p3
            i4(j)=p4
            elseif (j.eq.2) then
            i2(j)=p2
            i3(j)=p4
            i4(j)=p3
            elseif (j.eq.3) then
            i2(j)=p4
            i3(j)=p2
            i4(j)=p3
            endif
C---PRD55 Eq(21)
      apppp(j)=hmass**4/(za(i1(j),i2(j))*za(i2(j),i3(j))
     .  *za(i3(j),i4(j))*za(i4(j),i1(j)))
      enddo

      return
      end


      subroutine makemppp(p1,p2,p3,p4,za,zb,amppp,apmpp,appmp,apppm)
      implicit none
C     Taken from Kauffman hep-ph/9903330
C     and (older formula)
C     %\cite{Kauffman:1996ix}
C     \bibitem{Kauffman:1996ix}
C     R.~P.~Kauffman, S.~V.~Desai and D.~Risal,
C     %``Production of a Higgs boson plus two jets in hadronic collisions,''
C     Phys.\ Rev.\ D {\bf 55}, 4005 (1997)
C     [Erratum-ibid.\ D {\bf 58}, 119901 (1998)]
C     [arXiv:hep-ph/9610541].
C     %%CITATION = HEP-PH 9610541;%%

      integer j,k,p1,p2,p3,p4,j1,j2,jk(4),
     . i1(4),i2(4),i3(4),i4(4),k1(4),k2(4),k3(4),k4(4)
      double precision s123,s124,s134,s234
      double complex z2,amppp(3),apmpp(3),appmp(3),apppm(3),temp
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      data k1/1,2,3,4/
      data k2/2,3,4,1/
      data k3/3,4,1,2/
      data k4/4,1,2,3/
C---statement function
      z2(j1,j2)=-za(j1,p1)*zb(p1,j2)-za(j1,p2)*zb(p2,j2)
     .          -za(j1,p3)*zb(p3,j2)-za(j1,p4)*zb(p4,j2)
C---statement function
      jk(1)=p1
      jk(2)=p2
      jk(3)=p3
      jk(4)=p4
      do k=1,4
      do j=1,3
      i1(j)=jk(k1(k))
            if (j.eq.1) then 
            i2(j)=jk(k2(k))
            i3(j)=jk(k3(k))
            i4(j)=jk(k4(k))
      elseif (j.eq.2) then
            i2(j)=jk(k2(k))
            i3(j)=jk(k4(k))
            i4(j)=jk(k3(k))
      elseif (j.eq.3) then
            i2(j)=jk(k4(k))
            i3(j)=jk(k2(k))
            i4(j)=jk(k3(k))
      endif
      s124=s(i1(j),i2(j))+s(i1(j),i4(j))+s(i2(j),i4(j))
      s123=s(i1(j),i2(j))+s(i1(j),i3(j))+s(i2(j),i3(j))
      s134=s(i1(j),i3(j))+s(i1(j),i4(j))+s(i3(j),i4(j))
      s234=s(i2(j),i3(j))+s(i2(j),i4(j))+s(i3(j),i4(j))
C---PRD55 Eq(22)
c      amppp=
c     . -(z2(p1,p3)*zb(p2,p4))**2
c     . /((s(p1,p2)+s(p1,p4)+s(p2,p4))*s(p1,p2)*s(p1,p4)) 
c     . -(z2(p1,p4)*zb(p2,p3))**2/((s(p1,p2)+s(p1,p3)+s(p2,p3))*s(p1,p2)*s(p2,p3)) 
c     . -(z2(p1,p2)*zb(p3,p4))**2/((s(p1,p3)+s(p1,p4)+s(p3,p4))*s(p1,p4)*s(p3,p4)) 
c     . +zb(p2,p4)/(zb(p1,p2)*za(p2,p3)*za(p3,p4)*zb(p4,p1))
c     . *(+s(p2,p3)*z2(p1,p2)/za(p4,p1)
c     .   +s(p3,p4)*z2(p1,p4)/za(p1,p2)-zb(p2,p4)*(s(p2,p3)+s(p2,p4)+s(p3,p4)))
C---PRD55 Eq(A8+erratum)
c      amppp=
c     . -(z2(p1,p3)*zb(p2,p4))**2/((s(p1,p2)+s(p1,p4)+s(p2,p4))*s(p1,p2)*s(p1,p4)) 
c     . -(z2(p1,p4)*zb(p2,p3))**2/((s(p1,p2)+s(p1,p3)+s(p2,p3))*s(p1,p2)*s(p2,p3)) 
c     . -(z2(p1,p2)*zb(p3,p4))**2/((s(p1,p3)+s(p1,p4)+s(p3,p4))*s(p1,p4)*s(p3,p4)) 
c     . -zb(p2,p4)/(zb(p1,p2)*zb(p1,p4)*za(p1,p3))
c     . *(z2(p1,p2)**2/(za(p1,p4)*za(p3,p4))
c     .  +z2(p1,p4)**2/(za(p1,p2)*za(p2,p3)))

C---hep-ph/9903330 Eq(11)
      temp=
     . -(z2(i1(j),i3(j))*zb(i2(j),i4(j)))**2
     . /(s124*s(i1(j),i2(j))*s(i1(j),i4(j)))
     . -(z2(i1(j),i4(j))*zb(i2(j),i3(j)))**2
     . /(s123*s(i1(j),i2(j))*s(i2(j),i3(j))) 
     . -(z2(i1(j),i2(j))*zb(i3(j),i4(j)))**2
     . /(s134*s(i1(j),i4(j))*s(i3(j),i4(j)))
     . +zb(i2(j),i4(j))/
     . (zb(i1(j),i2(j))*za(i2(j),i3(j))*za(i3(j),i4(j))*zb(i4(j),i1(j)))
     . *(s(i2(j),i3(j))*z2(i1(j),i2(j))/za(i4(j),i1(j))
     . +s(i3(j),i4(j))*z2(i1(j),i4(j))/za(i1(j),i2(j))
     . -zb(i2(j),i4(j))*s234)
      if (k.eq.1) amppp(j)=temp
      if (k.eq.2) apmpp(j)=temp
      if (k.eq.3) appmp(j)=temp
      if (k.eq.4) apppm(j)=temp
      enddo
      enddo

      return
      end



      subroutine makemmpp(p1,p2,p3,p4,za,zb,
     . ammpp,ampmp,amppm,apmmp,apmpm,appmm)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
C     Taken from Kauffman hep-ph/9903330
C     and (older formula)
C     %\cite{Kauffman:1996ix}
C     \bibitem{Kauffman:1996ix}
C     R.~P.~Kauffman, S.~V.~Desai and D.~Risal,
C     %``Production of a Higgs boson plus two jets in hadronic collisions,''
C     Phys.\ Rev.\ D {\bf 55}, 4005 (1997)
C     [Erratum-ibid.\ D {\bf 58}, 119901 (1998)]
C     [arXiv:hep-ph/9610541].
C     %%CITATION = HEP-PH 9610541;%%
      integer j,k,p1,p2,p3,p4,jk(4),
     . i1(6),i2(6),i3(6),i4(6),k1(6),k2(6),k3(6),k4(6)
      double complex temp,
     . ammpp(3),ampmp(3),amppm(3),apmmp(3),apmpm(3),appmm(3)
      data k1/1,1,1,3,3,3/
      data k2/2,3,3,1,1,4/
      data k3/3,2,4,2,4,1/
      data k4/4,4,2,4,2,2/
      jk(1)=p1
      jk(2)=p2
      jk(3)=p3
      jk(4)=p4
      do k=1,6
            do j=1,3
            if (j.eq.1) then 
            i1(j)=jk(k1(k))
            i2(j)=jk(k2(k))
            i3(j)=jk(k3(k))
            i4(j)=jk(k4(k))
            elseif (j.eq.2) then
            i1(j)=jk(k2(k))
            i2(j)=jk(k1(k))
            i3(j)=jk(k3(k))
            i4(j)=jk(k4(k))
            elseif (j.eq.3) then
            i1(j)=jk(k2(k))
            i2(j)=jk(k3(k))
            i3(j)=jk(k1(k))
            i4(j)=jk(k4(k))
            endif
            temp=
     .      -za(jk(k1(k)),jk(k2(k)))**4/(za(i1(j),i2(j))*za(i2(j),i3(j))
     .      *za(i3(j),i4(j))*za(i4(j),i1(j)))
     .      -zb(jk(k3(k)),jk(k4(k)))**4/(zb(i1(j),i2(j))*zb(i2(j),i3(j))
     .      *zb(i3(j),i4(j))*zb(i4(j),i1(j)))
            if (k.eq.1) ammpp(j)=temp
            if (k.eq.2) ampmp(j)=temp
            if (k.eq.3) amppm(j)=temp
            if (k.eq.4) apmmp(j)=temp
            if (k.eq.5) apmpm(j)=temp
            if (k.eq.6) appmm(j)=temp
      enddo
      enddo

      return
      end
