      subroutine gg_hgg_v(p,msq)
C-----
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
      include 'scale.f'
      include 'epinv.f'
      include 'first_time.f'
      integer j,k,al
      double precision p(mxpart,4),msq(fn:nf,fn:nf),s34,rts34
      double precision hdecay,Asq,fac,q(4,5)
      double precision Hqarbvsq,qrqr,qarb,aqrb,aqbr,abab,qbra,bqra,bqar
      double precision Hqaqavsq,qaqa,aqqa,aqaq,qqqq,aaaa
      double precision Hqaggvsq,
     . qagg,aqgg,qgqg,gqqg,agag,gaag,ggqa
      double precision Hggggvsq,gggg
      double precision Hqarbvsqnum
      double precision th,phi,ga,de
      double precision sth,sphi,sga,sde
      double precision cth,cphi,cga,cde
      double precision e1,e2,e3,e4,Hqarbsq,Hqaqasq
      common/GZmom/q
C*************************************************** 
C----Think about what scheme we are actually in?
      scheme='tH-V'
C*************************************************** 

C      write(6,*) 'Entering gg_hgg_v'
C      write(6,*) 'epinv',epinv
C      write(6,*) 'scale',scale
c      epinv=0d0
      new_event = .true. 

C-----debug

      call dotem(6,p,s)
      Asq=(as/(3d0*pi))**2/vevsq

C   Deal with Higgs decay to b-bbar
      s34=s(3,4)+2d0*mb**2
c      s34=one 
      rts34=sqrt(s34)
      hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s34-4d0*mb**2)
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)


C---Translate momenta into gz notation to fill common block
      do al=1,4
      q(al,1)=p(1,al)/scale
      q(al,2)=p(2,al)/scale
      q(al,3)=p(5,al)/scale
      q(al,4)=p(6,al)/scale
      q(al,5)=-(q(al,1)+q(al,2)+q(al,3)+q(al,4))
      enddo


      fac=ason2pi*Asq*gsq**2*hdecay

C-----Basic process is q(-k6)+r(-k2)-->q(-k5)+r(-k1)
C-----Note that Hqarbvsq(1,2,3,4)=Hqarbvsq(4,3,2,1)
C---quark-quark
C     q(1)+r(2)->q(5)+r(6)
      qrqr=Hqarbvsq(6,2,5,1)

C----quark-antiquark annihilation (6-->5-->2-->6) wrt q(1)+r(2)->q(5)+r(6)
c     q(1)+a(2)->r(5)+b(6)
      qarb=Hqarbvsq(5,6,2,1)

C----antiquark-quark annihilation (1<-->2) wrt to the above
c     a(1)+q(2)->r(5)+b(6)
c      aqrb=Hqarbvsq(5,6,1,2)

c--- changed by JMC, 3/8/05 to agree with the LO routine
C----antiquark-quark annihilation (1<-->2, 5<-->6) wrt to the above
c     a(1)+q(2)->b(5)+r(6)
c      aqbr=Hqarbvsq(6,5,1,2)
      aqbr=qarb
            
C----quark-antiquark scattering (6<-->2) wrt q(1)+r(2)->q(5)+r(6)
c     q(1)+b(2)->q(5)+b(6)
      qbra=Hqarbvsq(2,6,5,1)

C----antiquark-quark scattering
c     q(2)+a(1)->r(5)+b(6) (1<-->2) wrt to the above
c      bqra=Hqarbvsq(1,6,5,2) 
c--- changed by JMC, 3/8/05 to agree with the LO routine
      bqar=qbra
      
C---antiquark-antiquark scattering (1<-->5,2<-->6) wrt q(1)+r(2)->q(5)+r(6)
C     a(1)+b(2)->a(5)+b(6)
      abab=Hqarbvsq(2,6,1,5)

C************Identical quarks
C---quark-antiquark
C     q(1)+q(2)->q(5)+q(6)
C      qqqq=Hqarbvsq(6,2,5,1)+Hqarbvsq(5,2,6,1)+Hqaqavsq(6,2,5,1)
      qqqq=qrqr+Hqarbvsq(5,2,6,1)+Hqaqavsq(6,2,5,1)

C     a(1)+a(2)->a(5)+a(6) (1<-->5,2<-->6) wrt q(1)+q(2)->q(5)+q(6)
C      aaaa=Hqarbvsq(2,6,1,5)+Hqarbvsq(2,5,1,6)+Hqaqavsq(2,6,1,5)
      aaaa=abab+Hqarbvsq(2,5,1,6)+Hqaqavsq(2,6,1,5)

C     q(1)+a(2)->q(5)+a(6) (2<-->6) wrt q(1)+q(2)->q(5)+q(6)
      qaqa=qbra+qarb+Hqaqavsq(2,6,5,1)

C     a(1)+q(2)->q(5)+a(6) (1<-->6) wrt q(1)+q(2)->q(5)+q(6)
C      aqqa=Hqarbvsq(1,2,5,6)+Hqarbvsq(5,2,1,6)+Hqaqavsq(1,2,5,6)
c      aqqa=aqrb+qbra+Hqaqavsq(1,2,5,6)

c--- changed by JMC, 3/8/05 to agree with the LO routine
      aqaq=qaqa
      
C      write(6,*) 'qrqr',qrqr
C      write(6,*) 'qarb',qarb
C      write(6,*) 'aqrb',aqrb
C      write(6,*) 'abab',abab
C      write(6,*) 'qbra',qbra
C      write(6,*) 'bqra',bqra

C      write(6,*) 'Identical'
C      write(6,*) 'qaqa',qaqa
C      write(6,*) 'aqqa',aqqa
C      write(6,*) 'qqqq',qqqq
C      write(6,*) 'aaaa',aaaa


C************quark gluon processes
C     q(1)+a(2)->g(3)+g(4)
      qagg=+Hqaggvsq(1,2,3,4)
C      write(*,*) 'qagg',qagg
C     a(1)+q(2)->g(3)+g(4)
      aqgg=+Hqaggvsq(2,1,3,4)
C      write(*,*) 'aqgg',aqgg
C     q(1)+g(2)->q(3)+g(4)

c      qgqg=-Hqaggvsq(1,3,2,4)
c--- Sign changed by JMC, 8/5/05
      qgqg=+Hqaggvsq(1,3,2,4)

C     g(1)+q(2)->q(3)+g(4)
c      gqqg=-Hqaggvsq(2,3,1,4)
c--- Sign changed by JMC, 8/5/05
      gqqg=+Hqaggvsq(2,3,1,4)

C     a(1)+g(2)->a(3)+g(4)
c      agag=-Hqaggvsq(3,1,2,4)
c--- Sign changed by JMC, 8/5/05
      agag=+Hqaggvsq(3,1,2,4)

C     g(1)+a(2)->a(3)+g(4)
c      gaag=-Hqaggvsq(3,2,1,4)
c--- Sign changed by JMC, 8/5/05
      gaag=+Hqaggvsq(3,2,1,4)

C     g(1)+g(2)->q(3)+a(4)
      ggqa=+Hqaggvsq(4,3,1,2)


C************gluon gluon process

C     g(1)+g(2)->g(3)+g(4)
      gggg=+Hggggvsq(1,2,3,4)
C      write(*,*) 'gggg',gggg

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0

      if ((j.eq.0).and.(k.eq.0)) then
C---gg - all poles cancelled
             msq(j,k)=fac*avegg*(half*gggg+dfloat(nf)*ggqa)

      elseif ((j.gt.0).and.(k.gt.0)) then
C---qq - all poles cancelled
             if (j.eq.k) then
             msq(j,k)=aveqq*fac*half*qqqq
             else
             msq(j,k)=aveqq*fac*qrqr
             endif

      elseif ((j.lt.0).and.(k.lt.0)) then
C---aa - all poles cancelled
             if (j.eq.k) then
             msq(j,k)=aveqq*fac*half*aaaa
             else
             msq(j,k)=aveqq*fac*abab
             endif
C----qa scattering - all poles cancelled
      elseif ((j.gt.0).and.(k.lt.0)) then
         if (j.eq.-k) then
         msq(j,k)=aveqq*fac*(dfloat(nf-1)*qarb+half*(qaqa+qagg))
c--- changed by JMC 8/5/05 - no factor of half for qaqa
         msq(j,k)=aveqq*fac*(dfloat(nf-1)*qarb+qaqa+half*qagg)
             else
         msq(j,k)=aveqq*fac*qbra
         endif
C----aq scattering - all poles cancelled
      elseif ((j.lt.0).and.(k.gt.0)) then
         if (j.eq.-k) then
         msq(j,k)=aveqq*fac*(dfloat(nf-1)*aqbr+half*(aqaq+aqgg))
c--- changed by JMC 8/5/05 - no factor of half for aqaq
         msq(j,k)=aveqq*fac*(dfloat(nf-1)*aqbr+aqaq+half*aqgg)
             else
         msq(j,k)=aveqq*fac*bqar
         endif

C----gq scattering - all poles cancelled
      elseif ((j.eq.0).and.(k.gt.0)) then
         msq(j,k)=aveqg*fac*gqqg

C----ga scattering - all poles cancelled
      elseif ((j.eq.0).and.(k.lt.0)) then
         msq(j,k)=aveqg*fac*gaag

C----qg scattering - all poles cancelled
      elseif ((j.gt.0).and.(k.eq.0)) then
         msq(j,k)=aveqg*fac*qgqg

C----ag scattering - all poles cancelled
      elseif ((j.lt.0).and.(k.eq.0)) then
         msq(j,k)=aveqg*fac*agag
      endif

      enddo
      enddo

      return
      end

      double precision function Hggggvsq(j1,j2,j3,j4)
      implicit none
      include 'epinv.f'
      include 'epinv2.f'
      include 'first_time.f' 
C---  matrix element squared for H--g(j1)+g(j2)+g(j3)+g(j4)
      integer al,j1,j2,j3,j4
      double precision q(4,5),qswap(4,5),sqres(-2:0)
      common/GZmom/q
      
      do al=1,4
      qswap(al,1)=q(al,j1)
      qswap(al,2)=q(al,j2)
      qswap(al,3)=q(al,j3)
      qswap(al,4)=q(al,j4)
      qswap(al,5)=-q(al,j1)-q(al,j2)-q(al,j3)-q(al,j4)
      enddo 

      call GZHggggvsqPoles(qswap,sqres) 
      Hggggvsq=epinv*epinv2*sqres(-2)+epinv*sqres(-1)+sqres(0)
      !call GZHggggvsq(qswap,sqres,first_time,new_event) 
      !if (first_time) first_time = .false. 
      !if (new_event) new_event = .false. 
      return
      end


      double precision function Hqaggvsq(j1,j2,j3,j4)
      implicit none
      include 'epinv.f'
      include 'epinv2.f'
      include 'first_time.f'
C---  matrix element squared for H--q(j1)+a(j2)+g(j3)+g(j4)
      integer al,j1,j2,j3,j4
      double precision q(4,5),qswap(4,5),sqres(-2:0)
      common/GZmom/q

      do al=1,4
      qswap(al,1)=q(al,j1)
      qswap(al,2)=q(al,j2)
      qswap(al,3)=q(al,j3)
      qswap(al,4)=q(al,j4)
      qswap(al,5)=-q(al,j1)-q(al,j2)-q(al,j3)-q(al,j4)
      enddo

      call GZHqaggvsqPoles(qswap,sqres)
      !call GZHqaggvsq(qswap,sqres,first_time,new_event) 
      Hqaggvsq=epinv*epinv2*sqres(-2)+epinv*sqres(-1)+sqres(0)
      !if (first_time) first_time = .false. 
      !if (new_event) new_event = .false. 
      return
      end

C      double precision function Hqarbvsqnum(j1,j2,j3,j4)
C      implicit none
C      include 'epinv.f'
C      include 'first_time.f'
CC---  matrix element squared for H--q(j1)+a(j2)+g(j3)+g(j4)
C      integer al,j1,j2,j3,j4
C      double precision q(4,5),qswap(4,5),sqres(-2:0)
C      common/GZmom/q
C      do al=1,4
C      qswap(al,1)=q(al,j1)
C      qswap(al,2)=q(al,j2)
C      qswap(al,3)=q(al,j3)
C      qswap(al,4)=q(al,j4)
C      qswap(al,5)=-q(al,j1)-q(al,j2)-q(al,j3)-q(al,j4)
C      enddo
C      call GZHqarbvsq(qswap,sqres,first_time,new_event) 
C      Hqarbvsqnum=epinv**2*sqres(-2)+epinv*sqres(-1)+sqres(0)
C      if (first_time) first_time = .false. 
C      if (new_event) new_event = .false. 
C      return
C      end


