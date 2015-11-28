      subroutine qqb_wz(p,msq)
      implicit none
C----Author R.K.Ellis December 1998
C----modified by JMC to include supplementary diagrams February 1999
c----Matrix element for WZ production
c----in the notation of DKS
C----averaged over initial oolours and spins
C For nwz=+1
c     u(-p1)+dbar(-p2)-->mu^-(p5)+mu^+(p6)+n(p3)+e^+(p4)
C For nwz=-1
c     d(-p1)+ubar(-p2)-->mu^-(p5)+mu^+(p6)+e^-(p3)+nbar(p4)
c---
c     Notation to allow room for p3 --- gluon emission.
      include 'constants.f'
      include 'masses.f'
      include 'ckm.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'dprodx.f'
      include 'sprodx.f'
      include 'zerowidth.f'
      include 'ewcharge.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),qdks(mxpart,4)     
      double complex AWZM,AWZP,propw,propz,props,cprop,a6treea,a6treeb
      double complex prop34,prop56,prop12
      double complex Fa123456,Fa126543,Fb123456
      double complex Fa123465,Fa125643,Fb123465
      double complex Fa213456,Fa216543,Fb213456
      double complex Fa213465,Fa215643,Fb213465
      double complex Fa346512,Fa342156,Fa652143
      double complex Fa345612,Fa342165,Fa653421
      double complex Fa346521,Fa341256,Fa651243
      double complex Fa345621,Fa341265,Fa653412
c      double complex Fa562143,Fa561243
      double complex ZgL(-nf:nf),ZgR(-nf:nf)
      double precision v2(2),cl1,cl2,en1,en2
      double precision ave,cotw

      common/nwz/nwz

      double precision FAC,FACM
      integer j,k,nwz

      parameter(ave=0.25d0/xn)


      FAC=-two*gwsq*esq
      if ((nwz.eq.1) .or. (nwz .eq. -1)) then
      FACM=nwz*FAC
      else
      write(6,*) 'nwz .ne. +1 or -1'
      stop
      endif 
      if     (nwz.eq.-1) then
        cl1=1d0
        cl2=0d0
        en1=le
        en2=ln
      elseif (nwz.eq.+1) then
        cl1=0d0
        cl2=1d0
        en1=ln
        en2=le
      endif
      v2(1)=l1
      v2(2)=r1
      cotw=sqrt((one-xw)/xw)

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
      enddo
      enddo

C----Change the momenta to DKS notation 
c   We have --- d(-p1)+ubar(-p2)-->nu(p3)+e^+(p4)+mu^-(p5)+mu^+(p6)
c   DKS have--- u( q2)+dbar( q1)-->nu(q3)+e^+(q4)+mu^-(q6)+mu^+(q5)

      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(3,j)
      qdks(4,j)=p(4,j)
      qdks(5,j)=p(6,j)
      qdks(6,j)=p(5,j)
      enddo

      call spinoru(6,qdks,za,zb)
      
c--   s returned from sprod (common block) is 2*dot product
c--   calculate propagators

      if     (zerowidth  .eqv. .true.) then
      prop12=s(1,2)/(s(1,2)-wmass**2+im*wmass*wwidth)
      prop34=s(3,4)/(s(3,4)-wmass**2+im*wmass*wwidth)
      prop56=s(5,6)/(s(5,6)-zmass**2+im*zmass*zwidth)
      cprop=dcmplx(1d0)
      elseif (zerowidth .neqv. .true.) then
      prop12=dcmplx(s(1,2)/(s(1,2)-wmass**2))
      prop34=dcmplx(s(3,4)/(s(3,4)-wmass**2))
      prop56=dcmplx(s(5,6)/(s(5,6)-zmass**2))
      props=(s(1,2)-wmass**2)/(s(1,2)-wmass**2+im*wmass*wwidth)
      propw=(s(3,4)-wmass**2)/(s(3,4)-wmass**2+im*wmass*wwidth)
      propz=(s(5,6)-zmass**2)/(s(5,6)-zmass**2+im*zmass*zwidth)
      cprop=props*propw*propz
      endif

c---case dbar-u or ubar-d
      Fb123456=A6treeb(1,2,3,4,5,6,za,zb)
      Fa123456=A6treea(1,2,3,4,5,6,za,zb)
      Fa126543=A6treea(1,2,6,5,4,3,za,zb)

      Fb123465=A6treeb(1,2,3,4,6,5,za,zb)
      Fa123465=A6treea(1,2,3,4,6,5,za,zb)
      Fa125643=A6treea(1,2,5,6,4,3,za,zb)

c---case u-dbar or d-ubar
      Fb213456=A6treeb(2,1,3,4,5,6,za,zb)
      Fa213456=A6treea(2,1,3,4,5,6,za,zb)
      Fa216543=A6treea(2,1,6,5,4,3,za,zb)

      Fb213465=A6treeb(2,1,3,4,6,5,za,zb)
      Fa213465=A6treea(2,1,3,4,6,5,za,zb)
      Fa215643=A6treea(2,1,5,6,4,3,za,zb)

      if (zerowidth .neqv. .true.) then
c---for supplementary diagrams.
      Fa346512=A6treea(3,4,6,5,1,2,za,zb)
      Fa342156=A6treea(3,4,2,1,5,6,za,zb)
      Fa652143=A6treea(6,5,2,1,4,3,za,zb)
      Fa345612=A6treea(3,4,5,6,1,2,za,zb)
      Fa342165=A6treea(3,4,2,1,6,5,za,zb)
      Fa346521=A6treea(3,4,6,5,2,1,za,zb)
      Fa341256=A6treea(3,4,1,2,5,6,za,zb)
      Fa651243=A6treea(6,5,1,2,4,3,za,zb)
      Fa345621=A6treea(3,4,5,6,2,1,za,zb)
      Fa341265=A6treea(3,4,1,2,6,5,za,zb)
      Fa653412=A6treea(6,5,3,4,1,2,za,zb)
      Fa653421=A6treea(6,5,3,4,2,1,za,zb)
C----removed because apparently not needed.
c      Fa562143=A6treea(5,6,2,1,4,3,za,zb)
c      Fa561243=A6treea(5,6,1,2,4,3,za,zb)
      endif

c---set up left/right handed couplings for both Z and gamma*
c---note that L/R labels the LEPTON coupling v2, NOT the quarks (all L)
      do j=-nf,nf
        ZgL(j)=L(j)*v2(1)*prop56+Q(j)*q1           
        ZgR(j)=L(j)*v2(2)*prop56+Q(j)*q1           
      enddo
      
      do j=-nf,nf
      do k=-nf,nf
c--no point in wasting time if it gives zero anyway
      if (Vsq(j,k) .ne. 0d0) then
          if ((j .gt. 0) .and. (k .lt. 0)) then
            AWZM=(FAC*(ZgL(+j)*Fa213456+ZgL(-k)*Fa216543)
     .           +FACM*(v2(1)*cotw*prop56+q1)*prop12*Fb213456)*prop34
            AWZP=(FAC*(ZgR(+j)*Fa213465+ZgR(-k)*Fa215643)
     .           +FACM*(v2(2)*cotw*prop56+q1)*prop12*Fb213465)*prop34
          elseif ((j .lt. 0) .and. (k. gt. 0)) then
            AWZM=(FAC*(ZgL(+k)*Fa123456+ZgL(-j)*Fa126543)
     .           +FACM*(v2(1)*cotw*prop56+q1)*prop12*Fb123456)*prop34
            AWZP=(FAC*(ZgR(+k)*Fa123465+ZgR(-j)*Fa125643)
     .           +FACM*(v2(2)*cotw*prop56+q1)*prop12*Fb123465)*prop34
          endif
          if (zerowidth .neqv. .true.) then
c---we need supplementary diagrams for gauge invariance.
c---now also assume that we have lepton decay products for W
c---so that v2(1)=le, v2(2)=re
c---1st term is diagram where Z couples to electron
c---2nd term is diagram where Z couples to neutrino
c---3rd term is diagram where gamma* couples to electron
c---4th term (l-h only) contains two W propagators
          if ((j .gt. 0) .and. (k .lt. 0)) then
             AWZM=AWZM+FAC*prop12*(
     .          (en1*Fa346512+en2*Fa342156)*v2(1)*prop56
     .          +q1**2*(cl1*Fa346512+cl2*Fa342156)
     .          +0.5d0/xw*prop34*(cl1*Fa652143+cl2*Fa653412))
            AWZP=AWZP+FAC*prop12*(
     .          (en1*Fa345612+en2*Fa342165)*v2(2)*prop56
     .          +q1**2*(cl1*Fa345612+cl2*Fa342165))
          elseif ((j .lt. 0) .and. (k. gt. 0)) then
            AWZM=AWZM+FAC*prop12*(
     .          (en1*Fa346521+en2*Fa341256)*v2(1)*prop56
     .          +q1**2*(cl1*Fa346521+cl2*Fa341256)
     .          +0.5d0/xw*prop34*(cl1*Fa651243+cl2*Fa653421))
            AWZP=AWZP+FAC*prop12*(
     .          (en1*Fa345621+en2*Fa341265)*v2(2)*prop56
     .          +q1**2*(cl1*Fa345621+cl2*Fa341265))
          endif
          endif

C-- Inclusion of width for W,Z a la Baur and Zeppenfeld
      AWZM=cprop*AWZM
      AWZP=cprop*AWZP

      msq(j,k)=Vsq(j,k)*ave*(abs(AWZM)**2+abs(AWZP)**2)
c      if (abs(j) .le. 2 .and. j .ne. 0) then
c        write(*,*) 'MCFM, j=',j
c        write(*,*) 'Minus amplitude = ',abs(AWZM)**2/3d0
c        write(*,*) ' Plus amplitude = ',abs(AWZP)**2/3d0
c      endif

      endif

      enddo
      enddo
      return
      end
