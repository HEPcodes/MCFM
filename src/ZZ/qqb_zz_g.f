      subroutine qqb_zz_g(P,msq)
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  q'(p4)+bar{q'}(p5) + n(p6)+ebar(p7)+ g(p3)
c   for the moment --- radiation only from initial line
      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zerowidth.f'
      include 'ewcharge.f'
      include 'srdiags.f'
      integer j,k,jk,hq,h34,h56,hg,jp,kp,h1,h2,h3,h4
      double precision fac,fac1,q34,q56,s127
      common/pchoice/j,k
      double precision P(mxpart,4),qdks(mxpart,4),msq(-nf:nf,-nf:nf)
      double precision ave,v34(2),v56(2)
      double complex zaa(2,2,2,2),z34(2,2,2,2),z56(2,2,2,2)
      double complex qqb(2,2,2,2),qbq(2,2,2,2)
      double complex qg(2,2,2,2),gq(2,2,2,2)
      double complex gqb(2,2,2,2),qbg(2,2,2,2)
      double complex aq12(2,2,2,2),aq34(2,2,2,2),aq56(2,2,2,2)
      double complex qa12(2,2,2,2),qa34(2,2,2,2),qa56(2,2,2,2)
      double complex gq12(2,2,2,2),gq34(2,2,2,2),gq56(2,2,2,2)
      double complex ag12(2,2,2,2),ag34(2,2,2,2),ag56(2,2,2,2)
      double complex ga12(2,2,2,2),ga34(2,2,2,2),ga56(2,2,2,2)
      double complex qg12(2,2,2,2),qg34(2,2,2,2),qg56(2,2,2,2)

      double complex cprop,propz1,propz2,propz3,amp
      double complex prop34,prop56,prop127

      do jp=-nf,nf
      do kp=-nf,nf
      msq(jp,kp)=0d0
      enddo
      enddo

      fac=-4d0*esq**2
      fac1=two*cf*xn*gsq

      v34(1)=l1
      v34(2)=r1
      q34=q1
      v56(1)=l2
      v56(2)=r2
      q56=q2

c--   s returned from sprodx (common block) is 2*dot product
      call spinoru(7,p,za,zb)

c--   calculate propagators
      if     (zerowidth  .eqv. .true.) then
      prop34=s(3,4)/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=s(5,6)/dcmplx(s(5,6)-zmass**2,zmass*zwidth)
      cprop=dcmplx(1d0)
c      s127=s(1,2)+s(1,7)+s(2,7)
c      prop127=s127/dcmplx(s127-zmass**2,zmass*zwidth)
      elseif (zerowidth .neqv. .true.) then
      s127=s(1,2)+s(1,7)+s(2,7)
      prop127=dcmplx(s127/(s127-zmass**2))
      prop34=dcmplx(s(3,4)/(s(3,4)-zmass**2))
      prop56=dcmplx(s(5,6)/(s(5,6)-zmass**2))
      propz1=(s(3,4)-zmass**2)/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      propz2=(s(5,6)-zmass**2)/dcmplx(s(5,6)-zmass**2,zmass*zwidth)
      propz3=(s127-zmass**2)/dcmplx(s127-zmass**2,zmass*zwidth)
      cprop=propz1*propz2*propz3
      endif

c--- Amplitude returned with arguments (hq,h34,h56,h7)
c---case qbar-q
      call zzgamp(1,2,3,4,5,6,7,za,zb,aq12,aq34,aq56)
c---case q-qbar
      call zzgamp(2,1,3,4,5,6,7,za,zb,qa12,qa34,qa56)
c---case qbar-g
      call zzgamp(1,7,3,4,5,6,2,za,zb,ag12,ag34,ag56)
c---case q-g
      call zzgamp(7,1,3,4,5,6,2,za,zb,qg12,qg34,qg56)
c---case g-q
      call zzgamp(7,2,3,4,5,6,1,za,zb,gq12,gq34,gq56)
c---case g-qbar
      call zzgamp(2,7,3,4,5,6,1,za,zb,ga12,ga34,ga56)


c--- Old code
cc---case qbar-q
c      call zzamps(1,2,3,4,5,6,7,za,zb,qbq)
cc---case q-qbar
c      call zzamps(2,1,3,4,5,6,7,za,zb,qqb)
cc---case qbar-g
c      call zzamps(1,7,3,4,5,6,2,za,zb,qbg)
cc---case q-g
c      call zzamps(7,1,3,4,5,6,2,za,zb,qg)
cc---case g-q
c      call zzamps(7,2,3,4,5,6,1,za,zb,gq)
cc---case g-qbar
c      call zzamps(2,7,3,4,5,6,1,za,zb,gqb)

c      do h1=1,2
c      do h2=1,2
c      do h3=1,2
c      do h4=1,2
cc      write(6,*) h1,h2,h3,h4,qbq(h1,h2,h3,h4)/aq12(h1,h2,h3,h4)
cc      write(6,*) h1,h2,h3,h4,qqb(h1,h2,h3,h4)/qa12(h1,h2,h3,h4)
cc      write(6,*) h1,h2,h3,h4,qg(h1,h2,h3,h4)/qg12(h1,h2,h3,h4)
cc      write(6,*) h1,h2,h3,h4,gq(h1,h2,h3,h4)/gq12(h1,h2,h3,h4)
cc      write(6,*) h1,h2,h3,h4,qbg(h1,h2,h3,h4)/ag12(h1,h2,h3,h4)
cc      write(6,*) h1,h2,h3,h4,gqb(h1,h2,h3,h4)/ga12(h1,h2,h3,h4)
c      enddo
c      enddo
c      enddo
c      enddo
cc      pause
       
      do j=-nf,nf
      do k=-nf,nf
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if ((j .eq. 0) .or. (k .eq. 0)) then
        jk=j+k
        ave=aveqg
      else
        jk=max(j,k)
        ave=aveqq
      endif

      if (jk .eq. 0) goto 19
       
      do hq=1,2
      do h34=1,2
      do h56=1,2
      do hg=1,2

      amp=0d0

c---case qbar-q
      if    ((j .lt. 0).and.(k .gt. 0)) then
      if (hq .eq. 1) then
      amp=(prop56*v56(h56)*l(k)+q56*q(k))
     &   *(prop34*v34(h34)*l(k)+q34*q(k))*aq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     &    +(prop56*v34(h34)*v56(h56)+q34*q56)
     &    *(prop127*v34(h34)*l(k)+q34*q(k))*aq56(hq,h34,h56,hg)
     &    +(prop34*v34(h34)*v56(h56)+q34*q56)
     &    *(prop127*v56(h56)*l(k)+q56*q(k))*aq34(hq,h34,h56,hg)
      endif
      elseif (hq .eq. 2) then
      amp=(prop56*v56(h56)*r(k)+q56*q(k))
     &   *(prop34*v34(h34)*r(k)+q34*q(k))*aq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     &   +(prop56*v34(h34)*v56(h56)+q34*q56)
     &   *(prop127*v34(h34)*r(k)+q34*q(k))*aq56(hq,h34,h56,hg)
     &   +(prop34*v34(h34)*v56(h56)+q34*q56)
     &   *(prop127*v56(h56)*r(k)+q56*q(k))*aq34(hq,h34,h56,hg)
      endif
      endif

c---case q-qbar
      elseif((j .gt. 0).and.(k .lt. 0)) then
      if (hq .eq. 1) then
      amp=(prop56*v56(h56)*l(j)+q56*q(j))
     &   *(prop34*v34(h34)*l(j)+q34*q(j))*qa12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)
     & *(prop127*v34(h34)*l(j)+q34*q(j))*qa56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)
     & *(prop127*v56(h56)*l(j)+q56*q(j))*qa34(hq,h34,h56,hg)
      endif
      elseif (hq .eq. 2) then
      amp=(prop56*v56(h56)*r(j)+q56*q(j))
     &   *(prop34*v34(h34)*r(j)+q34*q(j))*qa12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)
     & *(prop127*v34(h34)*r(j)+q34*q(j))*qa56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)
     & *(prop127*v56(h56)*r(j)+q56*q(j))*qa34(hq,h34,h56,hg)
      endif
      endif

c---case qbar-g
      elseif((j .lt. 0).and.(k .eq. 0)) then
      if (hq .eq. 1) then
      amp=(prop56*v56(h56)*l(-jk)+q56*q(-jk))
     &    *(prop34*v34(h34)*l(-jk)+q34*q(-jk))*ag12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     &+(prop56*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v34(h34)*l(-jk)+q34*q(-jk))*ag56(hq,h34,h56,hg)
     &+(prop34*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v56(h56)*l(-jk)+q56*q(-jk))*ag34(hq,h34,h56,hg)
      endif
      elseif (hq .eq. 2) then
      amp=(prop56*v56(h56)*r(-jk)+q56*q(-jk))
     &    *(prop34*v34(h34)*r(-jk)+q34*q(-jk))*ag12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     &+(prop56*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v34(h34)*r(-jk)+q34*q(-jk))*ag56(hq,h34,h56,hg)
     &+(prop34*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v56(h56)*r(-jk)+q56*q(-jk))*ag34(hq,h34,h56,hg)
      endif
      endif

c---case g-qbar
      elseif((k .lt. 0).and.(j .eq. 0)) then
      if (hq .eq. 1) then
      amp=(prop56*v56(h56)*l(-jk)+q56*q(-jk))
     &    *(prop34*v34(h34)*l(-jk)+q34*q(-jk))*ga12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     &+(prop56*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v34(h34)*l(-jk)+q34*q(-jk))*ga56(hq,h34,h56,hg)
     &+(prop34*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v56(h56)*l(-jk)+q56*q(-jk))*ga34(hq,h34,h56,hg)
      endif
      elseif (hq .eq. 2) then
      amp=(prop56*v56(h56)*r(-jk)+q56*q(-jk))
     &    *(prop34*v34(h34)*r(-jk)+q34*q(-jk))*ga12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     &+(prop56*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v34(h34)*r(-jk)+q34*q(-jk))*ga56(hq,h34,h56,hg)
     &+(prop34*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v56(h56)*r(-jk)+q56*q(-jk))*ga34(hq,h34,h56,hg)
      endif
      endif

c---case q-g
      elseif((j .gt. 0).and.(k .eq. 0)) then
      if (hq .eq. 1) then
      amp=(prop56*v56(h56)*l(+jk)+q56*q(+jk))
     &    *(prop34*v34(h34)*l(+jk)+q34*q(+jk))*qg12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     &+(prop56*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v34(h34)*l(+jk)+q34*q(+jk))*qg56(hq,h34,h56,hg)
     &+(prop34*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v56(h56)*l(+jk)+q56*q(+jk))*qg34(hq,h34,h56,hg)
      endif
      elseif (hq .eq. 2) then
      amp=(prop56*v56(h56)*r(+jk)+q56*q(+jk))
     &    *(prop34*v34(h34)*r(+jk)+q34*q(+jk))*qg12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     &+(prop56*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v34(h34)*r(+jk)+q34*q(+jk))*qg56(hq,h34,h56,hg)
     &+(prop34*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v56(h56)*r(+jk)+q56*q(+jk))*qg34(hq,h34,h56,hg)
      endif
      endif

      elseif((k .gt. 0).and.(j .eq. 0)) then
c---case g-q
      if (hq .eq. 1) then
      amp=(prop56*v56(h56)*l(+jk)+q56*q(+jk))
     &    *(prop34*v34(h34)*l(+jk)+q34*q(+jk))*gq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     &+(prop56*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v34(h34)*l(+jk)+q34*q(+jk))*gq56(hq,h34,h56,hg)
     &+(prop34*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v56(h56)*l(+jk)+q56*q(+jk))*gq34(hq,h34,h56,hg)
      endif
      elseif (hq .eq. 2) then
      amp=(prop56*v56(h56)*r(+jk)+q56*q(+jk))
     &    *(prop34*v34(h34)*r(+jk)+q34*q(+jk))*gq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     &+(prop56*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v34(h34)*r(+jk)+q34*q(+jk))*gq56(hq,h34,h56,hg)
     &+(prop34*v34(h34)*v56(h56)+q34*q56)
     &*(prop127*v56(h56)*r(+jk)+q56*q(+jk))*gq34(hq,h34,h56,hg)
      endif
      endif

      endif

C-- Inclusion of width a la Baur and Zeppenfeld
      amp=amp*fac*cprop
      msq(j,k)=msq(j,k)+fac1*ave*abs(amp)**2

      enddo
      enddo
      enddo
      enddo

      
   19 continue
      enddo
      enddo


      return
      end

      
      




