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
      integer j,k,jk,polq,pol1,pol2,pol3,jp,kp
      double precision fac,fac1
      common/pchoice/j,k
      double precision P(mxpart,4),qdks(mxpart,4),msq(-nf:nf,-nf:nf)
      double precision ave,v1(2),v2(2)
      double complex qqb(2,2,2,2),qbq(2,2,2,2)
      double complex qg(2,2,2,2),gq(2,2,2,2)
      double complex gqb(2,2,2,2),qbg(2,2,2,2)
      double complex cprop,propz1,propz2,amp
      double complex prop34,prop56

      do jp=-nf,nf
      do kp=-nf,nf
      msq(jp,kp)=0d0
      enddo
      enddo

      fac=-4D0*esq**2
      fac1=two*gsq*cf

      v1(1)=l1
      v1(2)=r1
      v2(1)=l2
      v2(2)=r2

C----Change the momenta to DKS notation 
c   We have --- q(p1)+qbar(p2)-->nu(p3)+e^+(p4)+b(p5)+bbar(p6)+g(p7)
c   DKS have--- q(q2)+qbar(q1)-->nu(q3)+e^+(q4)+b(q6)+bbar(q5)+g(q7)

      do jp=1,4
      qdks(1,jp)=p(1,jp)
      qdks(2,jp)=p(2,jp)
      qdks(3,jp)=p(3,jp)
      qdks(4,jp)=p(4,jp)
      qdks(5,jp)=p(6,jp)
      qdks(6,jp)=p(5,jp)
      qdks(7,jp)=p(7,jp)
      enddo

      call spinoru(7,qdks,za,zb)
c--   s returned from sprodx (common block) is 2*dot product

c--   calculate propagators
      if     (zerowidth  .eqv. .true.) then
      prop34=s(3,4)/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=s(5,6)/dcmplx(s(5,6)-zmass**2,zmass*zwidth)
      cprop=dcmplx(1d0)
      elseif (zerowidth .neqv. .true.) then
      prop34=dcmplx(s(3,4)/(s(3,4)-zmass**2))
      prop56=dcmplx(s(5,6)/(s(5,6)-zmass**2))
      propz1=(s(3,4)-zmass**2)/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      propz2=(s(5,6)-zmass**2)/dcmplx(s(5,6)-zmass**2,zmass*zwidth)
      cprop=propz1*propz2
      endif

c---case qbar-q
      call zzamps(1,2,3,4,5,6,7,za,zb,qbq)
c---case q-qbar
      call zzamps(2,1,3,4,5,6,7,za,zb,qqb)
c---case qbar-g
      call zzamps(1,7,3,4,5,6,2,za,zb,qbg)
c---case q-g
      call zzamps(7,1,3,4,5,6,2,za,zb,qg)
c---case g-q
      call zzamps(7,2,3,4,5,6,1,za,zb,gq)
c---case g-qbar
      call zzamps(2,7,3,4,5,6,1,za,zb,gqb)

      do j=-nf,nf
        do k=-nf,nf
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19
      jk=max(j,k)
        ave=xn*aveqq
        if (j .eq. 0 .or. k .eq. 0) then
          jk=j+k
        ave=xn*aveqg
        endif

c      if (abs(j) .le. 2 .and. j .ne. 0.and. k .ne. 0)
c     .write(*,*) 'MCFM, j=',j

        if (jk .eq. 0) goto 19
        
      do polq=1,2
      do pol1=1,2
      do pol2=1,2
      do pol3=1,2

      amp=0d0

      if    ((j .lt. 0).and.(k .gt. 0)) then
c---case qbar-q
        if (polq .eq. 1) then
          amp=(prop56*v2(pol1)*l(k)+q2*q(k))
     .       *(prop34*v1(pol2)*l(k)+q1*q(k))*qbq(polq,pol1,pol2,pol3)
        elseif (polq .eq. 2) then
          amp=(prop56*v2(pol1)*r(k)+q2*q(k))
     .       *(prop34*v1(pol2)*r(k)+q1*q(k))*qbq(polq,pol1,pol2,pol3)
        endif
      elseif((j .gt. 0).and.(k .lt. 0)) then
c---case q-qbar
        if (polq .eq. 1) then
          amp=(prop56*v2(pol1)*l(j)+q2*q(j))
     .       *(prop34*v1(pol2)*l(j)+q1*q(j))*qqb(polq,pol1,pol2,pol3)
        elseif (polq .eq. 2) then
          amp=(prop56*v2(pol1)*r(j)+q2*q(j))
     .       *(prop34*v1(pol2)*r(j)+q1*q(j))*qqb(polq,pol1,pol2,pol3)
        endif
      elseif((j .lt. 0).and.(k .eq. 0)) then
c---case qbar-g
        if (polq .eq. 1) then
          amp=(prop56*v2(pol1)*l(-jk)+q2*q(-jk))
     .      *(prop34*v1(pol2)*l(-jk)+q1*q(-jk))*qbg(polq,pol1,pol2,pol3)
        elseif (polq .eq. 2) then
          amp=(prop56*v2(pol1)*r(-jk)+q2*q(-jk))
     .      *(prop34*v1(pol2)*r(-jk)+q1*q(-jk))*qbg(polq,pol1,pol2,pol3)
        endif
      elseif((k .lt. 0).and.(j .eq. 0)) then
c---case g-qbar
        if (polq .eq. 1) then
          amp=(prop56*v2(pol1)*l(-jk)+q2*q(-jk))
     .      *(prop34*v1(pol2)*l(-jk)+q1*q(-jk))*gqb(polq,pol1,pol2,pol3)
        elseif (polq .eq. 2) then
          amp=(prop56*v2(pol1)*r(-jk)+q2*q(-jk))
     .      *(prop34*v1(pol2)*r(-jk)+q1*q(-jk))*gqb(polq,pol1,pol2,pol3)
        endif
      elseif((j .gt. 0).and.(k .eq. 0)) then
c---case q-g
        if (polq .eq. 1) then
          amp=(prop56*v2(pol1)*l(jk)+q2*q(jk))
     .       *(prop34*v1(pol2)*l(jk)+q1*q(jk))*qg(polq,pol1,pol2,pol3)
        elseif (polq .eq. 2) then
          amp=(prop56*v2(pol1)*r(jk)+q2*q(jk))
     .       *(prop34*v1(pol2)*r(jk)+q1*q(jk))*qg(polq,pol1,pol2,pol3)
        endif
      elseif((k .gt. 0).and.(j .eq. 0)) then
c---case g-q
        if (polq .eq. 1) then
          amp=(prop56*v2(pol1)*l(jk)+q2*q(jk))
     .       *(prop34*v1(pol2)*l(jk)+q1*q(jk))*gq(polq,pol1,pol2,pol3)
        elseif (polq .eq. 2) then
          amp=(prop56*v2(pol1)*r(jk)+q2*q(jk))
     .       *(prop34*v1(pol2)*r(jk)+q1*q(jk))*gq(polq,pol1,pol2,pol3)
        endif
      endif
      
C-- Inclusion of width a la Baur and Zeppenfeld
      amp=amp*FAC*cprop
      msq(j,k)=msq(j,k)+fac1*ave*abs(amp)**2

c      if (abs(j) .le. 2 .and. j .ne. 0.and. k .ne. 0) write(*,*)'(',
c     .polq,',',pol1,',',pol2,',',pol3,') = ',4d0*ave*fac1*abs(amp)**2

      enddo
      enddo
      enddo
      enddo

c      if ((msq(j,k).ne.0d0).and.(abs(s(1,7)).lt.10d0)
c     . .and.(abs(s(2,7)).lt.10d0)) then 
c      write(*,*) ' g(',j,',',k,') :',msq(j,k)
c      endif
        
   19 continue
      enddo
        enddo

      return
      end

      
      




