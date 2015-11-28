      subroutine qqb_wz_g(P,msq)
C---Author John Campbell Fri Feb 19 11:06:08 CST 1999
C---Modified to include supplementary diagrams by JC on Feb 24
c---Matrix element squared averaged over initial colors and spins
c---  averaged(summed) over initial(final) colours and spins
C For nwz=+1
c     u(-p1)+dbar(-p2)-->mu^-(p5)+mu^+(p6)+n(p3)+e^+(p4)+g(p7)
C For nwz=-1
c     d(-p1)+ubar(-p2)-->mu^-(p5)+mu^+(p6)+e^-(p3)+nbar(p4)+g(p7)
c---
c   for the moment --- radiation only from initial line
      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'dprodx.f'
      include 'sprodx.f'
      include 'zerowidth.f'
      include 'ewcharge.f'
      integer j,k,polg,polz,minus,mplus,jp,kp,nwz
      common/nwz/nwz
      double precision FAC,FACM,FAC1
      double complex prop12,prop34,prop56
      common/pchoice/j,k
      double precision P(mxpart,4),qdks(mxpart,4),msq(-nf:nf,-nf:nf)
      double precision ave,cotw,s127
      double complex 
     .  qu_qb(9,2,2),qu_gg(9,2,2),gg_qb(9,2,2),
     .  qb_qu(9,2,2),qb_gg(9,2,2),gg_qu(9,2,2),
     .  props,propw,propz,cprop,A(2,2)
      double precision v2(2),cl1,cl2,en1,en2
      double complex ZgLR(nf,2),c1(2),c2(2)
      data minus,mplus/1,2/
c      double precision fudge

      FAC=-2D0*gwsq*esq
      FAC1=two*gsq*cf
      if ((nwz.eq.1) .or. (nwz .eq. -1)) then
      FACM=nwz*FAC
      else
      write(6,*) 'nwz .ne. +1 or -1 in qqb_zz_g.f'
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

      do jp=-nf,nf
      do kp=-nf,nf
      msq(jp,kp)=0d0
      enddo
      enddo

C----Change the momenta to DKS notation 
c   We have --- d(-p1)+ubar(-p2)-->nu(p3)+e^+(p4)+mu^-(p5)+mu^+(p6)+g(p7)
c   DKS have--- u( q2)+dbar( q1)-->nu(q3)+e^+(q4)+mu^-(q6)+mu^+(q5)+g(p7)

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
      cotw=sqrt((one-xw)/xw)
      s127=s(1,2)+s(1,7)+s(2,7)
      if     (zerowidth  .eqv. .true.) then
      prop12=s127/(s127-wmass**2+im*wmass*wwidth)  
      prop34=s(3,4)/(s(3,4)-wmass**2+im*wmass*wwidth)
      prop56=s(5,6)/(s(5,6)-zmass**2+im*zmass*zwidth)
      cprop=dcmplx(1d0)
      elseif (zerowidth .neqv. .true.) then      
      prop12=dcmplx(s127/(s127-wmass**2))
      prop34=dcmplx(s(3,4)/(s(3,4)-wmass**2))
      prop56=dcmplx(s(5,6)/(s(5,6)-zmass**2))
      props=(s127-wmass**2)/(s127-wmass**2+im*wmass*wwidth)
      propw=(s(3,4)-wmass**2)/(s(3,4)-wmass**2+im*wmass*wwidth)
      propz=(s(5,6)-zmass**2)/(s(5,6)-zmass**2+im*zmass*zwidth)
      cprop=props*propw*propz
      endif

c---case dbar-u
      call wzamps(1,2,3,4,5,6,7,za,zb,qb_qu)
c---case u-dbar
      call wzamps(2,1,3,4,5,6,7,za,zb,qu_qb)

c---case g-u
      call wzamps(7,2,3,4,5,6,1,za,zb,gg_qu)
c---case u-g
      call wzamps(7,1,3,4,5,6,2,za,zb,qu_gg)

c---case dbar-g
      call wzamps(1,7,3,4,5,6,2,za,zb,qb_gg)
c---case g-dbar
      call wzamps(2,7,3,4,5,6,1,za,zb,gg_qb)

c---set up left/right handed couplings for both Z and gamma*
c---note that the second label corresponds to the helicity
c---of the LEPTON coupling v2, NOT the quarks (all L)
      do j=1,nf
        ZgLR(j,minus)=L(j)*v2(1)*prop56+Q(j)*q1           
        ZgLR(j,mplus)=L(j)*v2(2)*prop56+Q(j)*q1           
      enddo
      
      do polz=1,2
      if(nwz.eq.1) then
        c1(polz)=ZgLR(2,polz)
        c2(polz)=ZgLR(1,polz)
      else
        c1(polz)=ZgLR(1,polz)
        c2(polz)=ZgLR(2,polz)
      endif
      enddo
      
      do j=-nf,nf
      if (((j.eq.+1).or.(j.eq.+3).or.(j.eq.+5)
     . .or.(j.eq.-2).or.(j.eq.-4)) .and. (nwz .eq. +1))
     . go to 20
      if (((j.eq.-1).or.(j.eq.-3).or.(j.eq.-5)
     . .or.(j.eq.+2).or.(j.eq.+4)) .and. (nwz .eq. -1))
     . go to 20
      do k=-nf,nf
      if (((k.eq.+1).or.(k.eq.+3).or.(k.eq.+5)
     . .or.(k.eq.-2).or.(k.eq.-4)) .and. (nwz .eq. +1))
     . go to 19
      if (((k.eq.-1).or.(k.eq.-3).or.(k.eq.-5)
     . .or.(k.eq.+2).or.(k.eq.+4)) .and. (nwz .eq. -1))
     . go to 19
    
      if     ((j .gt. 0) .and. (k .lt. 0)) then

c---case u-db
          do polg=1,2
          do polz=1,2
          A(polg,polz)=((ZgLR(+j,polz)*qu_qb(2,polg,polz)
     .                  +ZgLR(-k,polz)*qu_qb(3,polg,polz))*FAC
     .                  +(cotw*v2(polz)*prop56+q1)
     .                    *prop12*qu_qb(1,polg,polz)*FACM)*prop34
     .            +FAC*((en1*v2(polz)*prop56+q1**2*cl1)
     .                   *prop12*qu_qb(5,polg,polz)
     .                 +(en2*v2(polz)*prop56+q1**2*cl2)
     .                   *prop12*qu_qb(4,polg,polz)
     .                 +0.5d0*prop34*prop12/xw*qu_qb(8,polg,polz)*cl1
     .                 +0.5d0*prop34*prop12/xw*qu_qb(9,polg,polz)*cl2)
     
c          A(polg,polz)=((L(+j)*qu_qb(2,polg,polz)
c     .                  +L(-k)*qu_qb(3,polg,polz))*FAC
c     .                  +cotw*prop12*qu_qb(1,polg,polz)*FACM)
c     .                 *prop34*prop56*v2(polz)
          enddo
          enddo
          ave=xn*aveqq*Vsq(j,k)

      elseif ((j .lt. 0) .and. (k .gt. 0)) then


c---case db-u
          do polg=1,2
          do polz=1,2
          A(polg,polz)=((ZgLR(+k,polz)*qb_qu(2,polg,polz)
     .                  +ZgLR(-j,polz)*qb_qu(3,polg,polz))*FAC
     .                  +(cotw*v2(polz)*prop56+q1)
     .                    *prop12*qb_qu(1,polg,polz)*FACM)*prop34
     .            +FAC*((en1*v2(polz)*prop56+q1**2*cl1)
     .                   *prop12*qb_qu(5,polg,polz)
     .                 +(en2*v2(polz)*prop56+q1**2*cl2)
     .                   *prop12*qb_qu(4,polg,polz)
     .                 +0.5d0*prop34*prop12/xw*qb_qu(6,polg,polz)*cl1
     .                 +0.5d0*prop34*prop12/xw*qb_qu(7,polg,polz)*cl2)

c          A(polg,polz)=((L(+k)*qb_qu(2,polg,polz)
c     .                  +L(-j)*qb_qu(3,polg,polz))*FAC
c     .                  +cotw*prop12*qb_qu(1,polg,polz)*FACM)
c     .                 *prop34*prop56*v2(polz)
          enddo
          enddo
          ave=xn*aveqq*Vsq(j,k)

      elseif ((j .gt. 0) .and. (k .eq. 0)) then
c---case u-g
          do polg=1,2
          do polz=1,2
          A(polg,polz)=((c1(polz)*qu_gg(2,polg,polz)
     .                  +c2(polz)*qu_gg(3,polg,polz))*FAC
     .                  +(cotw*v2(polz)*prop56+q1)
     .                    *prop12*qu_gg(1,polg,polz)*FACM)*prop34
     .            +FAC*((en1*v2(polz)*prop56+q1**2*cl1)
     .                   *prop12*qu_gg(5,polg,polz)
     .                 +(en2*v2(polz)*prop56+q1**2*cl2)
     .                   *prop12*qu_gg(4,polg,polz)
     .                 +0.5d0*prop34*prop12/xw*qu_gg(8,polg,polz)*cl1
     .                 +0.5d0*prop34*prop12/xw*qu_gg(9,polg,polz)*cl2)

c          A(polg,polz)=((c1*qu_gg(2,polg,polz)
c     .                  +c2*qu_gg(3,polg,polz))*FAC
c     .                  +cotw*prop12*qu_gg(1,polg,polz)*FACM)
c     .                 *prop34*prop56*v2(polz)
          enddo
          enddo
          ave=xn*aveqg*Vsum(j)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
c---case db-g
          do polg=1,2
          do polz=1,2
          A(polg,polz)=((c1(polz)*qb_gg(2,polg,polz)
     .                  +c2(polz)*qb_gg(3,polg,polz))*FAC
     .                  +(cotw*v2(polz)*prop56+q1)
     .                    *prop12*qb_gg(1,polg,polz)*FACM)*prop34
     .            +FAC*((en1*v2(polz)*prop56+q1**2*cl1)
     .                   *prop12*qb_gg(5,polg,polz)
     .                 +(en2*v2(polz)*prop56+q1**2*cl2)
     .                   *prop12*qb_gg(4,polg,polz)
     .                 +0.5d0*prop34*prop12/xw*qb_gg(6,polg,polz)*cl1
     .                 +0.5d0*prop34*prop12/xw*qb_gg(7,polg,polz)*cl2)

c          A(polg,polz)=((c1*qb_gg(2,polg,polz)
c     .                  +c2*qb_gg(3,polg,polz))*FAC
c     .                  +cotw*prop12*qb_gg(1,polg,polz)*FACM)
c     .                 *prop34*prop56*v2(polz)
          enddo
          enddo

          ave=xn*aveqg*Vsum(j)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
c---case g-u
          do polg=1,2
          do polz=1,2
          A(polg,polz)=((c1(polz)*gg_qu(2,polg,polz)
     .                  +c2(polz)*gg_qu(3,polg,polz))*FAC
     .                  +(cotw*v2(polz)*prop56+q1)
     .                    *prop12*gg_qu(1,polg,polz)*FACM)*prop34
     .            +FAC*((en1*v2(polz)*prop56+q1**2*cl1)
     .                   *prop12*gg_qu(5,polg,polz)
     .                 +(en2*v2(polz)*prop56+q1**2*cl2)
     .                   *prop12*gg_qu(4,polg,polz)
     .                 +0.5d0*prop34*prop12/xw*gg_qu(6,polg,polz)*cl1
     .                 +0.5d0*prop34*prop12/xw*gg_qu(7,polg,polz)*cl2)

c          A(polg,polz)=((c1(polz)*gg_qu(2,polg,polz)
c     .                  +c2(polz)*gg_qu(3,polg,polz))*FAC
c     .                  +cotw*prop12*gg_qu(1,polg,polz)*FACM)
c     .                 *prop34*prop56*v2(polz)
          enddo
          enddo
          ave=xn*aveqg*Vsum(k)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
c---case g-db
          do polg=1,2
          do polz=1,2
          A(polg,polz)=((c1(polz)*gg_qb(2,polg,polz)
     .                  +c2(polz)*gg_qb(3,polg,polz))*FAC
     .                  +(cotw*v2(polz)*prop56+q1)
     .                    *prop12*gg_qb(1,polg,polz)*FACM)*prop34
     .            +FAC*((en1*v2(polz)*prop56+q1**2*cl1)
     .                   *prop12*gg_qb(5,polg,polz)
     .                 +(en2*v2(polz)*prop56+q1**2*cl2)
     .                   *prop12*gg_qb(4,polg,polz)
     .                 +0.5d0*prop34*prop12/xw*gg_qb(8,polg,polz)*cl1
     .                 +0.5d0*prop34*prop12/xw*gg_qb(9,polg,polz)*cl2)
     
c          A(polg,polz)=((c1*gg_qb(2,polg,polz)
c     .                  +c2*gg_qb(3,polg,polz))*FAC
c     .                  +cotw*prop12*gg_qb(1,polg,polz)*FACM)
c     .                 *prop34*prop56*v2(polz)

          enddo
          enddo

          ave=xn*aveqg*Vsum(k)

      else
          ave=0d0
      endif
      
      if (ave.gt.0d0) then
      msq(j,k)=FAC1*ave*cdabs(cprop)**2
     .          *(cdabs(A(mplus,minus))**2+cdabs(A(minus,minus))**2
     .           +cdabs(A(mplus,mplus))**2+cdabs(A(minus,mplus))**2)
      endif

c      if (abs(j) .le. 2 .and. abs(k) .le. 2 .and.Vsq(j,k).ne.0d0) then        
c        write(*,*) 'MCFM, j=',j,', k=',k
c        fudge=1.0000221d0/Vsq(j,k)
c      write(*,*) '(-,-) = ',4d0*fac1*ave*cdabs(A(1,1))**2*fudge
c      write(*,*) '(-,+) = ',4d0*fac1*ave*cdabs(A(1,2))**2*fudge
c      write(*,*) '(+,-) = ',4d0*fac1*ave*cdabs(A(2,1))**2*fudge
c      write(*,*) '(+,+) = ',4d0*fac1*ave*cdabs(A(2,2))**2*fudge
c      endif

 19   continue
      enddo
 20   continue
      enddo

      return
      end

      
      




