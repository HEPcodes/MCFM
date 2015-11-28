      subroutine qqb_w_g_AMP(p,msq)
      implicit none
c----Matrix element for W production
C----averaged over initial colours and spins
C For nwz=+1
c     u(-p1)+dbar(-p2)--> W^+(n(p3)+e^+(p4)) + g(p5) 
C For nwz=-1
c     d(-p1)+ubar(-p2)--> W^-(e^-(p3)+nbar(p4)) + g(p5) 
c---
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'cutoff.f'
      include 'sprodx.f'
      include 'dprodx.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),prop,FAC
      double precision qqbWg,qbqWg,qgWq,qbgWqb,gqbWqb,gqWq
      double complex qqbWgp,qbqWgp,qgWqp,qbgWqbp,gqbWqbp,gqWqp
      double complex qqbWgm,qbqWgm,qgWqm,qbgWqbm,gqbWqbm,gqWqm

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      call spinoru(5,p,za,zb)

c---protect from soft and collinear singularities
      if  ((-s(1,5) .lt. cutoff) .or. (-s(2,5) .lt. cutoff)) return

      call Amps(1,2,3,4,5,za,zb,qqbWgp,qqbWgm)
      call Amps(5,2,3,4,1,za,zb,gqbWqbp,gqbWqbm)
      call Amps(1,5,3,4,2,za,zb,qgWqp,qgWqm)
      call Amps(2,1,3,4,5,za,zb,qbqWgp,qbqWgm)
      call Amps(5,1,3,4,2,za,zb,qbgWqbp,qbgWqbm)
      call Amps(2,5,3,4,1,za,zb,gqWqp,gqWqm)
      
      prop=((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
      FAC=CF*XN/four*gw**4*gsq/prop

      qqbWg=aveqq*fac*(abs(qqbWgp)**2+abs(qqbWgm)**2)
      qbqWg=aveqq*fac*(abs(qbqWgp)**2+abs(qbqWgm)**2)
      gqWq=aveqg*fac*(abs(gqWqp)**2+abs(gqWqm)**2)
      qgWq=aveqg*fac*(abs(qgWqp)**2+abs(qgWqm)**2)
      gqbWqb=aveqg*fac*(abs(gqbWqbp)**2+abs(gqbWqbm)**2)
      qbgWqb=aveqg*fac*(abs(qbgWqbp)**2+abs(qbgWqbm)**2)
      
      do j=-nf,nf
      do k=-nf,nf

      if     ((j .gt. 0) .and. (k .lt. 0)) then
          msq(j,k)=Vsq(j,k)*qqbWg
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=Vsq(j,k)*qbqWg
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &   (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWq
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqb
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWq
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqb
      endif

      enddo
      enddo

      return
      end
 
      subroutine Amps(j1,j2,j3,j4,j5,za,zb,Ampp,Ampm)
C---Amplitudes for the two polarizations of the process
c--- q(p1)+qbar(p2) --> l(p3)+a(p4)+g(p5) with s(3,4) removed
      implicit none 
      include 'constants.f'
      include 'sprodx.f'
      integer j1,j2,j3,j4,j5
      double complex Ampp,Ampm
      Ampp=twort2*za(j2,j3)*(zb(j4,j1)*za(j1,j2)+zb(j4,j5)*za(j5,j2))
     . /(za(j2,j5)*za(j1,j5))
      Ampm=twort2*zb(j4,j1)*(zb(j1,j2)*za(j2,j3)+zb(j1,j5)*za(j5,j3))
     . /(zb(j2,j5)*zb(j1,j5))
      return
      end




