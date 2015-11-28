      subroutine qqb_z2jet(p,msq)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) --> Z +g(p5) +g(p6)
c                          |
c                          --> l(p3)+a(p4)
c                            
c--all momenta incoming
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprodx.f'
      include 'dprodx.f'
      integer j,k,pq,pl,nquark,swap(2)
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),prop,fac,ggtemp,
     . qqbgg(2,2),qbqgg(2,2),qgqg(2,2),gqqg(2,2),
     . qbgqbg(2,2),gqbqbg(2,2),ggqbq(2,2)
      data swap/2,1/
      save swap
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo


      call spinorU(6,p,za,zb)

C---exclude the photon pole, 4*mbsq choosen as a scale approx above upsilon 
      if (s(3,4) .lt. 4d0*mbsq) return

      call z2jetsq(1,2,3,4,5,6,za,zb,qqbgg)
      call z2jetsq(1,5,3,4,2,6,za,zb,qgqg)
      call z2jetsq(2,5,3,4,1,6,za,zb,gqqg)
      call z2jetsq(5,6,3,4,1,2,za,zb,ggqbq)

c      call z2jetsq(2,1,3,4,5,6,za,zb,qbqgg)
c      call z2jetsq(5,1,3,4,2,6,za,zb,qbgqbg)
c      call z2jetsq(5,2,3,4,1,6,za,zb,gqbqbg)

      



      prop=s(3,4)/sqrt((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      fac=v*xn/four*(esq*gsq)**2

      do pq=1,2
      do pl=1,2
      qqbgg(pq,pl) =aveqq*fac*qqbgg(pq,pl)
      gqqg(pq,pl)  =aveqg*fac*gqqg(pq,pl)
      qgqg(pq,pl)  =aveqg*fac*qgqg(pq,pl)
c      qbqgg(pq,pl) =aveqq*fac*qbqgg(pq,pl)
c      gqbqbg(pq,pl)=aveqg*fac*gqbqbg(pq,pl)
c      qbgqbg(pq,pl)=aveqg*fac*qbgqbg(pq,pl)
      ggqbq(pq,pl) =avegg*fac*ggqbq(pq,pl)
      enddo      
      enddo      

      do pq=1,2
      do pl=1,2
      qbqgg(pq,pl)=qqbgg(swap(pq),pl)
      qbgqbg(pq,pl)=qgqg(swap(pq),pl)
      gqbqbg(pq,pl)=gqqg(swap(pq),pl)
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j .eq. 0) .and. (k .eq. 0)) then
           ggtemp=0d0
           do nquark=1,nf
           ggtemp=ggtemp
     .             +(Q(nquark)*q1+L(nquark)*l1*prop)**2*ggqbq(1,1)
     .             +(Q(nquark)*q1+R(nquark)*r1*prop)**2*ggqbq(2,2)
     .             +(Q(nquark)*q1+L(nquark)*r1*prop)**2*ggqbq(1,2)
     .             +(Q(nquark)*q1+R(nquark)*l1*prop)**2*ggqbq(2,1)
           enddo
          msq(j,k)=ggtemp
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          msq(j,k)=+(Q(j)*q1+L(j)*l1*prop)**2*qqbgg(1,1)
     .             +(Q(j)*q1+R(j)*r1*prop)**2*qqbgg(2,2)
     .             +(Q(j)*q1+L(j)*r1*prop)**2*qqbgg(1,2)
     .             +(Q(j)*q1+R(j)*l1*prop)**2*qqbgg(2,1)
c---Statistical factor
          msq(j,k)=half*msq(j,k)
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=+(Q(k)*q1+L(k)*l1*prop)**2*qbqgg(1,1)
     .             +(Q(k)*q1+R(k)*r1*prop)**2*qbqgg(2,2)
     .             +(Q(k)*q1+L(k)*r1*prop)**2*qbqgg(1,2)
     .             +(Q(k)*q1+R(k)*l1*prop)**2*qbqgg(2,1)
c---Statistical factor
          msq(j,k)=half*msq(j,k)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=+(Q(j)*q1+L(j)*l1*prop)**2*qgqg(1,1)
     .             +(Q(j)*q1+R(j)*r1*prop)**2*qgqg(2,2)
     .             +(Q(j)*q1+L(j)*r1*prop)**2*qgqg(1,2)
     .             +(Q(j)*q1+R(j)*l1*prop)**2*qgqg(2,1)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=+(Q(-j)*q1+L(-j)*l1*prop)**2*qbgqbg(1,1)
     .             +(Q(-j)*q1+R(-j)*r1*prop)**2*qbgqbg(2,2)
     .             +(Q(-j)*q1+L(-j)*r1*prop)**2*qbgqbg(1,2)
     .             +(Q(-j)*q1+R(-j)*l1*prop)**2*qbgqbg(2,1)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=+(Q(k)*q1+L(k)*l1*prop)**2*gqqg(1,1)
     .             +(Q(k)*q1+R(k)*r1*prop)**2*gqqg(2,2)
     .             +(Q(k)*q1+L(k)*r1*prop)**2*gqqg(1,2)
     .             +(Q(k)*q1+R(k)*l1*prop)**2*gqqg(2,1)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=+(Q(-k)*q1+L(-k)*l1*prop)**2*gqbqbg(1,1)
     .             +(Q(-k)*q1+R(-k)*r1*prop)**2*gqbqbg(2,2)
     .             +(Q(-k)*q1+L(-k)*r1*prop)**2*gqbqbg(1,2)
     .             +(Q(-k)*q1+R(-k)*l1*prop)**2*gqbqbg(2,1)
      endif

   19 continue
      enddo
      enddo

      return
      end
 
          
    
