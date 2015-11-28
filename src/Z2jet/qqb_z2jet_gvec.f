      subroutine qqb_z2jet_gvec(p,n,in,msq)
      implicit none
C-----Nov 10 99 --- checked that gives the right matrix element
C     when summed over polarizations.

c----Matrix element for Z+2jet production
C----averaged over initial colours and spins
c    line 6 contracted with the vector n(mu)
c     q(-p1)+qbar(-p2)--> g(p5)+ g(p6)+Z(f(p3)+af(p4))

c---but this routine contains no factor of 1/2 for identical gluons
c   in the final state.

      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprodx.f'
      include 'dprodx.f'

C ip is the label of the emitter
C in is the label of the contracted line
      integer j,k,pq,pl,nquark,in
      double precision fac,prop,n(4)
      double complex zab(mxpart,mxpart),zba(mxpart,mxpart)
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),ggtemp,
     . qqbgg(2,2),qbqgg(2,2),qgqg(2,2),gqqg(2,2),
     . qbgqbg(2,2),gqbqbg(2,2),ggqbq(2,2)


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call spinoru(6,p,za,zb)
C   zab=<i-|k|j-> zba=<i+|k|j+> where k is an arbitrary 4-vector 
c---Conventions of Bern, Dixon, Kosower, Weinzierl, 
c---ie, za(i,j)*zb(j,i)=s(i,j)
      call spinork(6,p,zab,zba,n)


C---exclude the photon pole, 4*mbsq choosen as a scale approx above upsilon 
      if (s(3,4) .lt. 4d0*mbsq) return

      do pq=1,2
      do pl=1,2
      qqbgg(pq,pl) =0d0
      gqqg(pq,pl)  =0d0
      qgqg(pq,pl)  =0d0
      qbqgg(pq,pl) =0d0
      gqbqbg(pq,pl)=0d0
      qbgqbg(pq,pl)=0d0
      ggqbq(pq,pl) =0d0
      enddo      
      enddo      

      write(6,*) 'in in qqb_z2jet_gvec',in
      if (in .eq. 1) then
Cargument 1-4 represent (1) incoming quark line
C                       (2) incoming quark line
C                       (3) outgoing gluon line
C                       (4) outgoing gluon line contracted with n
           call z2jetsqn(5,6,2,1,p,n,za,zb,zab,zba,ggqbq)
           call z2jetsqn(2,5,6,1,p,n,za,zb,zab,zba,gqqg)
           call z2jetsqn(5,2,6,1,p,n,za,zb,zab,zba,gqbqbg)
      elseif (in .eq. 2)  then
           call z2jetsqn(5,6,1,2,p,n,za,zb,zab,zba,ggqbq)
           call z2jetsqn(5,1,6,2,p,n,za,zb,zab,zba,qbgqbg)
           call z2jetsqn(1,5,6,2,p,n,za,zb,zab,zba,qgqg)
      elseif (in .eq. 5) then
          call z2jetsqn(1,2,6,5,p,n,za,zb,zab,zba,qqbgg)
          call z2jetsqn(2,1,6,5,p,n,za,zb,zab,zba,qbqgg)
          call z2jetsqn(1,6,2,5,p,n,za,zb,zab,zba,qgqg)
          call z2jetsqn(2,6,1,5,p,n,za,zb,zab,zba,gqqg)
          call z2jetsqn(6,1,2,5,p,n,za,zb,zab,zba,qbgqbg)
          call z2jetsqn(6,2,1,5,p,n,za,zb,zab,zba,gqbqbg)
      elseif (in .eq. 6) then
          call z2jetsqn(1,2,5,6,p,n,za,zb,zab,zba,qqbgg)
          call z2jetsqn(2,1,5,6,p,n,za,zb,zab,zba,qbqgg)
          call z2jetsqn(1,5,2,6,p,n,za,zb,zab,zba,qgqg)
          call z2jetsqn(2,5,1,6,p,n,za,zb,zab,zba,gqqg)
          call z2jetsqn(5,1,2,6,p,n,za,zb,zab,zba,qbgqbg)
          call z2jetsqn(5,2,1,6,p,n,za,zb,zab,zba,gqbqbg)
      endif

      prop=s(3,4)/sqrt((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      fac=v*xn/four*(esq*gsq)**2*two

      do pq=1,2
      do pl=1,2
      qqbgg(pq,pl) =aveqq*fac*qqbgg(pq,pl)
      gqqg(pq,pl)  =aveqg*fac*gqqg(pq,pl)
      qgqg(pq,pl)  =aveqg*fac*qgqg(pq,pl)

      qbqgg(pq,pl) =aveqq*fac*qbqgg(pq,pl)
      gqbqbg(pq,pl)=aveqg*fac*gqbqbg(pq,pl)
      qbgqbg(pq,pl)=aveqg*fac*qbgqbg(pq,pl)

      ggqbq(pq,pl) =avegg*fac*ggqbq(pq,pl)
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
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=+(Q(k)*q1+L(k)*l1*prop)**2*qbqgg(1,1)
     .             +(Q(k)*q1+R(k)*r1*prop)**2*qbqgg(2,2)
     .             +(Q(k)*q1+L(k)*r1*prop)**2*qbqgg(1,2)
     .             +(Q(k)*q1+R(k)*l1*prop)**2*qbqgg(2,1)
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
 
          
    
      subroutine z2jetsqn(i1,i2,i5,i6,p,n,za,zb,zab,zba,msq)
      implicit none
C-----Apart from overall factors returns the matrix element squared
C-----msq dependent on the helicities pq and pl of the quark and
C-----lepton lines for 
C-----q(-p1)+qbar(-p2)-->l(p3)+al(p4)+g(p5)+g(p6) where 
C-----where gluon 6 has been contracted with the vector n
Cargument 1-4 represent (i1) incoming quark line
C                       (i2) incoming quark line
C                       (i5) outgoing gluon line
C                       (i6) outgoing gluon line contracted with n
      include 'constants.f'
      include 'sprodx.f'
      double complex qcdabn(2,2,2),qcdban(2,2,2),qedn(2,2,2)
      double complex zab(mxpart,mxpart),zba(mxpart,mxpart)
      double precision msq(2,2),n(4),p(mxpart,4),nDp5,nDp6
      integer i1,i2,i3,i4,i5,i6,pg,pq,pl
      i3=3
      i4=4
      nDp6=n(4)*p(i6,4)-n(3)*p(i6,3)-n(2)*p(i6,2)-n(1)*p(i6,1)
      if (abs(nDp6) .gt. 1d-12) then 
      write(6,*) 'i6,nDp6',i6,nDp6
      pause
      endif
      nDp5=n(4)*p(i5,4)-n(3)*p(i5,3)-n(2)*p(i5,2)-n(1)*p(i5,1)
      call subqcdn(i1,i2,i3,i4,i5,i6,nDp5,za,zb,zab,zba,qcdabn,qcdban)
            
C--first argument is gluon line
C--second argument is polarization of i5 line pq
C--third argument is polarization of lepton line pl
C  1=L,2=R
      do pq=1,2
      do pl=1,2
      msq(pq,pl)=0d0
      enddo
      enddo

      do pg=1,2
      do pq=1,2
      do pl=1,2
      qedn(pg,pq,pl)=qcdabn(pg,pq,pl)+qcdban(pg,pq,pl)
      msq(pq,pl)=msq(pq,pl)-ninth*abs(qedn(pg,pq,pl))**2
     & +abs(qcdabn(pg,pq,pl))**2+abs(qcdban(pg,pq,pl))**2
      enddo
      enddo
      enddo
      return
      end



