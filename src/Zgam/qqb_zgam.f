      subroutine qqb_zgam(p,msq)
      implicit none
C-----Author Keith Ellis
C-----August 2002
c----Matrix element for Z/gamma+gamma production
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))+gamma(p5)
c---
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprodx.f'
      include 'dprodx.f'
      integer j,k,i1,i3,i5
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,s12,s34,
     . qqb(2),qbq(2)
      double complex prp12,prp34,ai(2,2,2),af(2,2,2),bi(2,2,2),bf(2,2,2)

      integer jj(-nf:nf),kk(-nf:nf)
      data jj/-1,-2,-1,-2,-1,0,1,2,1,2,1/
      data kk/-1,-2,-1,-2,-1,0,1,2,1,2,1/
      
c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call spinoru(5,p,za,zb)

      s34=s(3,4)
      s12=s(1,2)
      
      fac=aveqq*8d0*esq**3*xn
c      write(6,*) '4d0*pi/esq',4d0*pi/esq


c--   calculate propagator factors
      prp34=s34/Dcmplx((s34-zmass**2),zmass*zwidth)
      prp12=s12/Dcmplx((s12-zmass**2),zmass*zwidth)

      call zamps(1,2,3,4,5,za,zb,ai,af)
      call zamps(2,1,3,4,5,za,zb,bi,bf)

c      write(6,*) 'l1,r1',sqrt(esq)*l1,sqrt(esq)*r1
c      write(6,*) 'L(1),R(1)',sqrt(esq)*L(1),sqrt(esq)*R(1)
c      write(6,*) 'L(2),R(2)',sqrt(esq)*L(2),sqrt(esq)*R(2)
 
      do j=1,2
      qbq(j)=+cdabs(Q(j)*(Q(j)*q1+L(j)*l1*prp34)*ai(1,1,1)
     .               +q1*(Q(j)*q1+L(j)*l1*prp12)*af(1,1,1))**2
     .       +cdabs(Q(j)*(Q(j)*q1+L(j)*l1*prp34)*ai(1,1,2)
     .               +q1*(Q(j)*q1+L(j)*l1*prp12)*af(1,1,2))**2

     .       +cdabs(Q(j)*(Q(j)*q1+R(j)*r1*prp34)*ai(2,2,1)
     .               +q1*(Q(j)*q1+R(j)*r1*prp12)*af(2,2,1))**2
     .       +cdabs(Q(j)*(Q(j)*q1+R(j)*r1*prp34)*ai(2,2,2)
     .               +q1*(Q(j)*q1+R(j)*r1*prp12)*af(2,2,2))**2

     .       +cdabs(Q(j)*(Q(j)*q1+L(j)*r1*prp34)*ai(1,2,1)
     .               +q1*(Q(j)*q1+L(j)*r1*prp12)*af(1,2,1))**2
     .       +cdabs(Q(j)*(Q(j)*q1+L(j)*r1*prp34)*ai(1,2,2)
     .               +q1*(Q(j)*q1+L(j)*r1*prp12)*af(1,2,2))**2

     .       +cdabs(Q(j)*(Q(j)*q1+R(j)*l1*prp34)*ai(2,1,1)
     .               +q1*(Q(j)*q1+R(j)*l1*prp12)*af(2,1,1))**2
     .       +cdabs(Q(j)*(Q(j)*q1+R(j)*l1*prp34)*ai(2,1,2)
     .               +q1*(Q(j)*q1+R(j)*l1*prp12)*af(2,1,2))**2

      qqb(j)=+cdabs(Q(j)*(Q(j)*q1+L(j)*l1*prp34)*bi(1,1,1)
     .               +q1*(Q(j)*q1+L(j)*l1*prp12)*bf(1,1,1))**2
     .       +cdabs(Q(j)*(Q(j)*q1+L(j)*l1*prp34)*bi(1,1,2)
     .               +q1*(Q(j)*q1+L(j)*l1*prp12)*bf(1,1,2))**2

     .       +cdabs(Q(j)*(Q(j)*q1+R(j)*r1*prp34)*bi(2,2,1)
     .               +q1*(Q(j)*q1+R(j)*r1*prp12)*bf(2,2,1))**2
     .       +cdabs(Q(j)*(Q(j)*q1+R(j)*r1*prp34)*bi(2,2,2)
     .               +q1*(Q(j)*q1+R(j)*r1*prp12)*bf(2,2,2))**2

     .       +cdabs(Q(j)*(Q(j)*q1+L(j)*r1*prp34)*bi(1,2,1)
     .               +q1*(Q(j)*q1+L(j)*r1*prp12)*bf(1,2,1))**2
     .       +cdabs(Q(j)*(Q(j)*q1+L(j)*r1*prp34)*bi(1,2,2)
     .               +q1*(Q(j)*q1+L(j)*r1*prp12)*bf(1,2,2))**2

     .       +cdabs(Q(j)*(Q(j)*q1+R(j)*l1*prp34)*bi(2,1,1)
     .               +q1*(Q(j)*q1+R(j)*l1*prp12)*bf(2,1,1))**2
     .       +cdabs(Q(j)*(Q(j)*q1+R(j)*l1*prp34)*bi(2,1,2)
     .               +q1*(Q(j)*q1+R(j)*l1*prp12)*bf(2,1,2))**2
      qqb(j)=fac*qqb(j)
      qbq(j)=fac*qbq(j)
      enddo

      do j=-nf,nf
      k=-j
          if ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=0d0
          elseif ((j .gt. 0) .and. (k .lt. 0)) then
            msq(j,k)=qqb(jj(j))
          elseif ((j .lt. 0) .and. (k. gt. 0)) then
            msq(j,k)=qbq(kk(k))
          endif
      enddo

      return
      end


      subroutine zamps(p1,p2,p3,p4,p5,za,zb,ai,af)
      implicit none
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      integer p1,p2,p3,p4,p5
      double precision s12,s34
      double complex ai(2,2,2),af(2,2,2)
      s12=s(p1,p2)      
      s34=s(p3,p4)      
      ai(1,1,1)=-zb(p2,p4)**2*za(p3,p4)/(zb(p1,p5)*zb(p2,p5)*s34)
      af(1,1,1)=-zb(p2,p4)**2*za(p1,p2)/(zb(p3,p5)*zb(p4,p5)*s12)

c      ai(2,2,2)=+za(p2,p4)**2*zb(p3,p4)/(za(p1,p5)*za(p2,p5)*s34)
c      af(2,2,2)=+za(p2,p4)**2*zb(p1,p2)/(za(p3,p5)*za(p4,p5)*s12)
 
      ai(2,1,1)=+zb(p1,p4)**2*za(p3,p4)/(zb(p1,p5)*zb(p2,p5)*s34)
      af(2,1,1)=+zb(p1,p4)**2*za(p1,p2)/(zb(p3,p5)*zb(p4,p5)*s12)

c      ai(1,2,2)=-za(p1,p4)**2*zb(p3,p4)/(za(p1,p5)*za(p2,p5)*s34)
c      af(1,2,2)=-za(p1,p4)**2*zb(p1,p2)/(za(p3,p5)*za(p4,p5)*s12)

      ai(1,1,2)=+za(p1,p3)**2*zb(p3,p4)/(za(p1,p5)*za(p2,p5)*s34)
      af(1,1,2)=+za(p1,p3)**2*zb(p1,p2)/(za(p3,p5)*za(p4,p5)*s12)

c      ai(2,2,1)=-zb(p1,p3)**2*za(p3,p4)/(zb(p1,p5)*zb(p2,p5)*s34)
c      af(2,2,1)=-zb(p1,p3)**2*za(p1,p2)/(zb(p3,p5)*zb(p4,p5)*s12)


      ai(1,2,1)=+zb(p2,p3)**2*za(p3,p4)/(zb(p1,p5)*zb(p2,p5)*s34)
      af(1,2,1)=+zb(p2,p3)**2*za(p1,p2)/(zb(p3,p5)*zb(p4,p5)*s12)

c      ai(2,1,2)=-za(p2,p3)**2*zb(p3,p4)/(za(p1,p5)*za(p2,p5)*s34)
c      af(2,1,2)=-za(p2,p3)**2*zb(p1,p2)/(za(p3,p5)*za(p4,p5)*s12)

      ai(2,2,2)=Dconjg(ai(1,1,1))
      ai(1,2,2)=Dconjg(ai(2,1,1))
      ai(2,1,2)=Dconjg(ai(1,2,1))
      ai(2,2,1)=Dconjg(ai(1,1,2))

      af(2,2,2)=Dconjg(af(1,1,1))
      af(1,2,2)=Dconjg(af(2,1,1))
      af(2,1,2)=Dconjg(af(1,2,1))
      af(2,2,1)=Dconjg(af(1,1,2))

      return
      end

