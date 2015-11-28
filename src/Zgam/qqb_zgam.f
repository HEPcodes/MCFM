      subroutine qqb_zgam(p,msq)
      implicit none
C-----Author Keith Ellis
C-----October 2002
c----Matrix element for Z/gamma+gamma production
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))+gamma(p5)
c---
      include 'constants.f'
      include 'ewcouple.f'
      include 'zprods_decl.f'
      integer j,k,h12,h34,h5
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,qqb(2),qbq(2)
      double complex qbqamp(2,2,2,2),qqbamp(2,2,2,2)

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

      fac=aveqq*8d0*esq**3*xn

      call zamps(1,2,3,4,5,za,zb,qbqamp)
      call zamps(2,1,3,4,5,za,zb,qqbamp)

      do j=1,2
      qbq(j)=0d0
      qqb(j)=0d0
      do h12=1,2
      do h34=1,2
      do h5=1,2
      qbq(j)=qbq(j)+cdabs(qbqamp(j,h12,h34,h5))**2
      qqb(j)=qqb(j)+cdabs(qqbamp(j,h12,h34,h5))**2
      enddo
      enddo
      enddo
      qqb(j)=fac*qqb(j)
      qbq(j)=fac*qbq(j)
      enddo

      do j=-nf,nf
      k=-j
          if ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=0d0
          elseif ((j .gt. 0) .and. (k .lt. 0)) then
            msq(j,k)=qqb(jj(j))
          elseif ((j .lt. 0) .and. (k .gt. 0)) then
            msq(j,k)=qbq(kk(k))
          endif
      enddo

      return
      end


      subroutine zamps(p1,p2,p3,p4,p5,za,zb,qbqamp)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'zerowidth.f'
      integer p1,p2,p3,p4,p5,j
      double precision s12,s34
      double complex ai(2,2,2),af(2,2,2),prp34,prp12,qbqamp(2,2,2,2)
      s12=s(p1,p2)      
      s34=s(p3,p4)      

      if (zerowidth) then
C basic
      ai(1,1,2)=+za(p1,p3)**2*zb(p3,p4)/(za(p1,p5)*za(p2,p5)*s34)
C (3<-->4)
      ai(1,2,2)=+za(p1,p4)**2*zb(p4,p3)/(za(p1,p5)*za(p2,p5)*s34)
C Flip_2 (za<-->zb),(1<-->2),(3<-->4),
      ai(1,1,1)=+zb(p2,p4)**2*za(p4,p3)/(zb(p1,p5)*zb(p2,p5)*s34)
C Flip_2 (za<-->zb),(1<-->2)
      ai(1,2,1)=+zb(p2,p3)**2*za(p3,p4)/(zb(p1,p5)*zb(p2,p5)*s34)

      af(1,1,1)=czip
      af(2,1,1)=czip
      af(1,1,2)=czip
      af(1,2,1)=czip

      else
      ai(1,1,2)=+za(p1,p3)**2*zb(p3,p4)/(za(p1,p5)*za(p2,p5)*s34)
      ai(1,2,2)=+za(p1,p4)**2*zb(p4,p3)/(za(p1,p5)*za(p2,p5)*s34)
      ai(1,1,1)=+zb(p2,p4)**2*za(p4,p3)/(zb(p1,p5)*zb(p2,p5)*s34)
      ai(1,2,1)=+zb(p2,p3)**2*za(p3,p4)/(zb(p1,p5)*zb(p2,p5)*s34)

C basic
      af(1,1,2)=+za(p1,p3)**2*zb(p1,p2)/(za(p3,p5)*za(p4,p5)*s12)
C (1<-->2)
      af(2,1,2)=+za(p2,p3)**2*zb(p2,p1)/(za(p3,p5)*za(p4,p5)*s12)
C Flip_2 (za<-->zb),(1<-->2),(3<-->4)
      af(1,1,1)=+zb(p2,p4)**2*za(p2,p1)/(zb(p3,p5)*zb(p4,p5)*s12)
C Flip_2 (za<-->zb),(3<-->4)
      af(2,1,1)=+zb(p1,p4)**2*za(p1,p2)/(zb(p3,p5)*zb(p4,p5)*s12)

      endif

C     This is complex conjugation
      ai(2,2,1)=Dconjg(ai(1,1,2))
      ai(2,1,1)=Dconjg(ai(1,2,2))
      ai(2,1,2)=Dconjg(ai(1,2,1))
      ai(2,2,2)=Dconjg(ai(1,1,1))

      af(2,2,1)=Dconjg(af(1,1,2))
      af(1,2,1)=Dconjg(af(2,1,2))
      af(2,2,2)=Dconjg(af(1,1,1))
      af(1,2,2)=Dconjg(af(2,1,1))


c--   calculate propagator factors
      prp34=s34/Dcmplx((s34-zmass**2),zmass*zwidth)
      prp12=s12/Dcmplx((s12-zmass**2),zmass*zwidth)

      do j=1,2
      qbqamp(j,1,1,1)=Q(j)*(Q(j)*q1+L(j)*l1*prp34)*ai(1,1,1)
     .                 +q1*(Q(j)*q1+L(j)*l1*prp12)*af(1,1,1)
      qbqamp(j,1,1,2)=Q(j)*(Q(j)*q1+L(j)*l1*prp34)*ai(1,1,2)
     .                 +q1*(Q(j)*q1+L(j)*l1*prp12)*af(1,1,2)

      qbqamp(j,2,2,1)=Q(j)*(Q(j)*q1+R(j)*r1*prp34)*ai(2,2,1)
     .                 +q1*(Q(j)*q1+R(j)*r1*prp12)*af(2,2,1)
      qbqamp(j,2,2,2)=Q(j)*(Q(j)*q1+R(j)*r1*prp34)*ai(2,2,2)
     .                 +q1*(Q(j)*q1+R(j)*r1*prp12)*af(2,2,2)

      qbqamp(j,1,2,1)=Q(j)*(Q(j)*q1+L(j)*r1*prp34)*ai(1,2,1)
     .                 +q1*(Q(j)*q1+L(j)*r1*prp12)*af(1,2,1)
      qbqamp(j,1,2,2)=Q(j)*(Q(j)*q1+L(j)*r1*prp34)*ai(1,2,2)
     .                 +q1*(Q(j)*q1+L(j)*r1*prp12)*af(1,2,2)

      qbqamp(j,2,1,1)=Q(j)*(Q(j)*q1+R(j)*l1*prp34)*ai(2,1,1)
     .                 +q1*(Q(j)*q1+R(j)*l1*prp12)*af(2,1,1)
      qbqamp(j,2,1,2)=Q(j)*(Q(j)*q1+R(j)*l1*prp34)*ai(2,1,2)
     .                 +q1*(Q(j)*q1+R(j)*l1*prp12)*af(2,1,2)

      enddo
      return
      end

