      subroutine qqb_zbb_g(p,msq)
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  g*  + Z + g(p7)
c                           |     |
c                           |     --> e-(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+bb(p6)
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'prods.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'hardscale.f'
      integer j,k,nu,hq,Qh,hg,lh
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf),mmsq(2,2)
      double precision fac,pswap(mxpart,4),LRq(2),LRb(2),lr1(2)
      double precision msq_qqb(2,2,2,2,4),msq_qbq(2,2,2,2,4),
     .                 msq_qg(2,2,2,2,4),msq_qbg(2,2,2,2,4),
     .                 msq_gqb(2,2,2,2,4),msq_gq(2,2,2,2,4)
      double complex prop,czq,czb

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(1,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(7,nu)
      pswap(5,nu)=p(6,nu)
      pswap(6,nu)=p(3,nu)
      pswap(7,nu)=p(4,nu)
      enddo
      call spinoru(7,pswap,za,zb)
      call xzqqggg(mmsq)
      
C---call spinor routine and load common block twopij
      call spinoru(7,p,za,zb)
      prop=s(3,4)/(s(3,4)-zmass**2+im*zmass*zwidth)

      if (
     .      (s(5,6) .lt. four*hscalesq) 
     . .or. (s(1,5)*s(2,5)/s(1,2) .lt. hscalesq) 
     . .or. (s(1,6)*s(2,6)/s(1,2) .lt. hscalesq) ) return

c--- note that (2,1) and (4,3) are switched due to crossing from NT
      call msq_qqQQg(2,1,5,6,7,4,3,msq_qqb)
      call msq_qqQQg(1,2,5,6,7,4,3,msq_qbq)
      call msq_qqQQg(7,1,5,6,2,4,3,msq_qg)
      call msq_qqQQg(2,7,5,6,1,4,3,msq_gqb)
      call msq_qqQQg(7,2,5,6,1,4,3,msq_gq)
      call msq_qqQQg(1,7,5,6,2,4,3,msq_qbg)

c--- note the factor of 4d0*xw**2 relative to wbb
      fac=4d0*gsq**3*esq**2
c--- extra factor of 2**3=8 to compensate for Ta normalization
      fac=fac*8d0
       
      LRb(1)=L(1)
      LRb(2)=R(1)

      do j=-nf,nf
      do k=-nf,nf
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19
      
      msq(j,k)=0d0

      if     ((j .eq. 0) .and. (k .eq. 0)) then
C---L(1),R(1) for coupling to down quark
      msq(j,k)=abs(Q(1)*q1+prop*L(1)*l1)**2*mmsq(1,1)
     .        +abs(Q(1)*q1+prop*R(1)*l1)**2*mmsq(2,1)
     .        +abs(Q(1)*q1+prop*L(1)*r1)**2*mmsq(1,2)
     .        +abs(Q(1)*q1+prop*R(1)*r1)**2*mmsq(2,2)
c----debug
c          msq(j,k)=prop**2*(mmsq(1,1)+mmsq(2,1)+mmsq(1,2)+mmsq(2,2))
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
        LRq(1)=L(j)
        LRq(2)=R(j)
        lr1(1)=l1
        lr1(2)=r1
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
c--- couplings of Z to q (12) and Z to b (56)
c--- the 2nd of these is conjugated for convenience in interference
          czq=Q(j)*q1+prop*LRq(hq)*lr1(lh)
          czb=Q(1)*q1+Dconjg(prop)*LRb(Qh)*lr1(lh)
          msq(j,k)=msq(j,k)+aveqq*fac*(
     .      cdabs(czq)**2*msq_qqb(hq,Qh,hg,lh,1)
     .     +cdabs(czb)**2*msq_qqb(hq,Qh,hg,lh,2)
     .     +dreal(czq*czb*dcmplx(msq_qqb(hq,Qh,hg,lh,3),
     .                           msq_qqb(hq,Qh,hg,lh,4))))
        enddo
        enddo
        enddo
        enddo
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
        LRq(1)=L(k)
        LRq(2)=R(k)
        lr1(1)=l1
        lr1(2)=r1
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
          czq=Q(j)*q1+prop*LRq(hq)*lr1(lh)
          czb=Q(1)*q1+Dconjg(prop)*LRb(Qh)*lr1(lh)
          msq(j,k)=msq(j,k)+aveqq*fac*(
     .      cdabs(czq)**2*msq_qbq(hq,Qh,hg,lh,1)
     .     +cdabs(czb)**2*msq_qbq(hq,Qh,hg,lh,2)
     .     +dreal(czq*czb)*msq_qbq(hq,Qh,hg,lh,3))
        enddo
        enddo
        enddo
        enddo
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
        LRq(1)=L(j)
        LRq(2)=R(j)
        lr1(1)=l1
        lr1(2)=r1
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
          czq=Q(j)*q1+prop*LRq(hq)*lr1(lh)
          czb=Q(1)*q1+Dconjg(prop)*LRb(Qh)*lr1(lh)
          msq(j,k)=msq(j,k)+aveqg*fac*(
     .      cdabs(czq)**2*msq_qg(hq,Qh,hg,lh,1)
     .     +cdabs(czb)**2*msq_qg(hq,Qh,hg,lh,2)
     .     +dreal(czq*czb*dcmplx(msq_qg(hq,Qh,hg,lh,3),
     .                           msq_qg(hq,Qh,hg,lh,4))))
        enddo
        enddo
        enddo
        enddo
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
        LRq(1)=L(-j)
        LRq(2)=R(-j)
        lr1(1)=l1
        lr1(2)=r1
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
          czq=Q(j)*q1+prop*LRq(hq)*lr1(lh)
          czb=Q(1)*q1+Dconjg(prop)*LRb(Qh)*lr1(lh)
          msq(j,k)=msq(j,k)+aveqg*fac*(
     .      cdabs(czq)**2*msq_qbg(hq,Qh,hg,lh,1)
     .     +cdabs(czb)**2*msq_qbg(hq,Qh,hg,lh,2)
     .     +dreal(czq*czb*dcmplx(msq_qbg(hq,Qh,hg,lh,3),
     .                           msq_qbg(hq,Qh,hg,lh,4))))
        enddo
        enddo
        enddo
        enddo
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
        LRq(1)=L(k)
        LRq(2)=R(k)
        lr1(1)=l1
        lr1(2)=r1
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
          czq=Q(j)*q1+prop*LRq(hq)*lr1(lh)
          czb=Q(1)*q1+Dconjg(prop)*LRb(Qh)*lr1(lh)
          msq(j,k)=msq(j,k)+aveqg*fac*(
     .      cdabs(czq)**2*msq_gq(hq,Qh,hg,lh,1)
     .     +cdabs(czb)**2*msq_gq(hq,Qh,hg,lh,2)
     .     +dreal(czq*czb*dcmplx(msq_gq(hq,Qh,hg,lh,3),
     .                           msq_gq(hq,Qh,hg,lh,4))))
        enddo
        enddo
        enddo
        enddo
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
        LRq(1)=L(-k)
        LRq(2)=R(-k)
        lr1(1)=l1
        lr1(2)=r1
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
          czq=Q(j)*q1+prop*LRq(hq)*lr1(lh)
          czb=Q(1)*q1+Dconjg(prop)*LRb(Qh)*lr1(lh)
          msq(j,k)=msq(j,k)+aveqg*fac*(
     .      cdabs(czq)**2*msq_gqb(hq,Qh,hg,lh,1)
     .     +cdabs(czb)**2*msq_gqb(hq,Qh,hg,lh,2)
     .     +dreal(czq*czb*dcmplx(msq_gqb(hq,Qh,hg,lh,3),
     .                           msq_gqb(hq,Qh,hg,lh,4))))
        enddo
        enddo
        enddo
        enddo
      endif

   19 continue
      enddo
      enddo
      return
      end


      subroutine msq_qqQQg(i1,i2,i3,i4,i5,i6,i7,xsq)
c--- this will compute Z -> q qb Q Qb g amplitudes according to NT
c--- note that we simplify the situation by assuming that q != Q
c--- in which case all (1<-->3) and (2<-->4) amps vanish
c--- notation:  quark pairs (q,qb) are (1,2) and (3,4)
c---            lepton-antilepton (6,7) and gluon 5
c--- note that this is now updated so that the routine returns
c--- the matrix element separated according to coupling
c--- returned is msq(hq,Qh,hg,lh,j) where
c---   j=1 : Z couples to the (1,2) line in M and M*
c---   j=2 : Z couples to the (3,4) line in M and M*
c---   j=3 : Z couples to the (1,2) line in M and the (3,4) in M* (RE)
c---   j=4 : Z couples to the (1,2) line in M and the (3,4) in M* (IM)
      include 'constants.f'
      include 'qcdcouple.f'
      integer i1,i2,i3,i4,i5,i6,i7,j,hq,Qh,hg,lh
      double complex a1(2,2,2,2),a2(2,2,2,2),a3(2,2,2,2),a4(2,2,2,2),
     .               b1(2,2,2,2),b2(2,2,2,2),b3(2,2,2,2),b4(2,2,2,2)
      double complex mbar1,mbar2,mbar3,mbar4,temp
      double precision m0,mx,my,mz,mxx,mxy
      double precision mA,mB,mC,mD,mE,mF,mG,msq(2,2,2,2)
      double precision xA(4),xB(4),xC(4),xD(4),xE(4),xF(4),xG(4)
      double precision xy(4),xz(4),xxy(4),xsq(2,2,2,2,4)
      
c--- note that 'a' corresponds to the perm 1,2,3,4,5
      call nagyqqQQg(i1,i2,i3,i4,i5,i6,i7,a1,a2,a3,a4)
c---           'b' corresponds to the perm 3,4,1,2,5
      call nagyqqQQg(i3,i4,i1,i2,i5,i6,i7,b1,b2,b3,b4)
            
      do hq=1,2
      do Qh=1,2
      do hg=1,2
      do lh=1,2

c--- eqn (B57)
c        mbar1=a1(hq,Qh,hg,lh)+b3(Qh,hq,hg,lh)
c        mbar2=a2(hq,Qh,hg,lh)+b4(Qh,hq,hg,lh)  
c--- introduce new notation for eqn (B57) with (1<-->3, 2<-->4)
c--- mbar3 = m1(3,4,1,2) + m3(1,2,3,4)
c--- mbar4 = m2(3,4,1,2) + m4(1,2,3,4)
c        mbar3=b1(Qh,hq,hg,lh)+a3(hq,Qh,hg,lh)
c        mbar4=b2(Qh,hq,hg,lh)+a4(hq,Qh,hg,lh)
c--- eqn (B50)
c        mA=cdabs(mbar1)**2+cdabs(mbar2)**2
c     .    +cdabs(mbar3)**2+cdabs(mbar4)**2
        xA(1)=cdabs(a1(hq,Qh,hg,lh))**2+cdabs(a2(hq,Qh,hg,lh))**2
     .       +cdabs(a3(hq,Qh,hg,lh))**2+cdabs(a4(hq,Qh,hg,lh))**2
        xA(2)=cdabs(b1(Qh,hq,hg,lh))**2+cdabs(b2(Qh,hq,hg,lh))**2
     .       +cdabs(b3(Qh,hq,hg,lh))**2+cdabs(b4(Qh,hq,hg,lh))**2
        temp =2d0*(a1(hq,Qh,hg,lh)*Dconjg(b3(Qh,hq,hg,lh))
     .            +a2(hq,Qh,hg,lh)*Dconjg(b4(Qh,hq,hg,lh))
     .            +a3(hq,Qh,hg,lh)*Dconjg(b1(Qh,hq,hg,lh))
     .            +a4(hq,Qh,hg,lh)*Dconjg(b2(Qh,hq,hg,lh)))
        xA(3)=dreal(temp)
        xA(4)=dimag(temp)
c--- eqn (B51) - vanishes because different quarks
c        mB=0d0
c--- eqn (B52) - vanishes because different quarks
c        mC=0d0
c--- eqn (B53)
c        mD=2d0*dreal(mbar1*Dconjg(mbar2)+mbar3*Dconjg(mbar4))
        xD(1)=2d0*dreal(a1(hq,Qh,hg,lh)*Dconjg(a2(hq,Qh,hg,lh))
     .                 +a3(hq,Qh,hg,lh)*Dconjg(a4(hq,Qh,hg,lh)))
        xD(2)=2d0*dreal(b1(Qh,hq,hg,lh)*Dconjg(b2(Qh,hq,hg,lh))
     .                 +b3(Qh,hq,hg,lh)*Dconjg(b4(Qh,hq,hg,lh)))
        temp =2d0*(a1(hq,Qh,hg,lh)*Dconjg(b4(Qh,hq,hg,lh))
     .            +b3(Qh,hq,hg,lh)*Dconjg(a2(hq,Qh,hg,lh))
     .            +a3(hq,Qh,hg,lh)*Dconjg(b2(Qh,hq,hg,lh))
     .            +b1(Qh,hq,hg,lh)*Dconjg(a4(hq,Qh,hg,lh)))
        xD(3)=dreal(temp)
        xD(4)=dimag(temp)
c--- eqn (B54) - vanishes because different quarks
c        mE=0d0
c--- eqn (B55)
c        mF=2d0*dreal(mbar1*Dconjg(mbar3)+mbar2*Dconjg(mbar4))
        xF(1)=2d0*dreal(a1(hq,Qh,hg,lh)*Dconjg(a3(hq,Qh,hg,lh))
     .                 +a2(hq,Qh,hg,lh)*Dconjg(a4(hq,Qh,hg,lh)))
        xF(2)=2d0*dreal(b1(Qh,hq,hg,lh)*Dconjg(b3(Qh,hq,hg,lh))
     .                 +b2(Qh,hq,hg,lh)*Dconjg(b4(Qh,hq,hg,lh)))
        temp =2d0*(a1(hq,Qh,hg,lh)*Dconjg(b1(Qh,hq,hg,lh))
     .            +b3(Qh,hq,hg,lh)*Dconjg(a3(hq,Qh,hg,lh))
     .            +a2(hq,Qh,hg,lh)*Dconjg(b2(Qh,hq,hg,lh))
     .            +b4(Qh,hq,hg,lh)*Dconjg(a4(hq,Qh,hg,lh)))
        xF(3)=dreal(temp)
        xF(4)=dimag(temp)
c--- eqn (B56)
c        mG=2d0*dreal(mbar1*Dconjg(mbar4)+mbar2*Dconjg(mbar3))
        xG(1)=2d0*dreal(a1(hq,Qh,hg,lh)*Dconjg(a4(hq,Qh,hg,lh))
     .                 +a2(hq,Qh,hg,lh)*Dconjg(a3(hq,Qh,hg,lh)))
        xG(2)=2d0*dreal(b1(Qh,hq,hg,lh)*Dconjg(b4(Qh,hq,hg,lh))
     .                 +b2(Qh,hq,hg,lh)*Dconjg(b3(Qh,hq,hg,lh)))
        temp=2d0*(a1(hq,Qh,hg,lh)*Dconjg(b2(Qh,hq,hg,lh))
     .           +b3(Qh,hq,hg,lh)*Dconjg(a4(hq,Qh,hg,lh))
     .           +a2(hq,Qh,hg,lh)*Dconjg(b1(Qh,hq,hg,lh))
     .           +b4(Qh,hq,hg,lh)*Dconjg(a3(hq,Qh,hg,lh)))
        xG(3)=dreal(temp)
        xG(4)=dimag(temp)
c--- eqns (B44)-(B49)
c        m0=mB+mC+mE
c        mx=-0.5d0*(3d0*mC+2d0*mE+mB)
c        my=mA+mD
c        mz=mF+mG
c        mxx=0.25d0*(2d0*mC+mE)
c        mxy=-0.5d0*(mF+mD)
c--- eqn (B42) translated
c--- note that C3 = (N*CF**2-Tr*CF)/2 = CF/4*(N**2-2)
c--- in my notation, where Tr=1/2
c        msq(hq,Qh,hg,lh)=cf/2d0*(xn*cf*my+xn**2*mxy+(xn**2-2d0)/2d0*mz)
        do j=1,4
        xy(j)=xA(j)+xD(j)
        xz(j)=xF(j)+xG(j)
        xxy(j)=-0.5d0*(xF(j)+xD(j))
        xsq(hq,Qh,hg,lh,j)=cf/2d0*(xn*cf*xy(j)+xn**2*xxy(j)
     .                            +(xn**2-2d0)/2d0*xz(j))
        enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
