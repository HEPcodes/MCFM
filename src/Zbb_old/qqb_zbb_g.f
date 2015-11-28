      subroutine qqb_zbb_g(p,msq)
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  g*  + Z +g(p7)
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
      integer j,k,nu
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf),mmsq(2,2)
      double precision redmsqz,fac,prop,pswap(mxpart,4)
      double precision qqbZbbg1,qbqZbbg1,qgZbbq1,
     .                 gqZbbq1,gqbZbbqb1,qbgZbbqb1
      double precision qqbZbbg2,qbqZbbg2,qgZbbq2,
     .                 gqZbbq2,gqbZbbqb2,qbgZbbqb2


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
      prop=s(3,4)/sqrt(((s(3,4)-zmass**2)**2+(zmass*zwidth)**2))

      if (
     .      (s(5,6) .lt. four*hscalesq) 
     . .or. (s(1,5)*s(2,5)/s(1,2) .lt. hscalesq) 
     . .or. (s(1,6)*s(2,6)/s(1,2) .lt. hscalesq) ) return

c--- note the factor of 4d0*xw**2 relative to wbb
      fac=gsq**3*esq**2
      qqbZbbg1 =+redmsqz(1,2,7,5,6,3,4,4)*fac*aveqq
      qqbZbbg2 =+redmsqz(2,1,7,6,5,3,4,4)*fac*aveqq
      qbqZbbg1 =+redmsqz(2,1,7,5,6,3,4,4)*fac*aveqq
      qbqZbbg2 =+redmsqz(1,2,7,6,5,3,4,4)*fac*aveqq
      qgZbbq1  =+redmsqz(1,7,2,5,6,3,4,4)*fac*aveqg
      qgZbbq2  =+redmsqz(7,1,2,6,5,3,4,4)*fac*aveqg
      qbgZbbqb1=+redmsqz(7,1,2,5,6,3,4,4)*fac*aveqg
      qbgZbbqb2=+redmsqz(1,7,2,6,5,3,4,4)*fac*aveqg
      gqZbbq1  =+redmsqz(2,7,1,5,6,3,4,4)*fac*aveqg
      gqZbbq2  =+redmsqz(7,2,1,6,5,3,4,4)*fac*aveqg
      gqbZbbqb1=+redmsqz(7,2,1,5,6,3,4,4)*fac*aveqg
      gqbZbbqb2=+redmsqz(2,7,1,6,5,3,4,4)*fac*aveqg

      do j=-nf,nf
      do k=-nf,nf
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19
      
      msq(j,k)=0d0

      if     ((j .eq. 0) .and. (k .eq. 0)) then
C---L(1),R(1) for coupling to down quark
      msq(j,k)=(Q(1)*q1+prop*L(1)*l1)**2*mmsq(1,1)
     .        +(Q(1)*q1+prop*R(1)*l1)**2*mmsq(2,1)
     .        +(Q(1)*q1+prop*L(1)*r1)**2*mmsq(1,2)
     .        +(Q(1)*q1+prop*R(1)*r1)**2*mmsq(2,2)
c----debug
c          msq(j,k)=prop**2*(mmsq(1,1)+mmsq(2,1)+mmsq(1,2)+mmsq(2,2))
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
      msq(j,k)=((Q(j)*q1+prop*L(j)*l1)**2
     .         +(Q(j)*q1+prop*R(j)*r1)**2)*qqbZbbg1
     .        +((Q(j)*q1+prop*L(j)*r1)**2
     .         +(Q(j)*q1+prop*R(j)*l1)**2)*qqbZbbg2
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
      msq(j,k)=((Q(k)*q1+prop*L(k)*l1)**2
     .         +(Q(k)*q1+prop*R(k)*r1)**2)*qbqZbbg1
     .        +((Q(k)*q1+prop*L(k)*r1)**2
     .         +(Q(k)*q1+prop*R(k)*l1)**2)*qbqZbbg2    
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
      msq(j,k)=((Q(j)*q1+prop*L(j)*l1)**2
     .         +(Q(j)*q1+prop*R(j)*r1)**2)*qgZbbq1
     .        +((Q(j)*q1+prop*L(j)*r1)**2
     .         +(Q(j)*q1+prop*R(j)*l1)**2)*qgZbbq2
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
      msq(j,k)=((Q(-j)*q1+prop*L(-j)*l1)**2
     .         +(Q(-j)*q1+prop*R(-j)*r1)**2)*qbgZbbqb1
     .        +((Q(-j)*q1+prop*L(-j)*r1)**2
     .         +(Q(-j)*q1+prop*R(-j)*l1)**2)*qbgZbbqb2
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
      msq(j,k)=((Q(k)*q1+prop*L(k)*l1)**2
     .         +(Q(k)*q1+prop*R(k)*r1)**2)*gqZbbq1
     .        +((Q(k)*q1+prop*L(k)*r1)**2
     .         +(Q(k)*q1+prop*R(k)*l1)**2)*gqZbbq2
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
      msq(j,k)=((Q(-k)*q1+prop*L(-k)*l1)**2
     .         +(Q(-k)*q1+prop*R(-k)*r1)**2)*gqbZbbqb1
     .        +((Q(-k)*q1+prop*L(-k)*r1)**2
     .         +(Q(-k)*q1+prop*R(-k)*l1)**2)*gqbZbbqb2
      endif

   19 continue
      enddo
      enddo
      return
      end


      double precision function redmsqz(j1,j2,j7,j5,j6,j3,j4,jb)
      implicit none
c matrix element squared summed over colors and spins
      include 'constants.f'
      include 'masses.f'
      include 'prods.f'
      integer j1,j2,j7,j5,j6,j3,j4,jb
      double complex qedipp,qedimm,qedimp,qedipm,ipp,imm,imp,ipm
      double complex qedfpp,qedfmm,qedfmp,qedfpm,fpp,fmm,fmp,fpm
      double complex qcdapp,qcdamm,qcdamp,qcdapm,app,amm,amp,apm
      double complex qcdbpp,qcdbmm,qcdbmp,qcdbpm,bpp,bmm,bmp,bpm
      double precision pm,mm,pp,mp,prop
c---calculate the Z propagator
      prop=s(j3,j4)**2
      

      apm=qcdapm(j1,j2,j7,j5,j6,j3,j4,jb)
      bpm=qcdbpm(j1,j2,j7,j5,j6,j3,j4,jb)
      ipm=qedipm(j1,j2,j7,j5,j6,j3,j4,jb)
      fpm=qedfpm(j1,j2,j7,j5,j6,j3,j4,jb)

      pm=V*xn/eight*(abs(apm)**2+abs(bpm)**2)
     &+V/(eight*xn)*(abs(ipm)**2+abs(fpm)**2-two*(abs(ipm+fpm))**2)

c      pm=V*xn/eight*(abs(apm)**2+abs(bpm)**2)
c     &+0*V/(eight*xn)*(abs(ipm)**2+abs(fpm)**2-two*(abs(ipm+fpm))**2)
      
      amm=qcdamm(j1,j2,j7,j5,j6,j3,j4,jb)
      bmm=qcdbmm(j1,j2,j7,j5,j6,j3,j4,jb)
      imm=qedimm(j1,j2,j7,j5,j6,j3,j4,jb)
      fmm=qedfmm(j1,j2,j7,j5,j6,j3,j4,jb)

      mm=V*xn/eight*(abs(amm)**2+abs(bmm)**2)
     &+V/(eight*xn)*(abs(imm)**2+abs(fmm)**2-two*(abs(imm+fmm))**2)

c      mm=V*xn/eight*(abs(amm)**2+abs(bmm)**2)
c     &+0*V/(eight*xn)*(abs(imm)**2+abs(fmm)**2-two*(abs(imm+fmm))**2)

      amp=qcdamp(j1,j2,j7,j5,j6,j3,j4,jb)
      bmp=qcdbmp(j1,j2,j7,j5,j6,j3,j4,jb)
      imp=qedimp(j1,j2,j7,j5,j6,j3,j4,jb)
      fmp=qedfmp(j1,j2,j7,j5,j6,j3,j4,jb)

      mp=V*xn/eight*(abs(amp)**2+abs(bmp)**2)
     &+V/(eight*xn)*(abs(imp)**2+abs(fmp)**2-two*(abs(imp+fmp))**2)

c      mp=V*xn/eight*(abs(amp)**2+abs(bmp)**2)
c     &+0*V/(eight*xn)*(abs(imp)**2+abs(fmp)**2-two*(abs(imp+fmp))**2)

      app=qcdapp(j1,j2,j7,j5,j6,j3,j4,jb)
      bpp=qcdbpp(j1,j2,j7,j5,j6,j3,j4,jb)
      ipp=qedipp(j1,j2,j7,j5,j6,j3,j4,jb)
      fpp=qedfpp(j1,j2,j7,j5,j6,j3,j4,jb)

      pp=V*xn/eight*(abs(app)**2+abs(bpp)**2)
     &+V/(eight*xn)*(abs(ipp)**2+abs(fpp)**2-two*(abs(ipp+fpp))**2)

c      pp=V*xn/eight*(abs(app)**2+abs(bpp)**2)
c     &+0*V/(eight*xn)*(abs(ipp)**2+abs(fpp)**2-two*(abs(ipp+fpp))**2)

      redmsqz=(pm+mm+pp+mp)/prop
      return
      end
