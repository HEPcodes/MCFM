      subroutine qqb_wbb_g(p,msq)
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  g*  + W +g(p7)
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+bb(p6)
c   positively charged W only
      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'hardscale.f'
      integer j,k
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),redmsq,fac
      double precision qqbWbbg,qbqWbbg,qgWbbq,gqWbbq,gqbWbbqb,qbgWbbqb


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call spinoru(7,p,za,zb)

      if (
     .      (s(5,6) .lt. four*hscalesq) 
     . .or. (s(1,5)*s(2,5)/s(1,2) .lt. hscalesq) 
     . .or. (s(1,6)*s(2,6)/s(1,2) .lt. hscalesq) ) return

      fac=gsq**3*gw**4/4d0
      qqbWbbg =+redmsq(1,2,7,5,6,3,4,4)*fac*aveqq
      qbqWbbg =+redmsq(2,1,7,5,6,3,4,4)*fac*aveqq
      qgWbbq  =+redmsq(1,7,2,5,6,3,4,4)*fac*aveqg
      qbgWbbqb=+redmsq(7,1,2,5,6,3,4,4)*fac*aveqg
      gqWbbq  =+redmsq(2,7,1,5,6,3,4,4)*fac*aveqg
      gqbWbbqb=+redmsq(7,2,1,5,6,3,4,4)*fac*aveqg

      do j=-(nf-1),(nf-1)
      do k=-(nf-1),(nf-1)

      msq(j,k)=0d0

      if     ((j .gt. 0) .and. (k .lt. 0)) then
      msq(j,k)=Vsq(j,k)*qqbWbbg

      elseif ((j .lt. 0) .and. (k .gt. 0)) then
      msq(j,k)=Vsq(j,k)*qbqWbbg

      elseif ((j .gt. 0) .and. (k .eq. 0)) then
      msq(j,k)=Vsum(j)*qgWbbq 
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
      msq(j,k)=Vsum(j)*qbgWbbqb

      elseif ((j .eq. 0) .and. (k .gt. 0)) then
      msq(j,k)=Vsum(k)*gqWbbq

      elseif ((j .eq. 0) .and. (k .lt. 0)) then
      msq(j,k)=Vsum(k)*gqbWbbqb
      endif

      enddo
      enddo

      return
      end


      double precision function redmsq(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
c matrix element squared summed over colors and spins
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex qedipp,qedimm,qedimp,qedipm,ipp,imm,imp,ipm
      double complex qedfpp,qedfmm,qedfmp,qedfpm,fpp,fmm,fmp,fpm
      double complex qcdapp,qcdamm,qcdamp,qcdapm,app,amm,amp,apm
      double complex qcdbpp,qcdbmm,qcdbmp,qcdbpm,bpp,bmm,bmp,bpm
      double precision pm,mm,pp,mp,prop
c---calculate the W propagator
      prop=((s(j6,j7)-wmass**2)**2+(wmass*wwidth)**2)
      

      apm=qcdapm(j1,j2,j3,j4,j5,j6,j7,jb)
      bpm=qcdbpm(j1,j2,j3,j4,j5,j6,j7,jb)
      ipm=qedipm(j1,j2,j3,j4,j5,j6,j7,jb)
      fpm=qedfpm(j1,j2,j3,j4,j5,j6,j7,jb)
      
      pm=V*xn/eight*(cdabs(apm)**2+cdabs(bpm)**2)
     &+V/(eight*xn)*(cdabs(ipm)**2+cdabs(fpm)**2
     &-two*(cdabs(ipm+fpm))**2)
      
      amm=qcdamm(j1,j2,j3,j4,j5,j6,j7,jb)
      bmm=qcdbmm(j1,j2,j3,j4,j5,j6,j7,jb)
      imm=qedimm(j1,j2,j3,j4,j5,j6,j7,jb)
      fmm=qedfmm(j1,j2,j3,j4,j5,j6,j7,jb)

      mm=V*xn/eight*(cdabs(amm)**2+cdabs(bmm)**2)
     &+V/(eight*xn)*(cdabs(imm)**2+cdabs(fmm)**2
     &-two*(cdabs(imm+fmm))**2)

      amp=qcdamp(j1,j2,j3,j4,j5,j6,j7,jb)
      bmp=qcdbmp(j1,j2,j3,j4,j5,j6,j7,jb)
      imp=qedimp(j1,j2,j3,j4,j5,j6,j7,jb)
      fmp=qedfmp(j1,j2,j3,j4,j5,j6,j7,jb)

      mp=V*xn/eight*(cdabs(amp)**2+cdabs(bmp)**2)
     &+V/(eight*xn)*(cdabs(imp)**2+cdabs(fmp)**2
     &-two*(cdabs(imp+fmp))**2)

      app=qcdapp(j1,j2,j3,j4,j5,j6,j7,jb)
      bpp=qcdbpp(j1,j2,j3,j4,j5,j6,j7,jb)
      ipp=qedipp(j1,j2,j3,j4,j5,j6,j7,jb)
      fpp=qedfpp(j1,j2,j3,j4,j5,j6,j7,jb)

      pp=V*xn/eight*(cdabs(app)**2+cdabs(bpp)**2)
     &+V/(eight*xn)*(cdabs(ipp)**2+cdabs(fpp)**2
     &-two*(cdabs(ipp+fpp))**2)

      redmsq=(pm+mm+pp+mp)/prop
      return
      end
