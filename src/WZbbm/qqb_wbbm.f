      subroutine qqb_wbbm(p,msq)
c---  Matrix elements squared 
c     q(-p1)+qb(-p2) --> nu(p3)+e^+(p4)+b(p5)+bb(p6)
c---  averaged(summed) over initial(final) colours and spins
      implicit none
      include 'constants.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j,k
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf)
      double precision qqb,qbq,sumsq
      double precision faclo
      faclo=V*gsq**2*gwsq**2*aveqq

C---Initialize to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C---Fill dot-products
      call dotem(6,p,s)

      qqb=faclo*sumsq(1,2,3,4,6,5)
      qbq=faclo*sumsq(2,1,3,4,6,5)

      do j=-(nf-1),(nf-1)
      do k=-(nf-1),(nf-1)
      if     ((j .gt. 0) .and. (k .lt. 0)) then
               msq(j,k)=Vsq(j,k)*qqb
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
               msq(j,k)=Vsq(j,k)*qbq
      endif
      enddo
      enddo

      return
      end

      double precision function sumsq(p1,p2,p3,p4,p5,p6)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer p1,p2,p3,p4,p5,p6
      double precision s56,s134,s234,prop
      s134=s(p1,p3)+s(p1,p4)+s(p3,p4)
      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      s56=s134+s234+s(p1,p2)-s(p3,p4)
c---overall factor of 4 removed
      prop=s56**2*((s(p3,p4)-wmass**2)**2+(wmass*wwidth)**2)
      sumsq = 
     . +2d0*s(p1,p2)*s(p1,p2)*s(p3,p4)*mb**2/(s134*s234)
 
     . +2d0*s(p1,p2)*s(p1,p3)*s(p1,p4)*mb**2/s134**2
 
     . -2d0*s(p1,p2)*s(p1,p3)*s(p2,p4)*mb**2/(s134*s234)
 
     . +2d0*s(p1,p2)*s(p1,p3)*s(p3,p4)*mb**2/(s134*s234)
 
     . -2d0*s(p1,p2)*s(p1,p4)*s(p2,p3)*mb**2/(s134*s234)
 
     . +2d0*s(p1,p2)*s(p1,p4)*s(p3,p4)*mb**2/s134**2
 
     . +s(p1,p2)*s(p1,p5)*s(p2,p6)*s(p3,p4)/(s134*s234)
 
     . +s(p1,p2)*s(p1,p5)*s(p3,p4)*s(p3,p6)/(s134*s234)
 
     . +s(p1,p2)*s(p1,p6)*s(p2,p5)*s(p3,p4)/(s134*s234)
 
     . +s(p1,p2)*s(p1,p6)*s(p3,p4)*s(p3,p5)/(s134*s234)
 
     . +2d0*s(p1,p2)*s(p2,p3)*s(p2,p4)*mb**2/s234**2
 
     . +2d0*s(p1,p2)*s(p2,p3)*s(p3,p4)*mb**2/s234**2
 
     . +2d0*s(p1,p2)*s(p2,p4)*s(p3,p4)*mb**2/(s134*s234)
 
     . +s(p1,p2)*s(p2,p5)*s(p3,p4)*s(p4,p6)/(s134*s234)
 
     . +s(p1,p2)*s(p2,p6)*s(p3,p4)*s(p4,p5)/(s134*s234)
 
     . +2d0*s(p1,p2)*s(p3,p4)*s(p3,p4)*mb**2/(s134*s234)
 
     . +s(p1,p2)*s(p3,p4)*s(p3,p5)*s(p4,p6)/(s134*s234)
 
     . +s(p1,p2)*s(p3,p4)*s(p3,p6)*s(p4,p5)/(s134*s234)
 
     . -2d0*s(p1,p3)*s(p1,p3)*s(p2,p4)*mb**2/(s134*s234)
 
     . +s(p1,p3)*s(p1,p4)*s(p1,p5)*s(p2,p6)/s134**2
 
     . +s(p1,p3)*s(p1,p4)*s(p1,p6)*s(p2,p5)/s134**2
 
     . +2d0*s(p1,p3)*s(p1,p4)*s(p2,p3)*mb**2/(s134*s234)
 
     . +2d0*s(p1,p3)*s(p1,p4)*s(p2,p4)*mb**2/s134**2
 
     . +s(p1,p3)*s(p1,p4)*s(p2,p5)*s(p4,p6)/s134**2
 
     . +s(p1,p3)*s(p1,p4)*s(p2,p6)*s(p4,p5)/s134**2
 
     . +s(p1,p3)*s(p1,p5)*s(p2,p3)*s(p4,p6)/(s134*s234)
 
     . -s(p1,p3)*s(p1,p5)*s(p2,p4)*s(p2,p6)/(s134*s234)
 
     . -s(p1,p3)*s(p1,p5)*s(p2,p4)*s(p3,p6)/(s134*s234)
 
     . +s(p1,p3)*s(p1,p6)*s(p2,p3)*s(p4,p5)/(s134*s234)
 
     . -s(p1,p3)*s(p1,p6)*s(p2,p4)*s(p2,p5)/(s134*s234)
 
     . -s(p1,p3)*s(p1,p6)*s(p2,p4)*s(p3,p5)/(s134*s234)
 
     . +2d0*s(p1,p3)*s(p2,p3)*s(p2,p4)*mb**2/s234**2
 
     . +2d0*s(p1,p3)*s(p2,p3)*s(p3,p4)*mb**2/s234**2
 
     . +2d0*s(p1,p3)*s(p2,p3)*s(p4,p5)*s(p4,p6)/(s134*s234)
 
     . -2d0*s(p1,p3)*s(p2,p4)*s(p2,p4)*mb**2/(s134*s234)
 
     . -s(p1,p3)*s(p2,p4)*s(p2,p5)*s(p4,p6)/(s134*s234)
 
     . -s(p1,p3)*s(p2,p4)*s(p2,p6)*s(p4,p5)/(s134*s234)
 
     . -2d0*s(p1,p3)*s(p2,p4)*s(p3,p4)*mb**2/(s134*s234)
 
     . -s(p1,p3)*s(p2,p4)*s(p3,p5)*s(p4,p6)/(s134*s234)
 
     . -s(p1,p3)*s(p2,p4)*s(p3,p6)*s(p4,p5)/(s134*s234)
 
     . -2d0*s(p1,p4)*s(p1,p4)*s(p2,p3)*mb**2/s134**2
 
     . -s(p1,p4)*s(p1,p4)*s(p2,p5)*s(p3,p6)/s134**2
 
     . -s(p1,p4)*s(p1,p4)*s(p2,p6)*s(p3,p5)/s134**2
 
     . -s(p1,p4)*s(p1,p5)*s(p2,p3)*s(p2,p6)/(s134*s234)
 
     . +s(p1,p4)*s(p1,p5)*s(p2,p6)*s(p3,p4)/s134**2
 
     . -s(p1,p4)*s(p1,p6)*s(p2,p3)*s(p2,p5)/(s134*s234)
 
     . +s(p1,p4)*s(p1,p6)*s(p2,p5)*s(p3,p4)/s134**2
 
     . -2d0*s(p1,p4)*s(p2,p3)*s(p2,p3)*mb**2/s234**2
 
     . +2d0*s(p1,p4)*s(p2,p3)*s(p2,p4)*mb**2/(s134*s234)
 
     . +2d0*s(p1,p4)*s(p2,p3)*s(p3,p4)*mb**2/(s134*s234)
 
     . +2d0*s(p1,p4)*s(p2,p3)*s(p3,p4)*s(p5,p6)/(s134*s234)
 
     . -s(p1,p4)*s(p2,p3)*s(p3,p5)*s(p4,p6)/(s134*s234)
 
     . -s(p1,p4)*s(p2,p3)*s(p3,p6)*s(p4,p5)/(s134*s234)
 
     . +s(p1,p4)*s(p2,p4)*s(p2,p5)*s(p3,p6)/(s134*s234)
 
     . +s(p1,p4)*s(p2,p4)*s(p2,p6)*s(p3,p5)/(s134*s234)
 
     . +2d0*s(p1,p4)*s(p2,p4)*s(p3,p4)*mb**2/s134**2
 
     . +2d0*s(p1,p4)*s(p2,p4)*s(p3,p5)*s(p3,p6)/(s134*s234)
 
     . -2d0*s(p1,p4)*s(p2,p5)*s(p2,p6)*s(p3,p4)/(s134*s234)
 
     . -s(p1,p4)*s(p2,p5)*s(p3,p4)*s(p3,p6)/(s134*s234)
 
     . +s(p1,p4)*s(p2,p5)*s(p3,p4)*s(p4,p6)/s134**2
 
     . -s(p1,p4)*s(p2,p6)*s(p3,p4)*s(p3,p5)/(s134*s234)
 
     . +s(p1,p4)*s(p2,p6)*s(p3,p4)*s(p4,p5)/s134**2
 
     . -2d0*s(p1,p5)*s(p1,p6)*s(p2,p3)*s(p3,p4)/(s134*s234)
 
     . -s(p1,p5)*s(p2,p3)*s(p2,p3)*s(p4,p6)/s234**2
 
     . +s(p1,p5)*s(p2,p3)*s(p2,p4)*s(p2,p6)/s234**2
 
     . +s(p1,p5)*s(p2,p3)*s(p2,p4)*s(p3,p6)/s234**2
 
     . +s(p1,p5)*s(p2,p3)*s(p2,p6)*s(p3,p4)/s234**2
 
     . +s(p1,p5)*s(p2,p3)*s(p3,p4)*s(p3,p6)/s234**2
 
     . -s(p1,p5)*s(p2,p3)*s(p3,p4)*s(p4,p6)/(s134*s234)
 
     . -s(p1,p6)*s(p2,p3)*s(p2,p3)*s(p4,p5)/s234**2
 
     . +s(p1,p6)*s(p2,p3)*s(p2,p4)*s(p2,p5)/s234**2
 
     . +s(p1,p6)*s(p2,p3)*s(p2,p4)*s(p3,p5)/s234**2
 
     . +s(p1,p6)*s(p2,p3)*s(p2,p5)*s(p3,p4)/s234**2
 
     . +s(p1,p6)*s(p2,p3)*s(p3,p4)*s(p3,p5)/s234**2
 
     . -s(p1,p6)*s(p2,p3)*s(p3,p4)*s(p4,p5)/(s134*s234) 

c      write(6,*) 'sumsq',sumsq
      sumsq=sumsq/prop
      return
      end
