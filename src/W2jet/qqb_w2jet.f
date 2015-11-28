      subroutine qqb_w2jet(p,msq)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) --> W +g(p5) +g(p6)
c                          |
c                          --> nu(p3)+e^+(p4)
c                            
c   positively charged W only with p6 soft
c--all momenta incoming
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'sprodx.f'
      include 'dprodx.f'
      integer j,k,jj,kk,qp(6)
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision prop,FAC,ggWqbq2
      double precision qqbWgg2,qbqWgg2,qgWqg2,qbgWqbg2,gqbWqbg2,gqWqg2
      double precision qqWqq2(-nf:nf,-nf:nf)
      double precision xwqqqq,xwqqbqqb
      logical fourQ      
      data fourQ/.false./
      data qp/1,2,0,0,3,4/
      save qp,fourQ
      
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call spinoru(6,p,za,zb)

      call w2jetsq(1,2,3,4,5,6,za,zb,qqbWgg2)
      call w2jetsq(2,1,3,4,5,6,za,zb,qbqWgg2)
      call w2jetsq(1,5,3,4,2,6,za,zb,qgWqg2)
      call w2jetsq(2,5,3,4,1,6,za,zb,gqWqg2)
      call w2jetsq(5,1,3,4,2,6,za,zb,qbgWqbg2)
      call w2jetsq(5,2,3,4,1,6,za,zb,gqbWqbg2)
      call w2jetsq(5,6,3,4,1,2,za,zb,ggWqbq2)

      prop=s(3,4)**2/((s(3,4)-wmass**2)**2+wmass**2*wwidth**2)
      FAC=V*XN/four*(gwsq/2d0)**2*gsq**2*prop

c--- note statistical factor of one half for two gluons in the final state
      qqbWgg2 =half*aveqq*fac*qqbWgg2
      qbqWgg2 =half*aveqq*fac*qbqWgg2
      gqWqg2  =aveqg*fac*gqWqg2
      qgWqg2  =aveqg*fac*qgWqg2
      gqbWqbg2=aveqg*fac*gqbWqbg2
      qbgWqbg2=aveqg*fac*qbgWqbg2
      ggWqbq2 =avegg*fac*ggWqbq2
      
      call initqqqq(za,zb)
      
      do j=-nf,nf
      do k=-nf,nf      
          qqWqq2(j,k)=0d0
        if (fourQ) then
          if     ((j .gt. 0) .and. (k .lt. 0)) then
            qqWqq2(j,k)=xwqqbqqb(qp(1),qp(6),qp(5),qp(2),j,k)
          elseif ((j .lt. 0) .and. (k .gt. 0)) then
            qqWqq2(j,k)=xwqqbqqb(qp(6),qp(2),qp(5),qp(1),j,k)
          elseif ((j .gt. 0) .and. (k .gt. 0)) then
            qqWqq2(j,k)=xwqqqq(qp(1),qp(2),qp(5),qp(6),j,k)
          elseif ((j .lt. 0) .and. (k .lt. 0)) then
            qqWqq2(j,k)=xwqqqq(qp(5),qp(6),qp(1),qp(2),j,k)
          endif
        endif
      enddo
      enddo
      
      do j=-nf,nf
      do k=-nf,nf

      if     ((j .gt. 0) .and. (k .lt. 0)) then
          msq(j,k)=Vsq(j,k)*qqbWgg2+qqWqq2(j,k)
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=Vsq(j,k)*qbqWgg2+qqWqq2(j,k)
      elseif ((j .gt. 0) .and. (k .gt. 0)) then
          msq(j,k)=+qqWqq2(j,k)
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
          msq(j,k)=+qqWqq2(j,k)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &   (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWqg2
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqbg2
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWqg2
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqbg2
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
          msq(j,k)=0d0
          do jj=1,nf
            do kk=-nf,-1
              msq(j,k)=msq(j,k)+Vsq(jj,kk)
            enddo
          enddo
          msq(j,k)=msq(j,k)*ggWqbq2
      endif

      enddo
      enddo

      return
      end
 
      subroutine w2jetsq(i1,i2,i3,i4,i5,i6,za,zb,msq)
      implicit none
      include 'constants.f'
      include 'sprodx.f'
      double complex qcd1(-1:1,-1:1),qcd2(-1:1,-1:1),qed(-1:1,-1:1)
      double precision msq1,msq2,msqq,msq
      integer i1,i2,i3,i4,i5,i6
      
      call subqcd(i1,i2,i3,i4,i5,i6,za,zb,qcd1)
      call subqcd(i1,i2,i3,i4,i6,i5,za,zb,qcd2)
            
      qed(+1,+1)=qcd1(+1,+1)+qcd2(+1,+1)     
      qed(+1,-1)=qcd1(+1,-1)+qcd2(-1,+1)     
      qed(-1,+1)=qcd1(-1,+1)+qcd2(+1,-1)     
      qed(-1,-1)=qcd1(-1,-1)+qcd2(-1,-1)     
                   
      msq1= abs(qcd1(+1,+1))**2+abs(qcd1(+1,-1))**2
     .     +abs(qcd1(-1,+1))**2+abs(qcd1(-1,-1))**2 
     
      msq2= abs(qcd2(+1,+1))**2+abs(qcd2(+1,-1))**2
     .     +abs(qcd2(-1,+1))**2+abs(qcd2(-1,-1))**2 
     
      msqq= abs( qed(+1,+1))**2+abs( qed(+1,-1))**2
     .     +abs( qed(-1,+1))**2+abs( qed(-1,-1))**2 

      msq=msq1+msq2-ninth*msqq

      return
      end
          
    
