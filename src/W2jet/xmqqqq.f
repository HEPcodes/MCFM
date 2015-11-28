      real*8 function xmqqqq(n1,n2,nwz)
      implicit none
c*******************************************************************
c     the matrix elements of the  
C     helicity amplitudes for the QCD process
c     q(-q1)+q(-q2)+a(-q6) --> q(q3)+q(q4)+l(q5)
c     all squared
c     multiplied by (((a+l)^2-M**2)^2+M^2*Gam^2)/(a+l)^4/g^4
c     averaged over initial colour or spin
c*******************************************************************
      include 'constants.f'
      real*8 xzqqqq,xwqqqq,xwqqbqqb

      integer n1,n2,nwz,j


      logical first
      integer i(4),iq_q(4),iq_qb(4),iqb_q(4),iqb_qb(4)
      data iq_q  /1,2,3,4/
      data iq_qb /1,4,3,2/
      data iqb_q /4,2,3,1/
      data iqb_qb/3,4,1,2/
      data first/.true./

      
      if (first) then
      write(6,*) 'nwz',nwz
      first=.false.
      endif
      xmqqqq=0d0

      if ((n1 .eq. 0) .or. (n2 .eq. 0)) then
      xmqqqq=0d0
      return
      elseif ((n1 .gt. 0) .and. (n2 .gt. 0)) then
      do j=1,4
      i(j)=iq_q(j)
      enddo
      elseif ((n1 .gt. 0) .and. (n2 .lt. 0)) then
      do j=1,4
      i(j)=iq_qb(j)
      enddo
      elseif ((n1 .lt. 0) .and. (n2 .gt. 0)) then
      do j=1,4
      i(j)=iqb_q(j)
      enddo
      elseif ((n1 .lt. 0) .and. (n2. lt. 0)) then
      do j=1,4
      i(j)=iqb_qb(j)
      enddo
      endif


 
C************************************************************
Case Z0                                                     *
C************************************************************
      if (nwz .eq. 0) then
      xmqqqq=xzqqqq(i,n1,n2)

      elseif (nwz .ne. 0) then
      if ((n1 .gt. 0) .and. (n2 .gt. 0)) then
      xmqqqq=xwqqqq(i,n1,n2)
      elseif ((n1 .lt. 0) .and. (n2 .lt. 0)) then
      xmqqqq=xwqqqq(i,n1,n2)
      elseif ((n1 .lt. 0) .and. (n2 .gt. 0)) then
      xmqqqq=xwqqbqqb(i,n1,n2)
      elseif ((n1 .gt. 0) .and. (n2 .lt. 0)) then
      xmqqqq=xwqqbqqb(i,n1,n2)
      endif 
      endif 

      return
      end


