      real*8 function twojet(j,k,nwz)
      implicit none
      include 'constants.f'
      character*2 qq(-nf:nf),vb(-1:1)
      integer j,k,j1,k1,ipoint(-1:1,-1:1),nwz
      real*8 msqpc(9),msqpv(9)
      real*8 PC(-nf:nf,-nf:nf),PV(-nf:nf,-nf:nf)
      real*8 xmqqqq
      include 'ckm.f'
      common/qqggcouple/PC,PV     
      common/msq2/msqpc,msqpv     
      common/qq/qq
      common/vb/vb
case 1 ---  qj qk
case 2 ---  qj qj or qbarj qbarj
case 3 ---  qj qbk
case 4 ---  qbj qk
case 5 ---  qj g   
case 6 ---  qbarj g   
case 7 ---  g qbarj
case 8 ---  g qj
case 9 ---  g g

      data ipoint /2, 7, 3,         
     &             6, 9, 5,         
     &             4, 8, 1/
case 1 ---  q q
case 2 ---  qbar qbar
case 3 ---  q qb
case 4 ---  qb q
case 5 ---  q g   
case 6 ---  qbar g   
case 7 ---  g qbar
case 8 ---  g q
case 9 ---  g g
          
          
      twojet=0d0
      if (j.eq.0) then       
c---  particle one gluon
         j1 = 0     
      else
c     j1 set equal to +1 (quark) or -1 (antiquark) 
c sign(a,b)=abs(a)*sign(b)
         j1 = sign(1,j)      
         endif      
          
      if (k.eq.0) then       
c---  particle 2 gluon
         k1 = 0     
      else
         k1 = sign(1,k)      
      endif      

c      write(6,*) 'in twojet j,k',j,k
c      write(6,*) 'in twojet',vb(nwz),':',qq(j),' ',qq(k)
c      write(6,*) 'in twojet PC(j,k),PV(j,k)',PC(j,k),PV(j,k)

      twojet=+PC(j,k)*msqpc(ipoint(j1,k1))
     &       +PV(j,k)*msqpv(ipoint(j1,k1)) 
     &       +xmqqqq(j,k,nwz)   

      return        
      end 


