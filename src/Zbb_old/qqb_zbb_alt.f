      subroutine qqb_zbb_alt(p,msq)
c---  Matrix elements squared 
c     q(-p1)+qb(-p2) --> b(p5)+bb(p6)+e^-(p3)+e^+(p4)
c---  averaged(summed) over initial(final) colours and spins
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'prods.f'
      include 'zcouple.f'
c      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'hardscale.f'
      integer k,nu,iperm
      integer i1(2),i2(2),i3(2),i4(2),j,lh,h2,h3,hq,h(2:3)
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),msqzbb,mmsq(2,2),
     . prop,qqb,qbq,pswap(mxpart,4)
      double precision mqqb(2,2),m1,m0,x,fac
      double complex tempm0,m(2)
      double complex amp_qqgg,a6treeg,a6treeg1,tamp
      character*9 st(2,2)
      logical compare
      parameter(x=xn/cf)
      data i1/1,4/
      data i2/2,3/
      data i3/3,2/
      data i4/4,1/
      data st/'q+g-g-qb-','q+g-g+qb-','q+g+g-qb-','q+g+g+qb-'/
      
      compare=.false.
      
c--initialize to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(1,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(6,nu)
      pswap(5,nu)=p(3,nu)
      pswap(6,nu)=p(4,nu)
      enddo
      call spinoru(6,pswap,za,zb)

c--- USED TO COMPARE ONLY      
      if (compare) then
      do hq=1,2
      do lh=1,2
      mqqb(hq,lh)=0d0
      
      do h2=1,2
      do h3=1,2
        h(2)=h2
        h(3)=h3
C initialize loop sum to zero
        tempm0=czip
        do j=1,2
        m(j)=amp_qqgg(1,hq,i2(j),h(i2(j)),i3(j),h(i3(j)),4,lh)
c        write(*,*) hq,h2,h3,lh,i2(j),i3(j),m(j)
        tempm0=tempm0+m(j)
        enddo
      m1=Dble(Dconjg(m(1))*m(2))
      m0=abs(tempm0)**2
      mqqb(hq,lh)=mqqb(hq,lh)+fac*(m0-x*m1)
      enddo
      enddo
      write(*,*) 'old mqqb(',hq,',',lh,')',mqqb(hq,lh)
      enddo
      enddo
      endif
c--- USED TO COMPARE ONLY      
      
      do hq=1,2
      do lh=1,2
      mqqb(hq,lh)=0d0
c--- set iperm=1 for same-helicity quarks and leptons, 2 otherwise
      iperm=1+abs(hq-lh)

      do h2=1,2
      do h3=1,2
        h(2)=h2
        h(3)=h3
C initialize loop sum to zero
        tempm0=czip
        do j=1,2
        m(j)=a6treeg1(st(h(i2(j)),h(i3(j))),
     .                 i1(iperm),i2(j),i3(j),i4(iperm),5,6,za,zb)
c        write(*,*) hq,h2,h3,lh,i2(j),i3(j),m(j)
        tempm0=tempm0+m(j)
        enddo
      m1=Dble(Dconjg(m(1))*m(2))
      m0=abs(tempm0)**2
      mqqb(hq,lh)=mqqb(hq,lh)+fac*(m0-x*m1)
      enddo
      enddo
      if (compare) write(*,*) 'new mqqb(',hq,',',lh,')',mqqb(hq,lh)
      enddo
      enddo
      
c--- USED TO COMPARE ONLY      
      if (compare) then
      tamp=a6treeg('q+g+g+qb-',4,3,2,1,6,5,zb,za)
c      write(*,*) 'q+g-g-qb-',tamp
      tamp=a6treeg('q+g-g+qb-',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+g-g+qb-',tamp
      tamp=a6treeg('q+g+g-qb-',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+g+g-qb-',tamp
      tamp=a6treeg('q+g+g+qb-',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+g+g+qb-',tamp
      tamp=a6treeg('q+g+qb-g-',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+g+qb-g-',tamp
      tamp=a6treeg('q+g+qb-g+',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+g+qb-g+',tamp
      tamp=a6treeg('q+qb-g-g+',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+qb-g-g+',tamp
      tamp=a6treeg('q+qb-g+g-',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+qb-g+g-',tamp
      tamp=a6treeg('q+qb-g+g+',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+qb-g+g+',tamp
      pause
      endif
c--- USED TO COMPARE ONLY      
      
      return
      end
      
      double complex function a6treeg1(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
c----wrapper to a6treeg that also includes config st='q+g-g-qb-'
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      character*9 st
      double complex a6treeg

      if(st.eq.'q+g-g-qb-') then
        a6treeg1=a6treeg('q+g+g+qb-',j4,j3,j2,j1,j6,j5,zb,za)
      else
        a6treeg1=a6treeg(st,j1,j2,j3,j4,j5,j6,za,zb)
      endif
      
      return
      end
      
