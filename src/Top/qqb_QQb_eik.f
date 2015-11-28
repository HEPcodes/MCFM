      subroutine qqb_QQb_eik(p,msq) 
      implicit none

************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2002.                                                     *
*     calculate the element squared for the process                    *
c----My notation                                                       *
C      This is the four dimensional result for                         *
C      Quark antiquark annihilation in order alfa_s^3                  *
C                                                                      *
C      q(P1) + qbar(P2) --> Q(-P3) + Qbar(-P4)+g(-P5)                  *
c  Taken from                                                          *
C  %\cite{Ellis:1986ef}                                                *
c  \bibitem{Ellis:1986ef}                                              *
c  R.~K.~Ellis and J.~C.~Sexton,                                       *
c  %``Explicit Formulae For Heavy Flavor Production,''                 *
c  Nucl.\ Phys.\ B {\bf 282}, 642 (1987).                              *
c  %%CITATION = NUPHA,B282,642;%%                                      *
************************************************************************

      include 'constants.f'
      include 'qcdcouple.f'
      include 'ptilde.f'
      include 'masses.f'
      include 'msq_cs.f'
      include 'sprods_com.f'
      include 'qqgg.f'
      
      integer i,j,k,n2,n3
      logical first
      double precision msq(-nf:nf,-nf:nf),msqc(maxd,-nf:nf,-nf:nf),
     . p(mxpart,4)
      double precision gglo,qqb,wtqqb,wtqbq,wtgg,X,Y,
     . t1,t2,ro,mass2,width2,mass3,width3,
     . eik(4,4),masssq,feik1,feik2,feik3,feikqqb,feikqbq
      double precision 
     & msq15_2(fn:nf,fn:nf),msq25_1(fn:nf,fn:nf),
     & msq15_3(fn:nf,fn:nf),msq15_4(fn:nf,fn:nf),
     & msq35_1(fn:nf,fn:nf),msq45_1(fn:nf,fn:nf),
     & msq25_3(fn:nf,fn:nf),msq25_4(fn:nf,fn:nf),
     & msq35_2(fn:nf,fn:nf),msq45_2(fn:nf,fn:nf),
     & msq35_4(fn:nf,fn:nf),msq45_3(fn:nf,fn:nf)

      double precision
     & m15_2(0:2,fn:nf,fn:nf),m25_1(0:2,fn:nf,fn:nf),
     & m15_3(0:2,fn:nf,fn:nf),m35_1(0:2,fn:nf,fn:nf),
     & m15_4(0:2,fn:nf,fn:nf),m45_1(0:2,fn:nf,fn:nf),
     & m25_3(0:2,fn:nf,fn:nf),m35_2(0:2,fn:nf,fn:nf),
     & m25_4(0:2,fn:nf,fn:nf),m45_2(0:2,fn:nf,fn:nf),

     & m35_4(0:2,fn:nf,fn:nf),m45_3(0:2,fn:nf,fn:nf)

      double precision sub15_2(4),sub25_1(4),
     .                 sub35_1(4),sub15_3(4),
     .                 sub45_1(4),sub15_4(4),
     .                 sub35_2(4),sub25_3(4),
     .                 sub45_2(4),sub25_4(4),
     .                 sub35_4(4),sub45_3(4)
      common/breit/n2,n3,mass2,width2,mass3,width3
      data first/.true./

      
      if (first) then
      first=.false.
      write(6,*) 'mass2',mass2
      endif 
      masssq=mass2**2

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call dotem(5,p,s)

      eik(1,2)=+2d0*gsq*s(1,2)/(s(1,5)+s(2,5))/s(1,5)
      eik(2,1)=+2d0*gsq*s(1,2)/(s(1,5)+s(2,5))/s(2,5)

      eik(1,3)=+2d0*gsq*s(1,3)/(-s(1,5)+s(3,5))/s(1,5)
      eik(3,1)=-2d0*gsq*s(1,3)/(-s(1,5)+s(3,5))/s(3,5)

      eik(1,4)=+2d0*gsq*s(1,4)/(-s(1,5)+s(4,5))/s(1,5)
      eik(4,1)=-2d0*gsq*s(1,4)/(-s(1,5)+s(4,5))/s(4,5)

      eik(2,3)=+2d0*gsq*s(2,3)/(-s(2,5)+s(3,5))/s(2,5)

      eik(3,2)=-2d0*gsq*s(2,3)/(-s(2,5)+s(3,5))/s(3,5)

      eik(2,4)=+2d0*gsq*s(2,4)/(-s(2,5)+s(4,5))/s(2,5)
      eik(4,2)=-2d0*gsq*s(2,4)/(-s(2,5)+s(4,5))/s(4,5)

      eik(3,4)=+2d0*gsq*s(3,4)/(s(3,5)+s(4,5))/s(3,5)
      eik(4,3)=+2d0*gsq*s(3,4)/(s(3,5)+s(4,5))/s(4,5)

      eik(1,1)=0d0
      eik(2,2)=0d0
      eik(3,3)=gsq*masssq*(2d0/s(3,5))**2
      eik(4,4)=gsq*masssq*(2d0/s(4,5))**2

      t1=-s(1,3)/s(1,2)
      t2=-s(2,3)/s(1,2)
      ro=4d0*mass2**2/s(1,2)
      qqb=4d0/9d0*(t1**2+t2**2+ro/2d0)

      msq_cs(1,0,0)=V*xn*avegg*gsq**2
     . *(-2d0*(1d0+ro)*(1d0-1d0/t2)-4d0*t1**2-0.5d0*(ro/t2)**2)
      msq_cs(2,0,0)=V*xn*avegg*gsq**2
     . *(-2d0*(1d0+ro)*(1d0-1d0/t1)-4d0*t2**2-0.5d0*(ro/t1)**2)
      msq_cs(0,0,0)=-V*xn*avegg/xn**2*2d0*gsq**2
     . *(-2d0+(1d0+ro*(1d0-0.5d0*ro))/t1/t2
     .         -0.25d0*(ro/t1)**2-0.25d0*(ro/t2)**2)


      gglo=+msq_cs(1,0,0)+msq_cs(2,0,0)+msq_cs(0,0,0)

      qqb=4d0/9d0*gsq**2*(t1**2+t2**2+ro/2d0)



      msq15_2(1,-1)=qqb
      msq25_1(1,-1)=qqb

      msq15_3(1,-1)=qqb
      msq15_4(1,-1)=qqb
      msq25_3(1,-1)=qqb
      msq25_4(1,-1)=qqb

      msq35_1(1,-1)=qqb
      msq35_2(1,-1)=qqb
      msq45_1(1,-1)=qqb
      msq45_2(1,-1)=qqb

      msq45_3(1,-1)=qqb
      msq35_4(1,-1)=qqb


      msq15_2(-1,1)=qqb
      msq25_1(-1,1)=qqb

      msq15_3(-1,1)=qqb
      msq15_4(-1,1)=qqb
      msq25_3(-1,1)=qqb
      msq25_4(-1,1)=qqb

      msq35_1(-1,1)=qqb
      msq35_2(-1,1)=qqb
      msq45_1(-1,1)=qqb
      msq45_2(-1,1)=qqb

      msq45_3(-1,1)=qqb
      msq35_4(-1,1)=qqb

      msq15_2(0,0)=gg
      msq25_1(0,0)=gg

      msq15_3(0,0)=gg
      msq15_4(0,0)=gg
      msq25_3(0,0)=gg
      msq25_4(0,0)=gg

      msq35_1(0,0)=gg
      msq35_2(0,0)=gg
      msq45_1(0,0)=gg
      msq45_2(0,0)=gg

      msq45_3(0,0)=gg
      msq35_4(0,0)=gg

      do i=0,3
      m15_2(i,0,0)=msq_cs(i,0,0)
      m25_1(i,0,0)=msq_cs(i,0,0)
      m35_1(i,0,0)=msq_cs(i,0,0)
      m35_2(i,0,0)=msq_cs(i,0,0)
      m45_1(i,0,0)=msq_cs(i,0,0)
      m45_2(i,0,0)=msq_cs(i,0,0)
      m15_3(i,0,0)=msq_cs(i,0,0)
      m25_3(i,0,0)=msq_cs(i,0,0)
      m15_4(i,0,0)=msq_cs(i,0,0)
      m25_4(i,0,0)=msq_cs(i,0,0)
      m35_4(i,0,0)=msq_cs(i,0,0)
      m45_3(i,0,0)=msq_cs(i,0,0)
      enddo

      sub15_2(qq)=eik(1,2)
      sub25_1(qq)=eik(2,1)

      sub15_3(qq)=eik(1,3)
      sub35_1(qq)=eik(3,1)-0.5d0*eik(3,3)
      sub15_4(qq)=eik(1,4)
      sub45_1(qq)=eik(4,1)-0.5d0*eik(4,4)

      sub25_3(qq)=eik(2,3)
      sub35_2(qq)=eik(3,2)-0.5d0*eik(3,3)
      sub25_4(qq)=eik(2,4)
      sub45_2(qq)=eik(4,2)-0.5d0*eik(4,4)

      sub35_4(qq)=eik(3,4)-0.5d0*eik(3,3)
      sub45_3(qq)=eik(4,3)-0.5d0*eik(4,4)


      sub15_2(gg)=eik(1,2)
      sub25_1(gg)=eik(2,1)

      sub15_3(gg)=eik(1,3)
      sub35_1(qq)=eik(3,1)-0.5d0*eik(3,3)
      sub15_4(gg)=eik(1,4)
      sub45_1(qq)=eik(4,1)-0.5d0*eik(4,4)

      sub25_3(gg)=eik(2,3)
      sub35_2(qq)=eik(3,2)-0.5d0*eik(3,3)
      sub25_4(gg)=eik(2,4)
      sub45_2(qq)=eik(4,2)-0.5d0*eik(4,4)

      sub35_4(qq)=eik(3,4)-0.5d0*eik(3,3)
      sub45_3(qq)=eik(4,3)-0.5d0*eik(4,4)

      

c      write(6,*) 'eik:sub15_2(gg)',sub15_2(gg)
c      write(6,*) 'eik:sub25_1(gg)',sub25_1(gg)

c      write(6,*) 'eik:sub15_3(gg)',sub15_3(gg)
c      write(6,*) 'eik:sub15_4(gg)',sub15_4(gg)
c      write(6,*) 'eik:sub25_3(gg)',sub25_3(gg)
c      write(6,*) 'eik:sub25_4(gg)',sub25_4(gg)

      write(6,*) 'eik:sub15_2(qq)',sub15_2(qq)
      write(6,*) 'eik:sub25_1(qq)',sub25_1(qq)



c      write(6,*) 'eik:sub35_1(qq)',sub35_1(qq)
c      write(6,*) 'eik:sub35_2(qq)',sub35_2(qq)
c      write(6,*) 'eik:sub45_1(qq)',sub45_1(qq)
c      write(6,*) 'eik:sub45_2(qq)',sub45_2(qq)

c      write(6,*) 'eik:sub35_4(qq)',sub35_4(qq)
c      write(6,*) 'eik:sub45_3(qq)',sub45_3(qq)

      write(6,*) 'eik:(3,3)',eik(3,3)
      write(6,*) 'eik:(4,4)',eik(4,4)

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      k=-j
      if ((j .eq. 0) .and. (k.eq.0)) then
c      wtgg=
c     .  + nab1*xn* ( eik(1,2) + eik(1,4) + eik(2,3) + eik(3,4) )
c     .  + nab2*xn* ( eik(1,2) + eik(1,3) + eik(2,4) + eik(3,4) )
c     .  + gg*xn*( - half*eik(3,3) - half*eik(4,4) - eik(3,4) )
c     .  + gg/xn * ( half*eik(3,3) + half*eik(4,4) - eik(3,4) )
c     .  + ab*xn* ( eik(1,3) + eik(1,4) + eik(2,3) + eik(2,4) ) 

      msqc(1,j,k)= (+m15_2(1,j,k)+m15_2(2,j,k))*sub15_2(gg)*xn
      msqc(2,j,k)= (+m25_1(1,j,k)+m25_1(2,j,k))*sub25_1(gg)*xn

      msqc(3,j,k)= (+m15_3(2,j,k)+m15_3(0,j,k))*sub15_3(gg)*xn
      msqc(7,j,k)= (+m35_1(2,j,k)+m35_1(0,j,k))*sub35_1(qq)*xn
      msqc(5,j,k)= (+m15_4(1,j,k)+m15_4(0,j,k))*sub15_4(gg)*xn
      msqc(9,j,k)= (+m45_1(1,j,k)+m45_1(0,j,k))*sub45_1(qq)*xn

      msqc(4,j,k)= (+m25_3(1,j,k)+m25_3(0,j,k))*sub25_3(gg)*xn
      msqc(8,j,k)= (+m35_2(1,j,k)+m35_2(0,j,k))*sub35_2(qq)*xn
      msqc(6,j,k)= (+m25_4(2,j,k)+m25_4(0,j,k))*sub25_4(gg)*xn
      msqc(10,j,k)=(+m45_2(2,j,k)+m45_2(0,j,k))*sub45_2(qq)*xn

      msqc(11,j,k)=-(m35_4(2,j,k)+m35_4(0,j,k))*sub35_4(qq)
     . *(xn+1d0/xn)
      msqc(12,j,k)=-(m45_3(2,j,k)+m45_3(0,j,k))*sub45_3(qq)
     . *(xn+1d0/xn)

      msq(j,k)=
     . +msqc(1,j,k)+msqc(2,j,k)+msqc(3,j,k)+msqc(4,j,k)
     . +msqc(5,j,k)+msqc(6,j,k)+msqc(7,j,k)+msqc(8,j,k)
     . +msqc(9,j,k)+msqc(10,j,k)+msqc(11,j,k)+msqc(12,j,k)

      elseif ((j .gt. 0) .and. (k.lt.0)) then
      msqc(1,j,k)= -msq15_2(j,k)*sub15_2(qq)/xn
      msqc(2,j,k)= -msq25_1(j,k)*sub25_1(qq)/xn

      msqc(3,j,k)= +msq15_3(j,k)*sub15_3(qq)*(xn-2d0/xn)
      msqc(4,j,k)= +msq25_3(j,k)*sub25_3(qq)*2d0/xn
      msqc(5,j,k)= +msq15_4(j,k)*sub15_4(qq)*2d0/xn
      msqc(6,j,k)= +msq25_4(j,k)*sub25_4(qq)*(xn-2d0/xn)

      msqc(7,j,k)= +msq35_1(j,k)*sub35_1(qq)*(xn-2d0/xn)
      msqc(8,j,k)= +msq35_2(j,k)*sub35_2(qq)*2d0/xn
      msqc(9,j,k)= +msq45_1(j,k)*sub45_1(qq)*2d0/xn
      msqc(10,j,k)=+msq45_2(j,k)*sub45_2(qq)*(xn-2d0/xn)

      msqc(11,j,k)=-msq35_4(j,k)*sub35_4(qq)/xn
      msqc(12,j,k)=-msq45_3(j,k)*sub45_3(qq)/xn
      msq(j,k)=
     . +msqc(1,j,k)+msqc(2,j,k)+msqc(3,j,k)+msqc(4,j,k)
     . +msqc(5,j,k)+msqc(6,j,k)+msqc(7,j,k)+msqc(8,j,k)
     . +msqc(9,j,k)+msqc(10,j,k)+msqc(11,j,k)+msqc(12,j,k)

      write(6,*) 'eik',msq(1,-1)
      elseif ((j .lt. 0) .and. (k.gt.0)) then
      msqc(1,j,k)= -msq15_2(j,k)*sub15_2(qq)/xn
      msqc(2,j,k)= -msq25_1(j,k)*sub25_1(qq)/xn

      msqc(3,j,k)= +msq25_3(j,k)*sub25_3(qq)*(xn-two/xn)
      msqc(4,j,k)= +msq15_3(j,k)*sub15_3(qq)*two/xn
      msqc(5,j,k)= +msq25_4(j,k)*sub25_4(qq)*two/xn
      msqc(6,j,k)= +msq15_4(j,k)*sub15_4(qq)*(xn-two/xn)

      msqc(7,j,k)= +msq35_2(j,k)*sub35_2(qq)*(xn-two/xn)
      msqc(8,j,k)= +msq35_1(j,k)*sub35_1(qq)*two/xn
      msqc(9,j,k)= +msq45_2(j,k)*sub45_2(qq)*two/xn
      msqc(10,j,k)=+msq45_1(j,k)*sub45_1(qq)*(xn-two/xn)

      msqc(11,j,k)=-msq35_4(j,k)*sub35_4(qq)/xn
      msqc(12,j,k)=-msq45_3(j,k)*sub45_3(qq)/xn
      msq(j,k)=
     . +msqc(1,j,k)+msqc(2,j,k)+msqc(3,j,k)+msqc(4,j,k)
     . +msqc(5,j,k)+msqc(6,j,k)+msqc(7,j,k)+msqc(8,j,k)
     . +msqc(9,j,k)+msqc(10,j,k)+msqc(11,j,k)+msqc(12,j,k)

      endif
      enddo

      return
      end
