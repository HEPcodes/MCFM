      subroutine qqb_QQB_soft(p,msq) 
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
CCCCCCCCCC  UNchecked for the moment*********************
C           UNchecked for the moment*********************

      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      
      integer j,k,n2,n3
      logical first
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision qqb,gg,wtqqb,wtqbq,wtgg,
     . t1,t2,ro,mass2,width2,mass3,width3,nab1,nab2,ab,
     . eik(4,4),masssq,feik1,feik2,feik3,feikqqb,feikqbq
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


      eik(1,2)=+2d0*gsq*s(1,2)/(s(1,5)*s(2,5))
      eik(1,3)=+2d0*gsq*s(1,3)/(s(1,5)*s(3,5))
      eik(1,4)=+2d0*gsq*s(1,4)/(s(1,5)*s(4,5))
      eik(2,3)=+2d0*gsq*s(2,3)/(s(2,5)*s(3,5))
      eik(2,4)=+2d0*gsq*s(2,4)/(s(2,5)*s(4,5))
      eik(3,4)=+2d0*gsq*s(3,4)/(s(3,5)*s(4,5))
      eik(1,1)=0d0
      eik(2,2)=0d0
      eik(3,3)=gsq*masssq*(2d0/s(3,5))**2
      eik(4,4)=gsq*masssq*(2d0/s(4,5))**2

c      write(6,*) 'eik(1,2)',eik(1,2)
c      write(6,*) 'eik(1,3)',eik(1,3)
c      write(6,*) 'eik(2,3)',eik(2,3)
c      write(6,*) 'eik(1,4)',eik(1,4)
c      write(6,*) 'eik(2,4)',eik(2,4)
c      write(6,*) 'eik(3,4)',eik(3,4)



c      feik1=CF*(eik(1,3)+eik(1,4)+eik(2,3)+eik(2,4)
c     . -2d0*eik(1,2)-eik(3,3)-eik(4,4))+2d0*xn*eik(1,2)
c      feik2=(eik(1,4)+eik(2,3)-eik(1,3)-eik(2,4))
c      feik3=(2d0*(eik(1,2)+eik(3,4))
c     .  -eik(1,3)-eik(2,4)-eik(1,4)-eik(2,3))

      feikqqb=CF*(+2d0*eik(1,3)+2d0*eik(2,4)-eik(3,3)-eik(4,4))
     . +(2d0*eik(1,4)+2d0*eik(2,3)
     . -eik(1,3)-eik(2,4)-eik(1,2)-eik(3,4))/xn
      feikqbq=CF*(+2d0*eik(2,3)+2d0*eik(1,4)-eik(3,3)-eik(4,4))
     . +(2d0*eik(2,4)+2d0*eik(1,3)
     . -eik(2,3)-eik(1,4)-eik(1,2)-eik(3,4))/xn

c      X=gsq**2*xn**2/(4d0*V)*(
c     . -(s(1,2)+4d0*masssq)/s(2,3)+(s(1,2)+4d0*masssq)/s(1,3)
c     . +(2d0*masssq/s(1,3))**2-(2d0*masssq/s(2,3))**2
c     .  -2d0*(s(2,3)-s(1,3))/s(1,2))

c      Y=gsq**2/xn**2*((0.5d0*s(1,3))**2+(0.5d0*s(2,3))**2+masssq*s(1,2)
c     .  -(masssq*s(1,2))**2/(s(1,3)*s(2,3)))/(4d0*V)
c     . *(4d0/(s(1,3)*s(2,3))+8d0*xn**2/s(1,2)**2)
c      write(6,*) 'X',X
c      write(6,*) 'Y',Y

      t1=-s(1,3)/s(1,2)
      t2=-s(2,3)/s(1,2)
      ro=4d0*mass2**2/s(1,2)

c      X=+gsq**2*xn*(1d0/t2+ro/t2-1d0/t1-ro/t1
c     .  +(ro/2d0/t1)**2-(ro/2d0/t2)**2
c     .  +2d0*(t2-t1))
c
c      Y=gsq**2/xn*(t1**2+t2**2+ro-ro**2/(4d0*t1*t2))
c     . *(1d0/(xn**2*t1*t2)+2d0)

c      X=xn/(4d0*V)*X
c      Y=xn/(4d0*V)*Y
c      write(6,*) 'X',X
c      write(6,*) 'Y',Y


      nab1=V*xn*avegg*gsq**2
     . *(-2d0*(1d0+ro)*(1d0-1d0/t2)-4d0*t1**2-0.5d0*(ro/t2)**2)
      nab2=V*xn*avegg*gsq**2
     . *(-2d0*(1d0+ro)*(1d0-1d0/t1)-4d0*t2**2-0.5d0*(ro/t1)**2)
      ab=-V*xn*avegg/xn**2*2d0*gsq**2
     . *(-2d0+(1d0+ro*(1d0-0.5d0*ro))/t1/t2
     .         -0.25d0*(ro/t1)**2-0.25d0*(ro/t2)**2)
      wtgg=nab1+nab2+ab
c      X=0.5d0*(nab1-nab2)*xn
c      Y=0.5d0*(-nab1-nab2-(xn**2+1d0)*ab)/xn
c      write(6,*) 'wtgg',wtgg

c      write(6,*) 'xn*nab1',xn*nab1
c      write(6,*) 'xn*nab2',xn*nab2
c      write(6,*) 'X',X
c      write(6,*) 'Y',Y
c      pause

      qqb=4d0/9d0*gsq**2*(t1**2+t2**2+ro/2d0)
      gg=gsq**2*(1d0/(6d0*t1*t2)-3d0/8d0)
     . *(t1**2+t2**2+ro-0.25d0*ro**2/(t1*t2))
c      write(6,*) 'wtgg',wtgg

c      gg=(wtgg*feik1+X*feik2+Y*feik3)

c      write(6,*) '1,gg',gg

 
      gg=
     .  + nab1*xn* ( eik(1,2) + eik(1,4) + eik(2,3) + eik(3,4) )
     .  + nab2*xn* ( eik(1,2) + eik(1,3) + eik(2,4) + eik(3,4) )
     .  + wtgg*xn*( - half*eik(3,3) - half*eik(4,4) - eik(3,4) )
     .  + wtgg/xn * ( half*eik(3,3) + half*eik(4,4) - eik(3,4) )
     .  + ab*xn* ( eik(1,3) + eik(1,4) + eik(2,3) + eik(2,4) ) 

c      write(6,*) '2,gg',gg



      wtqqb=(qqb*feikqqb)
      wtqbq=(qqb*feikqbq)

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      k=-j
      if ((j .eq. 0) .and. (k.eq.0)) then
          msq(j,k)=gg
      elseif ((j .gt. 0) .and. (k.lt.0)) then
          msq(j,k)=wtqqb
      elseif ((j .lt. 0) .and. (k.gt.0)) then
          msq(j,k)=wtqbq
      endif
      enddo

      return
      end
