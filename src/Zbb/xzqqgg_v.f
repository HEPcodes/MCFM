      subroutine xzqqgg_v(mqqb,mqqb_vec,mqqb_ax)
      implicit none
C     Author J.M.Campbell, February 2000
C     Returns the interference of the tree and loop
C     amplitudes for the process
C     0---> q(p1)+g(p2)+g(p3)+qbar(p4)+l(p5)+a(p6)
C     mqqb(2,2) has two indices;the first for the helicity quark line;
C     the second for helicity of lepton line.
      include 'constants.f'
      include 'prods.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      integer k,nu,iperm
      integer i1(2),i2(2),i3(2),i4(2),i5(2),i6(2),j,lh,h2,h3,hq,h(2:3)
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),msqzbb,mmsq(2,2),
     . prop,qqb,qbq,pswap(mxpart,4)
      double precision mqqb(2,2),mqqb_vec(2,2),mqqb_ax(2,2),m1,m0,x,fac
      double complex tempm0,m(2),ml1(2),ml2(2),ml3,ml4(2),tamp
      double complex ml_vec(2),ml1_ax(2),ml2_ax(2)
      double complex amp_qqgg,a6treeg,a6treeg1,
     . a61g1lc,a61g1slc,a61g1nf,a63g1,a64v,a64ax,a65ax
      character*9 st1(2,2),st2(2,2),st3(2,2)
      logical compare
      parameter(x=xn/cf)
      data i1/1,4/
      data i2/2,3/
      data i3/3,2/
      data i4/4,1/
      data i5/6,5/
      data i6/5,6/
      data st1/'q+g-g-qb-','q+g-g+qb-','q+g+g-qb-','q+g+g+qb-'/
      data st2/'q+qb-g+g+','q+qb-g+g-','q+qb-g-g+','q+qb-g-g-'/
      data st3/'q+qb-g-g-','q+qb-g-g+','q+qb-g+g-','q+qb-g+g+'/
      
      compare=.false.
      
C ---final matrix element squared is needed as function of quark line helicity
C----and lepton line helicity
C----first argument is quark line helicity
C----second argument is lepton line helicity
      mqqb(1,1)=0d0
      mqqb(1,2)=0d0
      mqqb(2,1)=0d0
      mqqb(2,2)=0d0

      fac=avegg*8d0*gsq**2*esq**2*cf*xn**3*ason2pi
c--- no extra factor here since colour algebra is already done in (2.12)

      do hq=1,2
      do lh=1,2
      mqqb(hq,lh)=0d0
      mqqb_vec(hq,lh)=0d0
      mqqb_ax(hq,lh)=0d0
c--- set iperm=1 for same-helicity quarks and leptons, 2 otherwise
      iperm=1+abs(hq-lh)

      do h2=1,2
      do h3=1,2
        h(2)=h2
        h(3)=h3
        do j=1,2
        if (hq .eq. 1) then
        m(j)=  a6treeg1(st1(3-h(i2(j)),3-h(i3(j))),
     .     i1(1),i2(j),i3(j),i4(1),i6(lh),i5(lh),zb,za)
        ml1(j)= a61g1lc(st1(3-h(i2(j)),3-h(i3(j))),
     .     i1(1),i2(j),i3(j),i4(1),i6(lh),i5(lh),zb,za)
        ml2(j)=a61g1slc(st2(3-h(i2(j)),3-h(i3(j))),
     .     i1(1),i2(j),i3(j),i4(1),i6(lh),i5(lh),zb,za)
        ml4(j)=a61g1nf(st2(3-h(i2(j)),3-h(i3(j))),
     .     i1(1),i2(j),i3(j),i4(1),i6(lh),i5(lh),zb,za)
        ml_vec(j)=a64v(st2(3-h(i2(j)),3-h(i3(j))),
     .     i1(1),i4(1),i2(j),i3(j),i6(lh),i5(lh),zb,za)
        ml1_ax(j)=a64ax(st2(3-h(i2(j)),3-h(i3(j))),
     .     i1(1),i4(1),i2(j),i3(j),i6(lh),i5(lh),zb,za)
        ml2_ax(j)=a65ax(st2(3-h(i2(j)),3-h(i3(j))),
     .     i1(1),i4(1),i2(j),i3(j),i6(lh),i5(lh),zb,za)
        else
        m(j)=  a6treeg1(st1(h(i2(j)),h(i3(j))),
     .     i1(1),i2(j),i3(j),i4(1),i5(lh),i6(lh),za,zb)
        ml1(j)= a61g1lc(st1(h(i2(j)),h(i3(j))),
     .     i1(1),i2(j),i3(j),i4(1),i5(lh),i6(lh),za,zb)
        ml2(j)=a61g1slc(st2(h(i2(j)),h(i3(j))),
     .     i1(1),i2(j),i3(j),i4(1),i5(lh),i6(lh),za,zb)
        ml4(j)=a61g1nf(st2(h(i2(j)),h(i3(j))),
     .     i1(1),i2(j),i3(j),i4(1),i5(lh),i6(lh),za,zb)
        ml_vec(j)=a64v(st2(h(i2(j)),h(i3(j))),
     .     i1(1),i4(1),i2(j),i3(j),i5(lh),i6(lh),za,zb)
        ml1_ax(j)=a64ax(st2(h(i2(j)),h(i3(j))),
     .     i1(1),i4(1),i2(j),i3(j),i5(lh),i6(lh),za,zb)
        ml2_ax(j)=a65ax(st2(h(i2(j)),h(i3(j))),
     .     i1(1),i4(1),i2(j),i3(j),i5(lh),i6(lh),za,zb)
        endif
        enddo
        if (hq .eq. 1) then
        ml3=a63g1(st3(3-h2,3-h3),1,4,2,3,i6(lh),i5(lh),zb,za)
        else
        ml3=a63g1(st3(h2,h3),1,4,2,3,i5(lh),i6(lh),za,zb)
        endif

      mqqb(hq,lh)=mqqb(hq,lh)+fac*(
     .  Dble(Dconjg(m(1))*(
     .    ml1(1)
     .   -(ml1(1)+ml2(1)+ml1(2)-ml3)/xnsq
     .   +(ml2(1)+ml2(2))/xnsq**2
     .   +ml4(1)/xn
     .   -(ml4(1)+ml4(2))/xn**3))
     . +Dble(Dconjg(m(2))*(
     .    ml1(2)
     .   -(ml1(2)+ml2(2)+ml1(1)-ml3)/xnsq
     .   +(ml2(1)+ml2(2))/xnsq**2
     .   +ml4(2)/xn
     .   -(ml4(1)+ml4(2))/xn**3)))
     
      mqqb_vec(hq,lh)=mqqb_vec(hq,lh)+fac/xnsq*(
     .  Dble(Dconjg(m(1))*(
     .    (xn-4d0/xn)*ml_vec(1)))
     . +Dble(Dconjg(m(2))*(
     .    (xn-4d0/xn)*ml_vec(2))))

      mqqb_ax(hq,lh)=mqqb_ax(hq,lh)+fac/xnsq*(
     .  Dble(Dconjg(m(1))*(
     .    (xn-2d0/xn)*ml1_ax(1)-2d0/xn*ml1_ax(2)+one/xn*ml2_ax(1)))
     . +Dble(Dconjg(m(2))*(
     .    (xn-2d0/xn)*ml1_ax(2)-2d0/xn*ml1_ax(1)+one/xn*ml2_ax(2))))

      enddo
      enddo
      enddo
      enddo
         
      return
      end
      
      double complex function a61g1lc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
c----wrapper to a61g that also includes config st='q+g-g-qb-'
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      character*9 st
      double complex a61gcol

      if(st.eq.'q+g-g-qb-') then
        a61g1lc=a61gcol('q+g+g+qb-',j4,j3,j2,j1,j6,j5,zb,za,1)
      else
        a61g1lc=a61gcol(st,j1,j2,j3,j4,j5,j6,za,zb,1)
      endif
      
      return
      end

      double complex function a61g1slc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
c----wrapper to a61g that also includes config st='q+g-g-qb-'
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      character*9 st
      double complex a61gcol

      if(st.eq.'q+qb-g-g-') then
        a61g1slc=a61gcol('q+qb-g+g+',j1,j2,j3,j4,j5,j6,za,zb,2)
      else
        a61g1slc=a61gcol(st,j4,j3,j2,j1,j6,j5,zb,za,2)
      endif
      
      return
      end

      double complex function a61g1nf(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
c----wrapper to a61g that also includes config st='q+g-g-qb-'
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      character*9 st
      double complex a61gcol

      if(st.eq.'q+g-g-qb-') then
        a61g1nf=a61gcol('q+g+g+qb-',j4,j3,j2,j1,j6,j5,zb,za,3)
      else
        a61g1nf=a61gcol(st,j1,j2,j3,j4,j5,j6,za,zb,3)
      endif
      
      return
      end

      double complex function a61gcol(st,j1,j2,j3,j4,j5,j6,za,zb,ncol)
      implicit none
C---hep-ph/9708239, Eqn 2.13
      include 'constants.f'
      include 'sprodx.f'
      integer j1,j2,j3,j4,j5,j6,ncol
      character*9 st
      double complex a6g,a6sg,a6fg,a6tg

c--- use ncol=1 for leading colour piece, ncol=2 for subleading
      if     (ncol .eq. 1) then
c--- comes with natural colour factor (1)
      a61gcol=
     . +a6g(st,j1,j2,j3,j4,j5,j6,za,zb)
      elseif (ncol .eq. 2) then
c--- comes with natural colour factor (-1/xnsq)
      a61gcol=
     . +a6g(st,j1,j4,j3,j2,j5,j6,za,zb)
      elseif (ncol .eq. 3) then
c--- comes with natural colour factor (1/xn)
      a61gcol=
     . -dfloat(nf)*(a6sg(st,j1,j2,j3,j4,j5,j6,za,zb)
     .             +a6fg(st,j1,j2,j3,j4,j5,j6,za,zb))
     . +a6tg(st,j1,j2,j3,j4,j5,j6,za,zb)
      endif
      return
      end

      double complex function a63g1(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
c----wrapper to a63g that also includes config st='q+qb-g-g-'
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      character*9 st
      double complex a63g

      if(st.eq.'q+qb-g-g-') then
        a63g1=a63g('q+qb-g+g+',j4,j1,j2,j3,j6,j5,zb,za)
      else
        a63g1=a63g(st,j1,j4,j2,j3,j5,j6,za,zb)
      endif
      
      return
      end

      double complex function a64v(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
c----definition (2.13) of BDK, writes in terms of fvs and fvf
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      character*9 st
      double complex fvs,fvf

      if     (st.eq.'q+qb-g-g-') then
        a64v=-fvs('q+qb-g+g+',j4,j1,j2,j3,j6,j5,zb,za)
     .       -fvf('q+qb-g+g+',j4,j1,j2,j3,j6,j5,zb,za)
      elseif (st.eq.'q+qb-g-g+') then
        a64v=-fvs('q+qb-g+g-',j1,j4,j3,j2,j5,j6,za,zb)
     .       -fvf('q+qb-g+g-',j1,j4,j3,j2,j5,j6,za,zb)
      else
        a64v=-fvs(st,j1,j4,j2,j3,j5,j6,za,zb)
     .       -fvf(st,j1,j4,j2,j3,j5,j6,za,zb)
      endif
      
      return
      end

      double complex function a64ax(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
c----definition (2.13) of BDK
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      character*9 st
      double complex fax

      if (st.eq.'q+qb-g-g-') then
        a64ax=fax('q+qb-g+g+',j4,j1,j2,j3,j6,j5,zb,za)
      else
        a64ax=fax(st,j1,j4,j2,j3,j5,j6,za,zb)
      endif
      
      return
      end

      double complex function a65ax(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
c----definition (2.13) of BDK
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      character*9 st
      double complex faxsl

      if     (st.eq.'q+qb-g-g-') then
        a65ax=faxsl('q+qb-g+g+',j4,j1,j2,j3,j6,j5,zb,za)
      elseif (st.eq.'q+qb-g-g+') then
        a65ax=faxsl('q+qb-g+g-',j1,j4,j3,j2,j5,j6,za,zb)
      else
        a65ax=faxsl(st,j1,j4,j2,j3,j5,j6,za,zb)
      endif
      
      return
      end

