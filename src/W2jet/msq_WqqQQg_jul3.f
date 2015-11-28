      subroutine msq_WqqQQg(i1,i2,i3,i4,i5,i6,i7,MN,ofac)
************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2001.                                                     *
*     Return matrix elements squared as a function of f1,f2,f3,f4      *
*     summed over helicity, using the formulae of                      *
*     Nagy and Trocsanyi, PRD59 014020 (1999)                          *
************************************************************************
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'prods.f'
      include 'hardscale.f'
      integer Qh,hq,hg,f1,f2,f3,f4,i1,i2,i3,i4,i5,i6,i7,j
      double precision A(5,5,5,5),B(5,5,5,5),C(5,5,5,5),D(5,5,5,5),
     . E(5,5,5,5),F(5,5,5,5),G(5,5,5,5)
      double precision MN(5,5,5,5),
     . M0(5,5,5,5),Mx(5,5,5,5),My(5,5,5,5),
     . Mz(5,5,5,5),Mxx(5,5,5,5),Mxy(5,5,5,5),gq,gl,Qg,temp,ofac
      double precision x,y,z
      parameter(x=xn/cf,y=half/cf,z=0.25d0*(xn**2-two)/xn/cf**2)
      double complex prop
      double complex 
     .               mb1_1234(5,5,5,5,2,2,2),mb2_1234(5,5,5,5,2,2,2),
     .               mb1_3412(5,5,5,5,2,2,2),mb2_3412(5,5,5,5,2,2,2),
     .               mb1_3214(5,5,5,5,2,2,2),mb2_3214(5,5,5,5,2,2,2),
     .               mb1_1432(5,5,5,5,2,2,2),mb2_1432(5,5,5,5,2,2,2)
      

C---set everything to zero 
      do f1=1,5
      do f2=1,5
      do f3=1,5
      do f4=1,5
      MN(f1,f2,f3,f4)=zip
      A(f1,f2,f3,f4)=zip
      B(f1,f2,f3,f4)=zip
      C(f1,f2,f3,f4)=zip
      D(f1,f2,f3,f4)=zip
      E(f1,f2,f3,f4)=zip
      F(f1,f2,f3,f4)=zip
      G(f1,f2,f3,f4)=zip
      enddo
      enddo
      enddo
      enddo
      if (s(i6,i7) .lt. 4*mbsq) return
      
C---exclude the photon pole, 4*mbsq choosen as a scale approx above upsilon 

c---mb1_1234 etc, has 6 indices each with possible values 1 or 2
C---corresponding to f1,f3,hq,Qh,hg,lh
      call wmakemb(i1,i2,i3,i4,i5,i6,i7,mb1_1234,mb2_1234)
      call wmakemb(i3,i4,i1,i2,i5,i6,i7,mb1_3412,mb2_3412)
      call wmakemb(i3,i2,i1,i4,i5,i6,i7,mb1_3214,mb2_3214)
      call wmakemb(i1,i4,i3,i2,i5,i6,i7,mb1_1432,mb2_1432)

      do f1=1,5
      do f2=1,5
      do f3=1,5
      do f4=1,5
      
      temp=0d0
      
      do hq=1,2
      do Qh=1,2
      do hg=1,2
      
      A(f1,f2,f3,f4)=A(f1,f2,f3,f4)
     . +cdabs(mb1_1234(f1,f2,f3,f4,hq,Qh,hg))**2
     . +cdabs(mb2_1234(f1,f2,f3,f4,hq,Qh,hg))**2
     . +cdabs(mb1_3412(f3,f4,f1,f2,Qh,hq,hg))**2
     . +cdabs(mb2_3412(f3,f4,f1,f2,Qh,hq,hg))**2

      D(f1,f2,f3,f4)=D(f1,f2,f3,f4)+two*Dreal(
     .+mb1_1234(f1,f2,f3,f4,hq,Qh,hg)
     .*Dconjg(mb2_1234(f1,f2,f3,f4,hq,Qh,hg))
     .+mb1_3412(f3,f4,f1,f2,Qh,hq,hg)
     .*Dconjg(mb2_3412(f3,f4,f1,f2,Qh,hq,hg)))
      
      A(f1,f2,f3,f4)=A(f1,f2,f3,f4)
     . +cdabs(mb1_3214(f3,f2,f1,f4,Qh,hq,hg))**2
     . +cdabs(mb2_3214(f3,f2,f1,f4,Qh,hq,hg))**2
     . +cdabs(mb1_1432(f1,f4,f3,f2,hq,Qh,hg))**2
     . +cdabs(mb2_1432(f1,f4,f3,f2,hq,Qh,hg))**2

      D(f1,f2,f3,f4)=D(f1,f2,f3,f4)+two*Dreal(
     .+mb1_3214(f3,f2,f1,f4,Qh,hq,hg)
     .*Dconjg(mb2_3214(f3,f2,f1,f4,Qh,hq,hg))
     .+mb1_1432(f1,f4,f3,f2,hq,Qh,hg)
     .*Dconjg(mb2_1432(f1,f4,f3,f2,hq,Qh,hg)))

      if (hq .eq. Qh) then
      B(f1,f2,f3,f4)=B(f1,f2,f3,f4)-two*Dreal(
     .+mb1_1234(f1,f2,f3,f4,hq,Qh,hg)
     .*Dconjg(mb1_1432(f1,f4,f3,f2,hq,Qh,hg))
     .+mb2_1234(f1,f2,f3,f4,hq,Qh,hg)
     .*Dconjg(mb2_3214(f3,f2,f1,f4,Qh,hq,hg))
     .+mb1_3412(f3,f4,f1,f2,Qh,hq,hg)
     .*Dconjg(mb1_3214(f3,f2,f1,f4,Qh,hq,hg))
     .+mb2_3412(f3,f4,f1,f2,Qh,hq,hg)
     .*Dconjg(mb2_1432(f1,f4,f3,f2,hq,Qh,hg)))

      C(f1,f2,f3,f4)=C(f1,f2,f3,f4)-two*Dreal(
     .+mb1_1234(f1,f2,f3,f4,hq,Qh,hg)
     .*Dconjg(mb1_3214(f3,f2,f1,f4,Qh,hq,hg))
     .+mb2_1234(f1,f2,f3,f4,hq,Qh,hg)
     .*Dconjg(mb2_1432(f1,f4,f3,f2,hq,Qh,hg))
     .+mb1_3412(f3,f4,f1,f2,Qh,hq,hg)
     .*Dconjg(mb1_1432(f1,f4,f3,f2,hq,Qh,hg))
     .+mb2_3412(f3,f4,f1,f2,Qh,hq,hg)
     .*Dconjg(mb2_3214(f3,f2,f1,f4,Qh,hq,hg)))

      E(f1,f2,f3,f4)=E(f1,f2,f3,f4)-two*Dreal(
     .+(mb1_1234(f1,f2,f3,f4,hq,Qh,hg)
     . +mb1_3412(f3,f4,f1,f2,Qh,hq,hg))
     .*Dconjg(mb2_3214(f3,f2,f1,f4,Qh,hq,hg)
     .+mb2_1432(f1,f4,f3,f2,hq,Qh,hg))
     .+(mb2_1234(f1,f2,f3,f4,hq,Qh,hg)
     . +mb2_3412(f3,f4,f1,f2,Qh,hq,hg))
     .*Dconjg(mb1_3214(f3,f2,f1,f4,Qh,hq,hg)
     . +mb1_1432(f1,f4,f3,f2,hq,Qh,hg)))
      endif

       F(f1,f2,f3,f4)=F(f1,f2,f3,f4)+two*Dreal(
     . +mb1_3214(f3,f2,f1,f4,Qh,hq,hg)
     . *Dconjg(mb1_1432(f1,f4,f3,f2,hq,Qh,hg))
     . +mb2_3214(f3,f2,f1,f4,Qh,hq,hg)
     . *Dconjg(mb2_1432(f1,f4,f3,f2,hq,Qh,hg)))

       G(f1,f2,f3,f4)=G(f1,f2,f3,f4)+two*Dreal(
     . +mb1_3214(f3,f2,f1,f4,Qh,hq,hg)
     . *Dconjg(mb2_1432(f1,f4,f3,f2,hq,Qh,hg))
     . +mb2_3214(f3,f2,f1,f4,Qh,hq,hg)
     . *Dconjg(mb1_1432(f1,f4,f3,f2,hq,Qh,hg)))      

      F(f1,f2,f3,f4)=F(f1,f2,f3,f4)+two*Dreal(
     . +mb1_1234(f1,f2,f3,f4,hq,Qh,hg)
     . *Dconjg(mb1_3412(f3,f4,f1,f2,Qh,hq,hg))
     . +mb2_1234(f1,f2,f3,f4,hq,Qh,hg)
     . *Dconjg(mb2_3412(f3,f4,f1,f2,Qh,hq,hg)))

      G(f1,f2,f3,f4)=G(f1,f2,f3,f4)+two*Dreal(
     . +mb1_1234(f1,f2,f3,f4,hq,Qh,hg)
     . *Dconjg(mb2_3412(f3,f4,f1,f2,Qh,hq,hg))
     . +mb2_1234(f1,f2,f3,f4,hq,Qh,hg)
     . *Dconjg(mb1_3412(f3,f4,f1,f2,Qh,hq,hg)))

c--- DEBUG helicities
      M0(f1,f2,f3,f4)=
     . B(f1,f2,f3,f4)+C(f1,f2,f3,f4)+E(f1,f2,f3,f4)
      Mx(f1,f2,f3,f4)=-0.5d0
     . *(3d0*C(f1,f2,f3,f4)+2d0*E(f1,f2,f3,f4)+B(f1,f2,f3,f4))
      My(f1,f2,f3,f4)=A(f1,f2,f3,f4)+D(f1,f2,f3,f4)
      Mz(f1,f2,f3,f4)=F(f1,f2,f3,f4)+G(f1,f2,f3,f4)
      Mxx(f1,f2,f3,f4)=0.25d0*(2d0*C(f1,f2,f3,f4)+E(f1,f2,f3,f4))
      Mxy(f1,f2,f3,f4)=-0.5d0*(F(f1,f2,f3,f4)+D(f1,f2,f3,f4))

      MN(f1,f2,f3,f4)=
     . ofac*CF**3*xn*(M0(f1,f2,f3,f4)+x*Mx(f1,f2,f3,f4)
     . +y*My(f1,f2,f3,f4)+z*Mz(f1,f2,f3,f4)
     . +x**2*Mxx(f1,f2,f3,f4)+x*y*Mxy(f1,f2,f3,f4))
      
      if ((MN(f1,f2,f3,f4) .ne. 0d0) .and. (f1 .eq. 2) .and.
     . (f2 .eq. 1) .and. (f3 .eq. 2) .and. (f4 .eq. 2)) 
     . write(*,*) '(2,1,2,2)',hq,Qh,hg,MN(f1,f2,f3,f4)-temp
      if ((MN(f1,f2,f3,f4) .ne. 0d0) .and. (f1 .eq. 2) .and.
     . (f2 .eq. 2) .and. (f3 .eq. 2) .and. (f4 .eq. 1)) 
     . write(*,*) '(2,2,2,1)',hq,Qh,hg,MN(f1,f2,f3,f4)-temp
     
      temp=MN(f1,f2,f3,f4)
c--- DEBUG helicities
            
      enddo
      enddo
      enddo
      
      M0(f1,f2,f3,f4)=
     . B(f1,f2,f3,f4)+C(f1,f2,f3,f4)+E(f1,f2,f3,f4)
      Mx(f1,f2,f3,f4)=-0.5d0
     . *(3d0*C(f1,f2,f3,f4)+2d0*E(f1,f2,f3,f4)+B(f1,f2,f3,f4))
      My(f1,f2,f3,f4)=A(f1,f2,f3,f4)+D(f1,f2,f3,f4)
      Mz(f1,f2,f3,f4)=F(f1,f2,f3,f4)+G(f1,f2,f3,f4)
      Mxx(f1,f2,f3,f4)=0.25d0*(2d0*C(f1,f2,f3,f4)+E(f1,f2,f3,f4))
      Mxy(f1,f2,f3,f4)=-0.5d0*(F(f1,f2,f3,f4)+D(f1,f2,f3,f4))

      MN(f1,f2,f3,f4)=
     . CF**3*xn*(M0(f1,f2,f3,f4)+x*Mx(f1,f2,f3,f4)
     . +y*My(f1,f2,f3,f4)+z*Mz(f1,f2,f3,f4)
     . +x**2*Mxx(f1,f2,f3,f4)+x*y*Mxy(f1,f2,f3,f4))
         
      enddo
      enddo
      enddo
      enddo
 
      pause
 
      return 
      end

