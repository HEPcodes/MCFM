      subroutine qqb_ttb_g(p,msq)
      implicit none
C***********************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2002.                                                     *
*     calculate the element squared                                    *
*     for the process                                                  *
*----My notation                                                       *
*     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
*                       +t~(b~(p6)+e^-(p7)+nu(p8))+g(p9)               *
*                                                                      * 
*     Only five diagrams included, leading to 2 on-shell top quarks    *
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer b,j,k,h1,h2,h3,nu,j1,j2,j3
      double precision t(4),r(4),
     . msq(-nf:nf,-nf:nf),p(mxpart,4),ps(mxpart,4)
      double precision ttbqqbg_sq,fac,
     . wtgg,wtqqb,wtqbq,wtqg,wtqbarg,wtgq,wtgqbar
      double complex 
     . ttbgggppp,ttbgggmpp,ttbgggpmp,ttbgggppm,
     . ttbgggmmm,ttbgggpmm,ttbgggmpm,ttbgggmmp,      
     . a129(2,2,2),a192(2,2,2),a219(2,2,2),a291(2,2,2),
     . a912(2,2,2),a921(2,2,2),
     . a6sum,a3sum1a,a3sum1b,a3sum2a,a3sum2b,a3sum9a,a3sum9b
      double precision p3Dp5,p6Dp8,rDp7,tDp4,s34,s78
c--- these definitions are used for gauge check only
      double complex a,
     . ttbgggppp_full,ttbgggmpp_full,ttbgggpmp_full,ttbgggppm_full,
     . ttbgggmmm_full,ttbgggpmm_full,ttbgggmpm_full,ttbgggmmp_full,      
     . ttbqqbsqpp_full,ttbqqbsqpm_full,ttbqqbsqmp_full,ttbqqbsqmm_full,
     . ttbqqbtqpp_full,ttbqqbtqpm_full,ttbqqbtqmp_full,ttbqqbtqmm_full,
     . ttbqqbqqpp_full,ttbqqbqqpm_full,ttbqqbqqmp_full,ttbqqbqqmm_full,
     . ttbqqbrqpp_full,ttbqqbrqpm_full,ttbqqbrqmp_full,ttbqqbrqmm_full

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do nu=1,4
      do j=1,mxpart
      ps(j,nu)=p(j,nu)
      enddo
      enddo
      
      p3Dp5=p(3,4)*p(5,4)-p(3,3)*p(5,3)-p(3,2)*p(5,2)-p(3,1)*p(5,1) 
      p6Dp8=p(6,4)*p(8,4)-p(6,3)*p(8,3)-p(6,2)*p(8,2)-p(6,1)*p(8,1) 

      s34=2d0*(p(3,4)*p(4,4)-p(3,3)*p(4,3)-p(3,2)*p(4,2)-p(3,1)*p(4,1))
      s78=2d0*(p(7,4)*p(8,4)-p(7,3)*p(8,3)-p(7,2)*p(8,2)-p(7,1)*p(8,1))

c      we will have no further need for p3 and p5 
c      we will have no further need for p6 and p8 
      
      do nu=1,4
      t(nu)=p(3,nu)+p(4,nu)+p(5,nu)
      r(nu)=p(6,nu)+p(7,nu)+p(8,nu)
      enddo      
      tDp4=t(4)*p(4,4)-t(3)*p(4,3)-t(2)*p(4,2)-t(1)*p(4,1) 
      rDp7=r(4)*p(7,4)-r(3)*p(7,3)-r(2)*p(7,2)-r(1)*p(7,1)             
      do nu=1,4
c---t2
      ps(5,nu)=0.5d0*mt**2/tDp4*p(4,nu)
c---t1
      ps(3,nu)=t(nu)-ps(5,nu)

c---t2
      ps(8,nu)=0.5d0*mt**2/rDp7*p(7,nu)
c---t1
      ps(6,nu)=r(nu)-ps(8,nu)
      enddo

c      call writeout(ps)
c      pause

      call spinoru(9,ps,za,zb)

c--- gauge check of the ggg pieces (performed on 15/8/08)
c      j1=2
c      j2=1
c      write(6,*) 'ppp'
c      do j3=3,8
c      a=ttbgggppp_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,cdabs(a)
c      enddo
c      write(6,*) 'mmm'
c      do j3=3,8
c      a=ttbgggmmm_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,cdabs(a)
c      enddo
c      write(6,*) 'mpp'
c      do j3=3,8
c      a=ttbgggmpp_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,cdabs(a)
c      enddo
c      write(6,*) 'pmm'
c      do j3=3,8
c      a=ttbgggpmm_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,cdabs(a)
c      enddo
c      write(6,*) 'mpm'
c      do j3=3,8
c      a=ttbgggmpm_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,cdabs(a)
c      enddo
c      write(6,*) 'mmp'
c      do j3=3,8
c      a=ttbgggmmp_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,cdabs(a)
c      enddo	  
c      write(6,*) 'pmp'
c      do j3=3,8
c      a=ttbgggpmp_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,cdabs(a)
c      enddo
c      write(6,*) 'ppm'
c      do j3=3,8
c      a=ttbgggppm_full(1,2,9,5,3,6,8,j1,j2,j3)
c      write(6,*) j,a,cdabs(a)
c      enddo
c      pause

c--- gauge check of the qqbg pieces (performed on 15/8/08)
c      write(6,*) 'pp'
c      do j=3,8
c      write(6,*) j,cdabs(ttbqqbsqpp_full(1,2,9,5,3,6,8,j)),
c     .  	   cdabs(ttbqqbtqpp_full(1,2,9,5,3,6,8,j)),
c     .  	   cdabs(ttbqqbqqpp_full(1,2,9,5,3,6,8,j)),
c     .  	   cdabs(ttbqqbrqpp_full(1,2,9,5,3,6,8,j))
c      enddo
c      write(6,*) 'pm'
c      do j=3,8
c      write(6,*) j,cdabs(ttbqqbsqpm_full(1,2,9,5,3,6,8,j)),
c     .  	   cdabs(ttbqqbtqpm_full(1,2,9,5,3,6,8,j)),
c     .  	   cdabs(ttbqqbqqpm_full(1,2,9,5,3,6,8,j)),
c     .  	   cdabs(ttbqqbrqpm_full(1,2,9,5,3,6,8,j))
c      enddo
c      write(6,*) 'mp'
c      do j=3,8
c      write(6,*) j,cdabs(ttbqqbsqmp_full(1,2,9,5,3,6,8,j)),
c     .  	   cdabs(ttbqqbtqmp_full(1,2,9,5,3,6,8,j)),
c     .  	   cdabs(ttbqqbqqmp_full(1,2,9,5,3,6,8,j)),
c     .  	   cdabs(ttbqqbrqmp_full(1,2,9,5,3,6,8,j))
c      enddo
c      write(6,*) 'mm'
c      do j=3,8
c      write(6,*) j,cdabs(ttbqqbsqmm_full(1,2,9,5,3,6,8,j)),
c     .  	   cdabs(ttbqqbtqmm_full(1,2,9,5,3,6,8,j)),
c     .  	   cdabs(ttbqqbqqmm_full(1,2,9,5,3,6,8,j)),
c     .  	   cdabs(ttbqqbrqmm_full(1,2,9,5,3,6,8,j))
c      enddo
c      pause
      
c--- according to gfortify.frm and ttbggg.frm, the relationship between
c--- the indices of the amplitudes calls and the real momenta is:
c---  i1 -> p1                                         = ps(1)
c---  i2 -> p2                                         = ps(2)
c---  i3 -> p9                                         = ps(9)
c---  i4 -> pta                                        = ps(5)
c---  i5 -> ptb  with pta a rescaled version of p4     = ps(3)
c---  i6 -> pua                                        = ps(6)
c---  i7 -> pub  with pub a rescaled version of p7     = ps(8)

      a129(2,2,2)=ttbgggppp(1,2,9,5,3,6,8)
      a129(2,2,1)=ttbgggppm(1,2,9,5,3,6,8)
      a129(2,1,2)=ttbgggpmp(1,2,9,5,3,6,8)
      a129(2,1,1)=ttbgggpmm(1,2,9,5,3,6,8)
      a129(1,2,2)=ttbgggmpp(1,2,9,5,3,6,8)
      a129(1,2,1)=ttbgggmpm(1,2,9,5,3,6,8)
      a129(1,1,2)=ttbgggmmp(1,2,9,5,3,6,8)
      a129(1,1,1)=ttbgggmmm(1,2,9,5,3,6,8)

      a192(2,2,2)=ttbgggppp(1,9,2,5,3,6,8)
      a192(2,2,1)=ttbgggppm(1,9,2,5,3,6,8)
      a192(2,1,2)=ttbgggpmp(1,9,2,5,3,6,8)
      a192(2,1,1)=ttbgggpmm(1,9,2,5,3,6,8)
      a192(1,2,2)=ttbgggmpp(1,9,2,5,3,6,8)
      a192(1,2,1)=ttbgggmpm(1,9,2,5,3,6,8)
      a192(1,1,2)=ttbgggmmp(1,9,2,5,3,6,8)
      a192(1,1,1)=ttbgggmmm(1,9,2,5,3,6,8)

      a219(2,2,2)=ttbgggppp(2,1,9,5,3,6,8)
      a219(2,2,1)=ttbgggppm(2,1,9,5,3,6,8)
      a219(2,1,2)=ttbgggpmp(2,1,9,5,3,6,8)
      a219(2,1,1)=ttbgggpmm(2,1,9,5,3,6,8)
      a219(1,2,2)=ttbgggmpp(2,1,9,5,3,6,8)
      a219(1,2,1)=ttbgggmpm(2,1,9,5,3,6,8)
      a219(1,1,2)=ttbgggmmp(2,1,9,5,3,6,8)
      a219(1,1,1)=ttbgggmmm(2,1,9,5,3,6,8)

      a291(2,2,2)=ttbgggppp(2,9,1,5,3,6,8)
      a291(2,2,1)=ttbgggppm(2,9,1,5,3,6,8)
      a291(2,1,2)=ttbgggpmp(2,9,1,5,3,6,8)
      a291(2,1,1)=ttbgggpmm(2,9,1,5,3,6,8)
      a291(1,2,2)=ttbgggmpp(2,9,1,5,3,6,8)
      a291(1,2,1)=ttbgggmpm(2,9,1,5,3,6,8)
      a291(1,1,2)=ttbgggmmp(2,9,1,5,3,6,8)
      a291(1,1,1)=ttbgggmmm(2,9,1,5,3,6,8)

      a912(2,2,2)=ttbgggppp(9,1,2,5,3,6,8)
      a912(2,2,1)=ttbgggppm(9,1,2,5,3,6,8)
      a912(2,1,2)=ttbgggpmp(9,1,2,5,3,6,8)
      a912(2,1,1)=ttbgggpmm(9,1,2,5,3,6,8)
      a912(1,2,2)=ttbgggmpp(9,1,2,5,3,6,8)
      a912(1,2,1)=ttbgggmpm(9,1,2,5,3,6,8)
      a912(1,1,2)=ttbgggmmp(9,1,2,5,3,6,8)
      a912(1,1,1)=ttbgggmmm(9,1,2,5,3,6,8)

      a921(2,2,2)=ttbgggppp(9,2,1,5,3,6,8)
      a921(2,2,1)=ttbgggppm(9,2,1,5,3,6,8)
      a921(2,1,2)=ttbgggpmp(9,2,1,5,3,6,8)
      a921(2,1,1)=ttbgggpmm(9,2,1,5,3,6,8)
      a921(1,2,2)=ttbgggmpp(9,2,1,5,3,6,8)
      a921(1,2,1)=ttbgggmpm(9,2,1,5,3,6,8)
      a921(1,1,2)=ttbgggmmp(9,2,1,5,3,6,8)
      a921(1,1,1)=ttbgggmmm(9,2,1,5,3,6,8)

      wtgg=0d0
      do h1=1,2
      do h2=1,2
      do h3=1,2
c--- NB: make sure to permute helicity labels appropriately too
        a3sum9a=a129(h1,h2,h3)+a192(h1,h3,h2)+a912(h3,h1,h2)
        a3sum1a=a291(h1,h2,h3)+a219(h1,h3,h2)+a129(h3,h1,h2)
        a3sum2a=a912(h1,h2,h3)+a921(h1,h3,h2)+a291(h3,h1,h2)
        a3sum9b=a219(h2,h1,h3)+a291(h2,h3,h1)+a921(h3,h2,h1)
        a3sum1b=a921(h2,h1,h3)+a912(h2,h3,h1)+a192(h3,h2,h1)
        a3sum2b=a192(h2,h1,h3)+a129(h2,h3,h1)+a219(h3,h2,h1)
	a6sum=a3sum1a+a3sum1b
	wtgg=wtgg+xn**3*cf*(
     .   (cdabs(a129(h1,h2,h3))**2+cdabs(a192(h1,h3,h2))**2
     .   +cdabs(a219(h2,h1,h3))**2+cdabs(a291(h2,h3,h1))**2
     .   +cdabs(a912(h3,h1,h2))**2+cdabs(a921(h3,h2,h1))**2)
     .  -(cdabs(a3sum1a)**2+cdabs(a3sum2a)**2+cdabs(a3sum9a)**2
     .   +cdabs(a3sum1b)**2+cdabs(a3sum2b)**2+cdabs(a3sum9b)**2)/xn**2
     .  +(cdabs(a6sum)**2)*(xn**2+1d0)/xn**4
     .                     )
      enddo
      enddo
      enddo
            
      wtqqb=ttbqqbg_sq(1,2,9,5,3,6,8)
      wtqbq=ttbqqbg_sq(2,1,9,5,3,6,8)
      wtqg=ttbqqbg_sq(1,9,2,5,3,6,8)
      wtgq=ttbqqbg_sq(2,9,1,5,3,6,8)
      wtqbarg=ttbqqbg_sq(9,1,2,5,3,6,8)
      wtgqbar=ttbqqbg_sq(9,2,1,5,3,6,8)

c--- overall factor, starting with couplings
      fac=2d0*(gwsq/2d0)**4*gsq**3
c--- include top and anti-top propagators
      fac=fac/(mt*twidth)**4
c--- include W decays
      fac=fac
     . *8d0*p3Dp5/((s34-wmass**2)**2+(wmass*wwidth)**2)
     . *8d0*p6Dp8/((s78-wmass**2)**2+(wmass*wwidth)**2)
c--- correct normalization for p4 and p7     
      fac=fac
     . /(0.5d0*mt**2/tDp4)
     . /(0.5d0*mt**2/rDp7)
     
C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if     (j .gt. 0) then
          msq(j,-j)=aveqq*fac*wtqqb
	  msq(j, 0)=aveqg*fac*wtqg
	  msq(0, j)=aveqg*fac*wtgq
      elseif (j .eq. 0) then
          msq(j,j)=avegg*fac*wtgg
      elseif (j .lt. 0) then
          msq(j,-j)=aveqq*fac*wtqbq
	  msq(j, 0)=aveqg*fac*wtqbarg
	  msq(0, j)=aveqg*fac*wtgqbar
      endif
      enddo
      return
      end



      double precision function ttbqqbg_sq(i1,i2,i9,i5,i3,i6,i8)
c--- returns the summed squared helicity amplitudes for the
c--- ttbqqbg amplitudes
      implicit none
      include 'constants.f'
      integer i1,i2,i9,i5,i3,i6,i8
      double complex sq,tq,qq,rq,
     . ttbqqbsqpp,ttbqqbsqpm,ttbqqbsqmp,ttbqqbsqmm,
     . ttbqqbtqpp,ttbqqbtqpm,ttbqqbtqmp,ttbqqbtqmm,
     . ttbqqbqqpp,ttbqqbqqpm,ttbqqbqqmp,ttbqqbqqmm,
     . ttbqqbrqpp,ttbqqbrqpm,ttbqqbrqmp,ttbqqbrqmm
      double precision appsq,apmsq,ampsq,ammsq
      
c-- for q-qb , there are four colour amplitudes:
c---    sq	proportional to Ta(it,i1)*delta(i2,ib)
c---    tq	proportional to Ta(i2,ib)*delta(it,i1)
c---    qq	proportional to Ta(i2,i1)*delta(it,ib)
c---    rq	proportional to Ta(it,ib)*delta(i2,i1)
      sq=ttbqqbsqpp(i1,i2,i9,i5,i3,i6,i8)
      tq=ttbqqbtqpp(i1,i2,i9,i5,i3,i6,i8)
      qq=ttbqqbqqpp(i1,i2,i9,i5,i3,i6,i8)
      rq=ttbqqbrqpp(i1,i2,i9,i5,i3,i6,i8)
      appsq=xn*cf*(
     . +xn*cdabs(sq)**2+xn*cdabs(tq)**2
     . +xn*cdabs(qq)**2+xn*cdabs(rq)**2
     . +2d0*Dble((sq+tq)*dconjg(rq+qq)))
      
      sq=ttbqqbsqpm(i1,i2,i9,i5,i3,i6,i8)
      tq=ttbqqbtqpm(i1,i2,i9,i5,i3,i6,i8)
      qq=ttbqqbqqpm(i1,i2,i9,i5,i3,i6,i8)
      rq=ttbqqbrqpm(i1,i2,i9,i5,i3,i6,i8)
      apmsq=xn*cf*(
     . +xn*cdabs(sq)**2+xn*cdabs(tq)**2
     . +xn*cdabs(qq)**2+xn*cdabs(rq)**2
     . +2d0*Dble((sq+tq)*dconjg(rq+qq)))
      
      sq=ttbqqbsqmp(i1,i2,i9,i5,i3,i6,i8)
      tq=ttbqqbtqmp(i1,i2,i9,i5,i3,i6,i8)
      qq=ttbqqbqqmp(i1,i2,i9,i5,i3,i6,i8)
      rq=ttbqqbrqmp(i1,i2,i9,i5,i3,i6,i8)
      ampsq=xn*cf*(
     . +xn*cdabs(sq)**2+xn*cdabs(tq)**2
     . +xn*cdabs(qq)**2+xn*cdabs(rq)**2
     . +2d0*Dble((sq+tq)*dconjg(rq+qq)))
      
      sq=ttbqqbsqmm(i1,i2,i9,i5,i3,i6,i8)
      tq=ttbqqbtqmm(i1,i2,i9,i5,i3,i6,i8)
      qq=ttbqqbqqmm(i1,i2,i9,i5,i3,i6,i8)
      rq=ttbqqbrqmm(i1,i2,i9,i5,i3,i6,i8)
      ammsq=xn*cf*(
     . +xn*cdabs(sq)**2+xn*cdabs(tq)**2
     . +xn*cdabs(qq)**2+xn*cdabs(rq)**2
     . +2d0*Dble((sq+tq)*dconjg(rq+qq)))
      
      ttbqqbg_sq=appsq+apmsq+ampsq+ammsq

      return
      end
      
