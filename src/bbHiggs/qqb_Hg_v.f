      subroutine qqb_Hg_v(p,msq)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c     parton(-p1)+parton(-p2) --> H(p)+parton(p5)
c                                  |
c                                   --> b(p3)+bb(p4)
c                            
c--all momenta incoming
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'susycoup.f'
      include 'scheme.f'
      include 'runmb.f'
      include 'scale.f'
      integer j,k,j1,j2,j3
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),s
      double precision coupsq_eff,ghbb_eff
      double precision amz,alphas,mb_eff,mb_msbar,massfrun
      double precision fac,propsq,hdecay,coupsq,ghbb,bbghvirt
      common/couple/amz

c--susycoup is the deviation of Higgs coupling 
c-- from the standard model value

c--statement function 
      s(j,k)=2d0*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     .           -p(j,2)*p(k,2)-p(j,3)*p(k,3))
C---ur-amplitude is b(p1)+bbar(p2)+g(p3)+H(q)
c      amp(j1,j2,j3)=4d0*(s(j1,j2)**2+s(3,4)**2)/(s(j1,j3)*s(j2,j3))

      scheme='dred'


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      if (s(3,4) .lt. 4d0*mbsq) return

c--- calculate MS-bar mass from pole mass
c      mb_msbar=mb/(1d0+cf*alphas(mb,amz,2)/pi)
      mb_msbar=4.20d0
c--- run mb to appropriate scale, if necessary
      if (runmb) then
        mb_eff=massfrun(mb_msbar,scale,amz,2)
      else
        mb_eff=mb_msbar
      endif
      
c--- our definition
c      ghbb=sqrt(esq*mbsq/(xw*(1d0-xw)))/2d0/zmass
c--- definition according to Maltoni, Willenbrock
      ghbb=dsqrt(esq/xw)*mb/2d0/wmass
      coupsq=susycoup**2*ghbb**2
      hdecay=coupsq*2d0*(s(3,4)-4d0*mb**2)*xn 
      propsq=1d0/((s(3,4)-hmass**2)**2+(hmass*hwidth)**2)
c--- The _eff couplings include the running mass
c--- We need to separate these from the factors associated with the
c--- Higgs decay, because the Br. Ratio does not include running mb
      ghbb_eff=dsqrt(esq/xw)*mb_eff/2d0/wmass
      coupsq_eff=susycoup**2*ghbb_eff**2

      fac=coupsq_eff*propsq*hdecay
c      fac=CF*xn*gsq*coupsq*propsq*hdecay

      msq(0,+5)=+fac*aveqg*bbghvirt(2,5,1)
      msq(0,-5)=+fac*aveqg*bbghvirt(5,2,1)
      msq(+5,0)=+fac*aveqg*bbghvirt(1,5,2)
      msq(-5,0)=+fac*aveqg*bbghvirt(5,1,2)

      return
      end
          
    

      
     
