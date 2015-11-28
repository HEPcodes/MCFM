      subroutine qqb_higgs(p,msq)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c     f(-p1) + f(-p2) --> H + f(p5)
c                         |
c                         --> b(p3)+bbar(p4)
c                            
c--all momenta incoming

      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'dprodx.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),gg,qg,gq,qq,hbb
      double precision ehsvm3,ehsvm4,origmbsq

c--- matrix element for H -> bbbar
      hbb=xn*gw**2*0.5d0*mbsq/wmass**2*(s(3,4)-4d0*mb**2)/
     .     ((s(3,4)-hmass**2)**2+(hmass*hwidth)**2)

      origmbsq=mbsq
      mbsq=175d0**2
      gg=+avegg*ehsvm3(s(1,2),s(1,5),s(2,5))*hbb
      qq=+aveqq*ehsvm4(s(1,2),s(1,5),s(2,5))*hbb
      qg=-aveqg*ehsvm4(s(1,5),s(1,2),s(2,5))*hbb
      gq=-aveqg*ehsvm4(s(2,5),s(1,5),s(1,2))*hbb
      mbsq=origmbsq

      do j=-nf,nf    
      do k=-nf,nf
      msq(j,k)=0d0

      if ((j.eq. 0) .or. (k.eq.0)) then
           if ((j.eq. 0) .and. (k.eq.0)) then
                msq(j,k)=gg
           elseif ((j.eq.0).and.(k.ne.0)) then
                msq(j,k)=gq
           elseif ((j.ne.0).and.(k.eq.0)) then
                msq(j,k)=qg
           endif
      elseif ((j.eq.-k).and. (j.ne.0)) then
           msq(j,k)=qq
      endif

      enddo
      enddo

      return
      end
      
      
      double precision function ehsvm3(s,t,u)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
c---Matrix element squared Eqn 2.2 of EHSV
      double complex ehsva2,ehsva4
      double precision s,t,u,hmass2
      logical approx
      parameter(approx=.false.)
c--- approx TRUE uses the heavy fermion approximation to Msq


      if (approx) then
      ehsvm3=gwsq/pi*as**3*xn*V/9d0*(
     .        hmass**8+s**4+t**4+u**4)/s/t/u/wmass**2
      else
      ehsvm3=
     . abs(ehsva2(s,t,u))**2+abs(ehsva2(u,s,t))**2+abs(ehsva2(t,u,s))**2
     . +abs(ehsva4(s,t,u))**2 
      ehsvm3=gwsq/pi*as**3*xn*V*hmass**8/(s*t*u*wmass**2)*ehsvm3
      endif
      
      return
      end
      
      
      double precision function ehsvm4(s,t,u)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
c---Matrix element squared Eqn 2.6 of EHSV
      double complex ehsva5
      double precision s,t,u
      logical approx
      parameter(approx=.false.)
c--- approx TRUE uses the heavy fermion approximation to Msq


      if (approx) then
      ehsvm4=gwsq/pi*as**3*V/18d0*(u**2+t**2)/s/wmass**2
      else
      ehsvm4=abs(ehsva5(s,t,u))**2
      ehsvm4=gwsq/(4d0*pi)*as**3*V/2d0*(u**2+t**2)/(s*wmass**2)
     . *hmass**4/(u+t)**2*ehsvm4
      endif
      
      return
      end
      
      
 
