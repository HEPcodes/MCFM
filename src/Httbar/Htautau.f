      subroutine Htautau(p,msq)
      implicit none
c--- calculates the squared matrix elements for the process
c---  H -> tau  (-> vt (p3) + e-(p4) + ve~(p5) )
c---     + tau~ (-> vt~(p6) + ve(p7) + e+ (p8) )
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'dprodx.f'
      include 'sprodx.f'
      integer j,k
      double complex amp1,amp2
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      double precision s129,s345,s678,prop,msqhtautau
      double precision ehsvm3,ehsvm4,origmbsq,gg,qg,gq,qq
      
      s129=s(1,2)+s(1,9)+s(2,9)
      s345=s(3,4)+s(3,5)+s(4,5)
      s678=s(6,7)+s(6,8)+s(7,8)
      
c      fac=gw**10*0.5d0*mtau**2/wmass**2*(s129-4d0*mtau**2)/
c     .     ((s129-hmass**2)**2+(hmass*hwidth)**2)
      fac=gw**10*mtau**2/wmass**2/((s129-hmass**2)**2+(hmass*hwidth)**2)

c--- fill spinor products
      call spinoru(9,p,za,zb)

      amp1 = 4d0*za(3,4)*(zb(5,3)*za(3,6)*zb(6,8)
     .                   +zb(5,3)*za(3,7)*zb(7,8)
     .                   +zb(5,4)*za(4,6)*zb(6,8)
     .                   +zb(5,4)*za(4,7)*zb(7,8))*za(7,6)
     
      amp2 = 4d0*mtau**2*za(3,4)*zb(5,8)*za(7,6)
     
      prop=1d0/((s(4,5)-wmass**2)**2+(wmass*wwidth)**2)
     .        /((s(7,8)-wmass**2)**2+(wmass*wwidth)**2)
     .        /((s345-mtau**2)**2+(mtau*tauwidth)**2)
     .        /((s678-mtau**2)**2+(mtau*tauwidth)**2)
           
      msqhtautau=fac*cdabs(amp1+amp2)**2*prop
      
      origmbsq=mbsq
      mbsq=175d0**2
      gg=+avegg*ehsvm3(s(1,2),s(1,9),s(2,9))*msqhtautau
      qq=+aveqq*ehsvm4(s(1,2),s(1,9),s(2,9))*msqhtautau
      qg=-aveqg*ehsvm4(s(1,9),s(1,2),s(2,9))*msqhtautau
      gq=-aveqg*ehsvm4(s(2,9),s(1,9),s(1,2))*msqhtautau
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
      elseif ((j.eq.-k).and.(j.ne.0)) then
           msq(j,k)=qq
      endif

      enddo
      enddo

      return
      end
      
