      subroutine gg_hgg_gvec(p,n,in,msq)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'msq_struc.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
C  in is the label of the momentum contracted with n
      integer j,k,in
      double precision msq(-nf:nf,-nf:nf)
      double precision n(4),p(mxpart,4),hdecay,s34,fac,
     . qqgghn_ab,qqgghn_ba,qqgghn_sym,ab,ba,sym,
     . c1234,c1243,c1423,c1234x,c1243x,c1423x,p1p2(-1:1,-1:1),Asq

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C---fill dot products
c      call dotem(6,p,s)  
      call spinoru(6,p,za,zb)  

C   Deal with Higgs decay to b-bbar
      s34=s(3,4)+2d0*mb**2
      hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s34-4d0*mb**2) 
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      Asq=(as/(3d0*pi))**2/vevsq
      fac=gsq**2*Asq*hdecay

      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=0d0
      enddo
      enddo

c--- NOTE: there seems to be some redundancy in the function calls, e.g.
c--- the entries for (0,-1) = (0,+1). Need to check if this is true
c--- and eliminate where possible

C     The function qqgghn calculates q(p1)+qbar(p2)-->H+g(p3)+g(p4)
C     with p4 contracted with n
C     The function gggghn calculates g(p1)+g(p2)-->H+g(p3)+g(p4)
C     with p1 contracted with n

c--- Note that I have removed all references to p1p2(j,k) since
c--- the appropriate terms are actually used only in their 
c--- colour-separated forms

      if     (in .eq. 1) then
c--- check against squared 'n' amplitudes
        call qqgghn(2,5,6,1,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 1',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(2,5,6,1,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 1',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,0,-1)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_ba,0,-1)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_sym,0,-1)=-aveqg*fac*qqgghn_sym
c        p1p2(0,-1)=msq_strucv(igg_ab,0,-1)+msq_strucv(igg_ba,0,-1)
c     .            +msq_strucv(igg_sym,0,-1)
c        p1p2(0,+1)=-aveqg*fac*qqgghn_old(5,2,6,1,p,n)
c--- check against squared 'n' amplitudes
        call qqgghn(5,2,6,1,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 2',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(5,2,6,1,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 2',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,0,+1)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_ba,0,+1)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_sym,0,+1)=-aveqg*fac*qqgghn_sym
c        p1p2(0,+1)=msq_strucv(igg_ab,0,+1)+msq_strucv(igg_ba,0,+1)
c     .            +msq_strucv(igg_sym,0,+1)
c--- check against squared 'n' amplitudes
        call qqgghn(5,6,2,1,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 3',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(5,6,2,1,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 3',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,0,0)=avegg*fac*dfloat(nf)*qqgghn_ab
        msq_strucv(igg_ba,0,0)=avegg*fac*dfloat(nf)*qqgghn_ba
        msq_strucv(igg_sym,0,0)=avegg*fac*dfloat(nf)*qqgghn_sym
c        p1p2(0,0)=0.5d0*(msq_strucv(igggg_a,0,0)+msq_strucv(igggg_b,0,0)
c     .                  +msq_strucv(igggg_c,0,0))/2d0
c     .      +dfloat(nf)*(msq_strucv(igg_ab,0,0)+msq_strucv(igg_ba,0,0)
c     .                  +msq_strucv(igg_sym,0,0))

c--- check against squared 'n' amplitudes
        call gggghn_amp(1,2,5,6,p,n,c1234x,c1243x,c1423x)
        write(6,*) 'old a',c1234x,c1243x,c1423x
        call gggghn(1,2,5,6,p,n,c1234,c1243,c1423)
        write(6,*) 'new a',c1234,c1243,c1423
        write(6,*) 'new a',c1234/c1234x,c1243/c1243x,c1423/c1423x
        pause
        msq_strucv(igggg_a,0,0)=avegg*fac*0.5d0*c1234  
        msq_strucv(igggg_b,0,0)=avegg*fac*0.5d0*c1423 
        msq_strucv(igggg_c,0,0)=avegg*fac*0.5d0*c1243  
c        p1p2(0,0)=avegg*fac*(0.5d0*(c1234+c1243+c1423)
c     .                 +dfloat(nf)*qqgghn_old(5,6,2,1,p,n))

c        p1p2(0,-1)=-aveqg*fac*qqgghn_old(2,5,6,1,p,n)

      elseif (in .eq. 2) then
c        p1p2(+1,0)=-aveqg*fac*qqgghn_old(1,5,6,2,p,n)
c--- check against squared 'n' amplitudes
        call qqgghn(1,5,6,2,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 4',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(1,5,6,2,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 4',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,+1,0)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_ba,+1,0)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_sym,+1,0)=-aveqg*fac*qqgghn_sym
c        p1p2(+1,0)=msq_strucv(igg_ab,+1,0)+msq_strucv(igg_ba,+1,0)
c     .            +msq_strucv(igg_sym,+1,0)
c        p1p2(-1,0)=-aveqg*fac*qqgghn_old(5,1,6,2,p,n)
c--- check against squared 'n' amplitudes
        call qqgghn(5,1,6,2,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 5',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(5,1,6,2,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 5',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,-1,0)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_ba,-1,0)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_sym,-1,0)=-aveqg*fac*qqgghn_sym
c        p1p2(-1,0)=msq_strucv(igg_ab,-1,0)+msq_strucv(igg_ba,-1,0)
c     .            +msq_strucv(igg_sym,-1,0)
c--- check against squared 'n' amplitudes
        call qqgghn(5,6,1,2,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 6',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(5,6,1,2,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 6',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,0,0)=avegg*fac*dfloat(nf)*qqgghn_ab
        msq_strucv(igg_ba,0,0)=avegg*fac*dfloat(nf)*qqgghn_ba
        msq_strucv(igg_sym,0,0)=avegg*fac*dfloat(nf)*qqgghn_sym
c        p1p2(0,0)=0.5d0*(msq_strucv(igggg_a,0,0)+msq_strucv(igggg_b,0,0)
c     .                  +msq_strucv(igggg_c,0,0))/2d0
c     .      +dfloat(nf)*(msq_strucv(igg_ab,0,0)+msq_strucv(igg_ba,0,0)
c     .                  +msq_strucv(igg_sym,0,0))
c--- check against squared 'n' amplitudes
        call gggghn_amp(2,1,5,6,p,n,c1234x,c1243x,c1423x)
        write(6,*) 'old b',c1234x,c1243x,c1423x
        call gggghn(2,1,5,6,p,n,c1234,c1243,c1423)
        write(6,*) 'new b',c1234,c1243,c1423
        write(6,*) 'new b',c1234/c1234x,c1243/c1243x,c1423/c1423x
        pause
        msq_strucv(igggg_a,0,0)=avegg*fac*0.5d0*c1234     
        msq_strucv(igggg_b,0,0)=avegg*fac*0.5d0*c1423     
        msq_strucv(igggg_c,0,0)=avegg*fac*0.5d0*c1243     
c        p1p2(0,0)=+avegg*fac*(0.5d0*(c1234+c1243+c1324)
c     .                 +dfloat(nf)*qqgghn_old(5,6,1,2,p,n))

      elseif (in .eq. 5) then     
c        p1p2(1,-1)=+aveqq*fac*qqgghn_old(1,2,6,5,p,n)
c--- check against squared 'n' amplitudes
        call qqgghn(1,2,6,5,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 7',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(1,2,6,5,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 7',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,+1,-1)=aveqq*fac*0.5d0*qqgghn_ab
        msq_strucv(igg_ba,+1,-1)=aveqq*fac*0.5d0*qqgghn_ba
        msq_strucv(igg_sym,+1,-1)=aveqq*fac*0.5d0*qqgghn_sym
c        p1p2(+1,-1)=msq_strucv(igg_ab,+1,-1)+msq_strucv(igg_ba,+1,-1)
c     .             +msq_strucv(igg_sym,+1,-1)
c        p1p2(-1,1)=+aveqq*fac*qqgghn_old(2,1,6,5,p,n)
c--- check against squared 'n' amplitudes
        call qqgghn(2,1,6,5,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 8',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(2,1,6,5,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 8',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,-1,+1)=aveqq*fac*0.5d0*qqgghn_ab
        msq_strucv(igg_ba,-1,+1)=aveqq*fac*0.5d0*qqgghn_ba
        msq_strucv(igg_sym,-1,+1)=aveqq*fac*0.5d0*qqgghn_sym
c        p1p2(-1,+1)=msq_strucv(igg_ab,-1,+1)+msq_strucv(igg_ba,-1,+1)
c     .             +msq_strucv(igg_sym,-1,+1)
c--- for the qg, gq pieces, note that qbar-g and g-qbar are never used        
c--- check against squared 'n' amplitudes
        call qqgghn(1,6,2,5,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 9',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(1,6,2,5,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 9',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,+1,0)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_ba,+1,0)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_sym,+1,0)=-aveqg*fac*qqgghn_sym
c        p1p2(+1,0)=msq_strucv(igg_ab,+1,0)+msq_strucv(igg_ba,+1,0)
c     .            +msq_strucv(igg_sym,+1,0)
c--- check against squared 'n' amplitudes
        call qqgghn(2,6,1,5,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 10',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(2,6,1,5,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 10',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,0,+1)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_ba,0,+1)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_sym,0,+1)=-aveqg*fac*qqgghn_sym
c        p1p2(0,+1)=msq_strucv(igg_ab,0,+1)+msq_strucv(igg_ba,0,+1)
c     .            +msq_strucv(igg_sym,0,+1)
c--- check against squared 'n' amplitudes
        call gggghn_amp(5,6,1,2,p,n,c1234x,c1243x,c1423x)
        write(6,*) 'old c',c1234x,c1243x,c1423x
        call gggghn(5,6,1,2,p,n,c1234,c1243,c1423)
        write(6,*) 'new c',c1234,c1243,c1423
        write(6,*) 'new c',c1234/c1234x,c1243/c1243x,c1423/c1423x
        pause
        msq_strucv(igggg_a,0,0)=avegg*fac*0.5d0*c1234     
        msq_strucv(igggg_b,0,0)=avegg*fac*0.5d0*c1423     
        msq_strucv(igggg_c,0,0)=avegg*fac*0.5d0*c1243     
c        p1p2(0,0)=+0.5d0*avegg*fac*(c1234+c1243+c1423)

      elseif (in .eq. 6) then     
c        p1p2(1,-1)=+aveqq*fac*qqgghn_old(1,2,5,6,p,n)
c--- check against squared 'n' amplitudes
        call qqgghn(1,2,5,6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 11',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(1,2,5,6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 11',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,+1,-1)=aveqq*fac*0.5d0*qqgghn_ab
        msq_strucv(igg_ba,+1,-1)=aveqq*fac*0.5d0*qqgghn_ba
        msq_strucv(igg_sym,+1,-1)=aveqq*fac*0.5d0*qqgghn_sym
c        p1p2(+1,-1)=msq_strucv(igg_ab,+1,-1)+msq_strucv(igg_ba,+1,-1)
c     .             +msq_strucv(igg_sym,+1,-1)
c        p1p2(-1,1)=+aveqq*fac*qqgghn_old(2,1,5,6,p,n)
c--- check against squared 'n' amplitudes
        call qqgghn(2,1,5,6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 12',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(2,1,5,6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 12',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,-1,+1)=aveqq*fac*0.5d0*qqgghn_ab
        msq_strucv(igg_ba,-1,+1)=aveqq*fac*0.5d0*qqgghn_ba
        msq_strucv(igg_sym,-1,+1)=aveqq*fac*0.5d0*qqgghn_sym
c        p1p2(-1,+1)=msq_strucv(igg_ab,-1,+1)+msq_strucv(igg_ba,-1,+1)
c     .             +msq_strucv(igg_sym,-1,+1)
c--- for the qg, gq pieces, note that qbar-g and g-qbar are never used        
c--- check against squared 'n' amplitudes
        call qqgghn(1,5,2,6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 13',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(1,5,2,6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 13',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,+1,0)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_ba,+1,0)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_sym,+1,0)=-aveqg*fac*qqgghn_sym
c        p1p2(+1,0)=msq_strucv(igg_ab,+1,0)+msq_strucv(igg_ba,+1,0)
c     .            +msq_strucv(igg_sym,+1,0)
c--- check against squared 'n' amplitudes
        call qqgghn(2,5,1,6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'old 14',qqgghn_ab,qqgghn_ba,qqgghn_sym
        call qqgghn_dfm_amp(2,5,1,6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        write(6,*) 'new 14',qqgghn_ab,qqgghn_ba,qqgghn_sym
        msq_strucv(igg_ab,0,+1)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_ba,0,+1)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_sym,0,+1)=-aveqg*fac*qqgghn_sym
c        p1p2(0,+1)=msq_strucv(igg_ab,0,+1)+msq_strucv(igg_ba,0,+1)
c     .            +msq_strucv(igg_sym,0,+1)
c--- check against squared 'n' amplitudes
        call gggghn_amp(6,1,2,5,p,n,c1234x,c1243x,c1423x)
        write(6,*) 'old d',c1234x,c1243x,c1423x
        call gggghn(6,1,2,5,p,n,c1234,c1243,c1423)
        write(6,*) 'new d',c1234,c1243,c1423
        write(6,*) 'new d',c1234/c1234x,c1243/c1243x,c1423/c1423x
        pause
        msq_strucv(igggg_a,0,0)=avegg*fac*0.5d0*c1234     
        msq_strucv(igggg_b,0,0)=avegg*fac*0.5d0*c1423     
        msq_strucv(igggg_c,0,0)=avegg*fac*0.5d0*c1243     
c        p1p2(0,0)=+0.5d0*avegg*fac*(c1234+c1243+c1423)
      endif
 
      return
      end


