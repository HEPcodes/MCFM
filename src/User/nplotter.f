      subroutine bookplot(n,tag,titlex,var,wt,xmin,xmax,dx,llplot) 
      include 'nplot.f'
      integer n
      character titlex*8,llplot*3,tag*4
      double precision var,wt,xmin,xmax,dx
      if (tag.eq.'book') then
        call mbook(n,titlex,dx,xmin,xmax)
       elseif (tag .eq. 'plot') then
        call mfill(n,var,wt)
        linlog(n)=llplot
        titlearray(n)=titlex
      endif
      return
      end

      subroutine nplotter(vector,s,p,wt,switch)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'npart.f'
      include 'clustering.f'
      integer idum,n,switch,jets
      character tag*4,jetlabel(mxpart)*2
      double precision m56,m56_5,m56_10,m56_11,m56_12,m56_13,m56_15,
     . sigma,m34,m345,m346,m3456,m678,m47,etmiss,misset,
     . s(mxpart,mxpart),p(mxpart,4),eta,root,wt1
      double precision y3,y4,y5,y6,y7,y8,y34,y56
      double precision pt3,pt4,pt5,pt6,pt7,pt8,pt34,pt56,pt34a,pt34b
      double precision ptbbsq,ptbbpair,m56smw,gasdev,
     . pt,yrap,chi,cosphi,phi,var,vector(mxdim),wt,r,
     . kt12,kt14,kt15,kt24,kt25,kt56,costh,cosnew1,cosnew2,cosnew3
      double precision transm,transcm,dot,pttwo,yraptwo,sdot30m,mttbar
      double precision phill,thetall,fphi,ftheta,mtsqlet,mt1,mt2,mll
      double precision c4,cosnchi,nchi,m56psm20,m56psm40,smearp,swap
      double precision clustermass,m56clust,costhdd,cosllet,coslpairet
      double precision pjet(mxpart,4),bclustmass,rn,etvec(4),r56,r35,r36
      common/parts/jets,jetlabel
      double precision dsigdy,dsigdytmp
      double precision es17,es27,es56,es57,es67
      integer nproc,eventpart
      common/nproc/nproc
      logical first,bbproc
      data idum/34265765/
      data first/.true./
      save first,eta
      save es17,es27,es56,es57,es67
      if (first) then
        first=.false.
        tag='book'
c--- ensure we initialize all possible histograms
        eventpart=npart+3
      else
        tag='plot'

c--- set bbproc to TRUE if the process involves two b-jets
        if (
     .      (nproc .eq.  21)
     . .or. (nproc .eq.  26)
     . .or. (nproc .eq.  51)
     . .or. (nproc .eq.  52)
     . .or. (nproc .eq.  53)
     . .or. (nproc .eq.  73)
     . .or. (nproc .eq.  78)
     . .or. (nproc .eq.  84)
     . .or. (nproc .eq.  89)
     . .or. (nproc .eq.  91)
     . .or. (nproc .eq.  96)
     . .or. (nproc .eq.  101)
     . .or. (nproc .eq.  102)
     . .or. (nproc .eq.  151)
     . .or. (nproc .eq.  131)
     . .or. (nproc .eq.  152)
     . .or. (nproc .eq.  161)
     . .or. (nproc .eq.  171)
     . ) then
          bbproc=.true.
        else
          bbproc=.false.
        endif

c--- eventpart will contain the number of actual particles that have
c--- a defined momentum
c--- for lowest order and virtual terms switch=0 and eventpart=npart+2
c--- for real events switch=0 and eventpart=npart+2
c--- for real counter-events switch=1 and eventpart=npart+1
 
      m345=dsqrt(2d0*(dot(p,3,4)+dot(p,3,5)+dot(p,4,5)))
      m678=dsqrt(2d0*(dot(p,6,4)+dot(p,6,8)+dot(p,7,8)))
      m346=dsqrt(2d0*(dot(p,3,4)+dot(p,3,6)+dot(p,4,6)))
      m3456=dsqrt(2d0*(dot(p,3,4)+dot(p,3,5)+dot(p,3,6)
     .                +dot(p,4,5)+dot(p,4,6)+dot(p,5,6)))
      mttbar=2d0*(
     . +dot(p,3,4)+dot(p,3,5)+dot(p,3,6)+dot(p,3,7)+dot(p,3,8)
     .            +dot(p,4,5)+dot(p,4,6)+dot(p,4,7)+dot(p,4,8)
     .                       +dot(p,5,6)+dot(p,5,7)+dot(p,5,8)
     .                                  +dot(p,6,7)+dot(p,6,8)
     .                                             +dot(p,7,8))

      mttbar=dsqrt(mttbar)
c      write(6,*) 'm345',m345
c      write(6,*) 'm678',m678
c      write(6,*) 'mttbar',mttbar
c      pause
        eventpart=npart-switch+2
        if (jets .gt. 0) eventpart=4+jets
          
        if (switch .eq. 0) then
          es17=2d0*dot(p,1,7)  
          es27=2d0*dot(p,2,7)  
          es56=2d0*dot(p,5,6)  
          es57=2d0*dot(p,5,7)  
          es67=2d0*dot(p,6,7)  
        endif

        m34=dsqrt(2d0*dot(p,3,4))

        if (bbproc .and. clustering) then
c--- returns zero cluster mass if two b's are in one jet
          m56clust=dsqrt(bclustmass(jets,p,jetlabel))
        else
          m56clust=dsqrt(2d0*dot(p,5,6))
        endif  
        r56=r(p,5,6)
        r35=r(p,3,5)
        r36=r(p,3,6)
c--generate a gaussian kick only for event
        if (switch .eq. 0) eta=gasdev(idum)
        sigma=5d0
        m56_5=m56clust+sigma*eta
        sigma=10d0
        m56_10=m56clust+sigma*eta
        sigma=11d0
        m56_11=m56clust+sigma*eta
        sigma=12d0
        m56_12=m56clust+sigma*eta
        sigma=13d0
        m56_13=m56clust+sigma*eta
        sigma=15d0
        m56_15=m56clust+sigma*eta

        y3=yrap(3,p)
        pt3=pt(3,p)
        y4=yrap(4,p)
        pt4=pt(4,p)        
        y34=yraptwo(3,4,p)
        pt34=pttwo(3,4,p)
        if ((y34 .lt. -0.5d0) .or. (y34 .gt. 0.5d0)) then
          pt34a=-1d0
        else
          pt34a=pt34
        endif
        if ((y34 .lt. -0.1d0) .or. (y34 .gt. 0.1d0)) then
          pt34b=-1d0
        else
          pt34b=pt34*5d0
        endif
                
        if (eventpart .gt. 4) then        
        y5=yrap(5,p)
        pt5=pt(5,p)
        endif

        if (eventpart .gt. 5) then        
        y6=yrap(6,p)
        pt6=pt(6,p)
        y56=yraptwo(5,6,p)
        pt56=pttwo(5,6,p)
        transm=2d0*dsqrt(pttwo(4,5,p)**2+2d0*dot(p,4,5))
        transcm=dsqrt(pttwo(4,5,p)**2+2d0*dot(p,4,5))+pttwo(3,6,p)
        endif

        if (eventpart .gt. 6) then        
        y7=yrap(7,p)
        pt7=pt(7,p)
        endif

         if (eventpart .gt. 7) then        
         y8=yrap(8,p)
         pt8=pt(8,p)
         endif
          
         misset=etmiss(p,etvec)

c--- if we're doing W/Z+2 jets then make
c--- JET 5 = highest Et
c--- JET 6 = lowest Et
      if ((nproc .eq. 22) .or. (nproc .eq. 44)) then
        if (pt6 .gt. pt5) then
          swap=pt5
          pt5=pt6
          pt6=swap
          swap=y5
          y5=y6
          y6=swap
        endif
      endif

      endif

      n=1                  
      call bookplot(n,tag,'      y3',y3,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     pt3',pt3,wt,0d0,150d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'      y4',y4,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     pt4',pt4,wt,0d0,150d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'     y34',y34,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'    pt34',pt34,wt,0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'d/pt34^2',pt34,wt/2d0/pt34*pt34**4,
     .     20d0,200d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'pt34,y=0',pt34a,wt,10d0,150d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'pt34,y=0',pt34a,wt,20d0,480d0,40d0,'log')
      n=n+1
      call bookplot(n,tag,'     m34',m34,wt,20d0,140d0,5d0,'lin')
      n=n+1
      if (eventpart .gt. 4) then
      call bookplot(n,tag,'      y5',y5,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     pt5',pt5,wt,0d0,200d0,5d0,'log')
      n=n+1
      endif
      if (eventpart .gt. 5) then
      call bookplot(n,tag,'      y6',y6,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     pt6',pt6,wt,0d0,200d0,4d0,'log')
      n=n+1
      call bookplot(n,tag,'     y56',y56,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'    pt56',pt56,wt,10d0,150d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'     m56',m56clust,wt,20d0,200d0,4d0,'log')
      n=n+1
      call bookplot(n,tag,'     m56',m56clust,wt,84d0,117d0,3d0,'lin')
      n=n+1
      call bookplot(n,tag,'   m56_5',m56_5,wt,20d0,160d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'  m56_10',m56_10,wt,20d0,160d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'  m56_15',m56_15,wt,20d0,160d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'m56_10_a',m56_10,wt,86d0,114d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'m56_11_b',m56_11,wt,94d0,126d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'m56_12_c',m56_12,wt,103d0,137d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'m56_13_d',m56_13,wt,112d0,148d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'    m345',m345,wt,50d0,250d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'    m346',m346,wt,50d0,250d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'   m3456',m3456,wt,50d0,250d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'  misset',misset,wt,0d0,100d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'     r56',r56,wt,0d0,4d0,.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'     r35',r35,wt,0d0,4d0,.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'     r36',r36,wt,0d0,4d0,.1d0,'lin')
      n=n+1
c      call bookplot(n,tag,'  mttbar',mttbar,wt,300d0,1d3,20d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'  transm',transm,wt,20d0,200d0,20d0,'lin')
c      n=n+1
c      call bookplot(n,tag,' transcm',transcm,wt,20d0,200d0,20d0,'lin')
c      n=n+1
      endif

      if (eventpart .gt. 6) then
      call bookplot(n,tag,'      y7',y7,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     pt7',pt7,wt,0d0,100d0,5d0,'lin')
      n=n+1
      endif      

      if (eventpart .gt. 7) then
      call bookplot(n,tag,'      y8',y8,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     pt8',pt8,wt,0d0,100d0,5d0,'lin')
      endif      



c      if (rn(1) .le. 0.5d0) then
c        call mfill(n,m56,weight/2d0)
c        call mfill(n,m47,weight/2d0)
c      else
c        call mfill(n,m56,weight)
c      endif
c      n=n+1
c      call mfill(n,m345,weight)
c      n=n+1
c-15
c      call mfill(n,m346,weight)
c      n=n+1
c-15
c      if ((m56smw .gt. 84d0) .and. (m56smw .lt. 117d0))
c     . call mfill(n,costh,weight)
c      n=n+1

c      if ((m56smw .gt. 84d0) .and. (m56smw .lt. 117d0))
c     .  call mfill(n,cosnew1,weight)
c      n=n+1

c      if ((m56smw .gt. 84d0) .and. (m56smw .lt. 117d0))
c     .  call mfill(n,cosnew2,weight)
c      n=n+1
c-1
c      call mfill(n,cosnew3,weight)
c      n=n+1
c      phi=acos(cosphi)
c      call mfill(n,phi,weight)
c      write(6,*) 'phi',phi
c      pause
   

c      n=n+1
c      call mfill(n,ptbbpair,weight)
c      n=n+1
c      call mfill(n,pt5,weight)
c      n=n+1
c      call mfill(n,pt6,weight)
c      n=n+1
c      call mfill(n,pt7,weight)
c      n=n+1
c--- added by JMC
c      phill=fphi(4,7,p)
c      if (nnproc .eq. 61.or.nnproc .eq.111 .or. nnproc .eq. 181) then
c        call dittdrein(p,4,7,costhdd)
c        call mfill(n,180d0/pi*phill,weight)
c      elseif (nnproc .eq. 71) then
c        call mfill(n,180d0/pi*phill,weight/2d0)
c        phill=fphi(4,5,p)
c      call mfill(n,180d0/pi*phill,weight/2d0)c      
c      elseif (nnproc .eq. 82) then
c        call dittdrein(p,4,5,costhdd)
c        phill=fphi(4,5,p)
c        call mfill(n,180d0/pi*phill,weight)
c      endif
c      if(180d0/pi*phill .gt. 180d0) write(*,*) 'ERROR' 
c      n=n+1
c      thetall=ftheta(4,7,p)
c      if (nnproc .eq. 61.or.nnproc .eq.111 .or. nnproc .eq. 181) then
c        call mfill(n,180d0/pi*thetall,weight)
c      elseif (nnproc .eq. 71) then
c        call mfill(n,180d0/pi*thetall,weight/2d0)
c        thetall=ftheta(4,5,p)
c        call mfill(n,180d0/pi*thetall,weight/2d0)c      
c      elseif (nnproc .eq. 82) then
c        thetall=ftheta(4,5,p)
c        call mfill(n,180d0/pi*thetall,weight)
c      endif
c      if(180d0/pi*thetall .gt. 180d0) write(*,*) 'ERROR' 
c      n=n+1
c      if (nnproc .eq. 61.or.nnproc .eq.111) then
c        mt1=dsqrt(min(mtsqlet(4,5,6,p),mtsqlet(7,5,6,p)))
c        call mfill(n,mt1,weight)
c      elseif (nnproc .eq. 71) then
c        mt1=dsqrt(min(mtsqlet(4,6,0,p),mtsqlet(7,6,0,p)))
c        call mfill(n,mt1,weight/2d0)
c        mt1=dsqrt(min(mtsqlet(4,6,0,p),mtsqlet(5,6,0,p)))
c        call mfill(n,mt1,weight/2d0)
c      elseif (nnproc .eq. 82) then
c        mt1=dsqrt(min(mtsqlet(4,6,7,p),mtsqlet(5,6,7,p)))
c        call mfill(n,mt1,weight)
c      endif
c      n=n+1
c      if (nnproc .eq. 61.or.nnproc .eq.111) then
c        mll=dsqrt(s(4,7))
c        call mfill(n,mll,weight)
c      elseif (nnproc .eq. 71) then
c        call mfill(n,dsqrt(s(4,5)),weight/2d0)
c        call mfill(n,dsqrt(s(4,7)),weight/2d0)
c      elseif (nnproc .eq. 82) then
c        mll=dsqrt(s(4,5))
c        call mfill(n,mll,weight)
c      endif
c      n=n+1c      
c      call mfill(n,costhdd,weight)
c      n=n+1c      
c      if (nnproc .eq. 61.or.nnproc .eq.111) then
c        cosllet=coslpairet(4,7,5,6,p)
c        call mfill(n,cosllet,weight)
c      elseif (nnproc .eq. 82) then
c        cosllet=coslpairet(4,5,6,7,p)
c        call mfill(n,cosllet,weight)
c      endif
c      n=n+1c      

c      call mfill(n,pt5+pt6,weight)
c      n=n+1
c      call mfill(n,transcm,weight)
c      n=n+1
c      call mfill(n,nchi,weight)
c      n=n+1
c      call mfill(n,m56psm20,weight)
c      n=n+1
c      call mfill(n,m56psm40,weight)
c      n=n+1
c      call mfill(n,pt34,weight)
c      n=n+1
c      call mfill(n,pt34,weight)
c      n=n+1
c      if (switch .eq. 1 ) then
c      kt12=sqrt(s(1,3)*s(2,3)/s(1,2))
c      kt14=sqrt(s(1,3)*s(4,3)/s(1,4))
c      kt15=sqrt(s(1,3)*s(5,3)/s(1,5))
c      kt24=sqrt(s(2,3)*s(4,3)/s(2,4))
c      kt25=sqrt(s(2,3)*s(5,3)/s(2,5))
c      kt56=sqrt(s(4,3)*s(5,3)/s(4,5))

c      call mfill(n,-s(1,3),weight)
c      n=n+1
c      root=sqrt(-s(1,3))
c      wt1=root*weight
c      call mfill(n,-s(1,3),root*wt1)
c      n=n+1
c      call mfill(n,-s(2,3),weight)
c      n=n+1
c      root=sqrt(-s(2,3))
c      wt1=root*weight
c      call mfill(n,-s(2,3),root*wt1)
c      n=n+1
c      call mfill(n,+s(3,4),weight)
c      n=n+1
c      root=sqrt(s(3,4))
c      wt1=root*weight
c      call mfill(n,s(3,4),root*wt1)
c      n=n+1
c      call mfill(n,+s(3,5),weight)
c      n=n+1


      return 
      end

      double precision function pttwo(j,k,p)
      implicit none
      include 'constants.f'
      integer j,k
      double precision p(mxpart,4)
      pttwo=dsqrt((p(j,1)+p(k,1))**2+(p(j,2)+p(k,2))**2)
      return
      end

      double precision function yraptwo(j,k,p)
      implicit none
      include 'constants.f'
      integer j,k
      double precision p(mxpart,4)
      yraptwo=(p(j,4)+p(k,4)+p(j,3)+p(k,3))
     .       /(p(j,4)+p(k,4)-p(j,3)-p(k,3))
      if (yraptwo .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yraptwo=100d0
      else 
      yraptwo=0.5d0*dlog(yraptwo)
      endif
      return
      end

c--- this is trying to be decay angle, but need to unboost i back to (jk)
c--- rest frame before doing this
      double precision function sdot30m(i,j,k,p)
      implicit none
      include 'constants.f'
      integer i,j,k,n
      double precision p(mxpart,4),n1,n2
      
      sdot30m=0d0
      n1=0d0
      n2=0d0
      do n=1,3
      sdot30m=p(i,n)*(p(j,n)+p(k,n))
      n1=n1+p(i,n)**2
      n2=n2+(p(j,n)+p(k,n))**2
      enddo
      
      sdot30m=sdot30m/sqrt(n1*n2)
      
      return
      end
     
      double precision function cosnchi(i,j,k,l,p)
      implicit none
      include 'constants.f'
      integer i,j,k,l,n
      double precision p(mxpart,4),cr1(3),cr2(3),n1,n2
      
      call cross(p,i,j,cr1)
      call cross(p,k,l,cr2)
      
      n1=0d0
      n2=0d0
      cosnchi=0d0
      do n=1,3
        n1=n1+cr1(n)**2
        n2=n2+cr2(n)**2
        cosnchi=cosnchi+cr1(n)*cr2(n)
      enddo

      cosnchi=cosnchi/dsqrt(n1*n2)
      
      return
      end
      
      subroutine cross(p,i,j,r)
      implicit none
      include 'constants.f'
      integer i,j
      double precision p(mxpart,4),r(3)
      
      r(1)=p(i,2)*p(j,3)-p(j,2)*p(i,3)
      r(2)=p(i,3)*p(j,1)-p(j,3)*p(i,1)
      r(3)=p(i,1)*p(j,2)-p(j,1)*p(i,2)
      
      return
      end
            
      double precision function smearp(i,j,p,sd)
      implicit none
      include 'constants.f'
      include 'masses.f'     
      integer i,j,k,idum
      double precision p(mxpart,4),r1(4),r2(4),gasdev,sm1,sm2,sd
      data idum/56735345/

      sm1=1d0+gasdev(idum)/sd
      sm2=1d0+gasdev(idum)/sd

      do k=1,4
        r1(k)=p(i,k)*sm1
        r2(k)=p(j,k)*sm2
      enddo
      
      smearp=sqrt(2d0*(r1(4)*r2(4)-r1(1)*r2(1)-r1(2)*r2(2)-r1(3)*r2(3))
     . +mb**2)
      
      return
      end
      
      double precision function fphi(n1,n2,p)
      implicit none
      include 'constants.f'
      integer n1,n2
      double precision p(mxpart,4)
    
      fphi=p(n1,1)*p(n2,1)+p(n1,2)*p(n2,2)
      fphi=fphi/dsqrt(p(n1,1)**2+p(n1,2)**2)
      fphi=fphi/dsqrt(p(n2,1)**2+p(n2,2)**2)
      if (fphi .gt. 1d0) then
        fphi=0d0
      elseif (fphi .lt. -1d0) then
        fphi=pi
      else
        fphi=dacos(fphi)
      endif

      return
      end
          
      double precision function ftheta(n1,n2,p)
      implicit none
      include 'constants.f'
      integer n1,n2
      double precision p(mxpart,4)
    
      ftheta=p(n1,1)*p(n2,1)+p(n1,2)*p(n2,2)+p(n1,3)*p(n2,3)
      ftheta=ftheta/dsqrt(p(n1,1)**2+p(n1,2)**2+p(n1,3)**2)
      ftheta=ftheta/dsqrt(p(n2,1)**2+p(n2,2)**2+p(n2,3)**2)
      if (ftheta .gt. 1d0) then
        ftheta=0d0
      elseif (ftheta .lt. -1d0) then
        ftheta=pi
      else
        ftheta=dacos(ftheta)
      endif
   
      return
      end
      
      double precision function mtsqlet(n,nm1,nm2,p)
      implicit none
      include 'constants.f'
      integer n,nm1,nm2,i
      double precision p(mxpart,4),misset(4),etmiss,coslet,pt,pttwo
      
      if (nm2. eq. 0) then
        etmiss=pt(nm1,p)
        do i=1,4
          misset(i)=p(nm1,i)
        enddo
      else
        etmiss=pttwo(nm1,nm2,p)
        do i=1,4
          misset(i)=p(nm1,i)+p(nm2,i)
        enddo
      endif
      
      coslet=p(n,1)*misset(1)+p(n,2)*misset(2)
      coslet=coslet/dsqrt(p(n,1)**2+p(n,2)**2)
      coslet=coslet/dsqrt(misset(1)**2+misset(2)**2)
      mtsqlet=2d0*pt(n,p)*etmiss*(1d0-coslet)
      
      return
      end
            
      double precision function coslpairet(n1,n2,nm1,nm2,p)
      implicit none
      include 'constants.f'
      integer n1,n2,nm1,nm2,i
      double precision p(mxpart,4),misset(4),pp(4),coslet,pt,pttwo
      
      if (nm2. eq. 0) then
        do i=1,4
          misset(i)=p(nm1,i)
        enddo
      else
        do i=1,4
          misset(i)=p(nm1,i)+p(nm2,i)
        enddo
      endif
      
      do i=1,4
        pp(i)=p(n1,i)+p(n2,i)
      enddo
      
      coslpairet=pp(1)*misset(1)+pp(2)*misset(2)
      coslpairet=coslpairet/dsqrt(pp(1)**2+pp(2)**2)
      coslpairet=coslpairet/dsqrt(misset(1)**2+misset(2)**2)
      
      return
      end
            
      double precision function deltar(i,j,p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),phi1,phi2,yrap,dphi
      integer i,j
      
      phi1=atan2(p(i,1),p(i,2))
      phi2=atan2(p(j,1),p(j,2))
      dphi=phi1-phi2
      if (dphi .gt. pi) dphi=twopi-dphi
      if (dphi .lt. -pi) dphi=twopi+dphi
      deltar=(yrap(i,p)-yrap(j,p))**2+dphi**2
      deltar=dsqrt(deltar)
      
      return
      end
      
      double precision function deltarpjet(i,j,p,pjet)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pjet(mxpart,4),phi1,phi2,yrap,dphi
      integer i,j
      
      phi1=atan2(p(i,1),p(i,2))
      phi2=atan2(pjet(j,1),pjet(j,2))
      dphi=phi1-phi2
      if (dphi .gt. pi) dphi=twopi-dphi
      if (dphi .lt. -pi) dphi=twopi+dphi
      deltarpjet=(yrap(i,p)-yrap(j,p))**2+dphi**2
      deltarpjet=dsqrt(deltarpjet)
      
      return
      end
      
      double precision function ptqfour(q,j,k,l,m)
      implicit none
      include 'constants.f'
      integer j,k,l,m
      double precision q(mxpart,4)
      ptqfour=dsqrt((q(j,1)+q(k,1)+q(l,1)+q(m,1))**2
     .             +(q(j,2)+q(k,2)+q(l,2)+q(m,2))**2)
      return
      end

          
