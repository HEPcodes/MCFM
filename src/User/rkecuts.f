      logical function madezbbnncuts(p)
C--returns true if the event is to be passed
C--cuts appropriate for Z(nu nu) b bbar
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),etvec(4),ptminb,ymaxb,misset,etmiss,
     . missetmin,ptb,ptbb,yb,ybb,pt,ayrap,ptemin,yemax,ptjetmin,yjetmax
      logical jet7,jet8,lepton4
      integer nproc
      common/nproc/nproc
      parameter (ptminb=15d0,ymaxb=2d0,missetmin=35d0)
      parameter (ptemin=20d0,yemax=2.5d0)
      parameter (ptjetmin=8d0,yjetmax=2d0)
      
      madezbbnncuts=.false.
      misset=etmiss(p,etvec)
      ptb=pt(5,p)
      yb=ayrap(5,p)
      ptbb=pt(6,p)
      ybb=ayrap(6,p)

      lepton4=.false.
      jet7=.false.
      jet8=.false.
      if (nproc .eq. 151) then
C      ttbar> w^+b bbar   lepton^-1 n
      jet7=((pt(7,p) .gt. ptjetmin) .and. (ayrap(7,p).lt. yjetmax)) 
      lepton4=((pt(4,p) .gt. ptemin) .and. (ayrap(4,p) .lt. yemax)) 
      elseif (nproc .eq. 152) then
C      ttbar> w^+b bbar q qbar
      lepton4=((pt(4,p) .gt. ptemin) .and. (ayrap(4,p) .lt. yemax)) 
      jet7=((pt(7,p) .gt. ptjetmin) .and. (ayrap(7,p).lt. yjetmax)) 
      jet8=((pt(8,p) .gt. ptjetmin) .and. (ayrap(8,p).lt. yjetmax)) 
      elseif (nproc .eq. 161) then
C      qg> w^+b bbar q'
      lepton4=((pt(4,p) .gt. ptemin) .and. (ayrap(4,p) .lt. yemax)) 
      jet7=((pt(7,p) .gt. ptjetmin) .and. (ayrap(7,p).lt. yjetmax)) 
      elseif (nproc .eq. 171) then
C      q qbar> w^+b bbar
      lepton4=((pt(4,p) .gt. ptemin) .and. (ayrap(4,p) .lt. yemax)) 
      elseif (nproc .eq. 21 ) then
      lepton4=(pt(4,p) .gt. ptjetmin).and.(ayrap(4,p).lt. yjetmax) 
      elseif (nproc .eq. 73 ) then
      lepton4=(pt(4,p) .gt. ptjetmin).and.(ayrap(4,p).lt.yjetmax) 
      endif

C--- Cut on pt and rapidities of b-jets
      if (   (ptb .gt. ptminb) 
     . .and. (ptbb .gt. ptminb)
     . .and. (yb  .lt. ymaxb)
     . .and. (ybb .lt. ymaxb)
     . .and. (misset .gt. missetmin)
     . .and. (jet7 .eqv. .false.) 
     . .and. (jet8 .eqv. .false.) 
     . .and. (lepton4 .eqv. .false.) 
     . ) madezbbnncuts=.true.


      return
      end
