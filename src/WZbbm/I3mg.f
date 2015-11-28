      double precision function I3mg(s,msq)
C     expression for (-1)^n*C0(msq,msq,s;,0,msq,0) in s<0,msq>0 region
      include 'constants.f' 
      double precision s,msq,beta,lp,lm,rat,ddilog
      omro=1d0-4d0*msq/s
      if (omro .gt. 0d0) then
      beta=sqrt(omro)
      lp=0.5d0+0.5d0*beta
      lm=1d0-lp
      rat=-lm/lp
      write(6,*) 'beta,lp,lm',beta,lp,lm
      else
      write(6,*) 'beta imaginary in I3mg.f'
      stop
      endif
      I3mg=(-4d0*pisqo6-2d0*ddilog(rat)-0.5d0*log(rat)**2)/(s*beta)
      return
      end
