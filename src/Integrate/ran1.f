C  (C) Copr. 1986-92 Numerical Recipes Software ]2w.1,r1..

      double precision FUNCTION ran1(idum)
      implicit none
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      double precision AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1d0/IM,IQ=127773,
     .IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=3d-16,RNMX=1d0-EPS)
      INTEGER j,k,iv(NTAB),iy
      DATA iv /NTAB*0/, iy /0/
      SAVE iv,iy
      
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*dble(iy),RNMX)

      return
      END
