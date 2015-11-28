      double precision function qqgampsq(i1,i2,i3) 
      implicit none
      include 'constants.f'
      integer i1,i2,i3,j,h1,h3
      double complex samp(2,2),tamp(2,2),qamp(2,2),ramp(2,2)

      double complex ttbqqbsqpp,ttbqqbsqmp,ttbqqbsqpm,ttbqqbsqmm
      double complex ttbqqbrqpp,ttbqqbrqmp,ttbqqbrqpm,ttbqqbrqmm
      double complex ttbqqbtqpp,ttbqqbtqmp,ttbqqbtqpm,ttbqqbtqmm
      double complex ttbqqbqqpp,ttbqqbqqmp,ttbqqbqqpm,ttbqqbqqmm

      do h1=1,2
      do h3=1,2

      if ((h1.eq.1) .and. (h3.eq.1)) then
      samp(h1,h3)=ttbqqbsqmm(i1,i2,i3,4,5,6,7)
      ramp(h1,h3)=ttbqqbrqmm(i1,i2,i3,4,5,6,7)
      tamp(h1,h3)=ttbqqbtqmm(i1,i2,i3,4,5,6,7)
      qamp(h1,h3)=ttbqqbqqmm(i1,i2,i3,4,5,6,7)
      elseif ((h1.eq.1) .and. (h3.eq.2)) then
      samp(h1,h3)=ttbqqbsqmp(i1,i2,i3,4,5,6,7)
      ramp(h1,h3)=ttbqqbrqmp(i1,i2,i3,4,5,6,7)
      tamp(h1,h3)=ttbqqbtqmp(i1,i2,i3,4,5,6,7)
      qamp(h1,h3)=ttbqqbqqmp(i1,i2,i3,4,5,6,7)
      elseif ((h1.eq.2) .and. (h3.eq.1)) then
      samp(h1,h3)=ttbqqbsqpm(i1,i2,i3,4,5,6,7)
      ramp(h1,h3)=ttbqqbrqpm(i1,i2,i3,4,5,6,7)
      tamp(h1,h3)=ttbqqbtqpm(i1,i2,i3,4,5,6,7)
      qamp(h1,h3)=ttbqqbqqpm(i1,i2,i3,4,5,6,7)
      elseif ((h1.eq.2) .and. (h3.eq.2)) then
      samp(h1,h3)=ttbqqbsqpp(i1,i2,i3,4,5,6,7)
      ramp(h1,h3)=ttbqqbrqpp(i1,i2,i3,4,5,6,7)
      tamp(h1,h3)=ttbqqbtqpp(i1,i2,i3,4,5,6,7)
      qamp(h1,h3)=ttbqqbqqpp(i1,i2,i3,4,5,6,7)
      endif

      enddo
      enddo

      qqgampsq=0d0
      do h1=1,2
      do h3=1,2
      qqgampsq=qqgampsq
     . +xn*(
     . +cdabs(qamp(h1,h3))**2+cdabs(ramp(h1,h3))**2
     . +cdabs(samp(h1,h3))**2+cdabs(tamp(h1,h3))**2)
     . +(qamp(h1,h3)+ramp(h1,h3))*Dconjg(samp(h1,h3)+tamp(h1,h3))
     . +(samp(h1,h3)+tamp(h1,h3))*Dconjg(qamp(h1,h3)+ramp(h1,h3))

      enddo
      enddo
      qqgampsq=V*qqgampsq
      return
      end
