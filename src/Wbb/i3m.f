      double complex function I3m(s1,s2,s3)
      implicit none
      include 'constants.f'
      double precision s1,s2,s3,smax,smid,smin,del3,rtdel3
      double precision i3m1a,flag
      double complex i3m1b

      smax=max(s1,s2,s3)
      smin=min(s1,s2,s3)
      smid=s1+s2+s3-smax-smin
      del3=s1**2+s2**2+s3**2-two*(s1*s2+s2*s3+s3*s1)      

      if (del3 .gt. 0) then
      rtdel3=sqrt(del3)
         if (smax .lt. 0) then
c---case all negative
             flag=0
             i3m=i3m1b(smax,smid,smin,rtdel3,flag)
         elseif (smin .gt. 0) then
c---case all positive
             flag=0
             i3m=-i3m1b(-smin,-smid,-smax,rtdel3,flag)
         elseif ((smid .lt. 0) .and. (smin .lt. 0)) then
c---case two negative and one positive
             flag=+1
             i3m=i3m1b(smin,smid,smax,rtdel3,flag)
         elseif ((smax .gt. 0).and.(smid .gt. 0)) then
c---case two positive and one negative
             flag=-1
             i3m=-i3m1b(-smax,-smid,-smin,rtdel3,flag)
         endif
      elseif (del3 .lt. 0) then 
      rtdel3=sqrt(-del3)
         if (smax .lt. 0) then
c---case all negative
             i3m=+dcmplx(i3m1a(+s1,+s2,+s3,rtdel3))
         elseif (smin .gt. 0) then
c---case all positive
             i3m=-dcmplx(i3m1a(-s1,-s2,-s3,rtdel3))  
          endif
      endif

      return
      end     


      double precision function I3m1a(s1,s2,s3,rtmdel)
      implicit none
      include 'constants.f'
      double precision s1,s2,s3,d1,d2,d3,rtmdel,arg1,arg2,arg3,dclaus

      d1=s1-s2-s3
      d2=s2-s3-s1
      d3=s3-s1-s2
      
      arg1=two*datan(rtmdel/d1)
      arg2=two*datan(rtmdel/d2)
      arg3=two*datan(rtmdel/d3)
      i3m1a=two/rtmdel*(Dclaus(arg1)+Dclaus(arg2)+Dclaus(arg3))

      end

      
      double complex function I3m1b(s1,s2,s3,rtdel,flag)
      implicit none
      include 'constants.f'
      double precision s1,s2,s3,d3,temp,ddilog,xlog,ylog
      double precision x,y,rho,rtdel,argx,argy,argdlx,argdly,flag
      d3=s3-s1-s2
      x=s1/s3
      y=s2/s3
      rho=two*s3/(d3+rtdel)
      argx=rho*x
      argy=rho*y
      argdlx=-argx
      argdly=-argy

      if ((argdlx .gt. 1) .or. (argdly .gt. 1)) then
      write(6,*) 'problems with call of I3m1b'
      write(6,*) 'argdlx',argdlx
      write(6,*) 'argdly',argdly
      stop
      endif

      xlog=log(abs(argx))
      ylog=log(abs(argy))
      temp=xlog*ylog+pisq/3d0+(ylog-xlog)*log((one+argy)/(one+argx))
     & +two*(ddilog(argdlx)+ddilog(argdly))
      I3m1b=temp-abs(flag)*pisq+flag*impi*(xlog+ylog)
      I3m1b=-I3m1b/rtdel
      end





