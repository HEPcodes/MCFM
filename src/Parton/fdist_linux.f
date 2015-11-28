      subroutine fdist(pdlabel,ih,x,xmu,fx)
      implicit none
      character pdlabel*7
      double precision fx(-5:5),x,xmu
      double precision u_val,d_val,u_sea,d_sea,s_sea,c_sea,b_sea,gluon
      double precision Ctq3df,Ctq4Fn,Ctq5Pdf,Ctq6Pdf
      integer mode,Iprtn,ih,Irt
c---  ih1=+1 proton 
c---  ih1=-1 pbar 

C---set to zero if x out of range
      if (x .ge. 1d0) then
          do Iprtn=-5,5
             fx(Iprtn)=0d0
          enddo
          return
      endif
 
      if     ((pdlabel(1:3) .eq. 'mrs')
     .   .or. (pdlabel(2:4) .eq. 'mrs')) then
             if     (pdlabel .eq. 'mrs0119') then
             mode=1
             call mrst2001(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs0117') then
             mode=2
             call mrst2001(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs0121') then
             mode=3
             call mrst2001(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs01_j') then
             mode=4
             call mrst2001(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif     (pdlabel .eq. 'mrs99_1') then
             mode=1
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_2') then
             mode=2
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_3') then
             mode=3
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_4') then
             mode=4
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_5') then
             mode=5
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_6') then
             mode=6
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_7') then
             mode=7
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_8') then
             mode=8
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_9') then
             mode=9
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs9910') then
             mode=10
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs9911') then
             mode=11
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs9912') then
             mode=12
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98z1') then
             mode=1
             call mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98z2') then
             mode=2 
             call mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98z3') then
             mode=3
             call mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98z4') then
             mode=4
             call mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98z5') then
             mode=5
             call mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98ht') then
             mode=1
             call mrs98ht(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs96r1') then
             mode=1
             call mrs96(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs96r2') then
             mode=2 
             call mrs96(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs96r3') then
             mode=3
             call mrs96(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs96r4') then
             mode=4
             call mrs96(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'hmrs90e') then
             mode=1
             call mrsebh(x,xmu,mode,u_val,d_val,u_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             d_sea=u_sea
             elseif (pdlabel .eq. 'hmrs90b') then
             mode=2
             call mrsebh(x,xmu,mode,u_val,d_val,u_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             d_sea=u_sea
             elseif (pdlabel .eq. 'mrs95ap') then
             mode=20
             call mrseb(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs95_g') then
             mode=21
             call mrseb(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             endif
c-----assign mrs to standard grid
            fx(-5)=b_sea/x
            fx(-4)=c_sea/x
            fx(-3)=s_sea/x
            fx( 0)=gluon/x
            fx(+3)=fx(-3)
            fx(+4)=fx(-4)
            fx(+5)=fx(-5)
            if (ih.eq.1) then      
               fx(1)=(d_val+d_sea)/x
               fx(2)=(u_val+u_sea)/x
               fx(-1)=d_sea/x
               fx(-2)=u_sea/x
            elseif(ih.eq.-1) then      
               fx(-1)=(d_val+d_sea)/x
               fx(-2)=(u_val+u_sea)/x
               fx(+1)=d_sea/x
               fx(+2)=u_sea/x
            endif
      return

      elseif (pdlabel(1:5) .eq. 'cteq3') then
C   1      CTEQ3M   Standard MSbar scheme   0.116
C   3      CTEQ3L   Leading Order           0.116
C   2      CTEQ3D   Standard DIS scheme     0.116
          if (pdlabel .eq. 'cteq3_m') then
             mode=1
          elseif (pdlabel .eq. 'cteq3_l') then
             mode=2
          elseif (pdlabel .eq. 'cteq3_d') then
             mode=3
          endif
             fx(-5)=Ctq3df(mode,-5,x,xmu,Irt)/x
             fx(-4)=Ctq3df(mode,-4,x,xmu,Irt)/x
             fx(-3)=Ctq3df(mode,-3,x,xmu,Irt)/x

             fx(0)=Ctq3df(mode,0,x,xmu,Irt)/x

             fx(+3)=Ctq3df(mode,+3,x,xmu,Irt)/x
             fx(+4)=Ctq3df(mode,+4,x,xmu,Irt)/x
             fx(+5)=Ctq3df(mode,+5,x,xmu,Irt)/x
             if (ih.eq.1) then      
               fx(-1)=Ctq3df(mode,-2,x,xmu,Irt)/x
               fx(-2)=Ctq3df(mode,-1,x,xmu,Irt)/x
               fx(1)=Ctq3df(mode,+2,x,xmu,Irt)/x+fx(-1)
               fx(2)=Ctq3df(mode,+1,x,xmu,Irt)/x+fx(-2)
             elseif(ih.eq.-1) then      
               fx(1)=Ctq3df(mode,-2,x,xmu,Irt)/x
               fx(2)=Ctq3df(mode,-1,x,xmu,Irt)/x
               fx(-1)=Ctq3df(mode,+2,x,xmu,Irt)/x+fx(1)
               fx(-2)=Ctq3df(mode,+1,x,xmu,Irt)/x+fx(2)
             endif
             return


      elseif (pdlabel(1:5) .eq. 'cteq4') then
C   1      CTEQ4M   Standard MSbar scheme   0.116        1.6      cteq4m.tbl
C   2      CTEQ4D   Standard DIS scheme     0.116        1.6      cteq4d.tbl
C   3      CTEQ4L   Leading Order           0.116        1.6      cteq4l.tbl
C   4      CTEQ4A1  Alpha_s series          0.110        1.6      cteq4a1.tbl
C   5      CTEQ4A2  Alpha_s series          0.113        1.6      cteq4a2.tbl
C   6      CTEQ4A3  same as CTEQ4M          0.116        1.6      cteq4m.tbl
C   7      CTEQ4A4  Alpha_s series          0.119        1.6      cteq4a4.tbl
C   8      CTEQ4A5  Alpha_s series          0.122        1.6      cteq4a5.tbl
C   9      CTEQ4HJ  High Jet                0.116        1.6      cteq4hj.tbl
C   10     CTEQ4LQ  Low Q0                  0.114        0.7      cteq4lq.tbl

          if (pdlabel .eq. 'cteq4_m') then
             mode=1
          elseif (pdlabel .eq. 'cteq4_d') then
             mode=2
          elseif (pdlabel .eq. 'cteq4_l') then
             mode=3
          elseif (pdlabel .eq. 'cteq4a1') then
             mode=4
          elseif (pdlabel .eq. 'cteq4a2') then
             mode=5
          elseif (pdlabel .eq. 'cteq4a3') then
             mode=6
          elseif (pdlabel .eq. 'cteq4a4') then
             mode=7
          elseif (pdlabel .eq. 'cteq4a5') then
             mode=8
          elseif (pdlabel .eq. 'cteq4hj') then
             mode=9
          elseif (pdlabel .eq. 'cteq4lq') then
             mode=10
          endif

             fx(-5)=Ctq4Fn(mode,-5,x,xmu)
             fx(-4)=Ctq4Fn(mode,-4,x,xmu)
             fx(-3)=Ctq4Fn(mode,-3,x,xmu)

             fx(0)=Ctq4Fn(mode,0,x,xmu)

             fx(+3)=Ctq4Fn(mode,+3,x,xmu)
             fx(+4)=Ctq4Fn(mode,+4,x,xmu)
             fx(+5)=Ctq4Fn(mode,+5,x,xmu)
             if (ih.eq.1) then      
               fx(1)=Ctq4Fn(mode,+2,x,xmu)
               fx(2)=Ctq4Fn(mode,+1,x,xmu)
               fx(-1)=Ctq4Fn(mode,-2,x,xmu)
               fx(-2)=Ctq4Fn(mode,-1,x,xmu)
             elseif(ih.eq.-1) then      
               fx(1)=Ctq4Fn(mode,-2,x,xmu)
               fx(2)=Ctq4Fn(mode,-1,x,xmu)
               fx(-1)=Ctq4Fn(mode,+2,x,xmu)
               fx(-2)=Ctq4Fn(mode,+1,x,xmu)
             endif
             return

      elseif ((pdlabel(1:5) .eq. 'cteq5') .or. 
     .        (pdlabel(1:4) .eq. 'ctq5')) then
 
             fx(-5)=Ctq5Pdf(-5,x,xmu)
             fx(-4)=Ctq5Pdf(-4,x,xmu)
             fx(-3)=Ctq5Pdf(-3,x,xmu)
 
             fx(0)=Ctq5Pdf(0,x,xmu)

             fx(+3)=Ctq5Pdf(+3,x,xmu)
             fx(+4)=Ctq5Pdf(+4,x,xmu)
             fx(+5)=Ctq5Pdf(+5,x,xmu)

             if (ih.eq.1) then      
               fx(1)=Ctq5Pdf(+2,x,xmu)
               fx(2)=Ctq5Pdf(+1,x,xmu)
               fx(-1)=Ctq5Pdf(-2,x,xmu)
               fx(-2)=Ctq5Pdf(-1,x,xmu)
             elseif(ih.eq.-1) then      
               fx(1)=Ctq5Pdf(-2,x,xmu)
               fx(2)=Ctq5Pdf(-1,x,xmu)
               fx(-1)=Ctq5Pdf(+2,x,xmu)
               fx(-2)=Ctq5Pdf(+1,x,xmu)
             endif
             return


      elseif (pdlabel(1:5) .eq. 'cteq6') then

             fx(-5)=Ctq6Pdf(-5,x,xmu)
             fx(-4)=Ctq6Pdf(-4,x,xmu)
             fx(-3)=Ctq6Pdf(-3,x,xmu)
 
             fx(0)=Ctq6Pdf(0,x,xmu)

             fx(+3)=Ctq6Pdf(+3,x,xmu)
             fx(+4)=Ctq6Pdf(+4,x,xmu)
             fx(+5)=Ctq6Pdf(+5,x,xmu)

             if (ih.eq.1) then      
               fx(1)=Ctq6Pdf(+2,x,xmu)
               fx(2)=Ctq6Pdf(+1,x,xmu)
               fx(-1)=Ctq6Pdf(-2,x,xmu)
               fx(-2)=Ctq6Pdf(-1,x,xmu)
             elseif(ih.eq.-1) then      
               fx(1)=Ctq6Pdf(-2,x,xmu)
               fx(2)=Ctq6Pdf(-1,x,xmu)
               fx(-1)=Ctq6Pdf(+2,x,xmu)
               fx(-2)=Ctq6Pdf(+1,x,xmu)
             endif
             return

c--- NEW ATTEMPT
      elseif (pdlabel(1:5) .eq. 'mtung') then
            if     (pdlabel .eq. 'mtungs1') then
              mode=1
            elseif (pdlabel .eq. 'mtunge1') then
              mode=2
            elseif (pdlabel .eq. 'mtungb1') then
              mode=3
            elseif (pdlabel .eq. 'mtungb2') then
              mode=4
            elseif (pdlabel .eq. 'mtungn1') then
              mode=5
            endif
            call mt(x,xmu,mode,u_val,d_val,
     .               u_sea,s_sea,c_sea,b_sea,gluon) 
            d_sea=u_sea
c-----assign to standard grid
            fx(-5)=b_sea/x
            fx(-4)=c_sea/x
            fx(-3)=s_sea/x
            fx( 0)=gluon/x
            fx(+3)=fx(-3)
            fx(+4)=fx(-4)
            fx(+5)=fx(-5)
            if (ih.eq.1) then      
               fx(1)=(d_val+d_sea)/x
               fx(2)=(u_val+u_sea)/x
               fx(-1)=d_sea/x
               fx(-2)=u_sea/x
            elseif(ih.eq.-1) then      
               fx(-1)=(d_val+d_sea)/x
               fx(-2)=(u_val+u_sea)/x
               fx(+1)=d_sea/x
               fx(+2)=u_sea/x
            endif
 
      else
          write(6,*) 'Unimplemented mrs distribution' 
          write(6,*) 'pdlabel= ',pdlabel
          write(6,*) 'Implemented are: ',
     . 'mrs0119,mrs0177,mrs0121,mrs01_j,',
     . 'mrs99_1,mrs99_2,mrs99_3,mrs99_4,mrs99_5,mrs99_6,',
     . 'mrs99_7,mrs99_8,mrs99_9,mrs9910,mrs9911,mrs9912,',
     . 'mrs98z1,',
     . 'mrs98z2,',
     . 'mrs98z3,',
     . 'mrs98z4,',
     . 'mrs98z5,',
     . 'mrs96r1,',
     . 'mrs96r2,',
     . 'mrs96r3,',
     . 'mrs96r4,',
     . 'hmrs90e,',
     . 'hmrs90b,',
     . 'cteq4_m,',
     . 'cteq5_m,',
     . 'cteq5_d,',
     . 'cteq5_l,',
     . 'cteq5hj,',
     . 'cteq5hq,',
     . 'cteq5f3,',
     . 'cteq5f4,',
     . 'mtungs1,',
     . 'mtunge1,',
     . 'mtungb1,',
     . 'mtungb2,',
     . 'mtungn1'

         stop

      endif      

      end

  

