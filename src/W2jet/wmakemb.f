      subroutine wmakemb(i1,i2,i3,i4,i5,i6,i7,mb1,mb2)
      implicit none
C     Author: R.K. Ellis, March 2001
C     A subroutine calculating Nagy and Trocsanyi, PRD59 014020 (1999) 
C     Eq. B.56 with a factor of 2*i*e^2*g^3/s removed
      include 'constants.f'
      integer f1,f2,f3,f4,hq,Qh,hg,i1,i2,i3,i4,i5,i6,i7
      double complex mb1(5,5,5,5,2,2,2),mb2(5,5,5,5,2,2,2),
     .               m1_1234(5,5,5,5,2,2,2),m2_1234(5,5,5,5,2,2,2),
     .               m3_1234(5,5,5,5,2,2,2),m4_1234(5,5,5,5,2,2,2),
     .               m1_3412(5,5,5,5,2,2,2),m2_3412(5,5,5,5,2,2,2),
     .               m3_3412(5,5,5,5,2,2,2),m4_3412(5,5,5,5,2,2,2)

      call wmakem(i1,i2,i3,i4,i5,i6,i7,m1_1234,m2_1234,m3_3412,m4_3412)
      do f1=1,5
      do f2=1,5
      do f3=1,5
      do f4=1,5
      do hq=1,2
      do Qh=1,2
      do hg=1,2
      mb1(f1,f2,f3,f4,hq,Qh,hg)=
     . m1_1234(f1,f2,f3,f4,hq,Qh,hg)
     .+m3_3412(f3,f4,f1,f2,Qh,hq,hg)
      mb2(f1,f2,f3,f4,hq,Qh,hg)=
     . m2_1234(f1,f2,f3,f4,hq,Qh,hg)
     .+m4_3412(f3,f4,f1,f2,Qh,hq,hg)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
