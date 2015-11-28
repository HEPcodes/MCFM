      subroutine storedip(msq_dip,msq_dipv,dsub,dsubv,
     .                    sub_dip,sub_dipv,n)
c--- this routine transfers the information on the colour
c--- structure from a common block into separate arrays for
c--- each parton configuration
      implicit none
      include 'constants.f'
      integer i,j,k,n
      double precision msq_dip(6,0:2,-nf:nf,-nf:nf),dsub(4),sub_dip(6,4)
     .                ,msq_dipv(6,0:2,-nf:nf,-nf:nf),dsubv,sub_dipv(6)
      double precision msq_cs(0:2,-nf:nf,-nf:nf),mmsq_cs(0:2,2,2),
     .                msqv_cs(0:2,-nf:nf,-nf:nf),mmsqv_cs(0:2,2,2)
      common/msq_cols/msq_cs,mmsq_cs
      common/msqv_cols/msqv_cs,mmsqv_cs
      
      do i=0,2
        do j=-nf,nf
        do k=-nf,nf
          msq_dip(n,i,j,k)=msq_cs(i,j,k)
          msq_dipv(n,i,j,k)=msqv_cs(i,j,k)
        enddo
        enddo
      enddo
      
      do i=1,4
        sub_dip(n,i)=dsub(i)
      enddo
      sub_dipv(n)=dsubv
      
      return
      end
