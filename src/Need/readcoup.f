      subroutine readcoup
c--- reads in the anomalous couplings from the file anomcoup.DAT
      implicit none
      include 'anomcoup.f'
      logical newinput
      common/newinput/newinput
      
      if (newinput) return      
      
      open(unit=21,file='anomcoup.DAT',status='old',err=999)
      call checkversion(21,'anomcoup.DAT')
      read(21,*) delg1_z
      read(21,*) delk_z
      read(21,*) delk_g
      read(21,*) lambda_z
      read(21,*) lambda_g
      read(21,*) tevscale
      close(21)
c--- E-M gauge invariance requires that delg1_g=0
      delg1_g=0d0
      
      write(6,*)
      write(6,*)  '*************** Anomalous couplings ****************'
      write(6,*)  '*                                                  *'
      write(6,99) '*            Delta_g1(Z)  =  ',delg1_z,
     .                '                *'
      write(6,99) '*            Delta_g1(g)  =  ',0d0,
     .                '                *'
      write(6,99) '*            Delta_K(Z)   =  ',delk_z,
     .                '                *'
      write(6,99) '*            Delta_K(g)   =  ',delk_g,
     .                '                *'
      write(6,99) '*            Lambda(Z)    =  ',lambda_z,
     .                '                *'
      write(6,99) '*            Lambda(g)    =  ',lambda_g,
     .                '                *'
      write(6,99) '*            TeV-scale    =  ',tevscale,
     .                ' TeV            *'
      write(6,*)  '****************************************************'
      
      return
 
   99 format(1x,a29,f6.2,a17)
   
  999 write(6,*) 'Error reading anomcoup.DAT'
      call flush(6)
      stop
 
      end
