      subroutine fill_gnodkbasis(
     . gnodkbasis1,gnodkbasis2,gnodkbasis3,
     . gnodkbasis4,gnodkbasis5,gnodkbasis6,
     . gnodkbasis7,gnodkbasis8,gnodkbasis9,gnodkbasis10)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      double precision
     . gnodkbasis1,gnodkbasis2,gnodkbasis3,
     . gnodkbasis4,gnodkbasis5,gnodkbasis6,
     . gnodkbasis7,gnodkbasis8,gnodkbasis9,gnodkbasis10

      gnodkbasis1=  - 4.D0*t1**(-2)*ro + 8.D0*t1**(-1)

      gnodkbasis2=  - 8.D0 - 2.D0*t1**(-2)*ro**2 + 8.D0*t1**(-1)*ro + 2.
     & D0*t1**(-1)*ro**2 + 48.D0*t1 - 32.D0*t1**2 - 4.D0*t2**(-1)*ro + 
     & 2.D0*t2**(-1)*ro**2 - 16.D0*ro

      gnodkbasis3= 8.D0 - 4.D0*t1**(-2)*ro + 2.D0*t1**(-2)*ro**2 - 4.D0
     & *t1**(-1)*ro - 2.D0*t1**(-1)*ro**2 + 8.D0*t1**(-1) - 48.D0*t1 + 
     & 32.D0*t1**2 - 2.D0*t2**(-1)*ro**2 + 16.D0*ro

      gnodkbasis4= 2.D0*t1**(-2)*ro - 2.D0*t1**(-2)*ro**2 + 4.D0*
     & t1**(-1)*ro - 2.D0*t1**(-1)*ro**2 + 16.D0*t1 - 4.D0*t2**(-1)*ro
     &  - 2.D0*t2**(-1)*ro**2

      gnodkbasis5= 2.D0*t1**(-2)*ro**2 - 6.D0*t1**(-1)*ro + 2.D0*
     & t1**(-1)*ro**2 - 16.D0*t1 + 2.D0*t2**(-1)*ro + 2.D0*t2**(-1)*
     & ro**2

      gnodkbasis6=  - 4.D0*t1**(-1)*ro**2 - 32.D0*t1 - 4.D0*t2**(-2)*
     & ro**2 + 16.D0*t2**(-1)*ro - 4.D0*t2**(-1)*ro**2 + 8.D0*t2**(-1)

      gnodkbasis7=  - 4.D0*t1**(-2)*ro - 4.D0*t1**(-1)*ro + 16.D0*
     & t1**(-1) - 4.D0*t2**(-1)*ro

      gnodkbasis8=  - 32.D0 + 2.D0*t1**(-2)*ro - 4.D0*t1**(-2)*ro**2 + 
     & 26.D0*t1**(-1)*ro - 4.D0*t1**(-1)*ro**2 - 8.D0*t1**(-1) + 32.D0*
     & t1 + 10.D0*t2**(-1)*ro - 4.D0*t2**(-1)*ro**2 - 8.D0*t2**(-1)

      gnodkbasis9= 32.D0 - 2.D0*t1**(-2)*ro + 4.D0*t1**(-2)*ro**2 - 18.D
     & 0*t1**(-1)*ro + 4.D0*t1**(-1)*ro**2 - 32.D0*t1 - 2.D0*t2**(-1)*
     & ro + 4.D0*t2**(-1)*ro**2

      gnodkbasis10=  - 4.D0*t1**(-1)*ro + 8.D0*t1**(-1)*ro**2 + 64.D0*
     & t1 + 8.D0*t2**(-2)*ro**2 - 36.D0*t2**(-1)*ro + 8.D0*t2**(-1)*
     & ro**2 - 16.D0*t2**(-1)

      return
      end
