      subroutine fill_gnodkbasise(
     . gnodkbasisae0,gnodkbasisae1,gnodkbasisae2,
     . gnodkbasisbe0,gnodkbasisbe1,gnodkbasisbe2,
     . gnodkbasisce0,gnodkbasisce1,gnodkbasisce2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      double precision
     . gnodkbasisae0,gnodkbasisae1,gnodkbasisae2,
     . gnodkbasisbe0,gnodkbasisbe1,gnodkbasisbe2,
     . gnodkbasisce0,gnodkbasisce1,gnodkbasisce2

      gnodkbasisae0=  + ro*xnsq * (  - 8.D0 + 4.D0*t1**(-1) + 4.D0*
     &    t2**(-1) )
      gnodkbasisae0 = gnodkbasisae0 + ro * (  - 4.D0*t1**(-1) - 4.D0*
     &    t2**(-1) )
      gnodkbasisae0 = gnodkbasisae0 + ro**2*xnsq * (  - t1**(-2) - 
     &    t2**(-2) )
      gnodkbasisae0 = gnodkbasisae0 + ro**2 * ( t1**(-2) + 2.D0*
     &    t1**(-1) + t2**(-2) + 2.D0*t2**(-1) )
      gnodkbasisae0 = gnodkbasisae0 + xnsq * (  - 16.D0 + 4.D0*t1**(-1)
     &     + 16.D0*t1 - 16.D0*t1**2 + 4.D0*t2**(-1) )
      gnodkbasisae0 = gnodkbasisae0 + 8.D0 - 4.D0*t1**(-1) - 4.D0*
     &    t2**(-1)

      gnodkbasisae1=  + xnsq * ( 24.D0 - 8.D0*t1**(-1) - 16.D0*t1 + 16.D
     &    0*t1**2 - 8.D0*t2**(-1) )
      gnodkbasisae1 = gnodkbasisae1 - 8.D0 + 8.D0*t1**(-1) + 8.D0*
     &    t2**(-1)

      gnodkbasisae2=  + xnsq * (  - 8.D0 + 4.D0*t1**(-1) + 4.D0*
     &    t2**(-1) )
      gnodkbasisae2 = gnodkbasisae2 - 4.D0*t1**(-1) - 4.D0*t2**(-1)

      gnodkbasisbe0=  + ro * ( 4.D0*t1**(-1) + 4.D0*t2**(-1) )
      gnodkbasisbe0 = gnodkbasisbe0 + ro**2 * (  - t1**(-2) - 2.D0*
     &    t1**(-1) - t2**(-2) - 2.D0*t2**(-1) )
      gnodkbasisbe0 = gnodkbasisbe0 - 8.D0 + 4.D0*t1**(-1) + 4.D0*
     &    t2**(-1)

      gnodkbasisbe1=  + 8.D0 - 8.D0*t1**(-1) - 8.D0*t2**(-1)

      gnodkbasisbe2=  + 4.D0*t1**(-1) + 4.D0*t2**(-1)

      gnodkbasisce0=  + ro*xnsq * ( 4.D0*t1**(-1) - 4.D0*t2**(-1) )
      gnodkbasisce0 = gnodkbasisce0 + ro**2*xnsq * (  - t1**(-2) + 
     &    t2**(-2) )
      gnodkbasisce0 = gnodkbasisce0 + xnsq * (  - 8.D0 + 4.D0*t1**(-1)
     &     + 16.D0*t1 - 4.D0*t2**(-1) )

      gnodkbasisce1=  + xnsq * ( 8.D0 - 8.D0*t1**(-1) - 16.D0*t1 + 8.D0
     &    *t2**(-1) )

      gnodkbasisce2=  + xnsq * ( 4.D0*t1**(-1) - 4.D0*t2**(-1) )

      return
      end
