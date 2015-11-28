      double precision function msqzbb(i1,i2,i5,i6)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer i1,i2,i5,i6
      double complex aqqb_zbb
      double precision faclo
      faclo=4d0*V*gsq**2*esq**2*aveqq
      msqzbb=faclo
     &*(abs(aqqb_zbb(i1,i2,3,4,i5,i6))**2
     & +abs(aqqb_zbb(i1,i2,3,4,i6,i5))**2)
      return 
      end
