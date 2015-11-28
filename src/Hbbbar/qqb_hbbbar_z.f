      subroutine qqb_hbbbar_z(p,z)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR.f'
      double precision z,omz,lgomz,p(mxpart,4),xl12,dot

      omz=one-z
      lgomz=log(omz)
      xl12=log(two*dot(p,1,2)/musq)

      Rgg_g=ason2pi*2d0*xn*((2d0*lgomz+xl12-log(z))
     . *(one/z-2d0+z*omz)-log(z)/omz)
      Pgg_g=ason2pi*2d0*xn*(2d0*lgomz/omz+xl12/omz)

      Rg_gg=Rgg_g
      Pg_gg=Pgg_g

      Rqg_g=ason2pi*cf*((2d0*lgomz+xl12-log(z))*(1d0+omz**2)/z+z**2)
      Rg_qg=Rqg_g
      
      return
      end
