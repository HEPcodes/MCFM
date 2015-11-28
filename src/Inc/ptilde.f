      integer maxd,ndmax
C----maxd=The maximum possible number of dipoles
C----ndmax=The maximum number of dipoles for the problem at hand
      parameter (maxd=16)
      double precision ptilde(maxd,mxpart,4)
      common/ptildes/ptilde,ndmax
