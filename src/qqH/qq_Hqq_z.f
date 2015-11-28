      subroutine qq_Hqq_z(p,z)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer is
      double precision z,xl15,xl26,tempqq1,tempqq2,tempgq,tempqg,
     . p(mxpart,4),dot,if_qq,fi_qq,if_qg

      xl15=dlog(-two*dot(p,1,5)/musq)
      xl26=dlog(-two*dot(p,2,6)/musq)
c----contributions for one leg

      do is=1,3
      tempqq1=+ason2pi*cf*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
      tempqq2=+ason2pi*cf*(if_qq(z,xl26,is)+fi_qq(z,xl26,is))
      tempgq =+ason2pi*tr*if_qg(z,xl15,is)
      tempqg =+ason2pi*tr*if_qg(z,xl26,is)

      Q1(q,q,q,is)=tempqq1
      Q2(q,q,q,is)=tempqq2
      Q1(a,a,a,is)=tempqq1
      Q2(a,a,a,is)=tempqq2

      Q1(q,q,a,is)=tempqq1
      Q2(q,q,a,is)=tempqq2
      Q1(a,a,q,is)=tempqq1
      Q2(a,a,q,is)=tempqq2

      Q1(q,g,q,is)=tempgq
      Q1(a,g,a,is)=tempgq
      Q2(q,g,q,is)=tempqg
      Q2(a,g,a,is)=tempqg

      Q1(q,g,a,is)=tempgq
      Q1(a,g,q,is)=tempgq
      Q2(q,g,a,is)=tempqg
      Q2(a,g,q,is)=tempqg

      enddo

      return
      end
