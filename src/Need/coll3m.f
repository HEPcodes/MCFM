      subroutine coll3m(p,i,jpt)
      implicit none
      include 'constants.f'
      integer i,jpt,ii1(2),ii2(2),ii3(2),ii4(2),ii5(2),
     .              if1(1),if2(1),if3(1),if4(1),if5(1),
     .        iorf(5),pntr(5),ip
      double precision p(mxpart,4),dot,
     . pi1(4),pi2(4),pi3(4),pi4(4),pi5(4),pi6(4),
     . qi1(4),qi2(4),qi3(4),qi4(4),qi5(4),qi6(4),
     . ri1(4),ri2(4),ri3(4),ri4(4),ri5(4),ri6(4),
     . si1(4),si2(4),si3(4),si4(4),si5(4),si6(4),
     . ti1(4),ti2(4),ti3(4),ti4(4),ti5(4),ti6(4),
     . pf1(4),pf2(4),pf3(4),pf4(4),pf5(4),pf6(4),
     . qf1(4),qf2(4),qf3(4),qf4(4),qf5(4),qf6(4),
     . rf1(4),rf2(4),rf3(4),rf4(4),rf5(4),rf6(4),
     . sf1(4),sf2(4),sf3(4),sf4(4),sf5(4),sf6(4),
     . tf1(4),tf2(4),tf3(4),tf4(4),tf5(4),tf6(4)

c--- momenta for (1,6) collinear
      data pi1/     0.00000000000000d0,    0.00000000000000d0,
     .   -176.47075086043955d0, -176.47075086043955d0/
      data pi2/     0.00000000000000d0,    0.00000000000000d0,
     .    647.90174275436141d0, -647.90174275436141d0/
      data pi3/   -12.53683598250216d0,  -58.45127225917848d0,
     .   -283.61342217717402d0,  338.57834650818518d0/
      data pi4/    12.59197664857761d0,   58.58646028953971d0,
     .   -311.40844264418411d0,  362.20318794186812d0/
      data pi5/    -0.05514066607545d0,   -0.13518803036123d0,
     .    123.59087292743629d0,  123.59095916474769d0/

      data qi1/     0.00000000000000d0,    0.00000000000000d0,
     .   -234.74293846602535d0, -234.74293846602535d0/
      data qi2/     0.00000000000000d0,    0.00000000000000d0,
     .    386.05800910905845d0, -386.05800910905845d0/
      data qi3/   -36.89819887624374d0,   10.34106910214875d0,
     .   -143.49928802765712d0,  229.53313585391271d0/
      data qi4/    37.07646362808875d0,  -10.37171317104125d0,
     .   -159.31552775455921d0,  239.76795860369054d0/
      data qi5/    -0.17826475184500d0,    0.03064406889250d0,
     .    151.49974513918323d0,  151.49985311748054d0/

      data ri1/     0.00000000000000d0,    0.00000000000000d0,
     .   -167.27377053763155d0, -167.27377053763155d0/
      data ri2/     0.00000000000000d0,    0.00000000000000d0,
     .    303.63473902174195d0, -303.63473902174195d0/
      data ri3/     9.30416754233175d0,    6.23630168325287d0,
     .    -95.02527120588940d0,  199.44989636521038d0/
      data ri4/    -9.26380284401759d0,   -6.26069938428939d0,
     .   -107.24341377051445d0,  205.55087982557271d0/
      data ri5/    -0.04036469831416d0,    0.02439770103652d0,
     .     65.90771649229345d0,   65.90773336859041d0/

      data si1/     0.00000000000000d0,    0.00000000000000d0,
     .    -81.86687875267464d0,  -81.86687875267464d0/
      data si2/     0.00000000000000d0,    0.00000000000000d0,
     .    380.77001894611789d0, -380.77001894611789d0/
      data si3/     2.15800182215809d0,    3.80852000072804d0,
     .   -144.75677374757220d0,  227.15344008459610d0/
      data si4/    -2.16351447110542d0,   -3.84568431967243d0,
     .   -155.48948797918791d0,  234.13981070032341d0/
      data si5/     0.00551264894732d0,    0.03716431894439d0,
     .      1.34312153331683d0,    1.34364691387299d0/

      data ti1/     0.00000000000000d0,    0.00000000000000d0,
     .    -67.48072175627782d0,  -67.48072175627782d0/
      data ti2/     0.00000000000000d0,    0.00000000000000d0,
     .    767.79008133326272d0, -767.79008133326272d0/
      data ti3/    13.39574032062824d0,   19.31186182723337d0,
     .   -397.55892509875559d0,  435.00631350780395d0/
      data ti4/   -13.52025687472191d0,  -19.31502613408514d0,
     .   -329.33077372389141d0,  373.68385849819072d0/
      data ti5/     0.12451655409367d0,    0.00316430685177d0,
     .     26.58033924566212d0,   26.58063108354588d0/

c--- momenta for (5,6) collinear
    
c--- jpt picks the momentum configuration (1..5)
c--- k picks the collinear pair:
c---  1   2   3   4   5 
c---  15  25  35  45  34

c--- handle permutations (initial,final)
      data ii1/1,2/
      data ii2/2,1/
      data ii3/3,3/
      data ii4/4,4/
      data ii5/5,5/
c--- handle permutations (final,final)
      data if1/1/
      data if2/2/
      data if3/3/
      data if4/4/
      data if5/5/

c-- translate i to (i,f)/(f,f) and perm
      data iorf/1,1,2,1,1/
      data pntr/1,2,1,3,4/

      ip=pntr(i)

      if (iorf(i) .eq. 1) then
      
c--- initial-final case
      if (jpt. eq. 1) then
      p(ii1(ip),1)=pi1(1)
      p(ii1(ip),2)=pi1(2)
      p(ii1(ip),3)=pi1(3)
      p(ii1(ip),4)=pi1(4)

      p(ii2(ip),1)=pi2(1)
      p(ii2(ip),2)=pi2(2)
      p(ii2(ip),3)=pi2(3)
      p(ii2(ip),4)=pi2(4)

      p(ii3(ip),1)=pi3(1)
      p(ii3(ip),2)=pi3(2)
      p(ii3(ip),3)=pi3(3)
      p(ii3(ip),4)=pi3(4)

      p(ii4(ip),1)=pi4(1)
      p(ii4(ip),2)=pi4(2)
      p(ii4(ip),3)=pi4(3)
      p(ii4(ip),4)=pi4(4)

      p(ii5(ip),1)=pi5(1)
      p(ii5(ip),2)=pi5(2)
      p(ii5(ip),3)=pi5(3)
      p(ii5(ip),4)=pi5(4)

      elseif (jpt .eq. 2) then
      p(ii1(ip),1)=qi1(1)
      p(ii1(ip),2)=qi1(2)
      p(ii1(ip),3)=qi1(3)
      p(ii1(ip),4)=qi1(4)

      p(ii2(ip),1)=qi2(1)
      p(ii2(ip),2)=qi2(2)
      p(ii2(ip),3)=qi2(3)
      p(ii2(ip),4)=qi2(4)

      p(ii3(ip),1)=qi3(1)
      p(ii3(ip),2)=qi3(2)
      p(ii3(ip),3)=qi3(3)
      p(ii3(ip),4)=qi3(4)

      p(ii4(ip),1)=qi4(1)
      p(ii4(ip),2)=qi4(2)
      p(ii4(ip),3)=qi4(3)
      p(ii4(ip),4)=qi4(4)

      p(ii5(ip),1)=qi5(1)
      p(ii5(ip),2)=qi5(2)
      p(ii5(ip),3)=qi5(3)
      p(ii5(ip),4)=qi5(4)

      elseif (jpt .eq. 3) then
      p(ii1(ip),1)=ri1(1)
      p(ii1(ip),2)=ri1(2)
      p(ii1(ip),3)=ri1(3)
      p(ii1(ip),4)=ri1(4)

      p(ii2(ip),1)=ri2(1)
      p(ii2(ip),2)=ri2(2)
      p(ii2(ip),3)=ri2(3)
      p(ii2(ip),4)=ri2(4)

      p(ii3(ip),1)=ri3(1)
      p(ii3(ip),2)=ri3(2)
      p(ii3(ip),3)=ri3(3)
      p(ii3(ip),4)=ri3(4)

      p(ii4(ip),1)=ri4(1)
      p(ii4(ip),2)=ri4(2)
      p(ii4(ip),3)=ri4(3)
      p(ii4(ip),4)=ri4(4)

      p(ii5(ip),1)=ri5(1)
      p(ii5(ip),2)=ri5(2)
      p(ii5(ip),3)=ri5(3)
      p(ii5(ip),4)=ri5(4)

      elseif (jpt .eq. 4) then
      p(ii1(ip),1)=si1(1)
      p(ii1(ip),2)=si1(2)
      p(ii1(ip),3)=si1(3)
      p(ii1(ip),4)=si1(4)

      p(ii2(ip),1)=si2(1)
      p(ii2(ip),2)=si2(2)
      p(ii2(ip),3)=si2(3)
      p(ii2(ip),4)=si2(4)

      p(ii3(ip),1)=si3(1)
      p(ii3(ip),2)=si3(2)
      p(ii3(ip),3)=si3(3)
      p(ii3(ip),4)=si3(4)

      p(ii4(ip),1)=si4(1)
      p(ii4(ip),2)=si4(2)
      p(ii4(ip),3)=si4(3)
      p(ii4(ip),4)=si4(4)

      p(ii5(ip),1)=si5(1)
      p(ii5(ip),2)=si5(2)
      p(ii5(ip),3)=si5(3)
      p(ii5(ip),4)=si5(4)

      elseif (jpt .eq. 5) then
      p(ii1(ip),1)=ti1(1)
      p(ii1(ip),2)=ti1(2)
      p(ii1(ip),3)=ti1(3)
      p(ii1(ip),4)=ti1(4)

      p(ii2(ip),1)=ti2(1)
      p(ii2(ip),2)=ti2(2)
      p(ii2(ip),3)=ti2(3)
      p(ii2(ip),4)=ti2(4)

      p(ii3(ip),1)=ti3(1)
      p(ii3(ip),2)=ti3(2)
      p(ii3(ip),3)=ti3(3)
      p(ii3(ip),4)=ti3(4)

      p(ii4(ip),1)=ti4(1)
      p(ii4(ip),2)=ti4(2)
      p(ii4(ip),3)=ti4(3)
      p(ii4(ip),4)=ti4(4)

      p(ii5(ip),1)=ti5(1)
      p(ii5(ip),2)=ti5(2)
      p(ii5(ip),3)=ti5(3)
      p(ii5(ip),4)=ti5(4)
      endif
      write(*,*) 'Small: s(',ii1(ip),',',ii5(ip),') =',
     . 2d0*dot(p,ii1(ip),ii5(ip))
      else
c--- final-final case
      if (jpt. eq. 1) then
      p(if1(ip),1)=pf1(1)
      p(if1(ip),2)=pf1(2)
      p(if1(ip),3)=pf1(3)
      p(if1(ip),4)=pf1(4)

      p(if2(ip),1)=pf2(1)
      p(if2(ip),2)=pf2(2)
      p(if2(ip),3)=pf2(3)
      p(if2(ip),4)=pf2(4)

      p(if3(ip),1)=pf3(1)
      p(if3(ip),2)=pf3(2)
      p(if3(ip),3)=pf3(3)
      p(if3(ip),4)=pf3(4)

      p(if4(ip),1)=pf4(1)
      p(if4(ip),2)=pf4(2)
      p(if4(ip),3)=pf4(3)
      p(if4(ip),4)=pf4(4)

      p(if5(ip),1)=pf5(1)
      p(if5(ip),2)=pf5(2)
      p(if5(ip),3)=pf5(3)
      p(if5(ip),4)=pf5(4)

      elseif (jpt .eq. 2) then
      p(if1(ip),1)=qf1(1)
      p(if1(ip),2)=qf1(2)
      p(if1(ip),3)=qf1(3)
      p(if1(ip),4)=qf1(4)

      p(if2(ip),1)=qf2(1)
      p(if2(ip),2)=qf2(2)
      p(if2(ip),3)=qf2(3)
      p(if2(ip),4)=qf2(4)

      p(if3(ip),1)=qf3(1)
      p(if3(ip),2)=qf3(2)
      p(if3(ip),3)=qf3(3)
      p(if3(ip),4)=qf3(4)

      p(if4(ip),1)=qf4(1)
      p(if4(ip),2)=qf4(2)
      p(if4(ip),3)=qf4(3)
      p(if4(ip),4)=qf4(4)

      p(if5(ip),1)=qf5(1)
      p(if5(ip),2)=qf5(2)
      p(if5(ip),3)=qf5(3)
      p(if5(ip),4)=qf5(4)

      elseif (jpt .eq. 3) then
      p(if1(ip),1)=rf1(1)
      p(if1(ip),2)=rf1(2)
      p(if1(ip),3)=rf1(3)
      p(if1(ip),4)=rf1(4)

      p(if2(ip),1)=rf2(1)
      p(if2(ip),2)=rf2(2)
      p(if2(ip),3)=rf2(3)
      p(if2(ip),4)=rf2(4)

      p(if3(ip),1)=rf3(1)
      p(if3(ip),2)=rf3(2)
      p(if3(ip),3)=rf3(3)
      p(if3(ip),4)=rf3(4)

      p(if4(ip),1)=rf4(1)
      p(if4(ip),2)=rf4(2)
      p(if4(ip),3)=rf4(3)
      p(if4(ip),4)=rf4(4)

      p(if5(ip),1)=rf5(1)
      p(if5(ip),2)=rf5(2)
      p(if5(ip),3)=rf5(3)
      p(if5(ip),4)=rf5(4)

      elseif (jpt .eq. 4) then
      p(if1(ip),1)=sf1(1)
      p(if1(ip),2)=sf1(2)
      p(if1(ip),3)=sf1(3)
      p(if1(ip),4)=sf1(4)

      p(if2(ip),1)=sf2(1)
      p(if2(ip),2)=sf2(2)
      p(if2(ip),3)=sf2(3)
      p(if2(ip),4)=sf2(4)

      p(if3(ip),1)=sf3(1)
      p(if3(ip),2)=sf3(2)
      p(if3(ip),3)=sf3(3)
      p(if3(ip),4)=sf3(4)

      p(if4(ip),1)=sf4(1)
      p(if4(ip),2)=sf4(2)
      p(if4(ip),3)=sf4(3)
      p(if4(ip),4)=sf4(4)

      p(if5(ip),1)=sf5(1)
      p(if5(ip),2)=sf5(2)
      p(if5(ip),3)=sf5(3)
      p(if5(ip),4)=sf5(4)

      elseif (jpt .eq. 5) then
      p(if1(ip),1)=tf1(1)
      p(if1(ip),2)=tf1(2)
      p(if1(ip),3)=tf1(3)
      p(if1(ip),4)=tf1(4)

      p(if2(ip),1)=tf2(1)
      p(if2(ip),2)=tf2(2)
      p(if2(ip),3)=tf2(3)
      p(if2(ip),4)=tf2(4)

      p(if3(ip),1)=tf3(1)
      p(if3(ip),2)=tf3(2)
      p(if3(ip),3)=tf3(3)
      p(if3(ip),4)=tf3(4)

      p(if4(ip),1)=tf4(1)
      p(if4(ip),2)=tf4(2)
      p(if4(ip),3)=tf4(3)
      p(if4(ip),4)=tf4(4)

      p(if5(ip),1)=tf5(1)
      p(if5(ip),2)=tf5(2)
      p(if5(ip),3)=tf5(3)
      p(if5(ip),4)=tf5(4)

      endif

      write(*,*) 'Small: s(',if4(ip),',',if5(ip),') =',
     . 2d0*dot(p,if5(ip),if5(ip))
      endif      

      return
      end
















