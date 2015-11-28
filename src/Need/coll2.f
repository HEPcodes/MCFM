      subroutine coll2(p,i,jpt)
      implicit none
      include 'constants.f'
      integer i,jpt,ii1(6),ii2(6),ii3(6),ii4(6),ii5(6),ii6(6),ii7(6),
     .              if1(3),if2(3),if3(3),if4(3),if5(3),if6(3),if7(3),
     .        iorf(9),pntr(9),ip
      double precision p(mxpart,4),
     . pi1(4),pi2(4),pi3(4),pi4(4),pi5(4),pi6(4),pi7(4),
     . qi1(4),qi2(4),qi3(4),qi4(4),qi5(4),qi6(4),qi7(4),
     . ri1(4),ri2(4),ri3(4),ri4(4),ri5(4),ri6(4),ri7(4),
     . si1(4),si2(4),si3(4),si4(4),si5(4),si6(4),si7(4),
     . ti1(4),ti2(4),ti3(4),ti4(4),ti5(4),ti6(4),ti7(4),
     . pf1(4),pf2(4),pf3(4),pf4(4),pf5(4),pf6(4),pf7(4),
     . qf1(4),qf2(4),qf3(4),qf4(4),qf5(4),qf6(4),qf7(4),
     . rf1(4),rf2(4),rf3(4),rf4(4),rf5(4),rf6(4),rf7(4),
     . sf1(4),sf2(4),sf3(4),sf4(4),sf5(4),sf6(4),sf7(4),
     . tf1(4),tf2(4),tf3(4),tf4(4),tf5(4),tf6(4),tf7(4)

c--- momenta for (1,7) collinear
      data pi1/     0.00000000000000d0,    0.00000000000000d0,
     .   -349.60873219024251d0, -349.60873219024251d0/
      data pi2/     0.00000000000000d0,    0.00000000000000d0,
     .     22.07966336976217d0,  -22.07966336976217d0/
      data pi3/    -2.27335746957215d0,  -25.56841068321375d0,
     .    119.08257820094694d0,  121.81778281529583d0/
      data pi4/   -34.66161710346650d0,   58.69264048350193d0,
     .    102.95665816354408d0,  123.47601875407952d0/
      data pi5/    13.58587946656553d0,  -33.07394875089035d0,
     .     48.34561358455070d0,   60.13119456424645d0/
      data pi6/    23.36604178755529d0,   -0.04179680090137d0,
     .     25.37588706622342d0,   34.49506196817264d0/
      data pi7/    -0.01694668108218d0,   -0.00848424849647d0,
     .     31.76833180521521d0,   31.76833745821023d0/

      data qi1/     0.00000000000000d0,    0.00000000000000d0,
     .   -272.79599308146294d0, -272.79599308146294d0/
      data qi2/     0.00000000000000d0,    0.00000000000000d0,
     .     71.38433088744262d0,  -71.38433088744262d0/
      data qi3/     2.37579503601866d0,  -32.69467345146797d0,
     .     23.32110093016135d0,   40.23008603950615d0/
      data qi4/   -31.25970256614727d0,   58.12630270500311d0,
     .     66.31689061312551d0,   93.56156289455035d0/
      data qi5/    36.56429759162680d0,    3.30390102201961d0,
     .     37.84636564730677d0,   52.72770631311389d0/
      data qi6/    -7.73200265366031d0,  -28.73170274262586d0,
     .    -36.58045293386289d0,   47.15319866534383d0/
      data qi7/     0.05161259216212d0,   -0.00382753292889d0,
     .    110.50775793728958d0,  110.50777005639134d0/

      data ri1/     0.00000000000000d0,    0.00000000000000d0,
     .   -541.97569341200381d0, -541.97569341200381d0/
      data ri2/     0.00000000000000d0,    0.00000000000000d0,
     .     79.18086906887218d0,  -79.18086906887218d0/
      data ri3/    68.22451830549478d0,  -36.94673965220740d0,
     .    159.45152121191538d0,  177.32578516882052d0/
      data ri4/   -28.28025144318884d0,  -56.58109853027554d0,
     .    107.63934424838399d0,  124.84959656644810d0/
      data ri5/  -104.90477179673505d0,  105.57435716787907d0,
     .     57.30186174030460d0,  159.48184660343301d0/
      data ri6/    64.97931800805264d0,  -12.03175439386813d0,
     .     92.95015524303990d0,  114.04738595115285d0/
      data ri7/    -0.01881307362354d0,   -0.01476459152800d0,
     .     45.45194189948774d0,   45.45194819102153d0/

      data si1/     0.00000000000000d0,    0.00000000000000d0,
     .   -247.59848301528530d0, -247.59848301528530d0/
      data si2/     0.00000000000000d0,    0.00000000000000d0,
     .     29.66921098438979d0,  -29.66921098438979d0/
      data si3/    -8.58569124225996d0,  -11.17235858099033d0,
     .     54.20012205172404d0,   56.00168676737719d0/
      data si4/    -6.63838339592560d0,    3.86699769515934d0,
     .     -0.46302414902541d0,    7.69650548289606d0/
      data si5/    -7.46505802757188d0,  -29.89554381890737d0,
     .     84.34703762396867d0,   89.79918366846333d0/
      data si6/    22.58419466030091d0,   37.05471513111246d0,
     .     -0.52710458105796d0,   43.39787553719226d0/
      data si7/     0.10493800545654d0,    0.14618957362590d0,
     .     80.37224108528616d0,   80.37244254374625d0/

      data ti1/     0.00000000000000d0,    0.00000000000000d0,
     .   -242.65652903890091d0, -242.65652903890091d0/
      data ti2/     0.00000000000000d0,    0.00000000000000d0,
     .     66.68715589609039d0,  -66.68715589609039d0/
      data ti3/    -5.33292831741917d0,   -2.77791036483735d0,
     .      6.90284204071415d0,    9.15456927839059d0/
      data ti4/    -4.21865997501724d0,    1.50139963609006d0,
     .     -2.98823220757910d0,    5.38338411953642d0/
      data ti5/    46.75904002216957d0,   88.79591531683016d0,
     .      7.38272284956274d0,  100.62617451449692d0/
      data ti6/   -37.15340285454518d0,  -87.50159475794538d0,
     .    138.37544669231872d0,  167.88290167898145d0/
      data ti7/    -0.05404887518798d0,   -0.01780983013749d0,
     .     26.29659376779404d0,   26.29665534358588d0/

c--- momenta for (5,7) collinear
      data pf1/     0.00000000000000d0,    0.00000000000000d0,
     .   -406.58528138368723d0, -406.58528138368723d0/
      data pf2/     0.00000000000000d0,    0.00000000000000d0,
     .    102.91739913023866d0, -102.91739913023866d0/
      data pf3/   -23.35464687008511d0,  -37.88868825251215d0,
     .     62.89240276152825d0,   77.04833906737196d0/
      data pf4/   -80.35009104097168d0,  -11.67623014480668d0,
     .     -1.80433927253352d0,   81.21408203567492d0/
      data pf5/    22.04650797837901d0,   41.39387135725615d0,
     .     89.38183883877004d0,  100.93866560531529d0/
      data pf6/    43.88145797083443d0,  -62.97409837057458d0,
     .     -0.51967479623178d0,   76.75669013921035d0/
      data pf7/    37.77677196184334d0,   71.14514541063727d0,
     .    153.71765472191558d0,  173.54490366635338d0/

      data qf1/     0.00000000000000d0,    0.00000000000000d0,
     .    -69.93851802535644d0,  -69.93851802535644d0/
      data qf2/     0.00000000000000d0,    0.00000000000000d0,
     .     17.59840230937757d0,  -17.59840230937757d0/
      data qf3/    -4.04943635651898d0,    4.17672586416538d0,
     .     20.47376698170064d0,   21.28422204758401d0/
      data qf4/    -8.85075041071755d0,   -4.20189510942481d0,
     .      4.05844313371689d0,   10.60484163073824d0/
      data qf5/     4.69090439632005d0,   -7.17646483016576d0,
     .      8.21329285132285d0,   11.87284342420179d0/
      data qf6/    -0.35069665669931d0,   20.89041959655122d0,
     .      3.89099605369763d0,   21.25258735671687d0/
      data qf7/     8.55997902761579d0,  -13.68878552112603d0,
     .     15.70361669554087d0,   22.52242587549311d0/

      data rf1/     0.00000000000000d0,    0.00000000000000d0,
     .   -248.03990340131170d0, -248.03990340131170d0/
      data rf2/     0.00000000000000d0,    0.00000000000000d0,
     .     61.87170982661398d0,  -61.87170982661398d0/
      data rf3/    73.02267433694956d0,  -10.03648605278721d0,
     .     32.61509740060927d0,   80.60264634650164d0/
      data rf4/     0.72851350898851d0,   57.36228630304377d0,
     .     20.19000284058075d0,   60.81610671974818d0/
      data rf5/   -17.47766256468838d0,  -21.89561584380416d0,
     .     39.86079796052697d0,   48.72134948823265d0/
      data rf6/   -28.87195235432152d0,    9.02317081487299d0,
     .     30.90422806910142d0,   43.24440491961501d0/
      data rf7/   -27.40157292692819d0,  -34.45335522132539d0,
     .     62.59806730387933d0,   76.52710575382818d0/

      data sf1/     0.00000000000000d0,    0.00000000000000d0,
     .    -28.97247835777579d0,  -28.97247835777579d0/
      data sf2/     0.00000000000000d0,    0.00000000000000d0,
     .    276.98023263367696d0, -276.98023263367696d0/
      data sf3/    16.91618177486238d0,  -19.44983715135695d0,
     .    -12.79251521776157d0,   28.77675827210480d0/
      data sf4/   -28.45786079401319d0,   24.91853350246542d0,
     .   -105.93304990306726d0,  112.48375089163174d0/
      data sf5/     0.54017400656316d0,  -40.68438417132341d0,
     .    -60.04920732657512d0,   72.53563402843599d0/
      data sf6/    10.82806945202589d0,   47.25459190947976d0,
     .    -51.56004513847879d0,   70.77204108453660d0/
      data sf7/     0.17343556056176d0,  -12.03890408926482d0,
     .    -17.67293669001842d0,   21.38452671474363d0/

      data tf1/     0.00000000000000d0,    0.00000000000000d0,
     .   -240.75372938650861d0, -240.75372938650861d0/
      data tf2/     0.00000000000000d0,    0.00000000000000d0,
     .     18.00380642082524d0,  -18.00380642082524d0/
      data tf3/   -24.55369105372509d0,   20.06082113762909d0,
     .     42.77784699085886d0,   53.24720163774862d0/
      data tf4/    29.01769011330057d0,  -18.29322144283041d0,
     .     81.48171252674526d0,   88.40779245382799d0/
      data tf5/   -19.49737859048144d0,   -5.73390382973564d0,
     .     77.35082370905191d0,   79.97609238702556d0/
      data tf6/    19.08778659029235d0,    5.23311186495417d0,
     .      4.96769373241632d0,   20.40605394795271d0/
      data tf7/    -4.05440705938640d0,   -1.26680773001722d0,
     .     16.17184600661102d0,   16.72039538077897d0/

c--- jpt picks the momentum configuration (1..5)
c--- k picks the collinear pair:
c---  1   2   3   4   5   6   7   8   9
c---  17  27  57  67  16  26  56  15  25

c--- handle permutations (initial,final)
      data ii1/1,2,1,2,1,2/
      data ii2/2,1,2,1,2,1/
      data ii3/3,3,3,3,3,3/
      data ii4/4,4,4,4,4,4/
      data ii5/5,5,5,5,7,7/
      data ii6/6,6,7,7,6,6/
      data ii7/7,7,6,6,5,5/
c--- handle permutations (final,final)
      data if1/1,1,1/
      data if2/2,2,2/
      data if3/3,3,3/
      data if4/4,4,4/
      data if5/5,6,5/
      data if6/6,5,7/
      data if7/7,7,6/

c-- translate i to (i,f)/(f,f) and perm
      data iorf/1,1,2,2,1,1,2,1,1/
      data pntr/1,2,1,2,3,4,3,5,6/

      ip=pntr(i)

      if (iorf(i) .eq. 1) then
c--- initial-final case
      write(*,*) 'Small: s(',ii1(ip),',',ii7(ip),')'
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

      p(ii6(ip),1)=pi6(1)
      p(ii6(ip),2)=pi6(2)
      p(ii6(ip),3)=pi6(3)
      p(ii6(ip),4)=pi6(4)

      p(ii7(ip),1)=pi7(1)
      p(ii7(ip),2)=pi7(2)
      p(ii7(ip),3)=pi7(3)
      p(ii7(ip),4)=pi7(4)
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

      p(ii6(ip),1)=qi6(1)
      p(ii6(ip),2)=qi6(2)
      p(ii6(ip),3)=qi6(3)
      p(ii6(ip),4)=qi6(4)

      p(ii7(ip),1)=qi7(1)
      p(ii7(ip),2)=qi7(2)
      p(ii7(ip),3)=qi7(3)
      p(ii7(ip),4)=qi7(4)
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

      p(ii6(ip),1)=ri6(1)
      p(ii6(ip),2)=ri6(2)
      p(ii6(ip),3)=ri6(3)
      p(ii6(ip),4)=ri6(4)

      p(ii7(ip),1)=ri7(1)
      p(ii7(ip),2)=ri7(2)
      p(ii7(ip),3)=ri7(3)
      p(ii7(ip),4)=ri7(4)
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

      p(ii6(ip),1)=si6(1)
      p(ii6(ip),2)=si6(2)
      p(ii6(ip),3)=si6(3)
      p(ii6(ip),4)=si6(4)

      p(ii7(ip),1)=si7(1)
      p(ii7(ip),2)=si7(2)
      p(ii7(ip),3)=si7(3)
      p(ii7(ip),4)=si7(4)
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

      p(ii6(ip),1)=ti6(1)
      p(ii6(ip),2)=ti6(2)
      p(ii6(ip),3)=ti6(3)
      p(ii6(ip),4)=ti6(4)

      p(ii7(ip),1)=ti7(1)
      p(ii7(ip),2)=ti7(2)
      p(ii7(ip),3)=ti7(3)
      p(ii7(ip),4)=ti7(4)
      endif
      else
c--- final-final case
      write(*,*) 'Small: s(',if5(ip),',',if7(ip),')'
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

      p(if6(ip),1)=pf6(1)
      p(if6(ip),2)=pf6(2)
      p(if6(ip),3)=pf6(3)
      p(if6(ip),4)=pf6(4)

      p(if7(ip),1)=pf7(1)
      p(if7(ip),2)=pf7(2)
      p(if7(ip),3)=pf7(3)
      p(if7(ip),4)=pf7(4)
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

      p(if6(ip),1)=qf6(1)
      p(if6(ip),2)=qf6(2)
      p(if6(ip),3)=qf6(3)
      p(if6(ip),4)=qf6(4)

      p(if7(ip),1)=qf7(1)
      p(if7(ip),2)=qf7(2)
      p(if7(ip),3)=qf7(3)
      p(if7(ip),4)=qf7(4)
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

      p(if6(ip),1)=rf6(1)
      p(if6(ip),2)=rf6(2)
      p(if6(ip),3)=rf6(3)
      p(if6(ip),4)=rf6(4)

      p(if7(ip),1)=rf7(1)
      p(if7(ip),2)=rf7(2)
      p(if7(ip),3)=rf7(3)
      p(if7(ip),4)=rf7(4)
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

      p(if6(ip),1)=sf6(1)
      p(if6(ip),2)=sf6(2)
      p(if6(ip),3)=sf6(3)
      p(if6(ip),4)=sf6(4)

      p(if7(ip),1)=sf7(1)
      p(if7(ip),2)=sf7(2)
      p(if7(ip),3)=sf7(3)
      p(if7(ip),4)=sf7(4)
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

      p(if6(ip),1)=tf6(1)
      p(if6(ip),2)=tf6(2)
      p(if6(ip),3)=tf6(3)
      p(if6(ip),4)=tf6(4)

      p(if7(ip),1)=tf7(1)
      p(if7(ip),2)=tf7(2)
      p(if7(ip),3)=tf7(3)
      p(if7(ip),4)=tf7(4)
      endif

      endif      

      return
      end
















