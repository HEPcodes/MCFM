      SUBROUTINE GGTTWW(WT)
* THE CROSS SECTION FOR
*    G(1) G(2) ---> T  TBAR
* FOLLOWED BY THE DECAYS
*    TBAR ---> BBAR(3) E-(4) NUEBAR(5)
*    T    ---> B(6) NUE(7) E+(8)
*
*  NB This notation (which is correct) does not correspond to 
*  Eqs. (3,4) in Z. Phys C40 419, 1988. The latter is therefore wrong.
      IMPLICIT DOUBLE PRECISION (A-H,O-Y),COMPLEX*16 (Z)
      DOUBLE COMPLEX SPL(10,10),SMN(10,10)
*      COMMON/AA/PI,W,Q,QCDL,NF,MODE
      COMMON/PARS/RMT,RGT,RMW,RGW,RMB,RMTLO,RMTUP
      COMMON/MOM/PLAB(4,10)
      COMMON/COUPS/GW,GS
      COMMON/CSTD/SPL,SMN
      DATA INIT/0/
      SAVE RMW2,RMGW,RMB2,RMT2,RMGT,INIT
      IF(INIT.EQ.0) THEN
      INIT=1
       RMW2=RMW**2
       RMGW=RMW*RGW
       RMB2=RMB**2
       RMT2=RMT**2
       RMGT=RMT*RGT
      ENDIF
*
* THE AUXILIARY VECTORS FOR THE T AND TBAR MOMENTA
      D34=DOTKS(3,4)
      D35=DOTKS(3,5)
      D45=DOTKS(4,5)
      D67=DOTKS(6,7)
      D68=DOTKS(6,8)
      D78=DOTKS(7,8)
      C9 =(D34+D35+D45+0.5D0*RMB2)/(D34+D45)
      C10=(D67+D68+D78+0.5D0*RMB2)/(D68+D78)
      DO 1 K=1,4
      PLAB(K,9 )=PLAB(K,3)+(1.D0-C9 )*PLAB(K,4)+PLAB(K,5)
    1 PLAB(K,10)=PLAB(K,6)+PLAB(K,7)+(1.D0-C10)*PLAB(K,8)
*
* THE DENOMIANTORS
      XTBAR=2.D0*(D34+D35+D45)+RMB2
      XT   =2.D0*(D67+D68+D78)+RMB2
      ZD0=DCMPLX(2.D0*D45-RMW2,RMGW)*
     .    DCMPLX(2.D0*D78-RMW2,RMGW)*
     .    DCMPLX(XTBAR-RMT2,RMGT)*
     .    DCMPLX(XT   -RMT2,RMGT)
      D13=DOTKS(1,3)
      D14=DOTKS(1,4)
      D15=DOTKS(1,5)
      D1345=2.D0*(D34+D35+D45-D13-D14-D15)+RMB2
      D16=DOTKS(1,6)
      D17=DOTKS(1,7)
      D18=DOTKS(1,8)
      D1678=2.D0*(D67+D68+D78-D16-D17-D18)+RMB2
      ZDEN1=ZD0*DCMPLX(D1345-RMT2,RMGT)
      ZDEN2=ZD0*DCMPLX(D1678-RMT2,RMGT)
      D12=DOTKS(1,2)
      ZDEN3=ZD0*2.D0*D12
*
* THE SPINOR PARTS
      D1T=D16+D17+D18
      D2T=DOTKS(2,6)+DOTKS(2,7)+DOTKS(2,8)
      CALL STD
*
      Z81=SMN(8,10)*SPL(10,1)
      Z82=SMN(8,10)*SPL(10,2)
      Z14=SMN(1,9)*SPL(9,4)
      Z24=SMN(2,9)*SPL(9,4)
*
      ZG1PP=(-2.D0*(D1T-D12)*Z82*Z24
     .  +RMT2*(  Z82*SMN(1,2)*SPL(1,4)
     .         - SMN(8,1)*SPL(2,1)*Z24
     .         +2.D0*D2T*SMN(8,1)*SPL(1,4) ))/ZDEN1
c--- modified by JMC 1/24/00 to protect SQRT < 0
      SFAC=4.D0*D1T*D2T-2.D0*XT*D12
      IF (SFAC .LE. 0D0) THEN
        SFAC=0D0
      ELSE
        SFAC=DSQRT(SFAC)
      ENDIF
      ZG1PM=SFAC*( -Z81*Z24
     .             +RMT2*SMN(8,2)*SPL(1,4) )/ZDEN1
      ZG1MP=SFAC*( -Z82*Z14
     .             +RMT2*SMN(8,1)*SPL(2,4) )/ZDEN1
      ZG1MM=(-2.D0*D2T*Z81*Z14
     .  +RMT2*(  Z81*SMN(2,1)*SPL(2,4)
     .         - SMN(8,2)*SPL(1,2)*Z14
     .         + 2.D0*(D1T-D12)*SMN(8,2)*SPL(2,4) ))/ZDEN1
*
      ZG2PP=(-2.D0*(D2T-D12)*Z81*Z14
     .  +RMT2*(  Z81*SMN(2,1)*SPL(2,4)
     .         - SMN(8,2)*SPL(1,2)*Z14
     .         + 2.D0*D1T*SMN(8,2)*SPL(2,4) ))/ZDEN2
      ZG2PM=ZG1PM *ZDEN1/ZDEN2
      ZG2MP=ZG1MP *ZDEN1/ZDEN2
      ZG2MM=(-2.D0*D1T*Z82*Z24
     .  +RMT2*(  Z82*SMN(1,2)*SPL(1,4)
     .         - SMN(8,1)*SPL(2,1)*Z24
     .         + 2.D0*(D2T-D12)*SMN(8,1)*SPL(1,4) ))/ZDEN2
*
      ZG3PP=(-Z81*Z14+Z82*Z24
     .  +RMT2*( SMN(8,1)*SPL(1,4) - SMN(8,2)*SPL(2,4) ) )/ZDEN3
      ZG3PM=(0.D0,0.D0)
      ZG3MP=(0.D0,0.D0)
      ZG3MM=ZG3PP
*
      C=1.d0/(2.d0*D12)
      ZEVPP=C*( ZG1PP+ZG2PP)
      ZEVPM=C*( ZG1PM+ZG2PM)
      ZEVMP=C*( ZG1MP+ZG2MP)
      ZEVMM=C*( ZG1MM+ZG2MM)
      ZODPP=C*(-ZG1PP+ZG2PP)-ZG3PP
      ZODPM=C*(-ZG1PM+ZG2PM)-ZG3PM
      ZODMP=C*(-ZG1MP+ZG2MP)-ZG3MP
      ZODMM=C*(-ZG1MM+ZG2MM)-ZG3MM
*
      XEVEN=ABS(ZEVPP)**2+ABS(ZEVPM)**2
     .     +ABS(ZEVMP)**2+ABS(ZEVMM)**2
      XODD =ABS(ZODPP)**2+ABS(ZODPM)**2
     .     +ABS(ZODMP)**2+ABS(ZODMM)**2
* COLLECT ALL FACTORS
      WT = ( 28.d0/3.d0*XEVEN + 12.d0*XODD )
     .   * GW**8 * GS**4 * (128.d0)**2 *D67 *D35
     .             *4            /4.d0             /64.d0
      WT=WT/4.d0
* FACTORS FOR B,BAR SPINS, GG SPIN AVERAGE, GG COLOUR AVERAGE
      RETURN
      END
