!------------------------------------------------------------------
! ALEXI MODEL PARAMETERS
!------------------------------------------------------------------
      include "ALEXI_parm.inc"
!
!------------------------------------------------------------------
! ALEXI MODEL COMMON BLOCKS
!------------------------------------------------------------------
      common/absorption/aleafn,aleafv,aleafl,adeadv,adeadn,adeadl      
      common/aerodyn/xlog1,xlog2,a,uexp1,uexp2,expuxp,psima,psih     
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump,rsmin
      common/clumping/clumps1,clumps2,clump0,fveg
      common/cover2/perennial,iswater,fcbare
      common/cover3/iswater_inland,ai,aj
      common/constants/pi,cp,xk
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1, &
      & h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1, &
      & sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1, &
      & albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2, &
      & h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2, &
      & sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2, &
      & albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2
      common/data22/w2orig,eref2          
      common/emission/emleaf,emdead,emsoil,emcpy 
      common/emission2/esfc,aem,bem,eleaf,esoil
      common/flags/writeme,badinput,converged,stopiter            
      common/initial/hn0,psi0,fc0
      common/model/zta
      common/obs/hobs1,hobs2,xleobs1,xleobs2,gobs1,gobs2, &
      & rnobs1,rnobs2
      common/partition/difvis,difnir,dirvis,dirnir,fvis,fnir      
      common/partition_clear/difvisclr,difnirclr,dirvisclr, &
      & dirnirclr,fvisclr,fnirclr
      common/pbl/thrise,zpbli(mli),thpbli(mli),zpbl(ml),thpbl(ml), &
      & jzmax,jz1,nlev
      common/pbl2/tabtheta(mt),tabz2(mt),dtheta,itmax,thpblz1     
      common/reflection/rsoiln,rsoilv,rdcpyl
      common/reflection2/albv,albn,albvobs,albnobs,albobs
      common/site/xlat,xlong,stdlng   
      common/stability/zdlamx,zdlamn         
      common/timestamp/year,doy 
      common/transmission/taudn,taudv,taudl 
      common/view/theta,ftheta
      
      logical writeme,badinput,converged,stopiter,perennial
      logical iswater,iswater_inland
