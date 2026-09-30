c
c-------------------------------------------------------------------c
c     U S F L U X   I N P U T / O U T P U T   R O U T I N E S       c
c-------------------------------------------------------------------c
c
      subroutine extract_input(m,ibad,MDATE)
c     **************************************************************
c     *
c     *  Extracts input at grid cell IA JA.
c     *
c     *  Martha Anderson
c     *  Created:  10/5/01
c     *
c     **************************************************************
      include 'USflux_grids.inc'
      common/absorption/aleafn,aleafv,aleafl,adeadv,adeadn,adeadl
      common/availh2o/brz,awfrz,awcrz,bsfc,awfsfc,awcsfc
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump,rsmin
      common/clumping/clumps1,clumps2,clump0,fveg
      common/cover2/perennial,iswater,fcbare
      common/cover3/iswater_inland,ai,aj
      common/cover/iclass,beta,gvs,xndvi,fpar,dstom
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1,
     &       h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1,
     &       sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1,
     &       albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2
      common/data22/w2orig,eref2
      common/departure/xmc,xbc,xcc,xms,xbs,xcs
      common/emission/emleaf,emdead,emsoil,emcpy
      common/emission2/esfc,aem,bem,eleaf,esoil
      common/flags/writeme,badinput,converged,stopiter
      common/model/zta
      common/obs/hobs1,hobs2,xleobs1,xleobs2,gobs1,gobs2,
     &       rnobs1,rnobs2
      common/pbl/thrise,zpbli(mli),thpbli(mli),zpbl(ml),thpbl(ml),
     &       jzmax,jz1,nlev
      common/reflection/rsoiln,rsoilv,rdcpyl
      common/reflection2/albv,albn,albvobs,albnobs,albobs
      common/site/xlat,xlong,stdlng
      common/timestamp/year,doy
      common/usflags/clear
      common/constants/pi,cp,xk
      common/view/theta,ftheta
      common/lookup/lookup_i(ilg,jlg),lookup_j(ilg,jlg)
      common/slookup/insol_i(ilg,jlg),insol_j(ilg,jlg)
      common/hrdata_met/ctloc(kx,ky,kt),cta(kx,ky,kt),
     &       cea(kx,ky,kt),cwind(kx,ky,kt),csdn(kx,ky,kt),
     &       cxlwdn(kx,ky,kt),cpres(kx,ky,kt),
     &       clst(kx,ky,kt), clapse(ilg,jlg)

      logical writeme,badinput,converged,stopiter,clear
      logical present,WriteOutput
      logical iswater,perennial,haveSolar
      logical iswater_inland, isnan, investigate
      integer ip, jp
      real itime,dtime,elev_diff,elevation, wind_factor
      real clapse1, temp_lwdn, wt1, wt2, wt3, wt4
      real temp_lwdn2, new_ta
      real r15, r55
      real time(kt)
      integer tindex1, tindex2
      integer*4 m
      integer test
      real vars(14), mvars(10), svars(8), pvars(6), rvars(14)
      save


      read(102,*,iostat=test,end=200) i1,i2,i3,i4,(mvars(k),k=1,10)
      read(110,*,end=200) i1,i2,(svars(k),k=1,8)
      read(120,*,end=200) i1,i2,(pvars(k),k=1,6)
      read(130,*,end=200) i1,i2,(rvars(k),k=1,14)

      rsmin=pvars(6)
!      write(6,*) i1, i2, rsmin
!      write(6,*) "test = ", test, mvars(1), mvars(2)
!      write(6,*) mvars(:)
!      write(6,*) svars(:)
!      write(6,*) pvars(:)

      ai=i1
      aj=i2
      iflag=0
!  Initialize pixel flags
      writeme=.FALSE.
!      if (mod(ia*ja,50).eqv.0) writeme=.TRUE.
      converged=.FALSE.
      WriteOutput=.FALSE.
      badinput=.FALSE.
      investigate=.FALSE.
!                             P L A N T S
!                             -----------
!  Landcover information:
      iclass=5 !pvars(5)
      iswater=.FALSE.
      if (iclass.eq.0 .or. iclass.eq.WATER
     &                .or. iclass.eq.ICE) then
        iflag=3 ! iIFLAG=3 denotes water/ice
        iswater=.TRUE.
      endif
      iswater_inland=.FALSE.

!  LAI and cover fraction
      xlai=svars(8)/10.
      if (xlai.lt.0.0) then	! LSA SAF bad value is -999, and 0's typically indicate a bad solution
        xlai=BAD
        fc=BAD
        height=BAD
        z0=BAD
        disp=BAD
      else
        fc=1.0-exp(-0.5*xlai)
        if (fc.gt.1.0) fc=1.0
        if (fc.le.0.0) fc=0.01
        call cover_props(ia,ja)
        if (height.eq.0.0) then
         height=0.01
         disp=(2./3.)*height
        endif
!        if (z0t.ge.0.0) then
!         z0=z0+(1.0)*z0t*(z0/7+0.3)
!         height=8.*z0
!         disp=(2./3.)*height
!        endif
      endif
      fg=1.0 !svars(9)*1.2
      if (svars(9).eq.-9999.) then
       fg=1.0
      endif
!                               S O I L
!                               -------
!  Soil radiative properties:
      rsoilv=pvars(1)
      rsoiln=pvars(2)
      emsoil=0.94
      albedo2=BAD

!                              T I M E
!                              -------
!  Set two sunrise-relative times to local standard time as a
!  function of standard longitude (i.e. time zone).
      xlat=pvars(3)
      xlong=pvars(4)
      r15=svars(5)
      r55=svars(6)
      call getdgmt(xlong,dgmt,stdlng)
!      write(6,*) dgmt, r15, r55
      tloc1=r15+dgmt		! In local standard time
      tloc2=r55+dgmt
      if (tloc1.gt.24) then
       tloc1=tloc1-24.0
      endif
      if (tloc2.gt.24) then
       tloc2=tloc2-24.0
      endif

!			H O U R L Y  W E A T H E R
!                      ----------------------------
!  Build time series from WRF data.
!      call extract_hourly_input_nldas(ia,ja,dgmt)

!                     S U R F A C E   W E A T H E R
!                     -----------------------------
!  Weather at the 2 geostationary observation times:
       taobs1=mvars(1)-273.15
       taobs2=mvars(2)-273.15
       if (taobs1.gt.taobs2) then
        badinput=.TRUE.
       endif
       w1=(mvars(9)+mvars(10))/2.
       w2=(mvars(9)+mvars(10))/2.
       w1=abs(w1)
       w2=abs(w2)
       w2orig=abs(w2)
       ea1=mvars(3)
       ea2=mvars(4)
       pres1=mvars(5)
       pres2=mvars(6)
       xlwdn1=mvars(7)
       xlwdn2=mvars(8)
       sdn1=svars(3)
       sdn2=svars(4)
       call getxlwdn(ea1,taobs1,1.0,xlwdn1)
       call getxlwdn(ea2,taobs2,1.0,xlwdn2)
      if (r55-r15.le.0.0) then
       iflag=2
      endif

!      write(6,*) "after met"
!  Wind reference heights:
      refhtw = 30.0 		! Height of input wind data from WRF [m]
      zta    = 50.0		! Model blending height [m]


!                S U R F A C E   T E M P E R A T U R E
!                -------------------------------------
      theta=svars(7)
      trad1=svars(1)-273.15
      trad2=svars(2)-273.15
      dthr=(trad2-trad1)/(r55-r15)
      diff=dthr*rsmin-dthr
      offset=diff*(r55-r15)
      trad1=trad1-offset
      if (trad1.gt.trad2) then
       badinput=.TRUE.
      endif

!                     B O U N D A R Y   L A Y E R
!                     ----------------------------
!  PBL information:

      zpbli(1) = 0.0
      zpbli(2) = 100.0
      zpbli(3) = 300.0
      zpbli(4) = 500.0
      zpbli(5) = 700.0
      zpbli(6) = 1000.0
      zpbli(7) = 1400.0
      zpbli(8) = 1800.0
      zpbli(9) = 2200.0
      zpbli(10) = 2600.0
      zpbli(11) = 3000.0
      zpbli(12) = 3500.0
      zpbli(13) = 4000.0
      zpbli(14) = 4500.0

      do k = 1, 14
       thpbli(k) = rvars(k)
      enddo


!       thpbli(1) = 290.0
!       thpbli(2) = 293.0
!       thpbli(3) = 294.5
!       thpbli(4) = 295.5
!       thpbli(5) = 296.5
!       thpbli(6) = 297.5
!       thpbli(7) = 298.5
!       thpbli(8) = 299.5
!       thpbli(9) = 300.5


!                 SOIL AND LEAF RADIATIVE PROPERTIES
!                 ----------------------------------
!  Temporary surface emissivity calculation (for output)
      eleaf=0.97                ! Leaf emissivity
      b1=-1.97*cos(theta)+2.87
      b2=1.86*cos(theta)-2.62
      bem=b1*eleaf+b2
      aem=eleaf+0.025-bem-emsoil
      esfc=aem*fc*fc+bem*fc+emsoil

      if (investigate.eqv..TRUE.) then
       write(6,*) "-------------------------"
       write(6,*) "INPUTS......", xlat, xlong
       write(6,*) "LAI =       ", xlai, iflag, badinput
       write(6,*) "FC =        ", fc, height, disp
       write(6,*) "CLASS =     ", iclass
       write(6,*) "TLOC1 =     ", tloc1
       write(6,*) "TLOC2 =     ", tloc2
       write(6,*) "TRAD1 =     ", trad1
       write(6,*) "TRAD2 =     ", trad2
       write(6,*) "SDN1 =      ", sdn1
       write(6,*) "SDN2 =      ", sdn2
       write(6,*) "LWDN1 =     ", xlwdn1
       write(6,*) "LWDN2 =     ", xlwdn2
       write(6,*) "W1 =        ", w1
       write(6,*) "W2 =        ", w2
       write(6,*) "T1 =        ", taobs1
       write(6,*) "T2 =        ", taobs2
       write(6,*) "P1 =        ", pres1
       write(6,*) "P2 =        ", pres2
       write(6,*) "E1 =        ", ea1
       write(6,*) "E2 =        ", ea2
       write(6,*) "ANG =       ", theta
       write(6,*) "RSLV =      ", rsoilv
       write(6,*) "RSLN =      ", rsoiln
       write(6,*) "ALFV =      ", aleafv
       write(6,*) "ALFN =      ", aleafn
       write(6,*) "THETa =     ", theta
      endif

!                              F L A G S
!                              ---------
!  For the iflag array:
!    "0" denotes good input data
!    "1" is bad input data other than cloudy conditions
!    "2" denotes cloudy conditions
!    "3" indicates pixels classified as water or ice
!    "4" denotes that ALEXI was attempted but failed to converge

      if (iflag.eq.2) then
        iclear=0                   ! cloudy
        clear=.FALSE.
      else
        iclear=1                   ! clear
        clear=.TRUE.
      endif

!  Check variables for unreasonable values.
!  (Variables used only in ALEXI are checked in ALEXI.)
!      if (.not.iswater.and..not.badinput)
!     &         call checkUSinput(badinput,ibad,writeme)
      if (badinput) iflag=1

 200  if(test.lt.0) badinput=.TRUE.
      return
      end

      subroutine getdgmt(xnlon,dgmt,stdlng)
c     **************************************************************
c     *
c     *  Find standard longitude and offset of local time from
c     *  GMT.  XNLON is longitude in 0-360 deg.
c     *
c     *  Martha Anderson
c     *  Created:  04/08/02
c     *
c     ***************************************************************


      if (xnlon.ge.-172.5 .and. xnlon.lt.-157.5) then
        stdlng=-165.0
        dgmt=-11.0
      elseif (xnlon.ge.-157.5 .and. xnlon.lt.-142.5) then
        stdlng=-150.0
        dgmt=-10.0
      elseif (xnlon.ge.-142.5 .and. xnlon.lt.-127.5) then
        stdlng=-135.0
        dgmt=-9.0
      elseif (xnlon.ge.-127.5 .and. xnlon.lt.-112.5) then
        stdlng=-120.0
        dgmt=-8.0
      elseif (xnlon.ge.-112.5 .and. xnlon.lt.-97.5) then
        stdlng=-105.0
        dgmt=-7.0
      elseif (xnlon.ge.-97.5 .and. xnlon.lt.-82.5) then
        stdlng=-90.0
        dgmt=-6.0
      elseif (xnlon.ge.-82.5 .and. xnlon.lt.-67.5) then
        stdlng=-75.0
        dgmt=-5.0
      elseif (xnlon.ge.-67.5 .and. xnlon.lt.-52.5) then
        stdlng=-60.0
        dgmt=-4.0
      elseif (xnlon.ge.-52.5 .and. xnlon.lt.-37.5) then
        stdlng=-45.0
        dgmt=-3.0
      elseif (xnlon.ge.-37.5 .and. xnlon.lt.-22.5) then
        stdlng=-30.0
        dgmt=-2.0
      elseif (xnlon.ge.-22.5 .and. xnlon.lt.-7.5) then
        stdlng=-15.0
        dgmt=-1.0
      elseif (xnlon.ge.-7.5 .and. xnlon.lt.7.5) then
        stdlng=0.0
        dgmt=0.0
      elseif (xnlon.ge.7.5 .and. xnlon.lt.22.5) then
        stdlng=15.0
        dgmt=1.0
      elseif (xnlon.ge.22.5 .and. xnlon.lt.37.5) then
        stdlng=30.0
        dgmt=2.0
      elseif (xnlon.ge.37.5 .and. xnlon.lt.52.5) then
        stdlng=45.0
        dgmt=3.0
      elseif (xnlon.ge.52.5 .and. xnlon.lt.67.5) then
        stdlng=60.0
        dgmt=4.0
      elseif (xnlon.ge.67.5 .and. xnlon.lt.82.5) then
        stdlng=75.0
        dgmt=5.0
      elseif (xnlon.ge.82.5 .and. xnlon.lt.97.5) then
        stdlng=90.0
        dgmt=6.0
      elseif (xnlon.ge.97.5 .and. xnlon.lt.112.5) then
        stdlng=105.0
        dgmt=7.0
      elseif (xnlon.ge.112.5 .and. xnlon.lt.127.5) then
        stdlng=120.0
        dgmt=8.0
      elseif (xnlon.ge.127.5 .and. xnlon.lt.142.5) then
        stdlng=135.0
        dgmt=9.0
      elseif (xnlon.ge.142.5 .and. xnlon.lt.157.5) then
        stdlng=150.0
        dgmt=10.0
      elseif (xnlon.ge.157.5 .and. xnlon.lt.172.5) then
        stdlng=165.0
        dgmt=11.0
      endif

      return
      end



c-------------------------------------------------------------------c
c                 O U T P U T   R O U T I N E S
c-------------------------------------------------------------------c


c
      subroutine screen_output(ia,ja)
c     **************************************************************
c     *
c     *  Write pixel diagnostics to screen.
c     *
c     *  Martha Anderson
c     *  Created:  10/5/01
c     *
c     ***************************************************************
      include 'USflux_grids.inc'
      include 'date.inc'
      common/availh2o/brz,awfrz,awcrz,bsfc,awfsfc,awcsfc
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/cover/iclass,beta,gvs,xndvi,fpar,dstom
      common/cover3/iswater_inland,ai,aj
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1,
     &       h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1,
     &       sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1,
     &       albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2
      common/dayflux/aparday,acday,rnday,rnsday,rncday,gday,
     &       eday,esday,ecday,hday,hsday,hcday,sday
      common/dayflux2/erefday
      common/daypotflux/epotday,espotday,ecpotday
      common/flags/writeme,badinput,converged,stopiter
      common/lookup/itabclass(nclass),tablai(nbin,nclasm),
     &       tabfpar(nbin,nclasm),tabbeta(nclass),
     &       tabaleaf(nclass,3),tabadead(nclass,3),
     &       tabhmin(nclass),tabhmax(nclass),tabxl(nclass),
     &       tabgvs(nclass),tabdstom(nclass),tabperen(nclass),
     &       tabfcmin(nclass),tabdesc(nclass)
      common/initial/hn0,psi0,fc0
      common/site/xlat,xlong,stdlng
      common/stress/fpet,fsdn
      common/timestamp/year,doy
      common/usflags/clear
      logical writeme,badinput,converged,stopiter,clear
      character*50 tabdesc
      save
c  		Write out results to screen
c
      if (writeme.and.clear.and.converged) then
         write(6,1002) xlat,xlong
         write(6,*) iclass
         write(99,*)   ' LANDCOVER CLASS           '
         write(6,1003)' FCORIG (frac cover in)    ',fc0
         write(6,1003)' FC     (frac cover out)   ',fc
         write(6,1003)' XLAI   (LAI out)          ',xlai
         write(6,1003)' HEIGHT (canopy height)    ',height
         write(99,1003)' XL     (leaf size)	  ',xl
         write(6,1003)' Z0     (roughness)        ',z0
         write(6,1003)' TLOC1  (GOES obs1)        ',tloc1
         write(6,1003)' TLOC2  (GOES obs2)        ',tloc2
         write(6,1003)' TGMT1  (GOES obs1)        ',r15
         write(6,1003)' TGMT2  (GOES obs2)        ',r55
         write(6,1003)' W1     (wind speed1)      ',w1
         write(6,1003)' W2     (wind speed2)      ',w2
         write(99,1003)' TC1    (canopy temp1)     ',tc1
         write(99,1003)' TS1    (soil temp1)       ',ts1
         write(99,1003)' TC2    (canopy temp2)     ',tc2
         write(99,1003)' TS2    (soil temp2)       ',ts2
         write(6,1003)' XLEC2  (canopy latent)    ',xlec2
         write(6,1003)' XLES2  (soil latent)      ',xles2
         write(6,1003)' XLE2   (total latent)     ',xle2
         write(99,1003)' HC2    (canopy sensible)  ',hc2
         write(99,1003)' HS2    (soil sensible)    ',hs2
         write(6,1003)' H2     (total sensible)   ',h2
         write(6,1003)' G2     (ground cond)      ',g2
         write(99,1003)' ZEN1   (sun zenith T1)    ',zen1
         write(99,1003)' ZEN2   (sun zenith T2)    ',zen2
         write(99,1003)' RNET1  (net rad T1)       ',rnet1
         write(6,1003)' RNET2  (net rad T2)       ',rnet2
         write(99,1003)' TA1    (modelled air temp)',ta1
         write(99,1003)' TA2    (modelled air temp)',ta2
         write(99,1003)' TAOBS1 (observed air temp)',taobs1
         write(99,1003)' TAOBS2 (observed air temp)',taobs2
         write(6,1003)' TRAD1  (radiometric T1)   ',trad1
         write(6,1003)' TRAD2  (radiometric T2)   ',trad2
         write(6,1003)' Z2     (boundary layer ht)',z2
         write(99,1003)' TH1    (potential temp T1)',th1-273.15
         write(99,1003)' TH2    (potential temp T2)',th2-273.15
         write(6,1003)' RA2    (aerodynamic res)  ',ra2
         write(6,1003)' RS2    (soil res)         ',rs2
         write(6,1003)' RX2    (b.l. res)         ',rx2
         write(6,1003)' SWUP                      ',swup2
         write(6,1003)' SDN2                      ',sdn2
         write(6,1003)' LWUP                      ',xlwup2
         write(6,1003)' LWDN                      ',xlwdn2
!      elseif (writeme.and..not.clear) then
         write(99,1003)' RNDAY         ',rnday
         write(99,1003)' HDAY          ',hday
         write(99,1003)' GDAY          ',gday
         write(99,1003)' EDAY          ',eday
         write(99,1003)' ESDAY         ',esday
         write(99,1003)' ECDAY         ',ecday
         write(99,1003)' EREFDAY       ',erefday
         write(99,1003)' SDAY          ',sday
      endif

      write(103,1004) int(ai), int(aj), rnet2, xle2, h2, g2, xles2
!       write(103,1004) int(ai),int(aj),rnet2,xle2,h2,g2,xles2,xlecs,hs2,
!     &                hc2,z2,ra2,rs2,rx2,sdn2,swup2,xlwdn2,xlwup2

1002  format('Lat: ',f6.2,'  Lon: ',f7.2)
1003  format(a27,f12.3)
1004  format(2I6,5f10.1)
!1004  format(2I6,16f10.1)
      return
      end


      subroutine store_output_bin(i,j,ierr,ibad,iflag,iter)
c     **************************************************************
c     *
c     *  Store gridded output in binary files.
c     *
c     *  Martha Anderson
c     *  Created:  10/5/01
c     *
c     ***************************************************************
      include 'USflux_grids.inc'

      common/absorption/aleafn,aleafv,aleafl,adeadv,adeadn,adeadl
      common/availh2o/brz,awfrz,awcrz,bsfc,awfsfc,awcsfc
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/cover2/perennial,iswater,fcbare
      common/cover/iclass,beta,gvs,xndvi,fpar,dstom
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1,
     &       h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1,
     &       sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1,
     &       albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2
      common/dayflux/aparday,acday,rnday,rnsday,rncday,gday,
     &       eday,esday,ecday,hday,hsday,hcday,sday
      common/dayflux2/erefday
      common/dayflux3/xlwupday,xlwdnday,swupday
      common/daypotflux/epotday,espotday,ecpotday
      common/departure/xmc,xbc,xcc,xms,xbs,xcs
      common/emission2/esfc,aem,bem,eleaf,esoil
      common/flags/writeme,badinput,converged,stopiter
      common/hrdata/tloc(nohr),ta(nohr),ea(nohr),wind(nohr),
     &       sdn(nohr),xlwdn(nohr),rnet(nohr),rnsoil(nohr),
     &       g(nohr),pres(nohr),nohrin
      common/hrdata2/epot(nohr),ecpot(nohr),espot(nohr)
      common/reflection/rsoiln,rsoilv,rdcpyl
      common/reflection2/albv,albn,albvobs,albnobs,albobs
      common/site/xlat,xlong,stdlng
      common/stability/zdlamx,zdlamn
      common/stress/fpet,fsdn
      common/timestamp/year,doy
      common/usflags/clear
      common/view/theta,ftheta
      common/clumping/clumps1,clumps2,clump0,fveg
      common/partition/difvis,difnir,dirvis,dirnir,fvis,fnir

      logical writeme,badinput,converged,stopiter,clear
      logical perennial,iswater
      save

c  		Set values to bad: Non-convergence, bad inputs, cloudy

      xiter=iter
      if (.not.converged) then	! Didn't converge - set all instantaneous fluxes to bad & EDAY, HDAY
         fsdn=BAD
         h2=BAD
         xle2=BAD
         g2=BAD
         rnet2=BAD
         ts2=BAD
         tc2=BAD
         t2=BAD
         ta1=BAD
         ta2=BAD
         fpet=BAD
         z2=BAD
         xlec2=BAD
         xles2=BAD
         xiter=BAD
         hc2=BAD
         hs2=BAD
         rs2=BAD
         rx2=BAD
         ra2=BAD
         xlwup2=BAD
         swup2=BAD
         eday=BAD
         ecday=BAD
         esday=BAD
         hday=BAD
         hcday=BAD
         hsday=BAD
      endif

      if(clear.and..not.badinput.and..not.iswater)then	! Cloudy or bad inputs
         dtrad=trad2-trad1
      else
         dtrad=BAD
      endif

      tloc1o=tloc1
      tloc2o=tloc2

      if(badinput.or.iswater)then
        sday=BAD
        albedo2=BAD
        erefday=BAD
        xlwupday=BAD
        xlwdnday=BAD
        swupday=BAD
        rnday=BAD
        gday=BAD
        hday=BAD
        eday=BAD
        ecday=BAD
        esday=BAD
        xndvi=BAD
      endif

      if(iswater)then
        sdn1=BAD
        sdn2=BAD
        ierr=BAD
        trad2=BAD
        xlai=BAD
        height=BAD
        z0=BAD
        fc=BAD
        theta=BAD
        taobs1=BAD
        taobs2=BAD
        tloc1o=BAD
        tloc2o=BAD
        xlat=BAD
        xlong=BAD
        ibad=BAD
        iflag=3
        esfc=BAD
      endif

      if(badinput.or..not.clear.or.iswater)then
        w1=BAD		! These didn't get scaled ...
        w2=BAD
      endif

c  Some computations
      t1gmt=rise15(i,j)
      t2gmt=rise55(i,j)
      class=iclass
      cflag=iflag
      cierr=ierr
      cibad=ibad

c  Compute albedo2 assuming clear-sky partitioning
      if (albvobs.ne.BAD.and.albnobs.ne.BAD.and.
     &rsoilv.ne.BAD.and.rsoiln.ne.BAD.and.
     &aleafv.ne.BAD.and.aleafn.ne.BAD) then
        fvis=0.5
        fnir=0.5
        dirvis=0.8
        difvis=0.2
        dirnir=1.0
        difnir=0.0
        call getsunzen (xlat,xlong,stdlng,doy,year,tloc2,zen2)
        call getradprops (zen2,albedo2,taubtv,taubtn,clumps2)
      else
        albedo2=BAD
        albv=BAD
        albn=BAD
      endif

c  			Write data

c        call binwrite(i,j,ilg,jlg,130,cflag)
        call binwrite(i,j,ilg,jlg,136,fc)
        call binwrite(i,j,ilg,jlg,137,xlai)
c        call binwrite(i,j,ilg,jlg,139,height)
c        call binwrite(i,j,ilg,jlg,140,z0)

c        call binwrite(i,j,ilg,jlg,146,taobs1)
c        call binwrite(i,j,ilg,jlg,147,taobs2)
c        call binwrite(i,j,ilg,jlg,148,w1)
c        call binwrite(i,j,ilg,jlg,149,w2)
c        call binwrite(i,j,ilg,jlg,150,sdn1)
c        call binwrite(i,j,ilg,jlg,151,sdn2)
c        call binwrite(i,j,ilg,jlg,152,ea2)
c        call binwrite(i,j,ilg,jlg,153,pres2)
        call binwrite(i,j,ilg,jlg,154,xlwdn2)
c        call binwrite(i,j,ilg,jlg,155,zen2)

        call binwrite(i,j,ilg,jlg,156,trad1)
        call binwrite(i,j,ilg,jlg,157,trad2)
        call binwrite(i,j,ilg,jlg,158,dtrad)
c        call binwrite(i,j,ilg,jlg,159,theta)

        call binwrite(i,j,ilg,jlg,160,h2)
        call binwrite(i,j,ilg,jlg,161,xle2)
        call binwrite(i,j,ilg,jlg,162,g2)
        call binwrite(i,j,ilg,jlg,163,rnet2)
c        call binwrite(i,j,ilg,jlg,164,tc2)
c        call binwrite(i,j,ilg,jlg,165,ts2)
c        call binwrite(i,j,ilg,jlg,166,ta1)
c        call binwrite(i,j,ilg,jlg,167,ta2)
        call binwrite(i,j,ilg,jlg,168,z2)
        call binwrite(i,j,ilg,jlg,169,xlec2)
        call binwrite(i,j,ilg,jlg,170,xles2)
c        call binwrite(i,j,ilg,jlg,171,hc2)
c        call binwrite(i,j,ilg,jlg,172,hs2)
        call binwrite(i,j,ilg,jlg,173,xlwup2)
        call binwrite(i,j,ilg,jlg,174,swup2)

        call binwrite(i,j,ilg,jlg,175,ra2)
        call binwrite(i,j,ilg,jlg,176,rs2)
        call binwrite(i,j,ilg,jlg,177,rx2)

c        call binwrite(i,j,ilg,jlg,178,fpet)
c        call binwrite(i,j,ilg,jlg,207,fsdn)

c        call binwrite(i,j,ilg,jlg,179,rnday)
c        call binwrite(i,j,ilg,jlg,180,gday)
c        call binwrite(i,j,ilg,jlg,181,hday)
c        call binwrite(i,j,ilg,jlg,182,eday)
c        call binwrite(i,j,ilg,jlg,183,ecday)
c        call binwrite(i,j,ilg,jlg,184,esday)
c        call binwrite(i,j,ilg,jlg,185,sday)
c        call binwrite(i,j,ilg,jlg,186,swupday)
c        call binwrite(i,j,ilg,jlg,187,xlwupday)
c        call binwrite(i,j,ilg,jlg,188,xlwdnday)
c        call binwrite(i,j,ilg,jlg,189,erefday)

        call binwrite(i,j,ilg,jlg,191,rsoilv)
        call binwrite(i,j,ilg,jlg,192,rsoiln)
c        call binwrite(i,j,ilg,jlg,193,xndvi)

!        call binwrite(i,j,ilg,jlg,201,xmc)
!        call binwrite(i,j,ilg,jlg,202,xbc)
!        call binwrite(i,j,ilg,jlg,203,xcc)
!        call binwrite(i,j,ilg,jlg,204,xms)
!        call binwrite(i,j,ilg,jlg,205,xbs)
!        call binwrite(i,j,ilg,jlg,206,xcs)

      return
      end


      subroutine open_output(year,doy)
c     **************************************************************
c     *
c     *  Martha Anderson
c     *  Created:  10/5/01
c     *
c     ***************************************************************
      include 'USflux_dir.inc'
      character*256 infile
      character*256 odir, laidir, albdir

      odir='/data/data123/chain/4KM/GBIM/fluxes/'
      laidir='/data/data123/chain/4KM/GBIM/inputs/'
      albdir='/data/chain/CONUS/ALB/'

c     OPEN DIRECT ACCESS BINARY MATRIX FILES
c      call binopen("IFLG",odir,year,doy,130)

      call binopen("FCOV",odir,year,doy,136)
      call binopen("XLAI",odir,year,doy,137)
c      call binopen("VGHT",odir,year,doy,139)
c      call binopen("RUFF",odir,year,doy,140)

c      call binopen("TAI1",odir,year,doy,146)
c      call binopen("TAI2",odir,year,doy,147)
c      call binopen("WND1",odir,year,doy,148)
c      call binopen("WND2",odir,year,doy,149)
c      call binopen("SDN1",odir,year,doy,150)
c      call binopen("SDN2",odir,year,doy,151)
c      call binopen("EAI2",odir,year,doy,152)
c      call binopen("PRS2",odir,year,doy,153)
      call binopen("LWD2",odir,year,doy,154)
c      call binopen("ZEN2",odir,year,doy,155)

      call binopen("TRD1",odir,year,doy,156)
      call binopen("TRD2",odir,year,doy,157)
      call binopen("DTRD",odir,year,doy,158)
c      call binopen("THET",odir,year,doy,159)

      call binopen("SENS",odir,year,doy,160)
      call binopen("LATN",odir,year,doy,161)
      call binopen("GSOL",odir,year,doy,162)
      call binopen("RNET",odir,year,doy,163)
c      call binopen("TCAN",odir,year,doy,164)
c      call binopen("TSOL",odir,year,doy,165)
c      call binopen("TAO1",odir,year,doy,166)
c      call binopen("TAO2",odir,year,doy,167)
      call binopen("ZPBL",odir,year,doy,168)
      call binopen("XLEC",odir,year,doy,169)
      call binopen("XLES",odir,year,doy,170)
c      call binopen("HCAN",odir,year,doy,171)
c      call binopen("HSOL",odir,year,doy,172)
      call binopen("LWU2",odir,year,doy,173)
      call binopen("SWU2",odir,year,doy,174)

      call binopen("RA2_",odir,year,doy,175)
      call binopen("RS2_",odir,year,doy,176)
      call binopen("RX2_",odir,year,doy,177)

c      call binopen("FPET",odir,year,doy,178)
c      call binopen("FSDN",odir,year,doy,207)

c      call binopen("RDAY",odir,year,doy,179)
c      call binopen("GDAY",odir,year,doy,180)
c      call binopen("HDAY",odir,year,doy,181)
c      call binopen("EDAY",odir,year,doy,182)
c      call binopen("ECDY",odir,year,doy,183)
c      call binopen("ESDY",odir,year,doy,184)
c      call binopen("SDAY",odir,year,doy,185)
c      call binopen("SWUD",odir,year,doy,186)
c      call binopen("LWDD",odir,year,doy,187)
c      call binopen("LWUD",odir,year,doy,188)
c      call binopen("ERDY",odir,year,doy,189)

      call binopen("RSLV",odir,year,doy,191)
      call binopen("RSLN",odir,year,doy,192)
c      call binopen("FMAX",odir,year,doy,193)

c      call binopen("DXMC",year,doy,201)
c      call binopen("DXBC",year,doy,202)
c      call binopen("DXCC",year,doy,203)
c      call binopen("DXMS",year,doy,204)
c      call binopen("DXBS",year,doy,205)
c      call binopen("DXCS",year,doy,206)

      call binopen("MLAI",laidir,year,doy,301)
      call binopen2("ELEV",2012.,0,302)
      call binopen2("FMAX",2012.,0,303)
      call binopen2("CRUF",2012.,0,304)
      call binopen2("RSVA",2012.,0,305)
      call binopen2("RSNA",2012.,0,306)
      call binopen2("ALVA",2012.,0,307)
      call binopen2("ALNA",2012.,0,308)
      call binopen2("DIFF",2012.,0,309)

      return
      end

      subroutine binopen(cvar,dir,year,doy,iunit)
c     **************************************************************
c     *
c     *  Martha Anderson
c     *  Created:  10/5/01
c     *
c     ***************************************************************
      include 'USflux_dir.inc'
      include 'USflux_parm.inc'
      real*4 var(ilg,jlg)
      character*256 outfile,dir
      character*4   cvar
      character*7   cyyyyddd

c  Open binary file for Transform
      iyyyyddd=year*1000+doy
      write(cyyyyddd,'(i7)')iyyyyddd
      outfile=trim(dir)//trim(cvar)//cyyyyddd//'.dat'
      open (unit=iunit, file=outfile,
     &      form='unformatted',recl=4,access='direct')

      return
      end

      subroutine binopen2(cvar,year,doy,iunit)
c     **************************************************************
c     *
c     *  Martha Anderson
c     *  Created:  10/5/01
c     *
c     ***************************************************************
      include 'USflux_dir.inc'
      include 'USflux_parm.inc'
      real*4 var(ilg,jlg)
      character*256 outfile
      character*4   cvar
      character*7   cyyyyddd
c
c  Stored in TABDIR
      l1=index(TABDIR,' ')-1

c  Open binary file for Transform
      iyyyyddd=year*1000+doy
      write(cyyyyddd,'(i7)')iyyyyddd
      outfile=TABDIR(1:l1)//cvar//cyyyyddd//'.dat'
      l2=index(outfile,' ')-1
      open (unit=iunit, file=outfile,
     &      form='unformatted',recl=4,access='direct')

      return
      end


      subroutine binwrite(i,j,ilg,jlg,iunit,val)
c     **************************************************************
c     *
c     *  Martha Anderson
c     *  Created:  10/5/01
c     *
c     ***************************************************************
      real*4  rval

c  	Write one value to binary file
      rval=val
      jj=jlg-j+1
      irec=(jj-1)*ilg+i
      write (iunit,rec=irec) rval

      return
      end

      subroutine binread(i,j,ilg,jlg,iunit,val)
c     **************************************************************
c     *
c     *  Martha Anderson
c     *  Created:  9/11/04
c     *
c     ***************************************************************
      real*4  rval

c  	Read one value from binary file
      jj=jlg-j+1
      irec=(j-1)*ilg+i
      read (iunit,rec=irec) rval
      val=rval

      return
      end

           subroutine binread_swap(i,j,ilg,jlg,iunit,val)
c     **************************************************************
c     *
c     *  Martha Anderson
c     *  Created:  9/11/04
c     *
c     ***************************************************************
      real*4  rval

c  Read from binary file
      jj=jlg-j+1
      irec=(jj-1)*ilg+i
      read (iunit,rec=irec) rval
      call byteswapr4(rval)
      val=rval

      return
      end


c-------------------------------------------------------------------c
c            D E B U G G I N G    U T I L I T I E S
c-------------------------------------------------------------------c

      subroutine writepoint(ia,ja)
c     ***************************************************************
c     *
c     *  Writes out input parameters for a single point at IA,JA.
c     *  Useful for debugging.
c     *
c     ***************************************************************
      include 'USflux_grids.inc'
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1,
     &       h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1,
     &       sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1,
     &       albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2
      common/hrdata/tloc(nohr),ta(nohr),ea(nohr),wind(nohr),
     &       sdn(nohr),xlwdn(nohr),rnet(nohr),rnsoil(nohr),
     &       g(nohr),pres(nohr),nohrin
      common/hrobs/rnobs(nohr),hobs(nohr),xleobs(nohr),gobs(nohr),
     &       aobs(nohr),precp(nohr)
      common/obs/hobs1,hobs2,xleobs1,xleobs2,gobs1,gobs2,
     &       rnobs1,rnobs2
      common/pbl/thrise,zpbli(mli),thpbli(mli),zpbl(ml),thpbl(ml),
     &       jzmax,jz1,nlev
      common/reflection/rsoiln,rsoilv,rdcpyl
      common/site/xlat,xlong,stdlng
      common/timestamp/year,doy
      common/usflags/clear
      common/view/theta,ftheta

      logical clear
      dimension thpblk(mli)
c
      open(70,file="Alexi_Point.out")
      write(70,'(''Input fields for ALEXI for point (ia,ja) '',i3 ,i3)')
     & ia,ja
      write(70,*)'-----------------------------------------------------'
      write(70,*)'       navlat(ja)=',navlat(ja)
      write(70,*)'       navlon(ia)=',navlon(ia)
      write(70,*)'       radini(ia,ja)=',radini(ia,ja)
      write(70,*)'       radfin(ia,ja)=',radfin(ia,ja)
      write(70,*)'       swbeg(ia,ja)=',swbeg(ia,ja)
      write(70,*)'       swend(ia,ja)=',swend(ia,ja)
      write(70,*)'       xlwbeg(ia,ja)=',xlwbeg(ia,ja)
      write(70,*)'       xlwend(ia,ja)=',xlwend(ia,ja)
      write(70,*)'       satang(ia,ja)=',satang(ia,ja)
      write(70,*)'       psr15(ia,ja)=',psr15(ia,ja)
      write(70,*)'       psr55(ia,ja)=',psr55(ia,ja)
      write(70,*)'       wsr15(ia,ja)=',wsr15(ia,ja)
      write(70,*)'       wsr55(ia,ja)=',wsr55(ia,ja)
      write(70,*)'       tsr15(ia,ja)=',tsr15(ia,ja)
      write(70,*)'       tsr55(ia,ja)=',tsr55(ia,ja)
      write(70,*)'       vsr15(ia,ja)=',vsr15(ia,ja)
      write(70,*)'       vsr55(ia,ja)=',vsr55(ia,ja)
      write(70,*)'       pran(ia,ja)=',pran(ia,ja)
      write(70,*)'       htht(1,ia,ja)=',htht(1,ia,ja)
      write(70,*)'       htht(2,ia,ja)=',htht(2,ia,ja)
      write(70,*)'       htht(3,ia,ja)=',htht(3,ia,ja)
      write(70,*)'       htht(4,ia,ja)=',htht(4,ia,ja)
      write(70,*)'       htht(5,ia,ja)=',htht(5,ia,ja)
      write(70,*)'       htht(6,ia,ja)=',htht(6,ia,ja)
      write(70,*)'       htht(7,ia,ja)=',htht(7,ia,ja)
      write(70,*)'       htht(8,ia,ja)=',htht(8,ia,ja)
      write(70,*)'       htht(9,ia,ja)=',htht(9,ia,ja)
      write(70,*)'       htht(10,ia,ja)=',htht(10,ia,ja)
      write(70,*)'       htht(11,ia,ja)=',htht(11,ia,ja)
      write(70,*)'       htht(12,ia,ja)=',htht(12,ia,ja)
      write(70,*)'       htht(13,ia,ja)=',htht(13,ia,ja)
      write(70,*)'       htht(14,ia,ja)=',htht(14,ia,ja)
      write(70,*)'       htht(15,ia,ja)=',htht(15,ia,ja)
      write(70,*)'       htht(16,ia,ja)=',htht(16,ia,ja)
      write(70,*)'       htht(17,ia,ja)=',htht(17,ia,ja)
      write(70,*)'       htht(18,ia,ja)=',htht(18,ia,ja)
      write(70,*)'       htht(19,ia,ja)=',htht(19,ia,ja)
      write(70,*)'       htht(20,ia,ja)=',htht(20,ia,ja)
      write(70,*)'       htht(21,ia,ja)=',htht(21,ia,ja)
      write(70,*)'       htht(22,ia,ja)=',htht(22,ia,ja)
      write(70,*)'       htht(23,ia,ja)=',htht(23,ia,ja)
      write(70,*)'       htht(24,ia,ja)=',htht(24,ia,ja)
      write(70,*)'       htht(25,ia,ja)=',htht(25,ia,ja)
      write(70,*)'       htht(26,ia,ja)=',htht(26,ia,ja)
      write(70,*)'       htht(27,ia,ja)=',htht(27,ia,ja)
      write(70,*)'       htht(28,ia,ja)=',htht(28,ia,ja)
      write(70,*)'       htht(29,ia,ja)=',htht(29,ia,ja)
      write(70,*)'       htht(30,ia,ja)=',htht(30,ia,ja)
      write(70,*)'       htht(31,ia,ja)=',htht(31,ia,ja)
      write(70,*)'       htht(32,ia,ja)=',htht(32,ia,ja)
      write(70,*)'       htht(33,ia,ja)=',htht(33,ia,ja)
      write(70,*)'       htht(34,ia,ja)=',htht(34,ia,ja)
      write(70,*)'       htht(35,ia,ja)=',htht(35,ia,ja)
      write(70,*)'       htht(36,ia,ja)=',htht(36,ia,ja)
      write(70,*)'       htht(37,ia,ja)=',htht(37,ia,ja)
      write(70,*)'       htht(38,ia,ja)=',htht(38,ia,ja)
      write(70,*)'       htht(39,ia,ja)=',htht(39,ia,ja)
      write(70,*)'       htht(40,ia,ja)=',htht(40,ia,ja)
      write(70,*)'       htht(41,ia,ja)=',htht(41,ia,ja)
      write(70,*)'       lscls(ia,ja)=',lscls(ia,ja)
c      write(70,*)'       frcvr(ia,ja)=',frcvr(ia,ja)
      close(70)

c Write input file for single point run
      do jz=1,nlev
        thpblk=thpbli+273.15
      enddo
      iclear=0
      if (clear) then
        iclear=1
      endif
      open(70,file="us.input")
      write (70,*) xlat,xlong,stdlng
      write (70,*) refhtw
      write (70,*) xl,clump,fg
      write (70,*) xndvi,iclass
      write (70,*) xlai,height
      write (70,*) rsoilv,rsoiln,emsoil
      write (70,*) nlev
      write (70,*) (zpbli(jz) ,jz=1,nlev)
      write (70,*) (thpblk(jz),jz=1,nlev)
      write (70,*) year,doy,theta,iclear,nohrin
      write (70,*) tloc1,taobs1,ea1,w1,pres1,trad1,sdn1,xlwdn1,rnobs1,
     &            hobs1,xleobs1,gobs1
      write (70,*) tloc2,taobs2,ea2,w2,pres2,trad2,sdn2,xlwdn2,rnobs2,
     &            hobs2,xleobs2,gobs2
      do ihr=1,nohrin
         write (70,*) tloc(ihr),wind(ihr),ta(ihr),ea(ihr),sdn(ihr),
     &               pres(ihr),rnobs(ihr),hobs(ihr),xleobs(ihr),
     &               gobs(ihr),aobs(ihr),precp(ihr)
      enddo
      close(70)
c
      return
      end

      subroutine extract_hourly_input_nldas(ia,ja,dgmt)
c     **************************************************************
c     *
c     *  Extracts hourly weather/solar data at the current grid
c     *  point.
c     *
c     *  Martha Anderson
c     *  Created:  10/5/01
c     *
c     **************************************************************
      include 'USflux_grids.inc'
      common/hrdata/tloc(nohr),ta(nohr),ea(nohr),wind(nohr),
     &       sdn(nohr),xlwdn(nohr),rnet(nohr),rnsoil(nohr),
     &       g(nohr),pres(nohr),nohrin
      common/hrdata3/xlst(nohr)
      common/hrdata_met/ctloc(kx,ky,kt),cta(kx,ky,kt),
     &       cea(kx,ky,kt),cwind(kx,ky,kt),csdn(kx,ky,kt),
     &       cxlwdn(kx,ky,kt),cpres(kx,ky,kt),
     &       clst(kx,ky,kt), clapse(ilg,jlg)
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1,
     &       h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1,
     &       sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1,
     &       albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2

      do ihr=1,nohr
       sdn(ihr)=BAD
       ta(ihr)=BAD
       xlwdn(ihr)=BAD
       rnet(ihr)=BAD
       g(ihr)=BAD
       rnsoil(ihr)=BAD
       pres(ihr)=BAD
       tloc(ihr)=BAD
       ea(ihr)=BAD
       wind(ihr)=BAD
       xlst(ihr)=BAD
      enddo

      do ihr = 1, nohr
       tloc(ihr)=ctloc(ia,ja,ihr)+dgmt
       ea(ihr) = cea(ia,ja,ihr)*1000.
       ta(ihr) = cta(ia,ja,ihr)-273.15
       wind(ihr) = cwind(ia,ja,ihr)
       pres(ihr) = cpres(ia,ja,ihr)/100.
       sdn(ihr) = BAD
       xlst(ihr) = BAD
       xlwdn(ihr) = cxlwdn(ia,ja,ihr)
!       write(6,*) tloc(ihr), ta(ihr), sdn(ihr), xlst(ihr)
      enddo
      nohrin=nohr

      return
      end

      subroutine byteswapr4(r)

c     does a byteswap on real*4 number

      integer*1 ii(4), jj(4)
      real*4 r, s, t
      equivalence (s,ii)
      equivalence (t,jj)

      s = r

      jj(1) = ii(4)
      jj(2) = ii(3)
      jj(3) = ii(2)
      jj(4) = ii(1)

      r = t

      return
      end

