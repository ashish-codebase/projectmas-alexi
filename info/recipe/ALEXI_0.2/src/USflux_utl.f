c-------------------------------------------------------------------c
c            P O T E N T I A L   F L U X   R O U T I N E S          c
c-------------------------------------------------------------------c

      subroutine hourly_flux
c     **************************************************************
c     *
c     *  Computes hourly estimates of RN, RNSOIL, G and APAR from 
c     *  hourly values of XLWDN, SDN, TA and estimates of FPAR and FC.
c     *  The hourly fluxes are then integrated to produce daytime
c     *  total fluxes.
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     **************************************************************
      include 'USflux_parm.inc' 
      include 'ALEXI_parm.inc'
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump      
      common/cover/iclass,beta,gvs,xndvi,fpar,dstom         
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2          
      common/dayflux/aparday,acday,rnday,rnsday,rncday,gday,
     &       eday,esday,ecday,hday,hsday,hcday,sday  
      common/dayflux2/erefday 
      common/dayflux3/xlwupday,xlwdnday,swupday  
      common/dayobsflux/rnobsday,hobsday,eobsday,gobsday,
     &       aobsday,precpday 
      common/flags/writeme,badinput,converged,stopiter      
      common/hrdata/tloc(nohr),ta(nohr),ea(nohr),wind(nohr),
     &       sdn(nohr),xlwdn(nohr),rnet(nohr),rnsoil(nohr),
     &       g(nohr),pres(nohr),nohrin
      common/hrdata4/eref(nohr)
      common/hrdata5/xlwup(nohr),swup(nohr)      
      common/hrobs/rnobs(nohr),hobs(nohr),xleobs(nohr),gobs(nohr),
     &       aobs(nohr),precp(nohr)
      common/partition/difvis,difnir,dirvis,dirnir,fvis,fnir  
      common/departure/xmc,xbc,xcc,xms,xbs,xcs      
      common/site/xlat,xlong,stdlng 
      common/sun/strt,end 
      common/timestamp/year,doy             
      dimension apar(nohr)   
      logical flag
      logical writeme,badinput,converged,stopiter 
      save  
      
      if (converged) then
        call getTdepart2(xmc,xbc,xcc,tc2,taobs2,tloc2,strt,end)
        call getTdepart2(xms,xbs,xcs,ts2,taobs2,tloc2,strt,end)
      endif
           
      flag=.FALSE.
      
      do ihr=1,nohr
c        write(6,*)'TLOC SDN TA',tloc(ihr),sdn(ihr),ta(ihr)      
      enddo

      do ihr=1,nohr
c  	    Initialize to 0 for daily integration step
        rnet(ihr)  =0.
        rnsoil(ihr)=0.
        g(ihr)     =0.
        apar(ihr)  =0.
        eref(ihr)  =0.
        
        if (tloc(ihr).lt.strt.or.tloc(ihr).gt.end) sdn(ihr)=0.0
        call getsunzen (xlat,xlong,stdlng,doy,year,tloc(ihr),zen)
        if (tloc(ihr).ge.strt.and.tloc(ihr).le.end) then
          call checkvalue("SDN   ",40,sdn(ihr),-1.,1500.,flag,w,ibad)	! m/s
          call checkvalue("TA    ",38,ta(ihr),-50.,80.,flag,w,ibad)	! C
          if (flag) then
            badinput=.TRUE.
            ibad=49
            goto 500
          endif
          call getradcomps (sdn(ihr),zen,fclear)
          call getradprops (zen,albedo,taubtv,taubtn,clump)        ! CHANGE: clump should be clumps(ihr) !!!
c          if (cos(zen).lt.0.) goto 300   
            
c	    We don't know surface temperature (Ts and Tc) for all hours, so
c	    use air temperature as proxy baseline, applying dTs and dTc increments
          if (converged.or.xmc.ne.BAD) then
            dtc=xmc*tloc(ihr)*tloc(ihr)+xbc*tloc(ihr)+xcc
            dts=xms*tloc(ihr)*tloc(ihr)+xbs*tloc(ihr)+xcs
            if (dtc.lt.0) dtc=0.
            if (dts.lt.0) dts=0.         
c            tcproxy=ta(ihr)+dtc*fclear		! This causes artifacts in LWUP
c            tsproxy=ta(ihr)+dts*fclear
            tcproxy=ta(ihr)+dtc
            tsproxy=ta(ihr)+dts
           else
            tcproxy=ta(ihr)+2.0
            tsproxy=ta(ihr)+4.0 
          endif    
          call getnetrad (tsproxy,tcproxy,xlwdn(ihr),sdn(ihr),
     &                    albedo,taubtv,taubtn,rnet(ihr),
     &                    rnsoil(ihr),rndiv,swup(ihr),xlwup(ihr))
          call getsoilheat (rnsoil(ihr),g(ihr),tloc(ihr)) 
          call fao_PM(eref(ihr),tloc(ihr),sdn(ihr),zen,ta(ihr),
     &                    pres(ihr),ea(ihr),wind(ihr))

c          write(6,*) "rnet = ", rnet(ihr), tloc(ihr)  
 300	  continue 
         else  	! Nighttime
          if (sdn(ihr).eq.BAD) sdn(ihr)=0.0
          call checkvalue("SDN   ",40,sdn(ihr),-1.,1500.,flag,w,ibad)	! m/s
          call checkvalue("TA    ",38,ta(ihr),-50.,80.,flag,w,ibad)	! C
          if (flag) then
            badinput=.TRUE.
            ibad=49
            goto 500
          endif
          tcproxy=ta(ihr)
          tsproxy=ta(ihr) 
          call getradprops (zen,albedo,taubtv,taubtn,clump)        ! CHANGE: clump should be clumps(ihr) !!!
          call getnetradnight (tsproxy,tcproxy,xlwdn(ihr),sdn(ihr),
     &                    albedo,taubtv,taubtn,rnet(ihr),
     &                    dum,dum,swup(ihr),xlwup(ihr))
        endif                        
      enddo 	! IHR loop

c	Integrate 24-hr fluxes     
c      write(6,*) "RNET FINAL = ", rnet
      call integrate (rnet,rnday,strt,end,tloc,24)
      call integrate (rnsoil,rnsday,strt,end,tloc,24)
      call integrate (g,gday,strt,end,tloc,24)
      call integrate (sdn,sday,strt,end,tloc,24)       
      call integrate (swup,swupday,strt,end,tloc,24)      
      call integrate (xlwup,xlwupday,strt,end,tloc,24)      
      call integrate (xlwdn,xlwdnday,strt,end,tloc,24)      
      call integrate (eref,erefday,strt,end,tloc,24)      
      rncday=rnday-rnsday
 
 500  continue     
      return
      end

      subroutine fao_PM(eref,tloc,sdn,zen,ta,pres,ea,wind)
c     **************************************************************
c     *
c     *  Computes grass reference ET based primarily on FAO Penman-Monteith
c     *  (Equation 53 for hourly ETo)
c     *
c     *  Martha Anderson
c     *  Created:  08/31/11
c     *
c     **************************************************************
       
      parameter (cp=1010.e-6) 	! [MJ/kg-K]
      parameter (epsilon=0.622)  
      parameter (sigma=5.67e-8) ! [W/(m2-K4)]
      
      albedo=0.23 	! [ Hypothetical grass reference]
       
      xlam=(2.501-0.00237*ta)				! latent heat o' vaporization [MJ/kg]
      psych=0.1*pres*cp/(xlam*epsilon)			! psychometric constant [kPa/C]
      esat=0.6108*exp(17.2694*ta/(237.3+ta))		! [kPa]
      s=17.2694*237.3*esat/(237.3+ta)**2			! [kPa/C]	
      d=esat-0.1*ea					! [kPa]
      
      call getradcomps (sdn,zen,fclear)     						
      call getnetrad_simple (sdn,ta,fclear,zen,albedo,
     &                             rnet)	    
      rnetmj=rnet*3600./1.e6					! [MJ/m2-hr]
      gmj=0.1*rnetmj
            
      xnum1=0.408*s*(rnetmj-gmj)
      xnum2=psych*37.*wind*d/(ta+273.15)
      xnum=xnum1+xnum2
      xden=s+psych*(1.+0.24*wind)

      erefmmhr=xnum/xden						! [mm/hr]
      eref=erefmmhr*xlam*1.e6/3600.				! [W/m2]
       	     
      return
      end
      
      subroutine integrate(f,fint,tstrt,tend,t,nhr) 
c     **************************************************************
c     *
c     *
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     **************************************************************
      include 'ALEXI_parm.inc'  
      include 'USflux_parm.inc'  
                                 
      logical fbad	
      dimension f(nohr),t(nohr)		! CHANGE THESE DIMENSIONS!
      
      fbad=.FALSE.
      fint=0.      

      do i=1,nhr
c          write(6,*)t(i),f(i)
          if (f(i).eq.BAD) fbad=.TRUE.               
          fint=fint+f(i)          
      enddo	! I loop
      
      xn=24./(t(2)-t(1))
      
      if (fbad) then
         fint = BAD
      else
         fint=(fint/xn)*3600.*24/1.e6		! Integrated flux in MJ/day
      endif

      return
      end
           
      subroutine getTdepart2(xm,xb,xc,trad,ta,t2,tr,ts)
c     **************************************************************
c     *
c     *
c     *
c     *  Martha Anderson
c     *  Created:  02/22/11
c     *
c     **************************************************************
      
      a=trad-ta
      
      xb=a/(((t2*t2-tr*tr)*(tr-ts)/(ts*ts-tr*tr))+t2-tr)
      xm=xb*(tr-ts)/(ts*ts-tr*tr)
      xc=-xb*tr-xm*tr*tr

      return
      end
      
c-------------------------------------------------------------------c
c             I N P U T / O U T P U T   U T I L I T I E S           c
c-------------------------------------------------------------------c
                          
      subroutine checkUSinput(flag,ibad,w)
c     ***************************************************************
c     *
c     *  All input variables are checked for unreasonable values.
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************
      include 'USflux_parm.inc'
      include 'ALEXI_parm.inc'    
      common/availh2o/brz,awfrz,awcrz,bsfc,awfsfc,awcsfc        
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/cover/iclass,beta,gvs,xndvi,fpar,dstom
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
      common/pbl/thrise,zpbli(mli),thpbli(mli),zpbl(ml),thpbl(ml),
     &       jzmax,jz1,nlev      
      common/timestamp/year,doy 
      common/view/theta,ftheta      
      save
      
      logical flag,w,allbad

      flag=.FALSE.
      allbad=.TRUE.
      
c  Hourly quantities:
c  NOTE:  There are some CRAS files that are missing hours, meaning
c  not all nohrin (24) slots will be filled.  Check only those hour
c  slots that are actually used (sdn>0).
      do ihr=1,nohrin
        if (sdn(ihr).ne.BAD) allbad=.FALSE.
        if (sdn(ihr).gt.0.) then
c          call checkvalue("EA    ",37,ea(ihr),0.,80.,flag,w,ibad)	! mbar
c          call checkvalue("TA    ",38,ta(ihr),-50.,60.,flag,w,ibad)	! C
c          call checkvalue("WIND  ",39,wind(ihr),0.,100.,flag,w,ibad)	! m/s
c          call checkvalue("SDN   ",40,sdn(ihr),-1.,1500.,flag,w,ibad)	! m/s
c          call checkvalue("XLWDN ",41,xlwdn(ihr),0.,1000.,flag,w,ibad)	! m/s  
           call checkvalue("PRES  ",37,pres(ihr),0.,1500.,flag,w,ibad)	! mbar
           if (flag) goto 100
        endif
      enddo  
      
c  Check for missing hours
      do ihr=1,nohrin
        if (ihr.ne.1) then
          dtime=tloc(ihr)-tloc(ihr-1)
          if (dtime.gt.1) 
     &      write(6,*)'CRAS HR GAP BETWEEN TLOCS',tloc(ihr-1),tloc(ihr)
        endif      
      enddo 
      
c  Canopy parameters:
      xclass=iclass
      call checkvalue("ICLASS",42,xclass,1.,27.,flag,w,ibad)	
      call checkvalue("FG    ",44,fg,0.,1.,flag,w,ibad)
      call checkvalue("XLAI  ",22,xlai,0.,10.,flag,w,ibad)

c  Insolation data:
      if (allbad) then
        flag=.TRUE.
        ibad=47
      endif
      
 100  continue
 
      return
      end           

      subroutine setPBLheights
c     ***************************************************************
c     *
c     *  Set up the input PBL height array ZPBLI.  At present, PBL
c     *  potential temperatures are input at 200m spacings up to
c     *  a maximum height of 5km.
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************
      include 'ALEXI_parm.inc'             
      common/pbl/thrise,zpbli(mli),thpbli(mli),zpbl(ml),thpbl(ml),
     &           jzmax,jz1,nlev     
      save
      
      zinit=0.0				! Define lowest level to be at 0m
      zinc=200.0
      do j=1,mli
         zpbli(j)=zinit+(j-1.0)*zinc
      enddo
      
      return
      end 
       
