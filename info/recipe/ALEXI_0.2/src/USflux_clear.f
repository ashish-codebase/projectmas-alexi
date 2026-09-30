      subroutine clear_day_proc(ibad)
c     **************************************************************            
c     *
c     *  
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      include 'ALEXI_parm.inc'
      common/availh2o/brz,awfrz,awcrz,bsfc,awfsfc,awcsfc
      common/flags/writeme,badinput,converged,stopiter
      logical converged,writeme,badinput,stopiter
      
      call hourly_flux(ibad)    
      if (badinput) return
      call daily_flux_clear_SDN      

      return
      end

      subroutine daily_flux_clear_SDN
c     **************************************************************            
c     *
c     *  
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2 
      common/data22/w2orig,eref2               
      common/dayflux/aparday,acday,rnday,rnsday,rncday,gday,
     &       eday,esday,ecday,hday,hsday,hcday,sday
      common/daypotflux/epotday,espotday,ecpotday
      common/flags/writeme,badinput,converged,stopiter  
      common/moisture/faw,awf      
      common/moisture2/faw2
      common/obs/hobs1,hobs2,xleobs1,xleobs2,gobs1,gobs2,
     &       rnobs1,rnobs2                 
      common/stress/fpet,fsdn
      logical writeme,badinput,converged,stopiter
      save
      
c	System fluxes
      fsdn=sdn2/sday
      eday=xle2/fsdn
      gday=0.		! Assume soil heat flux integrates to 0.
      hday=rnday-eday-gday
c      write(6,*)'EDAY HDAY RNDAY GDAY',eday,hday,rnday,gday

c	Soil fluxes 
      if (xle2.ne.0.) then    
        esday=(xles2/xle2)*eday
      else
        esday=0.0
      endif
      hsday=rnsday-gday-esday
c      write(6,*)'XLES RNS G RNS-G ',xles2,rnsoil2,g2,
c     &           rnsoil2-g2

c	Canopy flux as residual      
      ecday=eday-esday
      hcday=rncday-ecday
            
c       Total evapotranspiration factor
      call fao_PM(eref2,tloc2,sdn2,zen2,taobs2,
     &                    pres2,ea2,w2orig)
      fpet=xle2/eref2
      if (fpet.lt.0.0) then
        if(writeme)write(6,*)'FPET = ',fpet,' - set to 0'
        fpet=0.0
      endif 

      fawsfc=BAD	! Currently not supported
      fawrz=BAD

      return
      end
            
      subroutine daily_flux_clear_EF
c     **************************************************************            
c     *
c     *  OLD - UNSUPPORTED 
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2          
      common/dayflux/aparday,acday,rnday,rnsday,rncday,gday,
     &       eday,esday,ecday,hday,hsday,hcday,sday
      common/daypotflux/epotday,espotday,ecpotday
      common/flags/writeme,badinput,converged,stopiter  
      common/moisture/faw,awf      
      common/moisture2/faw2
      common/obs/hobs1,hobs2,xleobs1,xleobs2,gobs1,gobs2,
     &       rnobs1,rnobs2                 
      common/stress/fpet,fsdn
      logical writeme,badinput,converged,stopiter
      save
      
c	System fluxes
      fevap=1.1*xle2/(rnet2-g2)
      eday=fevap*(rnday-gday)
      hday=rnday-gday-eday
      write(6,*)'EDAY HDAY RNDAY GDAY',eday,hday,rnday,gday

c	Soil fluxes     
      fevaps=1.1*xles2/(rnsoil2-g2)
      esday=fevaps*(rnsday-gday)
      hsday=rnsday-gday-esday
c      write(6,*)'XLES RNS G RNS-G FEVAPS ',xles2,rnsoil2,g2,
c     &           rnsoil2-g2,fevaps

c	Canopy flux as residual      
      ecday=eday-esday
      hcday=rncday-ecday

       xlam=(2.501-0.00237*ta2)*1e6			! latent heat o' vaporization [J/kg]
       esat=6.108*10**(7.5*ta2/(237.3+ta2))		! [mb]
       s=2.17e-3*xlam*esat/((273.15+ta2)*(273.15+ta2))	! [mb/K]	

c               Potential canopy transpiration (W/m2)
       ecpot2=rndiv2*1.3*fg*(s/(s+0.66)) 
                            	    
c		Potential soil evaporation (W/m2)
c		(following Jury & Tanner, '76, Agron J. 68, 239-242)            
       tauc=0.5
       alpha=1.3		! End-point PT coefficient
       tau=exp(-0.45*xlai/sqrt(2.0*cos(zen2)))
       if (tau.le.tauc) then
         pts=1.0
       else
         pts=alpha-((alpha-1.0)*(1.0-tau))/(1.0-tauc)
       endif    
       espot2=(rnsoil2-g2)*pts*(s/(s+0.66))

c		System ET as sum of soil and canopy components (W/m2)
       epot2=ecpot2+espot2 

c	Vegetation stress factors      
c      call getstress(taobs2,ea2,ft,fvpd)
c      write(6,*)ecpotday,ft,fvpd,ecday
      fawrz=ecday/(ecpotday*ft*fvpd)
c      fawrz=xlec2/(ecpot2*ft*fvpd)
c      write(6,*)'FAWRZ:',fawrz,ecday,ecpotday,ft,fvpd
c      write(6,*)'EC:',fevap,rnday,gday
      if (fawrz.gt.1.0) then
        if(writeme)write(6,*)'FAWRZ = ',fawrz,' - NOT set to 1'
c        fawrz=1.0
      endif 
      if (fawrz.lt.0.0) then
        if(writeme)write(6,*)'FAWRZ = ',fawrz,' - set to 0'
        fawrz=0.0
      endif  
      
c	Soil surface evaporation factor      
      fawsfc=esday/espotday
c      fawsfc=xles2/espot2
      if (fawsfc.gt.1.0) then
        if(writeme)write(6,*)'FAWSFC = ',fawsfc,' - NOT set to 1'
c        fawsfc=1.0
      endif       
      if (fawsfc.lt.0.0) then
        if(writeme)write(6,*)'FAWSFC = ',fawsfc,' - set to 0'
        fawsfc=0.0
      endif       
c      write(6,*)'ESDAY ESPOTDAY FAWSF',esday,espotday,fawsfc
            
c       Total evapotranspiration factor
      faw=eday/epotday
      faw2=xle2/epot2
c      write(6,*)'FAW:',faw,eday,epotday
      if (faw.gt.1.0) then
        if(writeme)write(6,*)'FAW = ',faw,' - NOT set to 1'
c        faw=1.0
      endif 
      if (faw.lt.0.0) then
        if(writeme)write(6,*)'FAW = ',faw,' - set to 0'
        faw=0.0
      endif 

      return
      end
      

       
       
