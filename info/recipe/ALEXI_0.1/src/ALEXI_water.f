      subroutine alexi_water(ia,ja,ierr,ibad,iter)
c     **************************************************************            
c     *
c     *  
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      include 'ALEXI.inc'
      save
      
      xlai0=xlai
      fc0=fc
      
c  Initialize run
      call runinit(badinput,ibad,writeme)
      if (badinput) then
         ierr=2
         go to 3000
      endif   
  
      delhmx=0.1 	! Convergence parameter for H [W/m2]
      fcnew=fc
      ftheta=fthetanew

1000  continue		! Return to here to try new value of FC

c  Initialize state variables      
      psima=psi0
      hn=hn0
      tc1=taobs1	
      tc2=taobs2
      ta1=taobs1
      ta2=taobs2 
      ts1=taobs1
      ts2=taobs2
      tac1=taobs1
      tac2=taobs2      
      
c  Initialize counters & flags
      iter=0 
      ierr=0     
      converged=.FALSE.
      stopiter=.FALSE.
      
c  BEGIN FLUX CONVERGENCE LOOP
c-----------------------------------------------------------------

      do while (.not.converged)
      
         call findflux_water (hn,hnnew,converged,stopiter,ierr) 
         if (stopiter) go to 2000
         delthn=hnnew-hn
         if (abs(delthn).gt.delhmx) then
            if(iter.lt.300)then
               f=0.25/(1+int(iter/20))
               dh=f*delthn
            else
               ierr=1
               go to 2000
            endif
            hn=hn+dh
            iter=iter+1
         else
            converged=.TRUE.
         end if     
        
      enddo	! while not converged

c  END FLUX CONVERGENCE LOOP      
c-----------------------------------------------------------------
2000  continue	! Escape from flux convergence loop


c  SOLUTION DID NOT CONVERGE: 
c ------------------------------------------------------------------
      if (.not.converged) then

c 		Move fcnew closer to 0.5 and try again.
c  		Bail if fcnew reaches 0.5 and it still 
c               doesn't converge.  

c          fcnew=fcnew+sign(0.01,0.5-fcnew)
c          if (abs(fcnew-0.5).le.1e-2) then
c             nfail=nfail+1
c          else
c             fc=fcnew
c             call cover_props
c             call updatefc             
c             go to 1000
c          endif
                    
c  SOLUTION CONVERGED:
c ------------------------------------------------------------------
      else
         nconv=nconv+1
         ierr=0
      
c  		Save these values as best first guess for next pixel  
        
c         hn0=hn
c         psi0=psima
      endif

3000  continue      
      call ALEXI_errorcode(writeme,ierr)
      
      xlai=xlai0
      ts2=trad2
      tc2=trad2
      
      return
      end
      
      subroutine findflux_water (hn,hnnew,converged,stopiter,ierr)
c     **************************************************************            
c     *
c     *  
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      include 'ALEXI_parm.inc' 
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1,
     &       h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1,
     &       sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1,
     &       albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2  
      common/pbl/thrise,zpbli(mli),thpbli(mli),zpbl(ml),thpbl(ml),
     &       jzmax,jz1,nlev  
      common/stability/zdlamx,zdlamn   
      common/view/theta,ftheta                
      logical converged,stopiter
      save
      
      h1=hn*t1/thrise
      h2=hn*t2/thrise
1000  continue      	! Return to here to try new value of FC

c  Compute net radiation and soil heat flux   

      call getnetrad_water(trad1,xlwdn1,sdn1,rnet1,swup1,lwup1)
      call getnetrad_water(trad2,xlwdn2,sdn2,rnet2,swup2,lwup2)
      if(rnet2.lt.0.0)then        
         ierr=9
         stopiter=.TRUE.
         go to 2000
      endif 
      call getsoilheat_water(rnet1,g1)           
      call getsoilheat_water(rnet2,g2)
      
c  Find air temperature and component fluxes at times t1 and t2  

1500  continue
      zdlamx=.4
      zdlamn=-.5  
      call findta_water (trad1,ta1,rnet1,h1,xle1,g1,ra1,rs1,rx1,w1,
     &                       rhocp1,stopiter,ierr)
      if(stopiter)go to 2000
      
c      zdlamx=-.2		! System should not be stable at time t2
      zdlamx=0.0
      zdlamn=-10.0
      call findta_water (trad2,ta2,rnet2,h2,xle2,g2,ra2,rs2,rx2,w2,
     &                       rhocp2,stopiter,ierr)
      if(stopiter)go to 2000
            
c  Trap bad solutions for ts and tc and alter fc accordingly.  
c  Appears to be a high fc-triggered condition

      if (ts2.lt.tc2.and.fc.gt.0.1) then
c         stopiter=.TRUE.
c         ierr=11
c         go to 2000
      endif 

c  Ta2 must be > Ta1 for PBL model.  
c  Ta2 < Ta1 appears to be a low fc, high hn-triggered condition.
c  Setting ta2=ta1+2 yields low hn guess for next iter.
            
      if(ta2.lt.ta1+0.5)then
         ta2=ta1+0.5
      endif
      
c  Now, are Ta1, Ta2 consistent with H1 -> H2 heat input into the PBL?
c  Find new value of Hn and return for convergence check.

      call growPBL(hnnew,stopiter,ierr)
      
c  Limit the next guess for H2 to be <= RNET2-G2.  
c  Otherwise we may wander into Never-never Land.

      h2new=hnnew*t2/thrise
       
      h2max=rnet2-g2
      if(h2new.gt.h2max)hnnew=h2max*thrise/t2
            
2000  continue		! Escape from flux computation
     
      return
      end 
                        
      subroutine findta_water (trad,ta,rn,h,xle,g,ra,rs,rx,wind,
     &                       rhocp,stopiter,ierr)
c     **************************************************************            
c     *
c     *  
c     *
c     *  Martha Anderson
c     *  Created:  01/28/02
c     *  
c     *************************************************************** 
      logical stopiter
      
      call getresistance (h,wind,ta,ta,ta,ra,rs,rx,rhocp)
      ta=trad-h*ra/rhocp
      xle=rn-h-g
      
      return 
      end   
             
      subroutine getsoilheat_water (rnet,g)
c     **************************************************************            
c     *
c     *  This ties G to DISP/HEIGHT.  At very low covers, wind speed 
c     *  near soil is very high and most of RN at soil is convected 
c     *  away from soil as H. DISP/HEIGHT exhibits this behavior too
c     *  (approaches 0 rapidly as fc->0, for fc<0.05)
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 

      g=0.55*rnet
               
      return
      end

      subroutine getnetrad_water (trad,xlwdn,sdn,rnet,swup,lwup)
c     **************************************************************            
c     *
c     *  Compute net radiation above the canopy and above
c     *  the soil surface
c     *
c     *      LW: using sky, canopy and soil thermal fluxes and
c     *          canopy thermal transmission coefficient.
c     *      SW: using SDN components, canopy transmission and
c     *          soil reflectivity.
c     *
c     *  Martha Anderson
c     *  Created:  08/10/01
c     *  
c     *************************************************************** 
      common/emission/emleaf,emdead,emsoil,emcpy 
      common/partition/difvis,difnir,dirvis,dirnir,fvis,fnir 
      common/reflection/rsoiln,rsoilv,rdcpyl      
      common/transmission/taudn,taudv,taudl 
      
      parameter(sigma=5.67e-8) 
      albedo=0.1
      emissivity=0.99
      
      rsky    = xlwdn
      rwater  = 0.99*sigma*(trad+273.15)**4
      
      rnetlw  = rsky-rwater
      rnetsw  = sdn*(1-albedo)
      
      rnet=rnetlw+rnetsw
      
      swup=sdn*albedo
      lwup=rwater

      return
      end
      
