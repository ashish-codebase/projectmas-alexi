      subroutine cloudy_day_proc(ibad)
c     **************************************************************            
c     *
c     *  
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      include 'ALEXI_parm.inc'
      common/flags/writeme,badinput,converged,stopiter
      logical converged,writeme,badinput,stopiter
      
      call hourly_flux(ibad)
      if (badinput) return      
      call daily_flux_cloudy
      
      return
      end
      
      subroutine daily_flux_cloudy
c     **************************************************************            
c     *
c     *  Do not compute daily fluxes for cloudy pixels - gap-fill 
c     *  in post-processing.
c     *  
c     ***************************************************************
      common/dayflux/aparday,acday,rnday,rnsday,rncday,gday,
     &       eday,esday,ecday,hday,hsday,hcday,sday
      
      eday=BAD
      hday=BAD
      ecday=BAD
      esday=BAD
      hcday=BAD
      hsday=BAD
           
      return
      end      
      

           
