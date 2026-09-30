c-------------------------------------------------------------------c
c         A E R O D Y N A M I C   R E S I S T A N C E   I I I       c
c								    c
c  This file contains code computing canopy architecture parameters c 
c  relating to aerodynamic resistances using equations developed by c  
c  Massman & Weil (1996).  		   			    c		   
c								    c
c  Variables and constants used:                                    c
c                                                                   c
c        cd: leaf drag coefficient (0.2)                            c
c      disp: displacement height                                    c
c        fc: fraction cover                                         c
c       LAI: leaf area index                                        c
c        ra: aerodynamic resistance                                 c
c        rs: soil surface resistance                                c
c         u: measured wind speed at height, z #INPUT#               c
c        uc: wind speed at canopy height, zc                        c
c        us: wind speed at soil surface                             c
c    refhtw: height at which u is measured  #INPUT#                 c
c        z0: canopy roughness length                                c                                                                   c
c                                                                   c
c  Subroutines: 						    c
c                                                                   c  
c       getresistance:  Calculates RS, RA, and RX.  Below FC=0.3,   c
c			RS is decreased toward 0.		    c
c                                                                   c
c          canopyarch:  Calculates values for Z0 and DISP for a     c
c                       given value of FC based on formalism        c
c                       of Massman & Weil (1996).  Also updates     c
c                       other variables dependent on FC.            c
c								    c
c              psimhn:  Calculates stability correction factors     c
c                       for aerodynamic resistances.                c
c                                                                   c  								    c
c-------------------------------------------------------------------c

      subroutine getresistance (h,wind,ta,ts,tc,ra,rs,rx,rhocp)
c     ***************************************************************            
c     *
c     *  Returns aerodynamic resistance (RA), soil resistance
c     *  (RS) and 1-sided leaf boundary layer resistance (RX).  All
c     *  resistances are returned in units of s m-1 and are
c     *  on a per unit ground area basis.
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************   
      common/aerodyn/xlog1,xlog2,a,uexp1,uexp2,expuxp,psima,psih      
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump,rsmin    
      common/cover2/perennial,iswater,fcbare 
      common/cover/iclass,beta,gvs,xndvi,fpar,dstom
      common/constants/pi,cp,xk 
      common/stability/zdlamx,zdlamn  
      logical perennial,iswater     
      save         
c 
c	Find Z/L and stability functions
c         z/L > 0 => stable atmosphere
c         z/L < 0 => unstable
c         z/L = 0 => neutral stability
     
c   	First, use stability from previous iteration

      ustarb=.4*wind/(xlog1-psima)      
      zdla=-(refhtw-disp)*(h)*9.8*0.4/
     &      (rhocp*(ta+273.)*ustarb*ustarb*ustarb)
      zdla=min(zdla,zdlamx)
      zdla=max(zdla,zdlamn)
      call psimhn_ky(zdla,psima,psih)  
      
c   	Now recalculate stability corrected ustar

      ustara=.4*wind/(xlog1-psima)
      zdla=zdla*ustarb*ustarb*ustarb/(ustara*ustara*ustara)
      zdla=min(zdla,zdlamx)
      zdla=max(zdla,zdlamn)  
      call psimhn_ky(zdla,psima,psih)
      ustar=.4*wind/(xlog1-psima)
      
c			BARE SOIL OR WATER      
c	Rs and Rx are set explicitly for bare soil or water case 
c 	Ra is computed using roughness length for heat in place
c	of Z0M.

      if (iswater .or.
     &   (.not.perennial .and.
     &               (fc.le.fcbare .or. height.eq.0.0))) then
        tas=273.15+(ta+ts)*0.5
        xmuo=1.8325e-5
        xmu=(296.16+393.16)/(tas+393.16)*(tas/296.16)**1.5*xmuo
        rho=rhocp/cp
        xnu=xmu/rho
        restar=ustar*z0/xnu
        z0h=z0*exp(-xk*(4.0*restar**(0.15)-5.))
        ra=(log((refhtw-disp)/z0h)-psih)/(ustar*xk)
        rs=0.0
        rx=10000.0
        return
      endif

c			ALL OTHER CASES      
c	Find Ra (Aerodynamic resistance)
      
      ra=(xlog1-psima)*(xlog1-psih)/(0.16*wind)

c   	Find Rs (Soil resistance)

c      tgrad=ts-tc
      tgrad=ts-ta	! TC gets crazy at low fc
      zdlamx2=0.4
      zdlamn2=-0.5       
      zdla2=zdla*(height-disp)/(refhtw-disp)
      zdla2=min(zdla2,zdlamx2)
      zdla2=max(zdla2,zdlamn2)  
      call psimhn_ky(zdla2,psima2,psih2) 
      uc=wind*((xlog2-psima2)/(xlog1-psima))       
      us=uc*expuxp
!      if (xndvi.gt.0.10.and.xndvi.le.0.40) then
!       val1=-0.0043*xndvi+0.0042
!       val2=-0.04*xndvi+0.028
!      endif
!      if (xndvi.gt.0.00) then
       val1=0.0025
       val2=0.012
!      endif
!      if (xndvi.le.0.10) then
!       val1=0.0038
!       val2=0.024
!      endif
      if(tgrad.gt.1.0)then
         rs=1./(val1*(tgrad)**0.33+val2*us)
      else
         rs=1./(val1+val2*us)
      endif

!      rs=10.
!       rs=rs*rsmin
!       write(6,*) "rsmin = ", rs, rsmin
c	Find Rx (Boundary layer resistance for 1 side of
c       a leaf, corresponding to the wind speed at z0+disp)

      udz=uc*exp(uexp2) 
      rx=(180.*sqrt(xl/udz))/xlai
!      rx=rsmin/xlai

      return
      end

      subroutine psimhn_ky(zdla,psima,psih)
c     ***************************************************************            
c     *
c     *  Calculates stability functions using forms based on the
c     *  work of Kader and Yaglom (1990) reported by Sugita et al.
c     *  (Flux determination over a smooth surface)
c     *
c     *  Martha Anderson
c     *  Created:  01/16/02
c     *  
c     *************************************************************** 

      y=-zdla
      y0=0.0
      if(y.lt.0.0059)then
        psih=0.0
        psima=0.0
      else
        psih=1.2*log((0.33+y**0.78)/0.33)
        y=min(y,15.025)
        psima=1.47*log((0.28+y**0.75)/(0.28+(0.0059+y0)**0.75))
     &        -1.29*(y**0.33-(0.0059+y0)**0.33)
      endif
      
      return
      end
      
      subroutine psimhn(zdla,psima,psih)
c     ***************************************************************            
c     *
c     *  Calculates stability functions using the simplest possible
c     *  scheme from the Cupid model
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 

      if(zdla.ge.0.0)then
        psima=-5.0*zdla
        psih=-5.0*zdla
      else
        psima=1.88+(1./(-.533+.790*zdla))  
        psih=exp(.598+.390*alog(-zdla)-
     &                              .09*(alog(-zdla))**2)
      endif  
      
      return
      end
      
      subroutine canopyarch
c     **************************************************************            
c     *
c     *  Updates canopy architecture factors derived from FC.
c     *  Formulae for dispdh and z0dh are from Massman & Weil (1996).
c     *
c     *  THIS VERSION USES WEIGHTED DISP AND Z0 COMPUTED OUTSIDE
c     *  OF ALEXI.
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      common/aerodyn/xlog1,xlog2,a,uexp1,uexp2,expuxp,psima,psih 
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/clumping/clumps1,clumps2,clump0,fveg
      common/cover2/perennial,iswater,fcbare
      common/constants/pi,cp,xk                  
      common/view/theta,ftheta 
      logical perennial,iswater              
      save  
      
      cd=0.20 
      z0s=0.005				! Roughness length for bare soil 
      z0w=0.00035			! Roughness length for water      
      
c	Cover fraction variables    
c      xlai=-2.0*alog(1.-fc)		! Compute LAI assuming random canopy
      ftheta=1.0-exp(-0.5*xlai*clump/cos(theta))
      if (ftheta.gt.0.8) ftheta=0.8

c	Rigid vegetation case - hold z0dh & dispdh fixed
c       --------------------------------------------------
      if (perennial) then
c         z0dh=0.1       
c         z0=z0dh*height
c         dispdh=0.7
c         disp=dispdh*height   

c	Surface water case - z0 is water roughness
c       --------------------------------------------------
      else if (iswater) then
         z0=z0w
         disp=0.0
         xlog1=alog(refhtw/z0)
         go to 500
                      
c	Bare soil case - z0 is soil roughness
c       --------------------------------------------------
      else if (fc.le.fcbare .or. height.eq.0.0) then
c         z0=z0s
c         disp=0.0
         xlog1=alog(refhtw/z0)
         go to 500
                  
c	All other cases - use Massman eqs for z0dh & dispdh 
c       --------------------------------------------------
      else
c         ufact=0.360-0.264*exp(-15.1*cd*xlai)	! THIS IS ALL BEING DONE IN LANDCOVER.F
c         xn=cd*xlai/(2.*(ufact)**2)
c         dispdh=0.7-(1./(5.*xn)*(1.-exp(-3.3*xn)))
c         if(dispdh.lt.0.)dispdh=0.       
c         disp=dispdh*height
c         z0dh=(1.-dispdh)*exp(-0.4/ufact)
c         z0=z0dh*height
         z0=max(z0,z0s)
c  FOR POSITVE XLOG2, NEED HEIGHT-DISP>Z0, otherwise negative UC
         height=max(height,disp+z0+0.001)
      endif
      
!	Wind extinction factor, A (Goudriaan '77, page 110)
      xld=xlai/height
      xlm=sqrt((4.*xl)/(pi*xld))
      a=sqrt(cd*clump*xlai*height/xlm)	! Assumes iw=0.5

c	Other factors (to save computation in oft-called GETRESISTANCE)      
      xlog1=alog((refhtw-disp)/z0) 
      xlog2=alog((height-disp)/z0)      
      if (height.gt.0.5) then
        uexp1=-a*(1.-0.05/height)
      else
        uexp1=-a*0.90
      endif
      expuxp=exp(uexp1)
      if(expuxp.gt.0.95)expuxp=0.95
      a=sqrt(cd*(xlai/fveg)*height/xlm)  ! Assumes iw=0.5
      uexp2=-a*(1.-(z0+disp)/height) 

 500  continue      
      return
      end
            

