       subroutine getsunzen (xlat,xlong,stdlng,doy,year,ftime,zen)
c     **************************************************************** 
c     *           
c     *  Computes sun zenith angle ZEN on day DAY, year YEAR at 
c     *  at coordinates XLAT, XLONG (from Cupid).
c     *
c     *  XLONG and STDLNG can be neg for W longitude, or 0-360
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ****************************************************************
      common/constants/pi,cp,xk
      common/sun/strt,end  
           
      pid180=pi/180.
      pid2=pi/2.

c		Latitude computations

      sinlat=sin(xlat*pid180)
      coslat=cos(xlat*pid180)
      
c 		Declination computations

      kday=(year-1977)*365+doy+28123
      xm=(-1.+.9856*kday)*pid180
      delnu=2.*.01674*sin(xm)+1.25*.01674*.01674*sin(2.*xm)
      slong=(-79.8280+.9856479*kday)*pid180+delnu      
      decmax=sin(23.44*pid180)
      decl=asin(decmax*sin(slong))
      sindec=sin(decl)
      cosdec=cos(decl)      
      hfday=12./pi*acos(-(sinlat*sindec)/(coslat*cosdec))
      eqtm=9.4564*sin(2.*slong)/cosdec-4.*delnu/pid180
      eqtm=eqtm/60.  
     
c		Longitude corrections  
      
      dlong=xlong-stdlng
      if (dlong.gt. 180) dlong=dlong-360.
      if (dlong.lt.-180) dlong=dlong+360.
      dlong=(dlong)/15.  		! DLONG positive E of STDLNG
      
      strt=12.-dlong-eqtm-hfday		! This is a true sunrise LT
      end =12.-dlong-eqtm+hfday

c		Get sun zenith angle
c      
      timsun=ftime+eqtm+dlong
      hrang=(timsun-12.)*pid2/6.
      zen=acos(sinlat*sindec+coslat*cosdec*cos(hrang))

      return
      end
           
      subroutine getradprops (zen,albedo,taubtv,taubtn,clumps)
c     **************************************************************            
c     *
c     *  Computes albedo based on analytical solutions from Goudriaan
c     *  1988 and summarized in Campbell & Norman 1998.  These assume
c     *  that leaf transmittivity equals leaf reflectivity.  All
c     *  we must specify is leaf absorptivity in NIR and VIS
c     *  wavebands (ALEAFN, ALEAFV).  Also soil reflectivities
c     *  RSOILN and RSOILV.
c     *
c     *  ZEN-dependent quantities (ALBEDO, TAUBTV and 
c     *  TAUBTN) are returned in argument list.
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      common/absorption/aleafn,aleafv,aleafl,adeadv,adeadn,adeadl
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump        
      common/emission/emleaf,emdead,emsoil,emcpy             
      common/partition/difvis,difnir,dirvis,dirnir,fvis,fnir                
      common/reflection/rsoiln,rsoilv,rdcpyl
      common/reflection2/albv,albn,albvobs,albnobs,albobs
      common/transmission/taudn,taudv,taudl 
      save
      
      coszen=cos(zen) 
c      write(6,*) "coszen zen =", coszen, zen     
c	Weighted live/dead leaf average properties
      ameanv = aleafv*fg + adeadv*(1-fg)
      ameann = aleafn*fg + adeadn*(1-fg)
      ameanl = aleafl*fg + adeadl*(1-fg)
c      write(6,*) "v = ", aleafv, adeadv
c      write(6,*) "n = ", aleafn, adeadn
c      write(6,*) "l = ", aleafl, adeadl

c		   D I F F U S E   C O M P O N E N T S
c                  -----------------------------------

c	Diffuse light canopy reflection coefficients 
c	for a deep canopy	
      akd=-0.0683*log(clumps*xlai)+0.804		! Fit to Fig 15.4 for x=1
      rcpyn=(1.0-sqrt(ameann))/(1.0+sqrt(ameann))   	! Eq 15.7   
      rcpyv=(1.0-sqrt(ameanv))/(1.0+sqrt(ameanv))
      rcpyl=(1.0-sqrt(ameanl))/(1.0+sqrt(ameanl))
      rdcpyn=2.0*akd*rcpyn/(akd+1.0)   			! Eq 15.8      
      rdcpyv=2.0*akd*rcpyv/(akd+1.0) 
      rdcpyl=2.0*akd*rcpyl/(akd+1.0) 
      
c	Diffuse canopy transmission coeff (visible) 				
      expfac = sqrt(ameanv)*akd*clumps*xlai
      xnum = (rdcpyv*rdcpyv-1.0)*exp(-expfac)
      xden = (rdcpyv*rsoilv-1.0)+rdcpyv*(rdcpyv-rsoilv)*exp(-2.0*expfac)
      taudv = xnum/xden					! Eq 15.11

c	Diffuse canopy transmission coeff (NIR) 				
      expfac = sqrt(ameann)*akd*clumps*xlai
      xnum = (rdcpyn*rdcpyn-1.0)*exp(-expfac)
      xden = (rdcpyn*rsoiln-1.0)+rdcpyn*(rdcpyn-rsoiln)*exp(-2.0*expfac)
      taudn = xnum/xden					! Eq 15.11
      
c	Diffuse canopy transmission coeff (longwave) 
c	Use deep canopy expression (ignore soil reflection effects)				
      taudl = exp(-sqrt(ameanl)*akd*clumps*xlai)	! Eq 15.6
c      emcpy=1.0-rdcpyl-taudl
c Bill K's emcpy:
      emcpy = 0.99*(1.-taudl)    
      
      
c	Diffuse radiation surface albedo
c	for a generic canopy
      fact=((rdcpyn-rsoiln)/(rdcpyn*rsoiln-1.0))*	! Eq 15.9
     &       exp(-2.0*sqrt(ameann)*akd*clumps*xlai)
      albdn=(rdcpyn+fact)/(1.0+rdcpyn*fact) 
      fact=((rdcpyv-rsoilv)/(rdcpyv*rsoilv-1.0))*	! Eq 15.9
     &       exp(-2.0*sqrt(ameanv)*akd*clumps*xlai)
      albdv=(rdcpyv+fact)/(1.0+rdcpyv*fact) 
      
c	Diffuse is all we can compute at low zenith angles...
!      if(coszen.lt.0.01)then
!        albedo=fvis*(difvis*albdv)+
!     &         fnir*(difnir*albdn)        
!        taubtv=0.0
!        taubtn=0.0     
!        goto 1000
!      endif   
          
c		     B E A M   C O M P O N E N T S
c                    -----------------------------
      
c     	Direct beam extinction coeff (spher. LAD)  
      akb=0.5/coszen	
c      write(6,*) "akb = ", akb
c	Direct beam canopy reflection coefficients 
c	for a deep canopy
      rcpyn=(1.0-sqrt(ameann))/(1.0+sqrt(ameann))   	! Eq 15.7   
      rcpyv=(1.0-sqrt(ameanv))/(1.0+sqrt(ameanv))
      rbcpyn=2.0*akb*rcpyn/(akb+1.0)   			! Eq 15.8      
      rbcpyv=2.0*akb*rcpyv/(akb+1.0) 
	
c	Direct beam radiation surface albedo 
c	for a generic canopy
      fact=((rbcpyn-rsoiln)/(rbcpyn*rsoiln-1.0))*	! Eq 15.9
     &       exp(-2.0*sqrt(ameann)*akb*clumps*xlai)
      albbn=(rbcpyn+fact)/(1.0+rbcpyn*fact) 
      fact=((rbcpyv-rsoilv)/(rbcpyv*rsoilv-1.0))*	! Eq 15.9
     &       exp(-2.0*sqrt(ameanv)*akb*clumps*xlai)
      albbv=(rbcpyv+fact)/(1.0+rbcpyv*fact)
c      write(6,*) "clumpts xlai =", clumps, xlai 
c	Weighted average albedo 
      albedo=fvis*(dirvis*albbv+difvis*albdv)+
     &       fnir*(dirnir*albbn+difnir*albdn)
      albv=dirvis*albbv+difvis*albdv
      albn=dirnir*albbn+difnir*albdn
c       write(6,*) "fvis fnir = ", fvis, fnir
c       write(6,*) "dirvis difvis =", dirvis, difvis
c       write(6,*) "albbv albdv = ", albbv, albdv
c       write(6,*) "dirnir difnir =", dirnir, difnir
c       write(6,*) "albbn albdn = ", albbn, albdn 
c       write(6,*) "albedo =", albedo
c	Direct beam+scattered canopy transmission coeff (visible) 				
      expfac = sqrt(ameanv)*akb*clumps*xlai
      xnum = (rbcpyv*rbcpyv-1.0)*exp(-expfac)
      xden = (rbcpyv*rsoilv-1.0)+rbcpyv*(rbcpyv-rsoilv)*exp(-2.0*expfac)
      taubtv = xnum/xden				! Eq 15.11
c      write(6,*) "taubtv (COVER) = ", taubtv

c	Direct beam+scattered canopy transmission coeff (NIR) 				
      expfac = sqrt(ameann)*akb*clumps*xlai
      xnum = (rbcpyn*rbcpyn-1.0)*exp(-expfac)
      xden = (rbcpyn*rsoiln-1.0)+rbcpyn*(rbcpyn-rsoiln)*exp(-2.0*expfac)
      taubtn = xnum/xden				! Eq 15.11

1000  continue

      return
      end
      
      subroutine getnetrad (ts,tc,xlwdn,sdn,albedo,taubtv,taubtn,
     &                      rnet,rnsoil,rndiv,swup,xlwup)
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
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      common/emission/emleaf,emdead,emsoil,emcpy 
      common/emission2/esfc,aem,bem,eleaf,esoil
      common/partition/difvis,difnir,dirvis,dirnir,fvis,fnir 
      common/reflection/rsoiln,rsoilv,rdcpyl      
      common/transmission/taudn,taudv,taudl 
      
      parameter(sigma=5.67e-8) 
      rsoill=1-emsoil
      tausolar=fvis*(difvis*taudv+dirvis*taubtv)+
     &         fnir*(difnir*taudn+dirnir*taubtn)
      tauthermal=taudl
      albedocanopy=albedo
      albedosoil=fvis*rsoilv+fnir*rsoiln
      
      rsky  = xlwdn
      rcpy  = 0.99*sigma*(tc+273.15)**4   ! BK's equations
      rsoil = emsoil*sigma*(ts+273.15)**4

c	RNET above soil
      rnsoil=tauthermal*rsky+(1-tauthermal)*rcpy-
     &         rsoil+tausolar*(1-albedosoil)*sdn
      
c	RNET divergence within the canopy  
      rndiv=((1.0-tauthermal)*(rsky+rsoil-(2.0*
     &            rcpy)))+((1.0-tausolar)*(1.0-
     &            albedocanopy)*sdn)  
       
c	RNET above canopy  
      rnet=rndiv+rnsoil
      
c	Upwelling components 
      rnsoillw=tauthermal*rsky+(1-tauthermal)*rcpy-rsoil
      rndivlw=(1.0-tauthermal)*(rsky+rsoil-(2.0*rcpy))
      rnlw=rnsoillw+rndivlw
      xlwup=xlwdn-rnlw
      
      rnsoilsw=tausolar*(1-albedosoil)*sdn
      rndivsw=((1.0-tausolar)*(1.0-albedocanopy)*sdn)
      rnsw=rnsoilsw+rndivsw
      swup=sdn-rnsw
      
      xlwup=sdn-swup+xlwdn-rnet
c      write(6,*)'TAUDL RNSLW RNCLW',taudl,rnsoillw,rndivlw
               
      return
      end

      subroutine getnetradnight (ts,tc,xlwdn,sdn,albedo,taubtv,
     &                      taubtn,rnet,rnsoil,rndiv,swup,xlwup)
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
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      common/emission/emleaf,emdead,emsoil,emcpy 
      common/emission2/esfc,aem,bem,eleaf,esoil
      common/partition/difvis,difnir,dirvis,dirnir,fvis,fnir 
      common/reflection/rsoiln,rsoilv,rdcpyl      
      common/transmission/taudn,taudv,taudl 
      
      parameter(sigma=5.67e-8) 
      rsoill=1-emsoil
      tausolar=fvis*(difvis*taudv+dirvis*taubtv)+
     &         fnir*(difnir*taudn+dirnir*taubtn)
      tauthermal=taudl
      albedocanopy=albedo
      albedosoil=fvis*rsoilv+fnir*rsoiln
      
      rsky  = xlwdn
      rcpy  = 0.99*sigma*(tc+273.15)**4   ! BK's equations
      rsoil = emsoil*sigma*(ts+273.15)**4

c	RNET above soil
      rnsoil=tauthermal*rsky+(1-tauthermal)*rcpy-
     &         rsoil
      
c	RNET divergence within the canopy  
      rndiv=((1.0-tauthermal)*(rsky+rsoil-(2.0*
     &            rcpy)))  
       
c	RNET above canopy  
      rnet=rndiv+rnsoil
      
c	Upwelling components 
      rnsoillw=tauthermal*rsky+(1-tauthermal)*rcpy-rsoil
      rndivlw=(1.0-tauthermal)*(rsky+rsoil-(2.0*rcpy))
      rnlw=rnsoillw+rndivlw
      xlwup=xlwdn-rnlw
      
      rnsoilsw=0.
      rndivsw=0.
      rnsw=rnsoilsw+rndivsw
      swup=sdn-rnsw
      
      xlwup=sdn-swup+xlwdn-rnet
c      write(6,*)'TAUDL RNSLW RNCLW',taudl,rnsoillw,rndivlw
               
      return
      end
      
      subroutine getxlwdn (ea,ta,fclear,xlwdn)
c     ***************************************************************
c     *
c     *  Calculates thermal radiation from sky with Brutsaert Eq.
c     *  Uses ratio to get weighted average of clear sky and 'clouds'.
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************
      parameter(sigma=5.67e-8)

      tak=ta+273.15
      esky=1.24*(ea/tak)**(1./7.)
      xlwdn=sigma*(tak**4)*(esky*fclear+1.-fclear)
      return
      end
      
      subroutine find_albedo_soil(tloc2)
c     ***************************************************************
c     *
c     *  
c     *
c     *  Martha Anderson
c     *  Created:  02/09/11
c     *
c     ***************************************************************
      common/absorption/aleafn,aleafv,aleafl,adeadv,adeadn,adeadl      
      common/clumping/clumps1,clumps2,clump0,fveg      
      common/partition/difvis,difnir,dirvis,dirnir,fvis,fnir                
      common/reflection/rsoiln,rsoilv,rdcpyl
      common/reflection2/albv,albn,albvobs,albnobs,albobs
      common/site/xlat,xlong,stdlng   
      common/timestamp/year,doy 
      logical done

      errmax=0.002

c  Assume clear-sky partitioning     
      fvis=0.5
      fnir=0.5
      dirvis=0.8
      difvis=0.2
      dirnir=1.0
      difnir=0.0
      
      call getsunzen (xlat,xlong,stdlng,doy,year,tloc2,zen2)
      
      done=.FALSE.
      iter=1
      alast=-9999
      do while (.not.done)
        call getradprops (zen2,albedo,taubtv,taubtn,clumps2)
        write(6,100)'RSOILV: ',iter,rsoilv,albv,albvobs       
        adiff=albvobs-albv
        if (abs(adiff).gt.errmax) then
          delta=albv-alast
          if(iter.eq.1) then
            deltr=0.02*abs(adiff)/adiff 	! change by 0.02 with sign of adiff
          else if(iter.lt.20)then
            dadr=delta/deltr
            dr=adiff/dadr
            deltr=0.75*dr			! nudge by 0.75*DR
          else
            rsoilv=BAD
            done=.TRUE.
          endif
          if (.not.done) then
            rsoilv=rsoilv+deltr
            alast=albv
            iter=iter+1
          endif
        else
          done=.TRUE.
        endif 
      enddo
     
      done=.FALSE.
      iter=1
      alast=-9999
      do while (.not.done)
        call getradprops (zen2,albedo,taubtv,taubtn,clumps2)
        write(6,100)'RSOILN: ',iter,rsoiln,albn,albnobs         
        adiff=albnobs-albn
        if (abs(adiff).gt.errmax) then
          delta=albn-alast
          if(iter.eq.1) then
            deltr=0.02*abs(adiff)/adiff 	! change by 0.02 with sign of adiff
          else if(iter.lt.20)then
            dadr=delta/deltr
            dr=adiff/dadr
            deltr=0.75*dr			! nudge by 0.75*DR
          else
            rsoiln=BAD
            done=.TRUE.
          endif
          if (.not.done) then
            rsoiln=rsoiln+deltr
            alast=albn
            iter=iter+1
          endif
        else
          done=.TRUE.
        endif
      enddo

 100  format(a10,i5,3f11.5)
       
      return
      end
      
      subroutine find_albedo_veg(tloc2)
c     ***************************************************************
c     *
c     *  
c     *
c     *  Martha Anderson
c     *  Created:  02/09/11
c     *
c     ***************************************************************
      common/absorption/aleafn,aleafv,aleafl,adeadv,adeadn,adeadl      
      common/clumping/clumps1,clumps2,clump0,fveg      
      common/partition/difvis,difnir,dirvis,dirnir,fvis,fnir                
      common/reflection/rsoiln,rsoilv,rdcpyl
      common/reflection2/albv,albn,albvobs,albnobs,albobs
      common/site/xlat,xlong,stdlng   
      common/timestamp/year,doy 
      logical done

      errmax=0.002

c  Assume clear-sky partitioning     
      fvis=0.5
      fnir=0.5
      dirvis=0.8
      difvis=0.2
      dirnir=1.0
      difnir=0.0
      
      call getsunzen (xlat,xlong,stdlng,doy,year,tloc2,zen2)
      
      done=.FALSE.
      iter=1
      alast=-9999
      do while (.not.done)
        call getradprops (zen2,albedo,taubtv,taubtn,clumps2)
        write(6,100)'ALEAFV: ',iter,aleafv,albv,albvobs 
        adiff=albvobs-albv
        if (abs(adiff).gt.errmax) then
          delta=albv-alast
          if(iter.eq.1) then
            deltr=-0.02*abs(adiff)/adiff 	! change by 0.02 with - sign of adiff
          else if(iter.lt.20)then
            dadr=delta/deltr
            dr=adiff/dadr
            deltr=0.75*dr			! nudge by 0.75*DR
          else
            aleafv=BAD
            done=.TRUE.
          endif
          if (.not.done) then
            aleafv=aleafv+deltr
            alast=albv
            iter=iter+1
          endif
        else
          done=.TRUE.
        endif 
      enddo
     
      done=.FALSE.
      iter=1
      alast=-9999
      do while (.not.done)
        call getradprops (zen2,albedo,taubtv,taubtn,clumps2)
        write(6,100)'ALEAFN: ',iter,aleafn,albn,albnobs          
        adiff=albnobs-albn
        if (abs(adiff).gt.errmax) then
          delta=albn-alast
          if(iter.eq.1) then
            deltr=-0.02*abs(adiff)/adiff 	! change by 0.02 with - sign of adiff
          else if(iter.lt.20)then
            dadr=delta/deltr
            dr=adiff/dadr
            deltr=0.75*dr			! nudge by 0.75*DR
          else
            aleafn=BAD
            done=.TRUE.
          endif
          if (.not.done) then
            aleafn=aleafn+deltr
            alast=albn
            iter=iter+1
          endif
        else
          done=.TRUE.
        endif
       enddo
      
 100  format(a10,i5,3f11.5)
      
      return
      end
