      subroutine getradcomps (sdn,zen,fclear)
c     ***************************************************************
c     *
c     *  Partitions SDN into VIS and NIR, diffuse and direct beam
c     *  components, following strategy in Cupid's RADIN4 subroutine.
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************
      common/constants/pi,cp,xk
      common/partition/difvis,difnir,dirvis,dirnir,fvis,fnir 
           
      pid2=pi/2.
      pid180=pi/180.
      ration=1.0	! Assume it is clear at night in absence
      			! of other information
      			           
      coszen=cos(zen)
      if(coszen.lt.0.01)goto 900	! Nighttime...
      
c		Calculate potential (clear-sky) visible and
c			NIR solar components
     
c  	Correct for curvature of atmos in airmas
      airmas=(sqrt(coszen**2+.0025)-coszen)/.00125      
c  	Correct for refraction(good to 89.5 deg.)
      airmas=airmas-2.8/(90.-zen/pid180)**2
      
      potbm1=600.*exp(-.160*airmas)
      potvis=(potbm1+(600.-potbm1)*.4)*coszen
      potdif=(600.-potbm1)*.4*coszen
      u=1.0/coszen
      axlog=alog10(u)
      a=10**(-1.195+.4459*axlog-.0345*axlog*axlog)
      watabs=1320.*a
      potbm2=720.*exp(-.05*airmas)-watabs
      if(potbm2.lt.0.)potbm2=0.
      eval=(720.-potbm2-watabs)*.54*coszen
      potnir=eval+potbm2*coszen      
      fclear=sdn/(potvis+potnir)
      fclear=min(1.0,fclear)
c
c		Partition SDN into VIS and NIR      
      
      fvis=potvis/(potvis+potnir)
      fnir=potnir/(potvis+potnir)
c
c		Estimate direct beam and diffuse fractions
c			in VIS and NIR wavebands
c      
      fb1=potbm1*coszen/potvis
      fb2=potbm2*coszen/potnir
      ratiox=fclear
      if(fclear.gt..9)ratiox=.9
      dirvis=fb1*(1.-((.9-ratiox)/.7)**.6667)
      if(fclear.gt.0.88)ratiox=.88
      dirnir=fb2*(1.-((.88-ratiox)/.68)**.6667)
      dirvis=max(0.0,dirvis)
      dirnir=max(0.0,dirnir)
      dirvis=min(fb1,dirvis)
      dirnir=min(fb2,dirnir)
      if(dirvis.lt..01.and.dirnir.gt..01)dirvis=.011
      if(dirnir.lt..01.and.dirvis.gt..01)dirnir=.011
      difvis=1.0-dirvis
      difnir=1.0-dirnir
      goto 1000
            
c		If it is nighttime... 
     
 900  continue
      fclear=ration
      fvis=0.5
      fnir=0.5
      difvis=1.0
      difnir=1.0
      dirvis=0.0
      dirnir=0.0
      if(sdn.le.1.0)sdn=0.0
      if(sdn.ne.0.0)then
!        write(6,*)"The sun is shining at night! ",
!     &                  sdn,zen,coszen
        sdn=0.0
      endif
      
1000  continue

      return
      end
      
      subroutine getnetrad_simple (sdn,ta,fclear,zen,albedo,
     &                             rnet)
c     **************************************************************            
c     *
c     *  Calculates net radiation using solar radiation and air temperature
c     *  estimates only.
c     *
c     *  Martha Anderson
c     *  Created:  09/04/01
c     *  
c     *************************************************************** 
      parameter(sigma=5.67e-8) 
          
      tak=ta+273.15
      eskyc=9.2e-6*tak*tak			! Clear: Swinbank '63 (C&N pg 163)
      esky=(1-0.84*fclear)*eskyc + 0.84*fclear	! Cloudy: Monteith & Unsworth '90 (C&N pg 164)
      esfc=0.98
      
      rnetl=esfc*(esky-1.0)*sigma*(tak**4)
      rnetl=esfc*esky*sigma*(tak**4)-esfc*sigma*((tak+4)**4)
      rnets=sdn*(1-albedo)
      rnet=rnetl+rnets
      
      return
      end     
            

