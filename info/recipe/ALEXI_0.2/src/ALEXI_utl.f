c-------------------------------------------------------------------c
c              A L E X I   U T I L I T Y   R O U T I N E S          c
c-------------------------------------------------------------------c

      subroutine set_constants
c     **************************************************************            
c     *
c     *  Set parameters held constant over whole model domain. 
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      include 'ALEXI_parm.inc'
      common/cover2/perennial,iswater,fcbare
      common/constants/pi,cp,xk
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1,
     &       h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1,
     &       sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1,
     &       albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2                
      common/initial/hn0,psi0,fc0
      common/partition_clear/difvisclr,difnirclr,dirvisclr,
     &       dirnirclr,fvisclr,fnirclr  
      common/pbl/thrise,zpbli(mli),thpbli(mli),zpbl(ml),thpbl(ml),
     &       jzmax,jz1,nlev
      logical perennial,iswater
      save                     
      
c  CONSTANTS:
      xk=0.4			! von Karman's constant
      pi=3.1415926537
      thrise=3.0*3600.0		! Flux rise time (Tennekes '73 suggests 3 hrs) [s]
      cp=1010.0			! Specific heat of air at const. pres. [J/kg-K]		
      z1=50.0			! Initial boundary layer height [m]
      fcbare=0.00		! fc.le.fcbare considered to be bare soil

c CLEAR-SKY RADIATION PARTITIONING:
      fvisclr=0.5
      fnirclr=0.5
      dirvisclr=0.8
      difvisclr=0.2
      dirnirclr=1.0
      difnirclr=0.0

c INITIALIZATIONS:
      hn0=80.       
      psi0=0.    

      return
      end
      
      subroutine getprofile(z1)
c     ***************************************************************
c     *
c     *  Interpolate the input potential temperature profile to a finer
c     *  vertical grid.  Find jz1, array index corresponding to height
c     *  z1.  For example:
c     *
c     *     Input arrays        Interpolated arrays
c     *     ---------------------------------------
c     *      THPBLI(41)	 -------->  THPBL(8000)
c     *      ZPBLI(41)	 -------->  ZPBL(8000)
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************     
      include 'ALEXI_parm.inc'
             
      common/pbl/thrise,zpbli(mli),thpbli(mli),zpbl(ml),thpbl(ml),
     &           jzmax,jz1,nlev
      save
c
c		Interpolate the profile to a fine-scale 
c			(1m spacings) array
c
      j0=1
      mxhtpbl=zpbli(nlev)
      zpbl(1)=zpbli(j0)
      thpbl(1)=thpbli(j0)
      do j=2,ml
         ht=zpbli(1)+j
         if (ht.gt.mxhtpbl) then
            jzmax=j-1
            go to 2000
         endif
         do while (ht.gt.zpbli(j0+1))
            j0=j0+1
         end do        
         f=(ht-zpbli(j0))/(zpbli(j0+1)-zpbli(j0))
         temp=thpbli(j0)+f*(thpbli(j0+1)-thpbli(j0))
         zpbl(j)=ht
         thpbl(j)=temp
         if(ht.eq.z1)jz1=j
      enddo	! j loop
      jzmax=j-1

 2000 continue
      return
      end
      
      subroutine getpbltable(z1,ierr)
c     **************************************************************            
c     *
c     *  Computes HNNEW associated with each possible TA2-TA1
c     *  and stores in lookup table HNTAB(I). 
c     *
c     *  Martha Anderson
c     *  Created:  04/02/02
c     *  
c     ***************************************************************
      include 'ALEXI_parm.inc'
      common/pbl/thrise,zpbli(mli),thpbli(mli),zpbl(ml),thpbl(ml),
     &       jzmax,jz1,nlev
      common/pbl2/tabtheta(mt),tabz2(mt),dtheta,itmax,thpblz1
      save
      
      dtheta=0.01		! DTa increment in table
c  Find potential temperature at height Z1
      j=1
      do while (z1.gt.zpbli(j+1))
        j=j+1
      enddo
      
      f=(z1-zpbli(j))/(zpbli(j+1)-zpbli(j))
      thpblz1=thpbli(j)+f*(thpbli(j+1)-thpbli(j))
      tlast=thpblz1
      zlast=z1 
c  Integrate pot. temp up from height Z1 to  Z in increments
c  of DTHETA.  Also find Z2 corresponding to each possible change
c  in temperature I*DTHETA
      tint=0.
      do it=1,mt
        t=tlast+dtheta		! DTAtot = I*DTHETA
        do while (t.gt.thpbli(j+1))
          j=j+1
          if (j.ge.mli) then
            itmax=it-1		! This is the biggest DTa we can 
            go to 100		! accommodate (itmax*dtheta)
          endif
        enddo
        f=(t-thpbli(j))/(thpbli(j+1)-thpbli(j))
        z=zpbli(j)+f*(zpbli(j+1)-zpbli(j))
        dz=z-zlast
        tint=tint+0.5*dz*(t+tlast)	! Trapezoid rule integration
        tabtheta(it)=tint
        tabz2(it)=z
        tlast=t
        zlast=z
      enddo  
      itmax=mt  

100   continue
      return
      end
            
      subroutine runinit(badinput,ibad,writeme)
c     **************************************************************            
c     *
c     *  Initializes variables that may vary from pixel to pixel
c     *  across domain.  Also perform necessary data corrrections. 
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     ***************************************************************
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/cover2/perennial,iswater,fcbare
      common/constants/pi,cp,xk      
      common/cover/iclass,beta,gvs,xndvi,fpar,dstom  
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1,
     &       h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1,
     &       sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1,
     &       albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2       
      common/emission/emleaf,emdead,emsoil,emcpy 
      common/emission2/esfc,aem,bem,eleaf,esoil            
      common/initial/hn0,psi0,fc0        
      common/model/zta
      common/partition/difvis,difnir,dirvis,dirnir,fvis,fnir      
      common/partition_clear/difvisclr,difnirclr,dirvisclr,
     &       dirnirclr,fvisclr,fnirclr         
      common/site/xlat,xlong,stdlng
      common/sun/strt,end        
      common/timestamp/year,doy       
      logical badinput,writeme,perennial,iswater    
      save 

c  Check validity of input before continuing with input computations 
         
      call checkinput(badinput,ibad,writeme)
      if(badinput)then
         go to 2000 
      endif 
      
c  Set partitioning factors to clear-sky values:

      fvis=fvisclr
      fnir=fnirclr
      dirvis=dirvisclr
      difvis=difvisclr
      dirnir=dirnirclr
      difnir=difnirclr
                        
c  Compute solar zenith angles at times T1 and T2

      call getsunzen (xlat,xlong,stdlng,doy,year,tloc1,zen1)
      call getsunzen (xlat,xlong,stdlng,doy,year,tloc2,zen2)

      coszen1=cos(zen1)
      coszen2=cos(zen2)
c      if (fc.le.0.5) then
c        dtloc=2*xndvi+0.75
c      else
c        dtloc=1.75
c      endif
c      if (fc.lt.0.10) then
c       dtloc=0.75
c      endif
      if (fc.lt.0.1) then
        dtloc=0.75
      else if (0.1.le.fc.and.fc.le.0.7) then
!        dtloc=2.0*fc+0.55
        dtloc=1.6667*fc+0.5833
      else
        dtloc=1.75
      endif
      dtloc=1.50
!      write(6,*) "DTLOC: ",dtloc, xndvi
      t1=(tloc1-strt-dtloc)*3600.  ! Referenced to DTLOC hours past sunrise [s]
      t2=(tloc2-strt-dtloc)*3600.  ! Referenced to DTLOC hours past sunrise [s]
!      t1=(tloc1-strt-1.75)*3600.  ! Referenced to DTLOC hours past sunrise [s]
!      t2=(tloc2-strt-1.75)*3600. 


            
c  Compute coeffs for directional thermal surface emissivity: 
c  (parameterized in terms of eleaf, esoil, fc and theta based on Cupid 
c  simulations - works well for theta < 70deg).
c  esfc=aem*fc*fc+bem*fc+esoil is computed in updatefc.
c  TEMPORARY: Replace with analytic expression for esfc taking
c             e.g. fg into account.

      eleaf=0.97		! Leaf emissivity
c      esoil=0.94		! Soil emissivity - USE EMSOIL AT GRID CELL     
      b1=-1.97*cos(theta)+2.87
      b2=1.86*cos(theta)-2.62
      bem=b1*eleaf+b2
      aem=eleaf+0.025-bem-emsoil   
      call updatefc
      
c  Correct radiometric temperatures for emissivity and sky emission:  
c  (Comment out - taken care of in atmos. correction routine...)
c
c      tsky1=(xlwdn1/5.67e-8)**.25-273.15
c      trad1=(((trad1+273.15)**4.
c     &     -(1.-esfc)*(tsky1+273.15)**4.)/esfc)**.25-273.15 
c      tsky2=(xlwdn2/5.67e-8)**.25-273.15
c      trad2=(((trad2+273.15)**4.
c     &     -(1.-esfc)*(tsky2+273.15)**4.)/esfc)**.25-273.15 
 
c  Scale wind measurements from REFHTW (height of measurement)
c  up to height ZTA (height of TA node) over grass.

      z0g=0.005                 ! Roughness length for grass [m]
      dispg=0.0                 ! Disp. height for grass [m]
      factw1=(alog(zta-dispg)-alog(z0g))/
     &                            (alog(refhtw-dispg)-alog(z0g))      
      w1=w1*factw1
      w2=w2*factw1 
      refhtw=zta 
      
      wmin=3.0
      wmax=20.0
      if(w1.lt.wmin)w1=wmin
      if(w2.lt.wmin)w2=wmin
      if(w1.gt.wmax)w1=wmax
      if(w2.gt.wmax)w2=wmax

c  Compute volumetric heat capacity of air [J/deg-m3]

      rho1=			! [kg/m3]
     &   pres1/(287.04*(taobs1+273.15))*(1.-.378*ea1/pres1)*100.
      rho2=			! [kg/m3]
     &   pres2/(287.04*(taobs2+273.15))*(1.-.378*ea2/pres2)*100.
      rhocp1=rho1*cp 		! [J/deg-m3] 
      rhocp2=rho2*cp  		! [J/deg-m3] 
       
c  Interpolate PBL profile to 1m levels

c      call getprofile(z1)   
      
c  Create lookup tables TABTHETA and TABZ2

      call getpbltable(z1,ierr)      

2000  continue
      return
      end

      subroutine updatefc
c     **************************************************************            
c     *
c     *  Computes factors and parameters involving fraction cover.  
c     *  Needs to be called whenever xlai or fc are changed.
c     *  Canopy architecture factors are updated in UPDATECANOPY.
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/clumping/clumps1,clumps2,clump0,fveg
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1,
     &       h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1,
     &       sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1,
     &       albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2               
      common/emission/emleaf,emdead,emsoil,emcpy 
      common/emission2/esfc,aem,bem,eleaf,esoil
      save

      call canopyarch
      
c	Thermal surface emissivity      
      esfc=aem*fc*fc+bem*fc+emsoil 
      
c	Radiation properties of system components
      call getradprops(zen1,albedo1,taubtv1,taubtn1,clumps1)
      call getradprops(zen2,albedo2,taubtv2,taubtn2,clumps2)
      
      return
      end
                        
      subroutine checkinput(flag,ibad,w)
c     ***************************************************************            
c     *
c     *  All input variables are checked for unreasonable values.
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      include 'ALEXI_parm.inc'
      common/absorption/aleafn,aleafv,aleafl,adeadv,adeadn,adeadl      
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
      common/reflection/rsoiln,rsoilv,rdcpyl
      common/site/xlat,xlong,stdlng  
      common/timestamp/year,doy 
      common/view/theta,ftheta      
      save
      
      logical flag,w    

      flag=.FALSE.
c
c  PBL parameters:
      do ka=1,nlev
c         call checkvalue("THPBLI",1,thpbli(ka),220.,450.,flag,w,ibad)	! K
c         call checkvalue("ZPBLI ",2,zpbli(ka),0.,10000.,flag,w,ibad)	! m
      enddo
c
c  Time Specific parameters:
      call checkvalue("EA1   ",3, ea1,0.,80.,flag,w,ibad)	! mbar
      call checkvalue("EA2   ",4, ea2,0.,80.,flag,w,ibad)
      call checkvalue("TAOBS1",5, taobs1,-50.,60.,flag,w,ibad)	! C
      call checkvalue("TAOBS2",6, taobs2,-50.,60.,flag,w,ibad)
      call checkvalue("TRAD1 ",7, trad1,-50.,60.,flag,w,ibad)	! C
      call checkvalue("TRAD2 ",8, trad2,-50.,60.,flag,w,ibad)
      call checkvalue("SDN1  ",9, sdn1,0.,1500.,flag,w,ibad)	! W/m-2
      call checkvalue("SDN2  ",10,sdn2,0.,1500.,flag,w,ibad)
      call checkvalue("PRES1 ",11,pres1,200.,1200.,flag,w,ibad)	! mbar
      call checkvalue("PRES2 ",12,pres2,200.,1200.,flag,w,ibad)
      call checkvalue("W1    ",13,w1,0.,100.,flag,w,ibad)	! m/s
      call checkvalue("W2    ",14,w2,0.,100.,flag,w,ibad)
      call checkvalue("XLWDN1",15,xlwdn1,0.,1000.,flag,w,ibad)	! W/m-2
      call checkvalue("XLWDN2",16,xlwdn2,0.,1000.,flag,w,ibad)
c	THIS IS TOO GENEROUS!  Need to vary t1, t2 with latitude....      
      call checkvalue("ZEN1  ",17,zen1,0.,3.14,flag,w,ibad)	! radians
      call checkvalue("ZEN2  ",18,zen2,0.,3.14,flag,w,ibad)
      call checkvalue("YEAR  ",19,year,1900.,2100.,flag,w,ibad) 
      call checkvalue("DOY   ",20,doy,0.,366.,flag,w,ibad)            
      call checkvalue("THETA ",21,theta,0.,80.,flag,w,ibad)	! degrees
c
c  Canopy parameters:
      call checkvalue("XLAI  ",22,xlai,0.,10.,flag,w,ibad)
      call checkvalue("FG    ",23,fg,0.,1.,flag,w,ibad)
      call checkvalue("HEIGHT",25,height,0.,35.,flag,w,ibad)	! m
      call checkvalue("CLUMP ",26,clump,0.,1.,flag,w,ibad)		
      call checkvalue("XL    ",27,xl,0.0,0.10,flag,w,ibad)	! m
c
c  Radiometric parameters:
      call checkvalue("ALEAFV",28,aleafv,0.,1.,flag,w,ibad)
      call checkvalue("ALEAFN",29,aleafn,0.,1.,flag,w,ibad)
      call checkvalue("ALEAFL",30,aleafl,0.,1.,flag,w,ibad)
      call checkvalue("RSOILV",31,rsoilv,0.,1.,flag,w,ibad)
      call checkvalue("RSOILN",32,rsoiln,0.,1.,flag,w,ibad)
      call checkvalue("EMSOIL",33,emsoil,0.,1.,flag,w,ibad)      
c
c  Checks on radiometric temperature
      dtrad=trad2-trad1
      dtair=taobs2-taobs1
      call checkvalue("DTRAD ",34,dtrad,0.,50.,flag,w,ibad)
c	DTRAD should be >~ DTAIR
      xmin=dtair-2.      
c      call checkvalue("DTRDTA",99,dtrad,xmin,50.,flag,w,ibad)
c	TRAD should not be too far below TAIR 
      xmin=taobs1-39.0     
      call checkvalue("TR:TA1",35,trad1,xmin,60.,flag,w,ibad)
      xmin=taobs2-29.0
      call checkvalue("TR:TA2",36,trad2,xmin,60.,flag,w,ibad)      
          
      return
      end          
c
c
      subroutine checkvalue(valname,itag,value,xmin,xmax,badinput,
     &                      writeme,ibad)
c     ***************************************************************            
c     *
c     *  The value of variable "valname" is scrutinized.  If it falls
c     *  outside of the prescribed range, an error message is printed
c     *  to the screen and the program pauses.
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     ***************************************************************            
      character*6 valname
      logical badinput,writeme
      
      if ((value.lt.xmin) .or. (value.gt.xmax)) then
         if(writeme)write(6,5000)valname,value,xmin,xmax
!         write(6,5000)valname,value,xmin,xmax
	 badinput=.TRUE.
         ibad=itag
      endif
 5000 format(a6,' = ',f9.2,' outside range ',f9.2,' to ',f9.2) 
      return
      end         
      
      
      subroutine ALEXI_errorcode(writeme,ierr)
c     **************************************************************            
c     *
c     *  
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     ***************************************************************
      logical writeme
      
      if(writeme)then
         if    (ierr.eq.0)then
            write(6,*)'*** CONVERGED ***'
         elseif(ierr.eq.1)then
            write(6,*)'*** HN did not converge. Bail ***'
         elseif(ierr.eq.2)then
            write(6,*)'*** Bad input.  Bail ***'
         elseif(ierr.eq.3)then
            write(6,*)'*** Exceeded max PBL profile layer. Bail ***'
         elseif(ierr.eq.4)then
            write(6,*)
         elseif(ierr.eq.5)then
            write(6,*)
         elseif(ierr.eq.6)then
            write(6,*)
         elseif(ierr.eq.7)then
            write(6,*)
         elseif(ierr.eq.8)then
            write(6,*)
         elseif(ierr.eq.9)then
c            write(6,*)'*** RNET < 0. Bail ***'
         elseif(ierr.eq.10)then
            write(6,*)
         elseif(ierr.eq.11)then
c            write(6,*)'*** TS < TC for all FC. Bail ***'  
         elseif(ierr.eq.12)then
c            write(6,*)'*** G/RNSOIL < 0.1. Bail ***'   
         elseif(ierr.eq.13)then
c            write(6,*)'*** TS-TRAD > 15. Bail ***'                      
         endif 
c         write(6,*)' '
      endif
         
      return
      end
