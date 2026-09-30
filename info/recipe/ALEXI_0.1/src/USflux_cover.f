      subroutine load_tables(MDATE)
c     **************************************************************            
c     *
c     *  Loads tables containing parameters associated with 
c     *  specific landcover classes. 
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      include 'USflux_parm.inc'
      include 'USflux_dir.inc'
      common/lookup/itabclass(nclass),tablai(nbin,nclasm),
     &       tabfpar(nbin,nclasm),tabbeta(nclass),
     &       tabaleaf(nclass,3),tabadead(nclass,3),
     &       tabhmin(nclass),tabhmax(nclass),tabxl(nclass),
     &       tabgvs(nclass),tabdstom(nclass),tabperen(nclass),
     &       tabfcmin(nclass),tabdesc(nclass)
      common/timestamp/year,doy
      character*50  tabdesc 
      character*256 dir,infile, vegfile   
      logical       present 
      save
      
c  Table files are contained in directory TABDIR      
      l1=index(TABDIR,' ')-1  
      
c     Assign properties by LDAS landcover class:
c     -------------------------------------------------         
c     NCLASS: Classes 1-14 are TMS classes.
c     TABALEAF, TABADEAD: The absorptivities here are 1-(rleaf+tleaf) 
c       listed for LDAS UMD classes (merged).
c     TABBETA: LUE values
c     TABHMIN, TABHMAX, TABXL: Min/max height and leaf size (m)
c     ITABCLASS: MODIS biome class corresponding to LDAS ICLASS
c     -------------------------------------------------         
      
      infile=TABDIR(1:l1)//'landcover.txt'
      inquire(file=infile,exist=present)
      if (.not.present) then
         write(6,*)'landcover.txt file does not exist'
         stop 'stopping in load_tables'         
      endif 
                
      open(unit=50,file=infile) 
      read(50,*)
      do while (1.eq.1)
        read(50,*,end=1000)iclass,tabaleaf(iclass,1),
     &  tabaleaf(iclass,2),tabaleaf(iclass,3),tabadead(iclass,1), 
     &  tabadead(iclass,2),tabadead(iclass,3),tabbeta(iclass),
     &  tabhmin(iclass),tabhmax(iclass),tabxl(iclass),
     &  itabclass(iclass),tabperen(iclass),tabfcmin(iclass)
      enddo 
1000  continue
      close(50)
            
      tabdesc(1)= 'Water (and Goodes Interrupted Space)'
      tabdesc(2)= 'Evergreen Needleleaf Forest'
      tabdesc(3)= 'Evergreen Broadleaf Forest'
      tabdesc(4)= 'Deciduous Needleleaf Forest'
      tabdesc(5)= 'Deciduous Broadleaf Forest'
      tabdesc(6)= 'Mixed Cover'
      tabdesc(7)= 'Woodland'
      tabdesc(8)= 'Wooded Grassland'
      tabdesc(9)=' Closed Shrubland'
      tabdesc(10)='Open Shrubland'
      tabdesc(11)='Grassland'
      tabdesc(12)='Cropland'
      tabdesc(13)='Bare Ground'
      tabdesc(14)='Urban and Built-Up'

c	Open file containing class-weighted canopy characteristics
c	Currently, these data overwrite values computed from the table paramaters.
c	Really only using the class names from the table (TABDESC)

!      l2=index(INDIR,' ')-1
!      vegfile=INDIR(1:l2)//'veg_us'
!      open (300,file=vegfile)

!      dir='/data/data123/chain/4KM/GBIM/veg/veg_us'
!      write(vegfile,'(a,I7)') trim(dir),MDATE
!      open (300,file=vegfile)

     
      return
      end
      
      subroutine cover_props
c     **************************************************************            
c     *
c     *  Computes vegetation cover properties based on FC and tabular
c     *  data.  Returns ALEAF, HEIGHT, FG, FC, XL, and CLUMPING
c     *
c     *  TEMPORARY FIX:  These values are now being computed offline
c     *  by landcover.f and read in from veg_us.  landcover.f also
c     *  outputs Z0 and DISP, accounting for subpixel variability in
c     *  canopy type and properties.  Eventually should read these in
c     *  as gridded binary files, but memory limits are excedded. For
c     *  now, just readining in pixel by pixel so whole grids are not
c     *  stored. 
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *  
c     *************************************************************** 
      include 'USflux_parm.inc'
      common/absorption/aleafn,aleafv,aleafl,adeadv,adeadn,adeadl   
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump,rsmin
      common/clumping/clumps1,clumps2,clump0,fveg
      common/cover2/perennial,iswater,fcbare      
      common/cover/iclass,beta,gvs,xndvi,fpar,dstom
      common/lookup/itabclass(nclass),tablai(nbin,nclasm),
     &       tabfpar(nbin,nclasm),tabbeta(nclass),
     &       tabaleaf(nclass,3),tabadead(nclass,3),
     &       tabhmin(nclass),tabhmax(nclass),tabxl(nclass),
     &       tabgvs(nclass),tabdstom(nclass),tabperen(nclass),
     &       tabfcmin(nclass),tabdesc(nclass)
      common/timestamp/year,doy
      common/view/theta,ftheta
      character*50 tabdesc
      logical perennial,iswater
      save      

c	Extract leaf absorptivities from table based on USGS class
      aleafv=tabaleaf(iclass,1)
      aleafn=tabaleaf(iclass,2)
      aleafl=tabaleaf(iclass,3)
      adeadv=tabadead(iclass,1)
      adeadn=tabadead(iclass,2)
      adeadl=tabadead(iclass,3)

c 	Distinguish between perennial cover, where stalks or trunks 
c	remain after fc(ndvi) -> 0, and annual cover (crops) which
c	approach bare soil conditions at 0 cover.  This will influence
c	surface aerodynamic properties.
      if (tabperen(iclass).eq.1) then
         perennial=.TRUE.
      else
         perennial=.FALSE.
      endif 
      perennial=.FALSE.     ! THIS FUNCTIONALITY IS NOT ROBUSTLY IMPLEMENTED YET
           
c	Compute canopy architecture parameters:
      xl=tabxl(iclass)			! [m]            
      hmax=tabhmax(iclass)		! [m]
      hmin=tabhmin(iclass)		! [m] 
      height=hmin+(hmax-hmin)*fc	! [m]
      
      if (iswater) then
         fg=1.0
         fc=fcbare
      else if (.not.perennial) then	! For annuals...
         fg=1.0				! - for now, always fully green
         if (fc.lt.fcbare) fc=fcbare
         xlai=(-2.0*alog(1.-fc))
         clump0=1
      else				! For perennials...
         fg=1.0
         fcmin=tabfcmin(iclass)
         if (fc.lt.fcmin) then		! - past peak greenness: adjust fc,fg
            fg=fc/fcmin
            fc=fcmin  
         endif 
         if (fc.lt.fcbare) fc=fcbare
      endif  

c	Set clumping factors - Really, this should be read in as a gridded input.
      if (iclass.eq.12) then              ! CROPS
        clump=0.9                         ! At GOES view angle
        clump0=0.9                        ! At nadir
        clumps1=1.0                       ! At sun angle T1
        clumps2=0.9                       ! At sun angle T2
      else
        clump=1.0                         ! At GOES view angle
        clump0=1.0                        ! At nadir
        clumps1=1.0                       ! At sun angle T1
        clumps2=1.0                       ! At sun angle T2
      endif

      call getfveg(clump0,xlai,fveg)
!      write(6,*)'FVEG CLUMP0 CLUMP CLASS XLAI:',
!     &             fveg,clump0,clump,iclass,xlai

c	Overwrite some values with class-averaged values read from veg_us
      iin=BAD
      jin=BAD
!      do while (iin.ne.ia.or.jin.ne.ja) 
        read(400,*)iin,jin,aleafv,aleafn,aleafl,adeadv,adeadn,adeadl,
     &           height,xl,z0,disp,rsmin2
!      enddo
      
      call canopyarch			! Fill z0, disp in case ALEXI not run       
       					! (for plotting purposes)
      return
      end  
 
      subroutine getfveg(clump0,xlai,fveg)

      gap=exp(-0.5*clump0*xlai)
      xmin=9999
      
      do i=1,100
        fv=i/100.
        rhs=fv*exp(-0.5*xlai/fv)+(1.-fv)
        diff=abs(gap-rhs)
        if (diff.lt.xmin) then
          fvmin=fv
          xmin=diff
        endif
      enddo
      fveg=fvmin
      fveg=min(fveg,1.0)
      fveg=max(fveg,0.1)

      return
      end
