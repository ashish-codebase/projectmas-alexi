!------------------------------------------------------------------
!           INCLUDE FILE FOR US FLUX PROJECT PARAMETERS
!------------------------------------------------------------------
!  Parameters describing grid
      parameter(kz=30,ilg=1456,jlg=625)            ! MENAE domain
      parameter(kx=1440,ky=600,kt=8)               ! GETD Input Grid
      parameter(sx=1440,sy=720)                   ! clearsky INSOL GRID
!  Maximum number of hours of data in which the solar radiation is
!  positive in CRAS model ("nohrc"), and how many are required for ALEXI
!  time series inputs ("nohr").
      parameter(nohrc =24)
      parameter(nohr  =24)
         
!  ICLASS values (SIB2 landcover classes)
      parameter(nclass=8)        
      parameter(ICE   =1)
      parameter(WATER =1)
      parameter(RIPARIAN=1)      
      
!  ICLASM values (MODIS lookup table)
      parameter(nclasm=6)
      parameter(nbin=20)	! Number of NDVI bins  
      
!  NCAR graphics parameters
      parameter(MXPL=50, MNCLR=43, MXCLR=250)  !NCAR Graphics settings
      parameter(MXPLTS=28)                !number of desired input fields.
      parameter(NOUTGRAPHICS= 8)          !number of desired output fields.
      

