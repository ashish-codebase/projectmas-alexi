C
c  John Mecikalski
c  Created:  04/25/96
c  Modifed:  08/23/01  For the US domain
c
c  =====================================================================
c  This program drives the processing steps required for execution of
c  the two-source, time-integrated boundary layer/canopy flux model.
c  We use the version of this model which requires knowledge of light-
c  use efficiency.
c
c  The routine's main purpose will be to collect all needed data for
c  use in the htsolve8.f routine.  This required reading the EXPACK file
c  PREP0005 as created in the preprocessing step (prepdrvr program).
c
c  Details of these data inputs are found in file "prepdrvr.f".
c
c
c  --------------------------------------------
c  The final products of the procedure include:
c
c    A) Total R(net), G, H, E and CO2 fluxes.
c    B) Soil R(net), G, H, and E.
c    C) Canopy delta R(net), H, E and CO2 fluxes.
c
c  =====================================================================
c
      program main 

c
      include "date.inc"
      include "USflux_grids.inc"
      include "USflux_dir.inc"
      common/timestamp/year,doy
      integer MONTH,DAY,YEARI
      character*2 cmonth,cday,trashc
      character*4 cyear
      character*9 months(12),datein
      character*17 datein2 
      character*256 datefile
      character*7 arg1
      character*10 arg2
      character*3 arg3
      character*16 t
      character*3 part

      data (months(k),k=1,12)/'January','February','March','April',
     +          'May','June','July','August','September','October',
     +          'November','December'/
c
c  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
c  See prepdrvr.f for a description of the block and arrays and constants.
c  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
c
c  Introductory text for niceness
      write(6,'(''################################################'')')
      write(6,'(''####                                        ####'')')
      write(6,'(''####   Begin Processing ...                 ####'')')
      write(6,'(''###    Atmospheric-Land Exchange Model -     ###'')')
      write(6,'(''##       Inverse (A L E X I) USA Version      ##'')')
      write(6,'(''###        [Anderson et al. 1997 &           ###'')')
      write(6,'(''####        Mecikalski et al. 1999]         ####'')')
      write(6,'(''####                                        ####'')')
      write(6,'(''################################################'')')
c
c      l1=index(INDIR,' ')-1
c      datefile=INDIR(1:l1)//'DATELIST'
c      open(19,file=datefile,status='old')
c      read(19,1009)MDATE,datein,CDATE
c 1009 format(i7,1x,a9,1x,a8)
c      YDATE=datein(1:2)//" "//datein(3:5)//" "//datein(6:9)
c      read(CDATE,'(i8)')JDATE
c      write(6,'(/''Beginning processing for: '',a11/)')YDATE      
c      write(6,*)'CDATE = ',CDATE
c      write(6,*)'MDATE = ',MDATE
cCC    MDATE=1999187              !TESTING
c      doy=mod(MDATE,1000)
c      year=MDATE/1000-mod(MDATE,1000)/1000
c      print*,year
c      close(19)
c
c
c  Perform processing by calling ALEXI; newst ALEXI model (09/2001)
c  ----------------------------------------------------------------------
c      MDATE=2004213
      call getarg(1,arg1)
      read(arg1,'(i7)') MDATE
      call getarg(2,t)
      call getarg(3,arg2)
      read(arg2,'(i10)') npoints
!     call getarg(4,arg3)
!     read(arg3,'(i3)') part
      call getarg(4,part)

      call landcover(MDATE,t,npoints, part)
      call USflux(MDATE,t,npoints, part)
c
      write(6,'(/''--->>  Processing is completed.'')')
      stop
      end
c
c  END ...
c
