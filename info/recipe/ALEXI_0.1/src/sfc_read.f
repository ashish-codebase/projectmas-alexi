!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Routine to read in gridded surface data.
!  
!  -----------------------------------------------------------------------

      subroutine sfc_read(MDATE)

c      include "date.inc"
      include "USflux_grids.inc"
      include "USflux_dir.inc"
      include "flag.inc"

      common/hrdata_met/ctloc(kx,ky,kt),cta(kx,ky,kt),
     &       cea(kx,ky,kt),cwind(kx,ky,kt),csdn(kx,ky,kt),
     &       cxlwdn(kx,ky,kt),cpres(kx,ky,kt),
     &       clst(kx,ky,kt), clapse(ilg,jlg)
      common/lookup/lookup_i(ilg,jlg),lookup_j(ilg,jlg) 
      common/slookup/insol_i(ilg,jlg), insol_j(ilg,jlg)
      real*4 temp(ilg,jlg)
      real*4 data42(ilg,jlg)
      character*256 dir, ldir
      character*256 cmd1, cmd2, cmd3
      character*256 cmd4, cmd5, cmd6, cmd7, cmd8, cmd9, cmd10
      character*256 omd1, omd2, omd3
      character*256 omd4, omd5, omd6, omd7, omd8, omd9, omd10
      character*256 cstr1, cstr2, cstr3, cstr4, cstr5
      character*256 ostr1, ostr2, ostr3, ostr4, ostr5
      character*256 f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11
      character*256 s1, s3, s5, s7, s9
      character*11 cmd
      character*9 zip
      character*19 fp1, fp2
      character*16 fp3
      integer IDATE, iyear, doy
      integer LDATE, lyear
      l1=index(INDIR,' ')-1  

! Read in NARR data

      cmd='gunzip -ff '
      zip='gzip -ff '
      dir='/work/waterforfood/bucricket/PROCESS_VIIRS/STATIC/CFSR'

      iyear=MDATE/1000-mod(MDATE,1000)/1000
      write(cstr1,'(I4,a,I7,a)') iyear,'/',MDATE,'/t2_series.bin.gz'
      write(cstr2,'(I4,a,I7,a)') iyear,'/',MDATE,'/q2_series.bin.gz'
      write(cstr3,'(I4,a,I7,a)') iyear,'/',MDATE,'/psfc_series.bin.gz'
      write(cstr4,'(I4,a,I7,a)') iyear,'/',MDATE,'/wind_surface.bin.gz'
      write(cstr5,'(I4,a,I7,a)') iyear,'/',MDATE,'/lwdn.bin.gz'
      write(cmd1,'(a,a,a)') cmd,trim(dir),trim(cstr1)
      write(cmd2,'(a,a,a)') cmd,trim(dir),trim(cstr2)
      write(cmd3,'(a,a,a)') cmd,trim(dir),trim(cstr3)
      write(cmd4,'(a,a,a)') cmd,trim(dir),trim(cstr4)
      write(cmd5,'(a,a,a)') cmd,trim(dir),trim(cstr5)

      write(ostr1,'(I4,a,I7,a)') iyear,'/',MDATE,'/t2_series.bin'
      write(ostr2,'(I4,a,I7,a)') iyear,'/',MDATE,'/q2_series.bin'
      write(ostr3,'(I4,a,I7,a)') iyear,'/',MDATE,'/psfc_series.bin'
      write(ostr4,'(I4,a,I7,a)') iyear,'/',MDATE,'/wind_surface.bin'
      write(ostr5,'(I4,a,I7,a)') iyear,'/',MDATE,'/lwdn.bin'
      write(omd1,'(a,a,a)') zip,trim(dir),trim(ostr1)
      write(omd2,'(a,a,a)') zip,trim(dir),trim(ostr2)
      write(omd3,'(a,a,a)') zip,trim(dir),trim(ostr3)
      write(omd4,'(a,a,a)') zip,trim(dir),trim(ostr4)
      write(omd5,'(a,a,a)') zip,trim(dir),trim(ostr5)

      call system(cmd1)
      call system(cmd2)
      call system(cmd3)
      call system(cmd4)
      call system(cmd5)


      s1='/t2_series.bin'
      s3='/q2_series.bin'
      s5='/psfc_series.bin'
      s7='/wind_surface.bin'
      s9='/lwdn.bin'
      write(f1,'(a,I4,a,I7,a)') trim(dir),iyear,'/',MDATE,trim(s1)
      write(f3,'(a,I4,a,I7,a)') trim(dir),iyear,'/',MDATE,trim(s3)
      write(f5,'(a,I4,a,I7,a)') trim(dir),iyear,'/',MDATE,trim(s5)
      write(f7,'(a,I4,a,I7,a)') trim(dir),iyear,'/',MDATE,trim(s7)
      write(f9,'(a,I4,a,I7,a)') trim(dir),iyear,'/',MDATE,trim(s9)
      f11='/data/data123/chain/4KM/GBIM/inputs/CLASS.dat'
      open(20,file=f1,form='unformatted',
     1 access='direct',recl=kx*ky*kt*4)
      open(22,file=f3,form='unformatted',
     1 access='direct',recl=kx*ky*kt*4)
      open(24,file=f5,form='unformatted',
     1 access='direct',recl=kx*ky*kt*4)
      open(26,file=f7,form='unformatted',
     1 access='direct',recl=kx*ky*kt*4)
      open(28,file=f9,form='unformatted',
     1 access='direct',recl=kx*ky*kt*4)
      open(30,file=f11,form='unformatted',
     1 access='direct',recl=ilg*jlg*4)

      read(20,rec=1) cta
      read(22,rec=1) cea
      read(24,rec=1) cpres
      read(26,rec=1) cwind
      read(28,rec=1) cxlwdn
      read(30,rec=1) data42 

      call system(omd1)
      call system(omd2)
      call system(omd3)
      call system(omd4)
      call system(omd5)

      end subroutine
c
c  John Mecikalski
c  Created:  96.04.30
c
c  Routine to setup domain lat/long arrays.
c
      subroutine arraynav
      include 'USflux_grids.inc' 
      CLAT = 37.30 
      CLON = -95.90
      dlat=0.04 !DY/DL                  ! increment of latitude per grid box
      hid=float(ilg/2)            ! half array dimension
      hjd=float(jlg/2)
      MINLAT=CLAT-(hjd*dlat)      ! southern latitude
      MAXLAT=CLAT+(hjd*dlat)      ! northern latitude
      MINLON=CLON-(hid*(dlat*cos(dlat*PI180)))    ! western longitude
      MAXLON=CLON+(hid*(dlat*cos(dlat*PI180)))    ! eastern longitude
      do ja=1,jlg
        navlat(ja)=MINLAT+((ja-1)*dlat)
        do ia=1,ilg
          navlon(ia)=MINLON+((ia-1)*dlat)
        enddo
      enddo
c
c      write(6,'(/'' Domain Boundaries:        '')')
c      write(6,'(''        northern latitude: '',f10.5)')MAXLAT
c      write(6,'(''        southern latitude: '',f10.5)')MINLAT
c      write(6,'(''        eastern longitude: '',f10.5)')MAXLON
c      write(6,'(''        western longitude: '',f10.5)')MINLON
c      write(6,'(/'' Analysis Grid Dimensions: '')')
c      write(6,'(''        east-west   (ilg): '',i4)')ilg
c      write(6,'(''        north-south (jlg): '',i4)')jlg
c      write(6,'(''        cellsize    (deg): '',f7.5)')dlat
c
      return
      end
      

c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  + END  END  END  END  END  END  END  END  END  END  END  END  END  END +
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
