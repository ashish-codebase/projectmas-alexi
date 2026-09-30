c
c  Include file for ALEXI input data fields, read into 2D arrays
c
c------------------------------------------------------------------
c  MODEL PARAMETERS
c------------------------------------------------------------------
      include "ALEXI_parm.inc"
      include "USflux_parm.inc"

c------------------------------------------------------------------
c  MODEL COMMON BLOCKS
c------------------------------------------------------------------
      real*4 MINLAT,MAXLAT,MINLON,MAXLON
      real*4 CLAT,CLON
      real*4 tsr15
     2      ,tsr55
     3      ,vsr15
     4      ,vsr55
     5      ,wsr15
     6      ,wsr55
     7      ,psr15
     8      ,psr55
     9      ,rise15
     +      ,rise55
     1      ,radini
     2      ,radfin
     3      ,raddts
     4      ,tc15
     5      ,tc55
     +      ,ztan
     1      ,zran
     3      ,zzan
     4      ,htht
     5      ,pran
     6      ,lscls
     9      ,satang
     +      ,swbeg
     1      ,swend
     2      ,xlwbeg
     3      ,xlwend
     4      ,navlat
     5      ,navlon
     7      ,mlai

c
c  ALEXI input and core data arrays
      common /canopy1/tsr15(ilg,jlg)
     2               ,tsr55(ilg,jlg)
     3               ,vsr15(ilg,jlg)
     4               ,vsr55(ilg,jlg)
     5               ,wsr15(ilg,jlg)
     6               ,wsr55(ilg,jlg)
     7               ,psr15(ilg,jlg)
     8               ,psr55(ilg,jlg)
     9               ,radini(ilg,jlg)
     +               ,radfin(ilg,jlg)
     8               ,ztan(kz,ilg,jlg)
     1               ,zzan(kz,ilg,jlg)
     2               ,htht(mli,ilg,jlg)
     4               ,pran(ilg,jlg)
     5               ,lscls(ilg,jlg)
     8               ,swbeg(ilg,jlg)
     9               ,swend(ilg,jlg)
     +               ,xlwbeg(ilg,jlg)
     1               ,xlwend(ilg,jlg)
      common /canopy2/
     3               satang(ilg,jlg)
     4               ,rise15(ilg,jlg)
     5               ,rise55(ilg,jlg)
     7               ,navlat(jlg)
     8               ,navlon(ilg)
     1               ,mlai(ilg,jlg)
      common /canopy3/CLAT,CLON,DX,DY,MINLAT,MINLON,MAXLAT,MAXLON
