!------------------------------------------------------------------
! US FLUX PROJECT PARAMETERS
!------------------------------------------------------------------
      include "USflux_parm.inc"
!
!------------------------------------------------------------------
! USFLUX PROJECT COMMON BLOCKS
!-----------------------------------------------------------------
      common/availh2o/brz,awfrz,awcrz,bsfc,awfsfc,awcsfc
      common/cover/iclass,beta,gvs,xndvi,fpar,dstom
      common/dayflux/aparday,acday,rnday,rnsday,rncday,gday,eday, &
      & esday,ecday,hday,hsday,hcday,sday
      common/dayflux2/erefday
      common/dayflux3/xlwupday,xlwdnday,swupday
      common/dayobsflux/rnobsday,hobsday,eobsday,gobsday,aobsday,precpday
      common/daypotflux/epotday,espotday,ecpotday
      common/departure/xmc,xbc,xcc,xms,xbs,xcs
      common/hrdata/tloc(nohr),ta(nohr),ea(nohr),wind(nohr),sdn(nohr), &
      & xlwdn(nohr),rnet(nohr),rnsoil(nohr),g(nohr),pres(nohr),nohrin
      common/hrdata2/epot(nohr),ecpot(nohr),espot(nohr)
      common/hrdata3/xlst(nohr)
      common/hrdata4/eref(nohr)
      common/hrdata5/xlwup(nohr),swup(nohr)
      common/hrdata_met/ctloc(kx,ky,kt),cta(kx,ky,kt),cea(kx,ky,kt), &
      & cwind(kx,ky,kt),csdn(kx,ky,kt),cxlwdn(kx,ky,kt),cpres(kx,ky,kt), &
      & clst(kx,ky,kt),clapse(ilg,jlg)
      common/lookup/lookup_i(ilg,jlg),lookup_j(ilg,jlg)
      common/slookup/insol_i(ilg,jlg),insol_j(ilg,jlg)
      common/hrobs/rnobs(nohr),hobs(nohr),xleobs(nohr),gobs(nohr), &
      & aobs(nohr),precp(nohr)
      common/lookup/itabclass(nclass),tablai(nbin,nclasm), &
      & tabfpar(nbin,nclasm),tabbeta(nclass),tabaleaf(nclass,3), &
      & tabadead(nclass,3),tabhmin(nclass),tabhmax(nclass), &
      & tabxl(nclass),tabgvs(nclass),tabdstom(nclass),tabperen(nclass), &
      & tabfcmin(nclass),tabdesc(nclass)
      common/moisture/faw,awf
      common/moisture2/faw2
      common/stress/fpet,fsdn
      common/sun/strt,end
      common/usflags/clear

      logical clear
      CHARACTER(LEN=50) tabdesc
