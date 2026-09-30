      subroutine landcover(MDATE,tile,npoints, part)

      include 'USflux_dir.inc'     
      include 'USflux_parm.inc'
      include 'ALEXI_parm.inc'
      real class_num(ilg,jlg)
      character*256 d,infile,cfile,outfile
      integer yr, dy, iunit1, iunit2
      character*256 classfile 
      dimension freq(nclass)
      common/lookup/alv(nclass),aln(nclass),all(nclass),adv(nclass),
     &              adn(nclass),adl(nclass),hmin(nclass),hmax(nclass),
     &              xl(nclass), rs(nclass)
      integer npoints
      real clai
      integer*4 OD
!     integer*10 npoints
      character*16 tile
      character*3 part

      OD=MDATE 
      write(6,*) "MDATE = ", MDATE
      call getdate(MDATE,year,doy)
      write(6,*) "YEAR DOY = ", year, doy
!      call laiopen("MLAI",year,doy,100)
      call loadclasstable

      d='./INPUTS/ALEXI_INPUT/'
      write(cfile,'(a,a,I7,a,a,a,a,a)') trim(d),'veg_',OD,'_',trim(tile),'_',trim(part),'.input'
      write(6,*) cfile
      open(unit=102,file=cfile)
      d='./INPUTS/ALEXI_INPUT/veg_us'
      write(outfile,'(a,I7,a,a,a,a,a)') trim(d),OD,'_',trim(tile),'_',trim(part),'.input'
      write(6,*) outfile
      open(unit=103,file=outfile)

!     npoints=14062500

      do ja=1,npoints

c opens cfile(Veg_class), reads it in        
        read(102,*,end=200)i,j,xntot,(freq(k),k=1,nclass),xlai
!        xntot=100
        xlai=xlai/10.
        if (xlai.lt.0) xlai=BAD
        if (xlai.gt.10) xlai=BAD
        xndvi=100
  
          if (xlai.ne.BAD) then
            fc=1.0-exp(-0.5*xlai)
            if (fc.lt.0) fc=0.0
            if (fc.gt.1) fc=1.0
            call weighted_avg(freq,height,xleaf,aleafv,aleafn,aleafl,
     &                        adeadv,adeadn,adeadl,z0eff,dispeff,fc,
     &                        class,xntot,rsmin)
          else
            xndvi=BAD
            fc=BAD
            aleafv=BAD
            aleafn=BAD
            aleafl=BAD
            adeadv=BAD
            adeadn=BAD
            adeadl=BAD
            height=BAD
            xleaf =BAD
            z0eff =BAD
            dispeff=BAD
            rsmin=BAD
          endif
          if (class.eq.1) then
            xndvi=BAD
            fc=BAD
            aleafv=BAD
            aleafn=BAD
            aleafl=BAD
            adeadv=BAD
            adeadn=BAD
            adeadl=BAD
            height=BAD
            xleaf =BAD
            z0eff =BAD
            dispeff=BAD
            rsmin=BAD
          endif

          
          write(103,1030)i,j,aleafv,aleafn,aleafl,adeadv,
     &                 adeadn,adeadl,height,xleaf,z0eff,dispeff,rsmin
          
      enddo
 200  continue
 
1020  format(4f9.2,f10.3,2f9.2)     
1030  format(2i5 ,8f13.5,3f13.5)

      close(100)
      close(102)
      close(103)   
 
      end
      
      
      subroutine weighted_avg(freq,height,xleaf,aleafv,aleafn,aleafl,
     &                        adeadv,adeadn,adeadl,z0eff,dispeff,fc,
     &                        class,xntot,rsmin)
c     **************************************************************            
c     *
c     *  Martha Anderson
c     *  Created:  10/5/04
c     *  
c     **************************************************************
      include 'USflux_parm.inc'
      dimension freq(nclass)
      common/lookup/alv(nclass),aln(nclass),all(nclass),adv(nclass),
     &              adn(nclass),adl(nclass),hmin(nclass),hmax(nclass),
     &              xl(nclass), rs(nclass)

      z0w=0.00035		! Roughness length for water   

      height=0.
      xleaf =0.
      aleafv=0.
      aleafn=0.
      aleafl=0.
      adeadv=0.
      adeadn=0.
      adeadl=0.
      zsum=0.
      dsum=0.
      xn  =0.
      xnv =0.
      rsmin =0. 

      if (class.eq.1) then   ! Pixel is predominantly water
        height=BAD
        xleaf =BAD
        aleafv=BAD
        aleafn=BAD
        aleafl=BAD
        adeadv=BAD
        adeadn=BAD
        adeadl=BAD
        z0eff =BAD
        dispeff=BAD
        rsmin=BAD
      else
        do ic=1,nclass
          if (freq(ic).gt.0) then
          if (ic.eq.1) then
            disp   =0.0
            z0     =z0w	
          else
            hc=hmin(ic)+fc*(hmax(ic)-hmin(ic))
            height =height +freq(ic)*hc
            xleaf  =xleaf  +freq(ic)*xl(ic)
            aleafv =aleafv +freq(ic)*alv(ic)
            aleafn =aleafn +freq(ic)*aln(ic)
            aleafl =aleafl +freq(ic)*all(ic)
            adeadv =adeadv +freq(ic)*adv(ic)
            adeadn =adeadn +freq(ic)*adn(ic)
            adeadl =adeadl +freq(ic)*adl(ic)
            rsmin  =rsmin  +freq(ic)*rs(ic)
            call canopyarch_lai(fc,hc,z0,disp)
            xnv=xnv+freq(ic)
          endif         
          zsum   =zsum   
     &            +freq(ic)/((log((50.-disp)/z0))*(log((50.-disp)/z0))) 
          dsum   =dsum   +freq(ic)*disp
          xn     =xn     +freq(ic)
          endif
        enddo

        xleaf =xleaf /xnv
        aleafv=aleafv/xnv
        aleafn=aleafn/xnv
        aleafl=aleafl/xnv
        adeadv=adeadv/xnv
        adeadn=adeadn/xnv
        adeadl=adeadl/xnv
        height=height/xn
        rsmin=rsmin/xn
        dispeff=dsum/xn
        zsum   =zsum/xn
        z0eff=(50.-dispeff)/(exp(1./sqrt(zsum))) 

      endif
      
      return
      end

      subroutine canopyarch_lai(fc,hc,z0,disp)
c     **************************************************************            
c     *
c     *  Martha Anderson
c     *  Created:  10/5/04
c     *  
c     **************************************************************
      cd=0.2 
      z0s=0.005			! Roughness length for bare soil 
      z0w=0.00035		! Roughness length for water   
      fcbare=0.0  		! fc.le.fcbare considered to be bare soil 
  
  
c	Bare soil case - z0 is soil roughness
c       --------------------------------------------------
      if (fc.le.fcbare .or. hc.eq.0.0) then
         z0=z0s
         disp=0.0
                  
c	All other cases - use Massman eqs for z0dh & dispdh 
c       --------------------------------------------------
      else 
c         $ufact=0.360-0.264*exp(-15.1*$cd*$xlai);
c         $xn=$cd*$xlai/(2.*($ufact*$ufact));
c         $dispdh=0.7-(1./(5.*$xn)*(1.-exp(-3.3*$xn)));
c         if ($dispdh<0) {$dispdh=0.0};       
c         $disp=$dispdh*$hc;
c         $z0dh=(1.-$dispdh)*exp(-0.4/$ufact);
c         $z0=$z0dh*$hc;
c         $tmp=1./8.;
c         print("ZODH: $z0dh $tmp\n");
         dispdh=2./3.
         z0dh=1./8.
         disp=dispdh*hc
         z0=z0dh*hc
         if (z0.lt.z0s) z0=z0s
      endif
      
      return
      end

 
      subroutine loadclasstable
c     **************************************************************            
c     *
c     *  Martha Anderson
c     *  Created:  10/5/04
c     *  
c     **************************************************************
      include 'USflux_dir.inc' 
      include 'USflux_parm.inc' 
      character*200 classfile
      common/lookup/alv(nclass),aln(nclass),all(nclass),adv(nclass),
     &              adn(nclass),adl(nclass),hmin(nclass),hmax(nclass),
     &              xl(nclass), rs(nclass)
      
      l1=index(HOMEDIR,' ')-1
      classfile=HOMEDIR(1:l1)//'store/landcover.txt'
      open(unit=30,file=classfile)
      read(30,*)
      do ic=1,nclass
        read(30,*)dum,alv(ic),aln(ic),all(ic),adv(ic),adn(ic),adl(ic),
     &            rs(ic),hmin(ic),hmax(ic),xl(ic)
      enddo
     
      write(6,*) "RS = ", rs 
      close(30)
      return
      end
 
      
      subroutine getdate(MDATE,year,doy)
c     **************************************************************            
c     *
c     *  Martha Anderson
c     *  Created:  10/5/04
c     *  
c     **************************************************************
      include 'USflux_dir.inc' 
      character*9 datein
      character*8 CDATE
      character*200 datefile
      
c      l1=index(INDIR,' ')-1
c      datefile=INDIR(1:l1)//'DATELIST'
c      open(19,file=datefile,status='old')
c      read(19,1009)MDATE,datein,CDATE
c 1009 format(i7,1x,a9,1x,a8)
c      read(CDATE,'(i8)')JDATE
c      write(6,*)'CDATE = ',CDATE
c      write(6,*)'MDATE = ',MDATE
      doy=mod(MDATE,1000)
      year=MDATE/1000-mod(MDATE,1000)/1000
!      close(19)
      
      return
      end
     
      
      subroutine laiopen(cvar,year,doy,iunit)
c     **************************************************************            
c     *
c     *  Martha Anderson
c     *  Created:  10/5/01
c     *  
c     *************************************************************** 
      include 'USflux_dir.inc'
      include 'USflux_parm.inc'
      real*4 var(ilg,jlg)
      character*256 dir,outfile
      character*4   cvar
      character*7   cyyyyddd
c
c  Output files are written to directory OUTDIR      
      
c  Open binary file for Transform
      iyyyyddd=year*1000+doy
      write(cyyyyddd,'(i7)')iyyyyddd
      
      dir='/data/data123/chain/4KM/GBIM/inputs/MLAI'
      write(outfile,'(a,I7,a)') trim(dir),iyyyyddd,'.dat'
!      outfile=OUTDIR(1:l1)//cvar//cyyyyddd//'.dat'
!      l2=index(outfile,' ')-1
c      write(6,*)outfile(1:l2)
      open (unit=iunit, file=outfile, 
     &      form='unformatted',recl=4,access='direct')
            
      return
      end  
          
      subroutine lairead(i,j,ilg,jlg,iunit,val)
c     **************************************************************
c     *
c     *  Martha Anderson
c     *  Created:  9/11/04
c     * 
c     ***************************************************************
      real*4  rval

c  Read from binary file
      jj=jlg-j+1
      jj=j
      irec=(jj-1)*ilg+i
      read (iunit,rec=irec) rval
c      call byteswapr4(rval)
      val=rval

      return
      end 

      subroutine laiswapr4(r)

c     does a byteswap on real*4 number

      integer*1 ii(4), jj(4)
      real*4 r, s, t
      equivalence (s,ii)
      equivalence (t,jj)

      s = r

      jj(1) = ii(4)
      jj(2) = ii(3)
      jj(3) = ii(2)
      jj(4) = ii(1)

      r = t

      return
      end
 
