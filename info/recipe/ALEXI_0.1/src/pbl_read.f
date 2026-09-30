      subroutine pbl_read(MDATE)
      include 'USflux_grids.inc' 
      include "USflux_dir.inc"

      common/lookup/lookup_i(ilg,jlg),lookup_j(ilg,jlg)
      character(len=11) :: cmd
      character(len=9) :: zip
      character(len=255) :: cmd1, cmd2, cmd3, omd1, omd2, omd3
      character(len=255) :: dir, f1, f2, f3, s1, s2
      character(len=255) :: cstr1, cstr2, ostr1, ostr2
      real :: theta(kx,ky,kz)
      real :: pres(kx,ky,kz)
      integer ip,jp,ti, flip
      real :: diff1, tdiff1
      real, dimension(30) :: above_hgt 
      real :: p0, e0
      integer iy

      call fill_hgt_domain(above_hgt)

      cmd='gunzip -ff '
      zip='gzip -ff '
      dir='/work/waterforfood/bucricket/PROCESS_VIIRS/STATIC/CFSR'

      s1='/temp_profile.bin'
      s2='/pressure_profile.bin'

      iy=MDATE/1000-mod(MDATE,1000)/1000
      write(cstr1,'(I4,a,I7,a)') iy,'/',MDATE,'/temp_profile.bin.gz'
      write(cstr2,'(I4,a,I7,a)') iy,'/',MDATE,'/pressure_profile.bin.gz'
      write(cmd1,'(a,a,a)') cmd,trim(dir),trim(cstr1)
      write(cmd2,'(a,a,a)') cmd,trim(dir),trim(cstr2)

      write(ostr1,'(I4,a,I7,a)') iy,'/',MDATE,'/temp_profile.bin'
      write(ostr2,'(I4,a,I7,a)') iy,'/',MDATE,'/pressure_profile.bin'
      write(omd1,'(a,a,a)') zip,trim(dir),trim(ostr1)
      write(omd2,'(a,a,a)') zip,trim(dir),trim(ostr2)

      call system(cmd1)
      call system(cmd2)

      write(f1,'(a,I4,a,I7,a)') trim(dir),iy,'/',MDATE,trim(s1)
      write(f2,'(a,I4,a,I7,a)') trim(dir),iy,'/',MDATE,trim(s2)

      open(36,file=f1,form='unformatted',
     &          access='direct',recl=kx*ky*kz*4)
      open(37,file=f2,form='unformatted',
     &          access='direct',recl=kx*ky*kz*4)

      read(36,rec=1) theta
      read(37,rec=1) pres

      call system(omd1)
      call system(omd2)

      p0=100000.
      e0=0.286 
      do ja = 1, jlg
      do ia = 1, ilg
      if (lookup_i(ia,ja).ne.-9999.and.lookup_j(ia,ja).ne.-9999.) then
!       ti=aint(rise15(ia,ja)/3)
       ip=lookup_i(ia,ja)
       jp=lookup_j(ia,ja)
       do kc = 1, kz
        ztan(kc,ia,ja)=theta(ip,jp,kc)*(p0/pres(ip,jp,kc))**(e0)
       enddo
      endif
      enddo
      enddo

      DZZ=200.0
      Rd=287.04
      Gr=9.81
c
      do kc=1,41
        do ja=1,jlg
          do ia=1,ilg
          if(lookup_i(ia,ja).ne.-9999.) then
          if(lookup_j(ia,ja).ne.-9999.) then
            kabove=0
            kbelow=0
            if (kc.eq.1) then
              hzht=0.0
              htht(kc,ia,ja)=ztan(kz,ia,ja)
              if (htht(kc,ia,ja).eq.0.0) htht(kc,ia,ja)=BADATA
              goto 55
            else
              hzht=DZZ*(kc-1)
            endif

c
c  Now do all levels above the surface level at 'kz'.
            do km=kz-1,1,-1
c           do km = 1, kz 
              hcz=above_hgt(km) !zzan(km,ia,ja)
              if (hcz.gt.hzht .and. kbelow.eq.0) then
                kbelow=km-1
              endif
            enddo
            kabove=kbelow+1                   !we go UP one level
            zzab=above_hgt(kabove) !zzan(kabove,ia,ja)
            zzbl=above_hgt(kbelow) !zzan(kbelow,ia,ja)
            if (zzab.eq.zzbl) then
              htht(kc,ia,ja)=BADATA
            else
              htht(kc,ia,ja)=((hzht-zzbl)/(zzab-zzbl))*
     1       (ztan(kabove,ia,ja)-ztan(kbelow,ia,ja))+ztan(kbelow,ia,ja)
            endif

  55        continue

          endif
          endif
          enddo
        enddo
      enddo
      
      return
      end


      subroutine fill_hgt_domain(above_hgt)

      real :: above_hgt(30)

      above_hgt(1) = 15000.
      above_hgt(2) = 14000.
      above_hgt(3) = 13000.
      above_hgt(4) = 12000.
      above_hgt(5) = 11000.
      above_hgt(6) = 10000.
      above_hgt(7) = 9500.
      above_hgt(8) = 9000.
      above_hgt(9) = 8500.
      above_hgt(10) = 8000. 
      above_hgt(11) = 7500.
      above_hgt(12) = 7000.
      above_hgt(13) = 6500.
      above_hgt(14) = 6000.
      above_hgt(15) = 5500.
      above_hgt(16) = 5000.
      above_hgt(17) = 4500.
      above_hgt(18) = 4000.
      above_hgt(19) = 3500.
      above_hgt(20) = 3000
      above_hgt(21) = 2600.
      above_hgt(22) = 2200.
      above_hgt(23) = 1800.
      above_hgt(24) = 1400.
      above_hgt(25) = 1000
      above_hgt(26) = 700.
      above_hgt(27) = 500.
      above_hgt(28) = 300.
      above_hgt(29) = 100.
      above_hgt(30) = 0.

      return 
      end

c
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  + END  END  END  END  END  END  END  END  END  END  END  END  END  END +
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
