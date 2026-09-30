      PROGRAM INTERP_LAI
      include "USflux_dir.inc"      
      parameter(ilgmx=2000, jlgmx=2000)
      parameter(BAD=-9999)
      character*256 filei,file1,file2,outfile

      l1=index(INPUT,' ')-1  
      filei=INPUT(1:l1)//'interp.in'
      open (unit=20,file=filei)
      read (20,*) file1,file2,outfile,f,ilg,jlg
      close (20)
      
      if (ilg.gt.ilgmx.or.jlg.gt.jlgmx) then
        write(6,*) 'Need to increase ILGMX JLGMX'
        stop
      endif
      
      open (unit=20, file=file1,status='old',
     &      form='unformatted',recl=4,access='direct') 
      open (unit=21, file=file2,status='old',
     &      form='unformatted',recl=4,access='direct') 
      open (unit=22, file=outfile,status='unknown',
     &      form='unformatted',recl=4,access='direct') 
      
      do ja=1,jlg
        do ia=1,ilg
          irec=(ja-1)*ilg+ia
          read (20,rec=irec) v1
          read (21,rec=irec) v2
          if (v1.ge.-.1.and.v2.ge.-.1) then
            val=v1+f*(v2-v1)
            if (val.le.0) val=0
          else
            val=BAD
          endif
          write (22,rec=irec) val 
        enddo
      enddo
      
      close (20)
      close (21)
      close (22)
                 
      stop
      end
      
      
 
