subroutine USflux(MDATE, t, npoints, part)
    !     **************************************************************
    !     *
    !     *  Main procedure for computing carbon, heat and water fluxes
    !     *  over the continental US.
    !     *
    !     *  Martha Anderson
    !     *  Created:  02/14/01
    !     *
    !     ***************************************************************
    INCLUDE "USflux.inc"
    INCLUDE "ALEXI_f90.inc"
    INCLUDE "flag.inc"
    dimension iflag(ilg, jlg)
    logical haveSolar
    CHARACTER(LEN=256) ofile, vegfile, cfile, d, logfile
    CHARACTER(LEN=256) rfile, mfile, pfile, sfile
    integer(kind = 4) m
    integer(kind = 4) OD
!    CHARACTER(LEN=16) t
    character*16 t
    CHARACTER(LEN=3) part
    !     integer*10 npoints
    save

    ! Initializations:
    write(6, *) "USflux MDATE = ", MDATE
    doy = mod(MDATE, 1000)
    year = MDATE / 1000 - mod(MDATE, 1000) / 1000
    d = './ALEXI_LOG/'
    write(logfile, '(a,a,I7,a,a,a,a,a)') trim(d), '/', MDATE, '_', trim(t), '_', trim(part), '.LOG'
    open(99, file = logfile)

    call set_constants

    ntot = 0
    nbad = 0
    nconv = 0
    nfail = 0
    ncloud = 0
    nwater = 0
    ! Loop over all pixels in model domain
    !------------------------------------

    d = './INPUTS/ALEXI_INPUT/'
    OD = MDATE
    write(mfile, '(a,a,I7,a,a,a,a,a)') trim(d), 'met_', OD, '_', trim(t), '_', trim(part), '.input'
    open(unit = 102, file = mfile)

    write(sfile, '(a,a,I7,a,a,a,a,a)') trim(d), 'sat_', OD, '_', trim(t), '_', trim(part), '.input'
    open(unit = 110, file = sfile)

    write(pfile, '(a,a,I7,a,a,a,a,a)') trim(d), 'nparm_', OD, '_', trim(t), '_', trim(part), '.input'
    open(unit = 120, file = pfile)

    write(rfile, '(a,a,I7,a,a,a,a,a)') trim(d), 'profile_', OD, '_', trim(t), '_', trim(part), '.input'
    open(unit = 130, file = rfile)

    d = './INPUTS/ALEXI_INPUT/veg_us'
    write(vegfile, '(a,I7,a,a,a,a,a)') trim(d), OD, '_', trim(t), '_', trim(part), '.input'
    open (400, file = vegfile)

    d = './OUTPUTS/ALEXI_OUTPUT/'
    write(ofile, '(a,a,I7,a,a,a,a,a)') trim(d), 'output_', OD, '_', trim(t), '_', trim(part), '.input'
    open(unit = 103, file = ofile)

    !      npoints=14062500
    do m = 1, npoints !179,179

        ntot = ntot + 1
        converged = .FALSE.
        ibad = BAD
        call extract_input(m, ibad, MDATE)
        fc0 = fc

        !		Check IFLAG array to determine how this pixel should
        !               be processed

        !     Skip point if flagged bad:
        if(badinput)then
            nbad = nbad + 1
            go to 1000
            !     Skip point if designated as water or ice:
        elseif (.not.iswater_inland.and.&
                (iclass.eq.WATER .or. iclass.eq.ICE .or.&
                        iflag(ia, ja).eq.3))then
            iswater = .TRUE.
            nwater = nwater + 1
            go to 1000
        elseif(writeme)then
        endif

        ! 		Execute clear- or cloudy-day flux procedure

        if (clear) then
            nclear = nclear + 1

            fc = fc0
            if (iswater_inland) then
                call ALEXI_water(ia, ja, ierr, ibad, iter)
            else
                call ALEXI(ia, ja, ierr, ibad, iter)
            endif

            if (badinput) then
                nbad = nbad + 1
                iflag(ia, ja) = 1
                go to 1000
            else if (.not.converged) then
                !                 iflag(ia,ja)=5
                go to 900
            endif
            900       if (converged) then
                nconv = nconv + 1
                call clear_day_proc(ibad)
            else
                nfail = nfail + 1
                go to 1000    ! Comment out if running cloudy proc
                call cloudy_day_proc(ibad)  ! If pixel doesn't converge, run cloud proc.
            endif

        else
            ncloud = ncloud + 1
            !cc           go to 1000  	! Comment out if running cloudy proc
            call cloudy_day_proc(ibad)
        endif

        ! 		Write diagnostics to screen
        1001    continue
        !        write(6,*) "before"
        call screen_output
        !        write(6,*) "after"

        ! 		Save any grid info to be retained
        1000    continue
        !        call store_output_bin(ia,ja,ierr,ibad,iter)

    enddo    ! JLG loop
    !-------------------------------------

    write(99, *)'          Bad input: ', nbad
    write(99, *)'       Water or ice: ', nwater
    write(99, *)'       Cloudy pixel: ', ncloud
    write(99, *)'          Converged: ', nconv
    write(99, *)'   Did not converge: ', nfail
    write(99, *)'       Total points: ', ntot

    !      do iunit=130,271
    !        close (iunit)
    !      enddo

    !      close (100)
    !      close (300)
    !      close (301)

    !      close(55)
    return
end
      