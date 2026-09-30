!
!  Include file for two-source time-integrated model over US domain
!
!------------------------------------------------------------------
!  MODEL PARAMETERS
!------------------------------------------------------------------
!  Parameters describing grid
      integer*8 filg, fjlg
      parameter(filg=1456,fjlg=625)
!
!  ALEXI input and core data arrays
      common /fflag/flag(filg,fjlg)

