c------------------------------------------------------------------
c           INCLUDE FILE FOR US FLUX PROJECT DIRECTORIES
c------------------------------------------------------------------
      character*500 HOMEDIR,INDIR,TABDIR,OUTDIR,INPUT

c Root directory
c HOMEDIR contains:
c - DATELIST
      data HOMEDIR/"./ALEXI/"/
      
c Inputs directory
c INDIR contains:
c - All gridded inputs 
      data INDIR/"./ALEXI/inputs/"/
      
c Tables directory
c TABDIR contains:
c - landcover.txt
c - Default ALEAF and RSOIL maps
      data TABDIR/"./ALEXI/store/"/

c Output directory
c OUTDIR contains:
c - binary matrices
      data OUTDIR/"./ALEXI/fluxes/"/

c Program command line input directory
c INPUT contains:
c - command line inputs for fortran programs
      data INPUT/"./ALEXI/inputs/"/


      
