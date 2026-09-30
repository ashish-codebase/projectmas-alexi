#!python

import os
import shutil
import platform
import subprocess


AddOption('--prefix',
          dest='prefix',
          type='string',
          nargs=1,
          action='store',
          metavar='DIR',
          help='installation prefix')
# env = Environment(PREFIX = GetOption('prefix'))
path = ["/Users/mschull/anaconda3/conda-bld/projectmas_alexi_1588608011610/_build_env/bin/", '/bin', '/usr/bin']
# env = Environment(PREFIX = GetOption('prefix'))
env = Environment(ENV={'PATH': path}, F90SUFFIXES=['.f', '.f90', '.f95'])
env["SHELL"]="/bin/bash"
# gfortran_path = os.path.join(os.path.dirname(os.environ.get('FC')), 'gfortran')
# shutil.copyfile(os.environ.get('FC'), gfortran_path)
# env = DefaultEnvironment(tools=[os.path.basename(os.environ.get('FC')), 'gnulink'],
#                          FC=gfortran_path)
#prefix = GetOption('prefix')
# prefix  = os.environ.get('PREFIX')
# build_prefix = os.environ.get('BUILD_PREFIX')
# src_dir = os.environ.get('SRC_DIR')
# base = os.path.abspath(os.path.join(prefix,os.pardir))
# base = os.path.join(base,'work')
# sourcePath = os.path.join(base,'source')
bin_path = "../bin/"
# prefix = os.environ.get('PREFIX')
# build_prefix = os.environ.get('BUILD_PREFIX')
# src_dir = os.environ.get('SRC_DIR')
# base = os.path.abspath(os.path.join(prefix, os.pardir))
# base = os.path.join(base, 'work')
# source_path = os.path.join(base, 'source')
# bin_path = os.path.join(prefix, 'bin')
lib_path = "./"
include_path = os.getcwd()
fortran_compiler = os.path.join(path[0], "x86_64-apple-darwin13.4.0-gfortran")
env["F90"]=fortran_compiler
env["FORTRAN"]=fortran_compiler
env.Append(F90FLAGS=['-g', '-fPIC', '-ffixed-line-length-132'])
env.Append(FORTRANFLAGS=['-g', '-fPIC', '-ffixed-line-length-132'])
# FLAGS  = -c -g -I$(SRCU) -I$(SRC) -fPIC -ffixed-line-length-132 # MAC
env.Append(FPATH=[include_path])
env.Append(LIBPATH = [lib_path])
# landcover = env.Program(target='landcover', source=['landcover.f'])
alexi = env.Program(target='alexi_proc', source=['USflux.f90', 'USflux_utl.f', 'USflux_rad.f',
                                                 'USflux_cover.f', 'USflux_run.f', 'USflux_clear.f',
                                                 'USflux_cloud.f', 'ALEXI.f', 'ALEXI_atmos.f', 'ALEXI_rad.f',
                                                 'ALEXI_utl.f', 'ALEXI_water.f', 'pbl_read.f', 'sfc_read.f','landcover.f',
                                                 'USflux_main.f'])
env.Install(bin_path, [alexi])
env.Alias('install', bin_path)
