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
path = [os.path.dirname(os.environ.get('FC')), '/bin', '/usr/bin']
# env = Environment(PREFIX = GetOption('prefix'))
env = Environment(PREFIX = GetOption('prefix'), ENV={'PATH': path}, F90SUFFIXES=['.f', '.f90', '.f95'])
# env["SHELL"]="/bin/bash" #MAC
# gfortran_path = os.path.join(os.path.dirname(os.environ.get('FC')), 'gfortran')
# shutil.copyfile(os.environ.get('FC'), gfortran_path)
# env = DefaultEnvironment(tools=[os.path.basename(os.environ.get('FC')), 'gnulink'],
#                          FC=gfortran_path)
# prefix = GetOption('prefix')
prefix = os.environ.get('PREFIX')
build_prefix = os.environ.get('BUILD_PREFIX')
src_dir = os.environ.get('SRC_DIR')
base = os.path.abspath(os.path.join(prefix, os.pardir))
base = os.path.join(base, 'work')
source_path = os.path.join(base, 'source')
bin_path = os.path.join(prefix, 'bin')
lib_path = os.path.join(build_prefix, 'lib')
include_path = os.path.join(src_dir, 'ALEXI_0.1', 'src')
fortran_compiler = os.path.join(path[0], os.path.basename(os.environ.get('FC')))
env["F90"]=fortran_compiler
env["FORTRAN"]=fortran_compiler
env.Append(F90FLAGS=['-g', '-fPIC', '-ffixed-line-length-132'])
env.Append(FORTRANFLAGS=['-g', '-fPIC', '-ffixed-line-length-132'])
# FLAGS  = -c -g -I$(SRCU) -I$(SRC) -fPIC -ffixed-line-length-132 # MAC
env.Append(FPATH=[include_path])
# env.Append(LIBPATH = [lib_path])
alexi = env.Program(target='alexi_proc', source=['USflux.f90', 'USflux_utl.f', 'USflux_rad.f',
                                                'USflux_cover.f','USflux_run.f', 'USflux_clear.f',
                                                 'USflux_cloud.f', 'ALEXI.f', 'ALEXI_atmos.f', 'ALEXI_rad.f',
                                                 'ALEXI_utl.f', 'ALEXI_water.f', 'pbl_read.f', 'sfc_read.f','landcover.f',
                                                 'USflux_main.f'])
# landcover = env.Program(target='landcover', source=['landcover.f'])
env.Install(bin_path, [alexi])
# build_executable = os.path.join(bin_path, "alexi_proc")
# prefix_bin_path = os.path.join(prefix, "bin")
# prefix_executable = os.path.join(prefix_bin_path, "alexi_proc")
# shutil.copyfile(build_executable,prefix_executable)
env.Alias('install', bin_path)
