#!python
import os
AddOption('--prefix',
          dest='prefix',
          type='string',
          nargs=1,
          action='store',
          metavar='DIR',
          help='installation prefix')

version = "0.2"
env = Environment(PREFIX = GetOption('prefix'))
prefix  = os.environ.get('PREFIX')
build_prefix = os.environ.get('BUILD_PREFIX')
src_dir = os.environ.get('SRC_DIR')
base = os.path.abspath(os.path.join(prefix,os.pardir))
base = os.path.join(base,'work')
# sourcePath = os.path.join(base,'source')
binPath = os.path.join(prefix,'bin')
lib_path = os.path.join(build_prefix, 'lib')
print(src_dir)
include_path = os.path.join(src_dir, 'ALEXI_{}'.format(version))
# include_path = os.path.join(src_dir, 'source', 'ALEXI_{}'.format(version))
env.Replace(FC = os.path.join(build_prefix,"bin", os.environ.get('FC')))
    
env.Append(FLAGS = ['-c', '-g', '-fPIC -ffixed-line-length-132'])
env.Append(FPATH = [include_path])
alexi = env.Program(target='alexi_proc', source=['USflux.f90', 'USflux_utl.f','USflux_rad.f',
                                                            'USflux_cover.f','USflux_run.f','USflux_clear.f',
                                                            'USflux_cloud.f','ALEXI.f', 'ALEXI_atmos.f','ALEXI_rad.f',
                                                            'ALEXI_utl.f', 'ALEXI_water.f','pbl_read.f', 'sfc_read.f',
                                                            'USflux_main.f'])
landcover = env.Program(target='landcover', source=['landcover.f'])
env.Install(binPath, [landcover, alexi])
env.Alias('install', binPath)
