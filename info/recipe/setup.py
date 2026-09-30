#!/usr/bin/env python
from setuptools import setup
import os
import shutil
import subprocess

version = "0.1"
try:
    from setuptools import setup

    setup_kwargs = {'entry_points': {'console_scripts': ['pyalexi=pyalexi.pyalexi:alexi']}}
except ImportError:
    from distutils.core import setup

    setup_kwargs = {'scripts': ['bin/pyalexi']}

base = os.getcwd()
prefix = os.environ.get('PREFIX')
src_dir = os.environ.get('SRC_DIR')
binDir = os.path.join(src_dir, 'bin')
process_dir = os.path.join(os.path.abspath(os.path.join(prefix, os.pardir)), 'work')

print("installing ALEXI...")
mkPath = os.path.join(process_dir, 'ALEXI_{}'.format(version), 'src')
os.chdir(mkPath)
subprocess.call(["scons", "-Q", "--prefix={}".format(prefix), "install"])
subprocess.call(["scons", "-c"])
os.chdir(base)

setup(
    name="alexi",
    version="{}".format(version),
    description="pythonic wrapper for alexi",
    author="Mitchell Schull",
    author_email="mitch.schull@noaa.gov",
    # packages=['pyalexi'],
    py_modules=['pyalexi.pyalexi'],
    platforms='Posix; MacOS X; Windows',
    license='BSD 3-Clause',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        # Uses dictionary comprehensions ==> 2.7 only
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: GIS',
    ],
)
