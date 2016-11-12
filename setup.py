#!/usr/bin/env python

from distutils.core import setup

#To prepare a new release
#python setup.py sdist upload

setup(name='vmap',
    version='0.1.0',
    description='Velocity map generation using the Ames Stereo Pipeline image correlator',
    author='David Shean',
    author_email='dshean@gmail.com',
    license='MIT',
    url='https://github.com/dshean/vmap',
    packages=['vmap'],
    long_description=open('README.md').read(),
    install_requires=['numpy','gdal','pygeotools'],
    #Note: this will write to /usr/local/bin
    scripts=['vmap/vmap.py', 'vmap/disp2v.py', 'vmap/disp_view.sh']
)

