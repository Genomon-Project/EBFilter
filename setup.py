#!/usr/bin/env python

from distutils.core import setup

setup(name='ebfilter',
      version='0.1.2',
      description='Python tools for filtering somatic mutations using beta-binomial sequencing error model.',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/friend1ws/EBFilter',
      package_dir = {'': 'lib'},
      packages=['ebfilter'],
      scripts=['EBFilter'],
      license='GPL-3'
     )

