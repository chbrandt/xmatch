from setuptools import setup
import glob
import os

import versioneer

scripts = [fname for fname in glob.glob(os.path.join('bin', '*'))]

setup(name='xmatch',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Astronimcal catalogs cross-matching package',
      url='http://github.com/chbrandt/xmatch',
      packages=['xmatch'],
      scripts=scripts
      )
