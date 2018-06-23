from setuptools import setup

import versioneer

scripts = ['bin/xmatch', 'bin/_xmatch.py']

setup(name='xmatch',
      description='Astronimcal catalogs cross-matching tool',
      packages=['xmatch'],
      scripts=scripts,
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      url='http://github.com/chbrandt/xmatch'
      )
