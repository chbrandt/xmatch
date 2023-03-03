from setuptools import setup

import versioneer

scripts = ['bin/xmatch', 'bin/_xmatch.py']

setup(
      scripts=scripts,
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
)