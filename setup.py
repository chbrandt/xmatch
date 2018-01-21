from setuptools import find_packages, setup

import versioneer

setup(name='xmatch',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Astronimcal catalogs cross-matching package',
      url='http://github.com/chbrandt/xmatch',
      packages=find_packages(),
      zip_safe=False,
      include_package_data=True)
