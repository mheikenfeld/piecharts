from setuptools import setup

setup(name='piecharts',
      use_scm_version=True,
      setup_requires=['setuptools_scm'],
      description='Plot an array of piecharts on matplotlib',
      url='http://github.com/mheikenfeld/piecharts',
      author='Max Heikenfeld',
      author_email='max.heikenfeld@physics.ox.ac.uk',
      license='BSD-3-Clause',
      packages=['piecharts'],
      zip_safe=False)
