from setuptools import setup

setup(name='piecharts',
      version='0.1',
      description='Plot an array of piecharts on matplotlib',
      url='http://github.com/mheikenfeld/piecharts',
      author='Max Heikenfeld',
      author_email='max.heikenfeld@physics.ox.ac.uk',
      license='GNU',
      packages=['piecharts'],
      install_requires=['numpy','matplotlib'],
      zip_safe=False)