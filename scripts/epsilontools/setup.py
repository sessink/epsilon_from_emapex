from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize


setup(
   name='epsilontools',
   version='1.0',
   description='A module containing tools to analyse microstructure temperature from EM-APEX floats.',
   author='Sebastian Essink',
   author_email='sebastianessink@gmail.com',
   packages=['epsilontools'],  #same as name
   # install_requires=['pandas', 'numpy'], #external packages as dependencies,
   zip_safe=False
)
