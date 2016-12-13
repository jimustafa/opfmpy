from __future__ import absolute_import, division, print_function

from setuptools import setup, find_packages


with open('README.rst', 'r') as f:
    long_description = f.read()

setup(
    name='opfmpy',
    version='0.0.0',
    description='a Python package implementing the optimized projection functions method',
    author='Jamal I. Mustafa',
    author_email='jimustafa@gmail.com',
    license='BSD',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    packages=find_packages(exclude=['docs', 'tests']),
    scripts=['opfmpy/scripts/opfm.py'],
    long_description=long_description,
    install_requires=[
        'numpy',
        'scipy',
        'wannier90-utils',
        'codiag',
        'f90nml',
        ],
    tests_require=['pytest'],
)
