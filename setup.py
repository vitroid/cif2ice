#!/usr/bin/env python2

from setuptools import setup
import os
import codecs
import re

#Copied from wheel package
here = os.path.abspath(os.path.dirname(__file__))
#README = codecs.open(os.path.join(here, 'README.txt'), encoding='utf8').read()
#CHANGES = codecs.open(os.path.join(here, 'CHANGES.txt'), encoding='utf8').read()

with codecs.open(os.path.join(os.path.dirname(__file__), 'cif2ice', '__init__.py'),
                 encoding='utf8') as version_file:
    metadata = dict(re.findall(r"""__([a-z]+)__ = "([^"]+)""", version_file.read()))

setup(name='cif2ice',
      version=metadata['version'],
      description='Prepare ice modules for GenIce from CIF file',
      #long_description=README + '\n\n' +  CHANGES,
      classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        ],
      author='Masakazu Matsumoto',
      author_email='vitroid@gmail.com',
      url='https://github.com/vitroid/cif2ice/',
      keywords=['cif2ice',],
      license='MIT',
      packages=['cif2ice'],
      install_requires=['requests', 'validators', 'pycifrw', 'numpy' ],
      entry_points = {
              'console_scripts': [
                  'cif2ice = cif2ice.__main__:main'
              ]
          }
      )

