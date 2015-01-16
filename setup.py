from distutils.core import setup
import os,re

setup(name='bioopt',
      version='1.0',
      description='ioopt is a library for creating, manipulating and converting and exporting BioOpt files',
      author='Sergej Andrejev',
      author_email='sandrejev@gmail.com',
      url='https://github.com/sandrejev/bioopt',
      py_modules=[p[:-3] for p in os.listdir(".") if re.match(r"^[^_].*\.py$", p) and p != "setup.py"],
      data_files=[('examples', ['toy.bioopt'])]
 )