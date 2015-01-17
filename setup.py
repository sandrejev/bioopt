from distutils.core import setup
import os
import re

version_file = "__init__.py"
version_str = open(version_file, "rt").read()
mo = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_str, re.M)

if mo:
    version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in {0}.".format(version_file))

setup(name='bioopt',
      version=version,
      description='ioopt is a library for creating, manipulating and converting and exporting BioOpt files',
      author='Sergej Andrejev',
      author_email='sandrejev@gmail.com',
      url='https://github.com/sandrejev/bioopt',
      py_modules=[p[:-3] for p in os.listdir(".") if re.match(r"^[^_].*\.py$", p) and p != "setup.py"],
      data_files=[('examples', ['toy.bioopt'])]
 )