__author__ = 'Sergej Andrejev'
__version__ = '0.2.0'

from .bioopt_parser import *
from .output_parser import *
from .converter import *
from .model import *


def toy_path():
    import os.path
    dir, filename = os.path.split(__file__)
    return os.path.abspath(os.path.join(dir, "toy.txt"))