from unittest import TestCase
from test_partse import *
__author__ = 'Aleksej'

def power(a, b):
    return a**b

class TestParse(TestCase):
    def test_parse3(self):
        model = """
-REACTIONS
R1: A + B -> 3 c
-CONSTRAINTS
R1[0, 100]
-EXTERNAL METABOLITES
C
-OBJ
R1 1 1
-DESIGNOBJ
R1 1 1
"""
        r = parse(model)

class Metabolite:
    def __init__(self, name, coefficient):
        if not isinstance(name, basestring):
            raise TypeError("Name is not a string")

        if not isinstance(coefficient, (int, float)):
            raise TypeError("Coefficient is not a number")

        coefficient = float(coefficient)

        if coefficient <= 0:
            raise ValueError("Coefficient value is zero or negative")

        self.__name = name
        self.__coefficient = coefficient

    def __eq__(self, other):
        return self.get_name() == other.get_name()

    def get_name(self):
        return self.__name
    def get_coefficient(self):
        return self.__coefficient





a = Metabolite("bla", 1)
b = Metabolite("bla", 1)

test = [a,b]

c = Metabolite("bla", 1)

if c in test:
    print True

print a==b
class TestMetabolite(TestCase):
    def test_new(self):
        m = Metabolite("H2O", 3)

        self.assertTrue(m.name == "H2O")
        self.assertTrue(m.coefficient == 3)

    def test_invalid_coef(self):
        with self.assertRaises(TypeError, message="Providing not string name should rise TypeError"):
            Metabolite(15, "a")

        with self.assertRaises(TypeError, message="Providing character coefficient should raise TypeError"):
            Metabolite("H2O", "a")

        with self.assertRaises(ValueError, message="Providing zero or negative coefficient coefficient should rise ValueError"):
            Metabolite("H2O", 0)

        with self.assertRaises(ValueError, message="Providing zero or negative coefficient coefficient should rise ValueError"):
            Metabolite("H2O", -1)


