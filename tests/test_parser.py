from unittest import TestCase
from src.model import *
from src.parser import *
class TestBiooptParser(TestCase):
    model = """-REACTIONS
R1: A + B -> 3 C
R2: B + C <-> 1 E
-CONSTRAINTS
R1[-100, 100]
R2[-100, 100]
-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 1 1
"""
    def test_strip_comments(self):
        parser = BiooptParser()
        self.assertEquals(parser.strip_comments("test"), "test")
        self.assertEquals(parser.strip_comments("    test"), "test")
        self.assertEquals(parser.strip_comments("test    "), "test")
        self.assertEquals(parser.strip_comments("test# test    "), "test")
        self.assertEquals(parser.strip_comments("#test# test    "), "")

    def test_parse_reaction_member(self):
        parser = BiooptParser()

        self.assertRaises(TypeError, parser.parse_reaction_member, None)
        self.assertRaises(TypeError, parser.parse_reaction_member, 123)
        self.assertRaises(ValueError, parser.parse_reaction_member, "")

        rm = parser.parse_reaction_member("A")
        self.assertEquals("A", rm.metabolite.name)
        self.assertEquals(1, rm.coefficient)

        rm = parser.parse_reaction_member("2 B")
        self.assertEquals("B", rm.metabolite.name)
        self.assertEquals(2, rm.coefficient)

        rm = parser.parse_reaction_member(".5 C")
        self.assertEquals("C", rm.metabolite.name)
        self.assertEquals(.5, rm.coefficient)

        rm = parser.parse_reaction_member("23e-2 D")
        self.assertEquals("D", rm.metabolite.name)
        self.assertEquals(.23, rm.coefficient)

        rm = parser.parse_reaction_member("23e+2 E")
        self.assertEquals("E", rm.metabolite.name)
        self.assertEquals(2300, rm.coefficient)

        rm = parser.parse_reaction_member("A B")
        self.assertEquals("A B", rm.metabolite.name)
        self.assertEquals(1, rm.coefficient)

        rm = parser.parse_reaction_member("A . B")
        self.assertEquals("A . B", rm.metabolite.name)
        self.assertEquals(1, rm.coefficient)

        rm = parser.parse_reaction_member("A -> B -> C")
        self.assertEquals("A -> B -> C", rm.metabolite.name)
        self.assertEquals(1, rm.coefficient)

        rm = parser.parse_reaction_member("5'end-something")
        self.assertEquals("5'end-something", rm.metabolite.name)
        self.assertEquals(1, rm.coefficient)

        rm = parser.parse_reaction_member("5\"end-something")
        self.assertEquals("5\"end-something", rm.metabolite.name)
        self.assertEquals(1, rm.coefficient)

        rm = parser.parse_reaction_member("A B")
        self.assertEquals("A B", rm.metabolite.name)
        self.assertEquals(1, rm.coefficient)

        rm = parser.parse_reaction_member("2 A . B")
        self.assertEquals("A . B", rm.metabolite.name)
        self.assertEquals(2, rm.coefficient)

        rm = parser.parse_reaction_member("2 A -> B -> C")
        self.assertEquals("A -> B -> C", rm.metabolite.name)
        self.assertEquals(2, rm.coefficient)

        rm = parser.parse_reaction_member("3 5'end-something")
        self.assertEquals("5'end-something", rm.metabolite.name)
        self.assertEquals(3, rm.coefficient)

        rm = parser.parse_reaction_member("4 5\"end-something")
        self.assertEquals("5\"end-something", rm.metabolite.name)
        self.assertEquals(4, rm.coefficient)

        rm = parser.parse_reaction_member("1e-1 A . B")
        self.assertEquals("A . B", rm.metabolite.name)
        self.assertEquals(.1, rm.coefficient)

        rm = parser.parse_reaction_member("2e-1 A -> B -> C")
        self.assertEquals("A -> B -> C", rm.metabolite.name)
        self.assertEquals(.2, rm.coefficient)

        rm = parser.parse_reaction_member("3e-1  5'end-something")
        self.assertEquals("5'end-something", rm.metabolite.name)
        self.assertEquals(.3, rm.coefficient)

        rm = parser.parse_reaction_member("4e-1  5\"end-something")
        self.assertEquals("5\"end-something", rm.metabolite.name)
        self.assertEquals(.4, rm.coefficient)

        rm = parser.parse_reaction_member("A B")
        self.assertEquals("A B", rm.metabolite.name)
        self.assertEquals(1, rm.coefficient)

        rm = parser.parse_reaction_member("(2) A . B")
        self.assertEquals("A . B", rm.metabolite.name)
        self.assertEquals(2, rm.coefficient)

        rm = parser.parse_reaction_member("(2) A -> B -> C")
        self.assertEquals("A -> B -> C", rm.metabolite.name)
        self.assertEquals(2, rm.coefficient)

        rm = parser.parse_reaction_member("(3e+1) 5'end-something")
        self.assertEquals("5'end-something", rm.metabolite.name)
        self.assertEquals(30, rm.coefficient)

        rm = parser.parse_reaction_member("(4e-1) 5\"end-something")
        self.assertEquals("5\"end-something", rm.metabolite.name)
        self.assertEquals(.4, rm.coefficient)

    def test_parse_reaction(self):
        parser = BiooptParser()
        self.assertEquals(parser.parse_reaction(None), None)
        self.assertEquals(parser.parse_reaction(""), None)
        self.assertRaises(TypeError, parser.parse_reaction, 123)

        r_rev = parser.parse_reaction("R1: A + 2.5 B <-> 3 C")
        self.assertEquals(r_rev.direction, Direction.reversible())

        r1 = parser.parse_reaction("R1: A + 2.5 B -> 3 C")
        self.assertEquals(r1.name, "R1")
        self.assertEquals(r1.direction, Direction.forward())
        self.assertEquals(len(r1.reactants), 2)
        self.assertEquals(len(r1.products), 1)
        self.assertEquals(r1.reactants[0].metabolite.name, "A")
        self.assertEquals(r1.reactants[0].coefficient, 1)
        self.assertEquals(r1.reactants[1].metabolite.name, "B")
        self.assertEquals(r1.reactants[1].coefficient, 2.5)
        self.assertEquals(r1.products[0].metabolite.name, "C")
        self.assertEquals(r1.products[0].coefficient, 3)

        r_bkw = parser.parse_reaction("R1: A + .5 B <- 3 C")
        self.assertEquals(r_bkw.direction, Direction.forward())
        self.assertEquals(len(r_bkw.reactants), 1)
        self.assertEquals(len(r_bkw.products), 2)
        self.assertEquals(r_bkw.reactants[0].metabolite.name, "C")
        self.assertEquals(r_bkw.reactants[0].coefficient, 3)
        self.assertEquals(r_bkw.products[0].metabolite.name, "A")
        self.assertEquals(r_bkw.products[0].coefficient, 1)
        self.assertEquals(r_bkw.products[1].metabolite.name, "B")
        self.assertEquals(r_bkw.products[1].coefficient, .5)

#    def test_parse_reaction_section(self):
#        reaction_section = """
#R1: A + B -> 3 C
#R2: B + C <-> 1 E
#"""
#
#        parser = BiooptParser()
#        rl = parser.parse_reactions_section(reaction_section)
#
#        print rl