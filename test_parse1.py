from unittest import TestCase
from bioopt import *
from bioopt import Metabolite as M

def power(a, b):
    return a**b

#class TestParse(TestCase):
#    def test_parse3(self):
#        model = """
#-REACTIONS
#R1: A + B -> 3 c
#-CONSTRAINTS
#R1[0, 100]
#-EXTERNAL METABOLITES
#C
#-OBJ
#R1 1 1
#-DESIGNOBJ
#R1 1 1
#"""
#        r = parse(model)

#class TestMetabolite(TestCase):
#    def test_new(self):
#        m = Metabolite("H2O", 3)
#
#        self.assertTrue(m.name == "H2O")
#        self.assertTrue(m.coefficient == 3)
#
#    def test_invalid_coef(self):
#        with self.assertRaises(TypeError, message="Providing not string name should rise TypeError"):
#            Metabolite(15, "a")
#
#        with self.assertRaises(TypeError, message="Providing character coefficient should raise TypeError"):
#            Metabolite("H2O", "a")
#
#        with self.assertRaises(ValueError, message="Providing zero or negative coefficient coefficient should rise ValueError"):
#            Metabolite("H2O", 0)
#
#        with self.assertRaises(ValueError, message="Providing zero or negative coefficient coefficient should rise ValueError"):
#            Metabolite("H2O", -1)

class TestBounds(TestCase):
    def test_new(self):
        self.assertRaises(TypeError, Bounds, "a", 1)
        self.assertRaises(TypeError, Bounds, 1, "a")
        self.assertRaises(ValueError, Bounds, 2, 1)

        try:
            Bounds(1, 2)
        except:
            self.fail("Creating Bounds object failed")

        self.assertEquals(Bounds(1, 2).lb, 1)
        self.assertEquals(Bounds(1, 2).ub, 2)

    def test_setters(self):
        b = Bounds(0, 0)
        self.assertRaises(TypeError, b.lb, "a")
        self.assertRaises(TypeError, b.ub, "a")
        self.assertRaises(ValueError, b.__setattr__, "ub", -1)
        self.assertRaises(ValueError, b.__setattr__, "lb", 1)
        b.lb = -1
        self.assertEquals(b.lb, -1)
        b.ub = 2
        self.assertEquals(b.ub, 2)

class TestMetabolite(TestCase):
    def test_new(self):
        self.assertRaises(TypeError, Metabolite, 1)
        self.assertRaises(TypeError, Metabolite, None)
        self.assertRaises(TypeError, Metabolite, "name", 2)
        self.assertRaises(ValueError, Metabolite, "")

        try:
            m = Metabolite("H2O", True)
        except:
            self.fail("Creating Metabolite object failed")

        self.assertEquals(m.name, "H2O")
        self.assertEquals(m.boundary, True)

    def test_arithmetics(self):
        try:
            -1*Metabolite("Na")
            self.fail("Multiplying metabolite by negative coefficient should rise ValueError")
        except ValueError:
            pass
        except Exception:
            self.fail("Multiplying metabolite by negative coefficient should rise ValueError")

        try:
            "test"*Metabolite("Na")
            self.fail("Multiplying metabolite by anything else than a number should rise TypeError")
        except TypeError:
            pass
        except Exception, ex:
            print ex
            self.fail("Multiplying metabolite by anything else than a number should rise TypeError")

        self.assertEquals(type(2*M("Na")), ReactionMember)
        self.assertEquals((2*M("Na")).coefficient, 2)


class TestReactionMember(TestCase):
    def test_new(self):
        m = Metabolite("H2O")
        self.assertRaises(TypeError, ReactionMember, m, m)
        self.assertRaises(TypeError, ReactionMember, m, None)
        self.assertRaises(ValueError, ReactionMember, m, 0)
        self.assertRaises(ValueError, ReactionMember, m, -1)

        try:
            rm = ReactionMember(m, 2)
        except:
            self.fail("Creating ReactionMember object failed")

        self.assertEquals(rm.coefficient, 2)

    def test_arithmetics(self):
        Na = 1*Metabolite("Na")
        H2O = 2*Metabolite("H2O")

        try:
            Na + 2
            self.fail("Adding number to ReactionMember should raise an error")
        except TypeError:
            pass
        except Exception:
            self.fail("Adding number to ReactionMember should raise an error")

        try:
            2 + Na
            self.fail("Adding number to ReactionMember should raise an error")
        except TypeError:
            pass
        except Exception:
            self.fail("Adding number to ReactionMember should raise an error")


        self.assertEquals(type(Na + H2O), ReactionMemberList)
        self.assertEquals((Na + H2O)[0], Na)
        self.assertEquals((Na + H2O)[1], H2O)

        reactants = Na + H2O
        self.assertEquals(len(reactants + Na), 3)
        self.assertEquals(len(reactants + H2O), 3)
        self.assertEquals(len(Na + reactants), 3)
        self.assertEquals(len(H2O + reactants), 3)

        self.assertTrue((H2O + reactants)[0] is H2O)
        self.assertTrue((reactants + Na)[0] is Na)

class TestReaction(TestCase):
    def test_new(self):
        name = "Burn"
        bounds = Bounds(-10, float("Inf"))
        direction = Direction.forward()

        m_2_Na = ReactionMember(Metabolite("Na"), 2)
        m_2_H2O = ReactionMember(Metabolite("H2O"), 2)
        m_2_NaOH = ReactionMember(Metabolite("NaOH"), 2)
        m_1_H2 = ReactionMember(Metabolite("H2"), 1)
        reactants = ReactionMemberList([m_2_Na, m_2_H2O])
        products = ReactionMemberList([m_2_NaOH, m_1_H2])

        self.assertRaises(TypeError, Reaction, 1, reactants, products, direction, bounds)
        self.assertRaises(TypeError, Reaction, name, 1, products, direction, bounds)
        self.assertRaises(TypeError, Reaction, name, reactants, 1, direction, bounds)
        self.assertRaises(TypeError, Reaction, name, reactants, products, 1, bounds)
        self.assertRaises(TypeError, Reaction, name, reactants, products, direction, 1)

        try:
            r = Reaction(name, reactants, products, direction, bounds)
        except Exception, ex:
            self.fail("Creating Reaction object failed {0}".format(ex))

        self.assertEquals(r.name, name)
        self.assertEquals(r.reactants, reactants)
        self.assertEquals(r.products, products)
        self.assertEquals(r.direction, direction)
        self.assertEquals(r.bounds, bounds)
        self.assertEquals(r.find_effective_bounds().lb, 0)
        self.assertEquals(r.find_effective_bounds().ub, float("Inf"))

class TestArithmetic(TestCase):
    def test_reaction_member(self):
        self.assertEquals((2*Metabolite("Na")).coefficient, 2)



class TestModel(TestCase):
    def test_new(self):
        Na = Metabolite("Na")
        H2O = Metabolite("H2O", True)
        NaOH = Metabolite("NaOH")
        H2 = Metabolite("H2")

        NH3_H4N = Metabolite("NH3_H4N")
        Homocysteine = Metabolite("Homocysteine")
        Oxobutyrate = Metabolite("2-Oxobutyrate")
        H2S = Metabolite("H2S")

        r1 = Reaction("Burn", 2*Na + 2*H2O, 2*NaOH + 1*H2, Direction.forward(), Bounds())
        r2 = Reaction("L-Methionine methanethiol-lyase", 1*H2O + 1*Homocysteine, 1*NH3_H4N + 1*Oxobutyrate + 1*H2S, Direction.reversible(), Bounds())

        model = Model()
        model.reactions.append(r1)
        model.reactions.append(r2)

        #print ReactionMember(Homocysteine, 1) + ReactionMember(Oxobutyrate, 1)
        #print model

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