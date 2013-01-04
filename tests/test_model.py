from unittest import TestCase
from bioopt_parser import *
from model import *
import itertools

class TestBounds(TestCase):
    def test_new(self):
        self.assertRaises(TypeError, Bounds, "a", 1)
        self.assertRaises(TypeError, Bounds, 1, "a")
        self.assertRaises(ValueError, Bounds, 2, 1)

        try:
            Bounds(1, 2)
        except Exception, ex:
            self.fail("Creating Bounds object failed: {0}".format(ex))

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

    def test_setters(self):
        m = Metabolite("H2O")

        self.assertRaises(TypeError, m.__setattr__, "name", 1)
        m.name = "Na"
        self.assertEquals(m.name, "Na")

        self.assertRaises(TypeError, m.__setattr__, "boundary", "true")
        new_boundary = not m.boundary
        self.assertNotEqual(new_boundary, m.boundary)
        m.boundary = new_boundary
        self.assertEquals(m.boundary, new_boundary)

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
            self.fail("Multiplying metabolite by anything else than a number should rise TypeError")

        self.assertEquals(type(2*Metabolite("Na")), ReactionMember)
        self.assertEquals((2*Metabolite("Na")).coefficient, 2)


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

    def test_setters(self):
        m = ReactionMember(Metabolite("Na"), 1)
        H2O = Metabolite("H2O")

        self.assertRaises(TypeError, m.__setattr__, "metabolite", 1)
        m.metabolite = H2O
        self.assertEquals(m.metabolite, H2O)

        self.assertRaises(TypeError, m.__setattr__, "coefficient", "test")
        self.assertRaises(ValueError, m.__setattr__, "coefficient", -1)
        m.coefficient = 2
        self.assertEquals(m.coefficient, 2)

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
        self.assertTrue(all(reactant in r.find_participants() for reactant in reactants))
        self.assertTrue(all(product in r.find_participants() for product in products))
        self.assertEquals(r.direction, direction)
        self.assertEquals(r.bounds, bounds)
        self.assertEquals(r.find_effective_bounds().lb, 0)
        self.assertEquals(r.find_effective_bounds().ub, float("Inf"))

    def test_setters(self):
        name = "Burn"
        bounds = Bounds(-10, float("Inf"))
        direction = Direction.forward()

        m_2_Na = ReactionMember(Metabolite("Na"), 2)
        m_2_H2O = ReactionMember(Metabolite("H2O"), 2)
        m_2_NaOH = ReactionMember(Metabolite("NaOH"), 2)
        m_1_H2 = ReactionMember(Metabolite("H2"), 1)
        reactants = ReactionMemberList([m_2_Na, m_2_H2O])
        products = ReactionMemberList([m_2_NaOH, m_1_H2])

        r = Reaction("")

        self.assertRaises(TypeError, r.__setattr__, "name", 1)
        r.name = name
        self.assertEquals(r.name, name)

        self.assertRaises(TypeError, r.__setattr__, "bounds", 1)
        r.bounds = bounds
        self.assertEquals(r.bounds, bounds)

        self.assertRaises(TypeError, r.__setattr__, "direction", 1)
        r.direction = direction
        self.assertEquals(r.direction, direction)

        self.assertRaises(TypeError, r.__setattr__, "reactants", 1)
        r.reactants = reactants
        self.assertEquals(r.reactants, reactants)

        self.assertRaises(TypeError, r.__setattr__, "products", 1)
        r.products = products
        self.assertEquals(r.products, products)

class TestOperation(TestCase):
    def test_equality(self):
        self.assertTrue(Operation.multiplication() is Operation.multiplication())
        self.assertTrue(Operation.addition() is Operation.addition())
        self.assertTrue(Operation.division() is Operation.division())
        self.assertTrue(Operation.subtraction() is Operation.subtraction())
        self.assertNotEquals(Operation.multiplication(), Operation.addition())
        self.assertNotEquals(Operation.multiplication(), Operation.division())
        self.assertNotEquals(Operation.multiplication(), Operation.subtraction())
        self.assertNotEquals(Operation.addition(), Operation.division())
        self.assertNotEquals(Operation.addition(), Operation.subtraction())
        self.assertNotEquals(Operation.division(), Operation.subtraction())

class TestArithmetic(TestCase):
    def test_reaction_member(self):
        self.assertEquals((2*Metabolite("Na")).coefficient, 2)

class TestMathExpression(TestCase):
    def test_new(self):
        self.assertRaises(TypeError, MathExpression, "a", [1])
        self.assertRaises(ValueError, MathExpression, Operation.addition(), [])
        self.assertRaises(ValueError, MathExpression, Operation.addition(), [1])
        self.assertRaises(ValueError, MathExpression, Operation.addition(), [])

        try:
            ex = MathExpression(Operation.addition(), [Reaction("R1"), 2])
            MathExpression(Operation.negation(), [True])
        except Exception, ex:
            self.fail("Creating <MathExpression> object failed object failed: {0}".format(ex))

        self.assertEquals(Operation.addition(), ex.operation)
        self.assertEquals([Reaction("R1"), 2], ex.operands)


    def test_setters(self):
        ex = MathExpression(Operation.addition(), [Reaction("R1"), 2])
        self.assertRaises(TypeError, ex.__setattr__, "operation", "a")
        ex.operation = Operation.division()
        self.assertEquals(Operation.division(), ex.operation)
        self.assertRaises(ValueError, ex.__setattr__, "operands", [10])


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
        r2 = Reaction("L-Methionine methanethiol-lyase", 1*H2O + 1*Homocysteine, 1*NH3_H4N + 1*Oxobutyrate + 1*H2S)

        model = Model()
        model.reactions.append(r1)
        model.reactions.append(r2)

        #print ReactionMember(Homocysteine, 1) + ReactionMember(Oxobutyrate, 1)
        #print model