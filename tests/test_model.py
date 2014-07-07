from unittest import TestCase
from bioopt_parser import *
from model import Metabolite as M
from model import Reaction as R
from model import Bounds as B
from model import ReactionMemberList as RML
from model import MathExpression as ME

warnings.simplefilter("ignore")

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
        bounds = Bounds(-10, Bounds.inf())
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
        self.assertEquals(r.find_effective_bounds().ub, Bounds.inf())

    def test_setters(self):
        name = "Burn"
        bounds = Bounds(-10, Bounds.inf())
        direction = Direction.forward()

        m_2_Na = ReactionMember(Metabolite("Na"), 2)
        m_2_H2O = ReactionMember(Metabolite("H2O"), 2)
        m_2_NaOH = ReactionMember(Metabolite("NaOH"), 2)
        m_1_H2 = ReactionMember(Metabolite("H2"), 1)
        reactants = ReactionMemberList([m_2_Na, m_2_H2O])
        products = ReactionMemberList([m_2_NaOH, m_1_H2])

        r = Reaction("whatsoever")

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

    def test_reverse(self):
        fwd = Direction.forward()
        rev = Direction.reversible()

        m_2_Na = ReactionMember(Metabolite("Na"), 2)
        m_2_H2O = ReactionMember(Metabolite("H2O"), 2)
        m_2_NaOH = ReactionMember(Metabolite("NaOH"), 2)
        m_1_H2 = ReactionMember(Metabolite("H2"), 1)

        r = R("r", m_2_Na + m_2_H2O, m_2_NaOH + m_1_H2, direction=fwd, bounds=Bounds(-10, Bounds.inf()))
        self.assertRaises(RuntimeError, r.reverse)

        r = R("r", m_2_Na + m_2_H2O, m_2_NaOH + m_1_H2, direction=fwd, bounds=Bounds(1, 10))
        self.assertRaises(RuntimeError, r.reverse)

        r = R("r", m_2_Na + m_2_H2O, m_2_NaOH + m_1_H2, direction=rev, bounds=Bounds(1, 10))
        self.assertRaises(RuntimeError, r.reverse)

        products = m_2_Na + m_2_H2O
        reactants = m_2_NaOH + m_1_H2
        r = R("r", reactants, products, direction=rev, bounds=Bounds(-7, 8))
        try:
            r.reverse()
        except Exception, ex:
            self.fail("Creating Bounds object failed: {0}".format(ex))
        self.assertEquals(-8, r.bounds.lb)
        self.assertEquals(7, r.bounds.ub)
        self.assertEquals(reactants, r.products)
        self.assertEquals(products, r.reactants)


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

    def test_variables(self):
        mult = Operation.multiplication()
        add = Operation.addition()

        r1 = R("R1")
        r2 = R("R2")
        ex = ME(add, [ME(mult, [r2, r1, 1]), ME(mult, [r2, 1]), ME(mult, [r1, 2])])
        self.assertEquals([r2, r1, 1, 2], ex.find_variables())

    def test_representation(self):
        # TODO: test representation
        #print ME(Operation.addition(), [ME(mult, [R("R2"), R("R1"), float(1)]), ME(mult, [R("R2"), R("R1"), float(1)]), ME(mult, [R("R2"), R("R1"), float(1)])]).__repr__(tree=True)
        #print ME(mult, [ME(add, [R("R2"), R("R1"), float(1)]), ME(add, [R("R2"), R("R1"), float(1)]), ME(add, [R("R2"), R("R1"), float(1)])]).__repr__()
        pass

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

        # TODO: finish this test
        #print ReactionMember(Homocysteine, 1) + ReactionMember(Oxobutyrate, 1)
        #print model

    def test_references(self):
        model = Model()
        r1 = R("R1", 1*M("A") + 1*M("B"), 3*M("C"), direction=Direction.forward(), bounds=Bounds(-100, 100))
        r2 = R("R2", 1*M("B") + 1*M("C"), 1*M("E", boundary=True), direction=Direction.reversible(), bounds=Bounds(-100, 100))
        model.reactions = [r1, r2]
        model.objective = ME(Operation.multiplication(), [R("R1"), 1, 1])
        model.design_objective = ME(Operation.multiplication(), [R("R2"), R("R1"), 1])

        self.assertFalse(model.objective.operands[0] is r1)
        self.assertFalse(model.design_objective.operands[0] is r2)
        self.assertFalse(model.design_objective.operands[1] is r1)
        self.assertFalse(r1.reactants[1].metabolite is r2.reactants[0].metabolite)

        model.unify_references()

        self.assertTrue(model.objective.operands[0] is r1)
        self.assertTrue(model.design_objective.operands[0] is r2)
        self.assertTrue(model.design_objective.operands[1] is r1)
        self.assertTrue(r1.reactants[1].metabolite is r2.reactants[0].metabolite)

        self.assertFalse(model.objective.operands[1] is r1)
        self.assertFalse(r1.reactants[1].metabolite is r2.reactants[1].metabolite)

    def test_commune(self):
        fwd = Direction.forward()
        rev = Direction.reversible()

        model1 = Model()
        r1 = R("R1", 1*M("A") + 1*M("B"), 3*M("C"), direction=fwd, bounds=Bounds(-100, 100))
        r2 = R("R2", 1*M("B") + 1*M("C"), 1*M("E", boundary=True), direction=Direction.reversible(), bounds=Bounds(-100, 100))
        model1.reactions = [r1, r2]
        model2 = Model()
        r3 = R("R1", 1*M("A"), 3*M("B"), direction=fwd, bounds=Bounds(-100, 100))
        r4 = R("R1", 1*M("A"), 3*M("B"), direction=fwd, bounds=Bounds(-100, 100))
        r5 = R("R1", 1*M("C"), 3*M("B"), direction=fwd, bounds=Bounds(-100, 100))
        r6 = R("R2", 1*M("B") + 1*M("C"), 1*M("E", boundary=True), direction=Direction.reversible(), bounds=Bounds(-100, 100))
        model2.reactions = [r3, r4, r5, r6]
        model3 = Model()
        r7 = R("R1", 1*M("A") + 1*M("B"), 3*M("C"), direction=fwd, bounds=Bounds(-100, 100))
        r8 = R("R2", 1*M("B") + 1*M("C"), 1*M("E", boundary=True), direction=Direction.reversible(), bounds=Bounds(-100, 100))
        model3.reactions = [r7, r8]

        com_model = Model.commune([model1, model2, model3])

        com_model_true = Model()
        com_model_true.reactions = [
            R("ML0000_R1", 1*M("ML0000_A") + 1*M("ML0000_B"), 3*M("ML0000_C"), direction=fwd, bounds=B(-100, 100)),
            R("ML0000_R2", 1*M("ML0000_B") + 1*M("ML0000_C"), 1*M("ML0000_E"), direction=rev, bounds=B(-100, 100)),
            R("ML0000_OUT_E", 1*M("ML0000_E"), 1*M("ENV_E"), direction=fwd, bounds=B(0, B.inf())),
            R("ML0000_IN_E", 1*M("ENV_E"), 1*M("ML0000_E"), direction=fwd, bounds=B(0, B.inf())),

            R("ML0001_R1", 1*M("ML0001_A"), 3*M("ML0001_B"), direction=fwd, bounds=Bounds(-100, 100)),
            R("ML0001_R1", 1*M("ML0001_A"), 3*M("ML0001_B"), direction=fwd, bounds=Bounds(-100, 100)),
            R("ML0001_R1", 1*M("ML0001_C"), 3*M("ML0001_B"), direction=fwd, bounds=Bounds(-100, 100)),
            R("ML0001_R2", 1*M("ML0001_B") + 1*M("ML0001_C"), 1*M("ML0001_E"), direction=rev, bounds=B(-100, 100)),
            R("ML0001_OUT_E", 1*M("ML0001_E"), 1*M("ENV_E"), direction=fwd, bounds=B(0, B.inf())),
            R("ML0001_IN_E", 1*M("ENV_E"), 1*M("ML0001_E"), direction=fwd, bounds=B(0, B.inf())),

            R("ML0002_R1", 1*M("ML0002_A") + 1*M("ML0002_B"), 3*M("ML0002_C"), direction=fwd, bounds=B(-100, 100)),
            R("ML0002_R2", 1*M("ML0002_B") + 1*M("ML0002_C"), 1*M("ML0002_E"), direction=rev, bounds=B(-100, 100)),
            R("ML0002_OUT_E", 1*M("ML0002_E"), 1*M("ENV_E"), direction=fwd, bounds=B(0, B.inf())),
            R("ML0002_IN_E", 1*M("ENV_E"), 1*M("ML0002_E"), direction=fwd, bounds=B(0, B.inf())),

            R("ENV_ExtO", 1*M("ENV_E"), 1*M("ENV_ExtX", boundary=True), direction=fwd, bounds=B(0, B.inf())),
            R("ENV_ExtI", 1*M("ENV_ExtX", boundary=True), 1*M("ENV_E"), direction=fwd, bounds=B(0, B.inf())),
        ]
        com_model_true.unify_references()

        self.assertEquals(com_model, com_model_true)


    def test_save(self):
        model = Model()
        r1 = R("R1", 1*M("A") + 1*M("B"), 3*M("C"), direction=Direction.forward(), bounds=B(-100, 100))
        r2 = R("R2", 1*M("B") + 1*M("C"), 1*M("E", boundary=True), direction=Direction.reversible())
        model.reactions = [r1, r2]
        model.objective = ME(Operation.multiplication(), [R("R1"), 1, 1])
        model.design_objective = ME(Operation.multiplication(), [R("R2"), R("R1"), 1])

        expected = """-REACTIONS
R1	:	A + B -> 3 C
R2	:	B + C <-> E

-CONSTRAINTS
R1	[-100, 100]
R2	[-1000, 1000]

-EXTERNAL METABOLITES
E

-OBJECTIVE
R1 1 1

-DESIGN OBJECTIVE
R2 R1 1

"""
        self.assertEquals(expected, model.save())

    def test_sbml(self):
        fwd = Direction.forward()
        rev = Direction.reversible()

        a = M("A_c")
        b = M("B_c")
        c = M("C_e")
        e = M("E_e", boundary=True)

        model = Model()
        model.reactions = [
            R("R1", 1*a, 3*b, direction=fwd, bounds=Bounds(-100, 100)),
            R("R2", 1*c, 3*b, direction=fwd, bounds=Bounds(-100, 100)),
            R("R3", 1*b + 1*c, 1*e, direction=rev, bounds=Bounds(-100, 100))]

        import libsbml
        sbml = model.sbml()
        sbml_export = libsbml.writeSBMLToString(sbml)

        doc = libsbml.readSBMLFromString(sbml_export)
        self.assertEquals(0, doc.getNumErrors())

        exp_export = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version3" level="2" version="3">
  <model>
    <listOfUnitDefinitions>
      <unitDefinition id="mmol_per_gDW_per_hr">
        <listOfUnits>
          <unit kind="mole" scale="-3"/>
          <unit kind="second" exponent="-1" multiplier="0.00027778"/>
        </listOfUnits>
      </unitDefinition>
    </listOfUnitDefinitions>
    <listOfCompartments>
      <compartment id="C_0000" name="c"/>
      <compartment id="C_0001" name="e"/>
    </listOfCompartments>
    <listOfSpecies>
      <species id="M_0001" name="A_c" compartment="C_0000" initialAmount="0" boundaryCondition="false"/>
      <species id="M_0002" name="B_c" compartment="C_0000" initialAmount="0" boundaryCondition="false"/>
      <species id="M_0003" name="C_e" compartment="C_0001" initialAmount="0" boundaryCondition="false"/>
      <species id="M_0004" name="E_e" compartment="C_0001" initialAmount="0" boundaryCondition="true"/>
    </listOfSpecies>
    <listOfReactions>
      <reaction id="R_0001" name="R1" reversible="false">
        <listOfReactants>
          <speciesReference species="M_0001" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="M_0002" stoichiometry="3"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <ci> FLUX_VALUE </ci>
          </math>
          <listOfParameters>
            <parameter id="R_0001_LB" name="LOWER_BOUND" value="-100" units="mmol_per_gDW_per_hr"/>
            <parameter id="R_0001_UB" name="UPPER_BOUND" value="100" units="mmol_per_gDW_per_hr"/>
            <parameter id="R_0001_OBJ" name="OBJECTIVE_COEFFICIENT" value="0" units="dimensionless"/>
            <parameter id="FLUX_VALUE" name="FLUX_VALUE" value="0" units="mmol_per_gDW_per_hr"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R_0002" name="R2" reversible="false">
        <listOfReactants>
          <speciesReference species="M_0003" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="M_0002" stoichiometry="3"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <ci> FLUX_VALUE </ci>
          </math>
          <listOfParameters>
            <parameter id="R_0002_LB" name="LOWER_BOUND" value="-100" units="mmol_per_gDW_per_hr"/>
            <parameter id="R_0002_UB" name="UPPER_BOUND" value="100" units="mmol_per_gDW_per_hr"/>
            <parameter id="R_0002_OBJ" name="OBJECTIVE_COEFFICIENT" value="0" units="dimensionless"/>
            <parameter id="FLUX_VALUE" name="FLUX_VALUE" value="0" units="mmol_per_gDW_per_hr"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R_0003" name="R3" reversible="true">
        <listOfReactants>
          <speciesReference species="M_0002" stoichiometry="1"/>
          <speciesReference species="M_0003" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="M_0004" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <ci> FLUX_VALUE </ci>
          </math>
          <listOfParameters>
            <parameter id="R_0003_LB" name="LOWER_BOUND" value="-100" units="mmol_per_gDW_per_hr"/>
            <parameter id="R_0003_UB" name="UPPER_BOUND" value="100" units="mmol_per_gDW_per_hr"/>
            <parameter id="R_0003_OBJ" name="OBJECTIVE_COEFFICIENT" value="0" units="dimensionless"/>
            <parameter id="FLUX_VALUE" name="FLUX_VALUE" value="0" units="mmol_per_gDW_per_hr"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>
"""


        self.assertEquals(exp_export, sbml_export)
