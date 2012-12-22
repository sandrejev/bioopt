from unittest import TestCase
from bioopt_parser import *
from model import Metabolite as M
from model import Reaction as R
from model import ReactionMemberList as RML

class TestBiooptParser(TestCase):
    def setUp(self):
        self.fwd = Direction.forward()
        self.rev = Direction.reversible()

#    model = """-REACTIONS
#R1: A + B -> 3 C
#R2: B + C <-> 1 E
#-CONSTRAINTS
#R1[-100, 100]
#R2[-100, 100]
#-EXTERNAL METABOLITES
#E
#-OBJ
#R2 1 1
#-DESIGNOBJ
#R2 1 1
#"""
    def test_strip_comments(self):
        parser = BiooptParser()
        self.assertEquals(parser.strip_comments("test"), "test")
        self.assertEquals(parser.strip_comments("    test"), "test")
        self.assertEquals(parser.strip_comments("test    "), "test")
        self.assertEquals(parser.strip_comments("test# test    "), "test")
        self.assertEquals(parser.strip_comments("test # test    "), "test")
        self.assertEquals(parser.strip_comments("#test# test    "), "")

    def test_parse_reaction_member(self):
        parser = BiooptParser()

        self.assertRaises(TypeError, parser.parse_reaction_member, None)
        self.assertRaises(TypeError, parser.parse_reaction_member, 123)
        self.assertRaises(ValueError, parser.parse_reaction_member, "")

        self.assertEquals(1*M("A"), parser.parse_reaction_member("A"))
        self.assertEquals(2*M("B"), parser.parse_reaction_member("2 B"))
        self.assertEquals(.5*M("C"), parser.parse_reaction_member(".5 C"))

        self.assertEquals(.23*M("D"), parser.parse_reaction_member("23e-2 D"))
        self.assertEquals(2300*M("E"), parser.parse_reaction_member("23e+2 E"))

        self.assertEquals(1*M("A B"), parser.parse_reaction_member("A B"))
        self.assertEquals(1*M("A . B"), parser.parse_reaction_member("A . B"))
        self.assertEquals(1*M("A -> B -> C"), parser.parse_reaction_member("A -> B -> C"))
        self.assertEquals(1*M("5'end-something"), parser.parse_reaction_member("5'end-something"))

        self.assertEquals(1*M("5\"end-something"), parser.parse_reaction_member("5\"end-something"))
        self.assertEquals(1*M("A B"), parser.parse_reaction_member("A B"))
        self.assertEquals(2*M("A . B"), parser.parse_reaction_member("2 A . B"))

        self.assertEquals(2*M("A -> B -> C"), parser.parse_reaction_member("2 A -> B -> C"))
        self.assertEquals(3*M("5'end-something"), parser.parse_reaction_member("3 5'end-something"))
        self.assertEquals(4*M("5\"end-something"), parser.parse_reaction_member("4 5\"end-something"))

        self.assertEquals(.1*M("A . B"), parser.parse_reaction_member("1e-1 A . B"))
        self.assertEquals(.2*M("A -> B -> C"), parser.parse_reaction_member("2e-1 A -> B -> C"))
        self.assertEquals(.3*M("5'end-something"), parser.parse_reaction_member("3e-1  5'end-something"))
        self.assertEquals(.4*M("5\"end-something"), parser.parse_reaction_member("4e-1  5\"end-something"))

        self.assertEquals(1*M("A B"), parser.parse_reaction_member("A B"))
        self.assertEquals(2*M("A . B"), parser.parse_reaction_member("(2) A . B"))
        self.assertEquals(2*M("A -> B -> C"), parser.parse_reaction_member("(2) A -> B -> C"))
        self.assertEquals(30*M("5'end-something"), parser.parse_reaction_member("(3e+1) 5'end-something"))
        self.assertEquals(.4*M("5\"end-something"), parser.parse_reaction_member("(4e-1) 5\"end-something"))
        self.assertEquals(1*M("3e- A"), parser.parse_reaction_member("3e- A"))

    def test_parse_reaction_member_list(self):
        parser = BiooptParser()

        self.assertRaises(TypeError, parser.parse_reaction_member_list, None)
        self.assertRaises(TypeError, parser.parse_reaction_member_list, 123)
        self.assertRaises(ValueError, parser.parse_reaction_member_list, "")

        self.assertEquals(1*M("A") + 2.5*M("B"), parser.parse_reaction_member_list("A + 2.5 B"))
        self.assertEquals(1*M("A") + 25*M("5'B") + .03*M("C,1,1"), parser.parse_reaction_member_list("A + 2.5e+1 5'B + 3e-2 C,1,1"))

    def test_parse_reaction(self):
        parser = BiooptParser()
        self.assertRaises(TypeError, parser.parse_reaction, None)
        self.assertRaises(ValueError, parser.parse_reaction, "")
        self.assertRaises(TypeError, parser.parse_reaction, 123)

        r_rev = parser.parse_reaction("R1: A + 2.5 B <-> 3 C")
        self.assertEquals(r_rev.direction, self.rev)

        r1 = parser.parse_reaction("R1: A + 2.5 B -> 3 C")
        self.assertEquals(r1.name, "R1")
        self.assertEquals(r1.direction, self.fwd)
        self.assertEquals(1*M("A") + 2.5*M("B"), r1.reactants)
        self.assertEquals(RML([3*M("C")]), r1.products)

        r_bkw = parser.parse_reaction("R1: A + .5 B <- 3 C")
        self.assertEquals(r_bkw.direction, self.fwd)
        self.assertEquals(RML([3*M("C")]), r_bkw.reactants)
        self.assertEquals(1*M("A") + .5*M("B"), r_bkw.products)

    def test_parse_reaction_section(self):
        reaction_section = """
            # Reaction 1
            R1: A + B -> 3 C
            # Reaction 2
            R2: B + C <-> 1 E # reaction 2 extra
            """

        parser = BiooptParser()
        rl = parser.parse_reactions_section(reaction_section)
        self.assertEquals(Reaction("R1", 1*M("A") + 1*M("B"), 3*M("C"), self.fwd), rl[0])
        self.assertEquals(Reaction("R2", 1*M("B") + 1*M("C"), 1*M("E"), self.rev), rl[1])

    def test_parse_constraint(self):
        parser = BiooptParser()
        self.assertEquals(R("A", bounds=Bounds(-1, 1)), parser.parse_constraint("A [-1,1] # comment"))
        self.assertEquals(R("A", bounds=Bounds(-1, 1)), parser.parse_constraint("A[-1,1]"))
        self.assertEquals(R("B-10", bounds=Bounds(-1, 1)), parser.parse_constraint("B-10 [-1, 1]"))
        self.assertEquals(R("10'C", bounds=Bounds(-1, 1)), parser.parse_constraint("10'C[-1, 1]"))
        self.assertEquals(R("A", bounds=Bounds(-1, 1)), parser.parse_constraint("A\t[\t-1\t, \t1\t]\t"))
        self.assertEquals(R("A", bounds=Bounds(-1, 1)), parser.parse_constraint("A\t[ -1, 1 ]"))
        self.assertEquals(R("3\"-something[e]", bounds=Bounds(1, 10)), parser.parse_constraint("3\"-something[e] [1, 10]"))
        self.assertEquals(R("H2O[e]", bounds=Bounds(-1.4, 2500)), parser.parse_constraint("H2O[e] [-14e-1, 25e+2]"))

    def test_parse_constraints_section(self):
        constraints_section = """
            # Constraint 1
            R1 [-1000, 1000]
            # Constraint 2
            R2 [-1e+2, 1e+2]
            """

        parser = BiooptParser()
        cons = parser.parse_constraints_section(constraints_section)
        self.assertEquals(Reaction("R1", bounds=Bounds(-1000, 1000)), cons[0])
        self.assertEquals(Reaction("R2", bounds=Bounds(-100, 100)), cons[1])

    def test_parse_external_metabolites_section(self):
        constraints_section = """
            A
            'A
            #comment
            A B  # test
            """

        parser = BiooptParser()
        self.assertEquals([M("A"), M("'A"), M("A B")], parser.parse_external_metabolites_section(constraints_section))

