from __future__ import unicode_literals, print_function
import pypeg2 as p
from pypeg2 import Whitespace as ws
from pypeg2 import optional as opt
from pypeg2 import attr
import re

re_comment_short = re.compile(r"\#.*")
re_comment_long = re.compile(r"(?m)/%.*?%/")
re_capital = re.compile(r"(?m)[a-z]")
re_number_str = r"(?:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)|(?:[-+]?(?:[0-9]*\.[0-9]+|[0-9]+))"

class BiooptWhitespace(p.Symbol):
    regex = re.compile("[ \t]+")
ows = p.omit(opt(BiooptWhitespace))

nl = p.omit(["\r", "\n", "\r\n"])


class Plus(p.Keyword):
    regex = re.compile(r"[ \t]*\+[ \t]*")

class Direction(p.Keyword):
    regex = re.compile(r"<?->")

class Restline(p.Symbol):
    regex = re.compile(r".*?(?=(\n|\r|$))")

class ReactionSeparator(p.Symbol):
    regex = re.compile(r"[ \t]*:[ \t]*")



class Metabolite(p.Symbol):
    regex = re.compile(r".+?(?=[ \t]*(\+|->|<->|\n|\r|$))")


class Number(p.Symbol):
    regex = re.compile(re_number_str)


class ReactionMember(p.List):
    grammar = opt(attr("coefficient", Number), ws), attr("metabolite", Metabolite), ows

    def __repr__(self):
        return "('{coef}', '{met}')".format(coef=self.__dict__.get("coefficient", 1), met=self.metabolite)


class ReactionMemberList(p.List):
    grammar = opt(ReactionMember), p.maybe_some(p.omit(Plus), ReactionMember)

    def __repr__(self):
        if len(self) > 0:
            return "[" + "] + [".join(re_capital.sub("", str(type(m).__name__))+":"+str(m) for m in self) + "]"
        else:
            return ""


class ReactionName(p.Symbol):
    regex = re.compile(".*?(?=([ \t]*[:\[])|$)")


class Reaction(p.List):
    grammar = attr("name", ReactionName), p.omit(ReactionSeparator), \
              attr("lhs", ReactionMemberList), ows, \
              attr("direction", Direction), ows, \
              attr("rhs", ReactionMemberList), p.omit(Restline)

    def __repr__(self):
        return "{r} : {lhs} {dir} {rhs}".format(r=self.name, lhs=self.lhs, rhs=self.rhs, dir=self.direction)





class ConstraintBoundaries(p.List):
    grammar = "[", ows, attr("lb", Number), ows, ",", ows, attr("ub", Number), ows, "]"

    def __repr__(self):
        return "[{0}, {1}]".format(self.lb, self.ub)

class Constraint(p.List):
    grammar = attr("reaction", ReactionName), ows, attr("boundaries", ConstraintBoundaries)

    def __repr__(self):
        return "{0} {1}".format(self.reaction, self.boundaries)


class Section(p.List):
    def __repr__(self):
        items = "\n".join(str(r) for r in self.__dict__.get("items", []))
        return "[{type}]{section}:\n{items}".format(type=type(self).__name__, section=self.name, items=items)


class ReactionsSection(Section):
    grammar = "-", p.attr("name", Restline), nl, attr("items", p.some([nl, Reaction]))


class ConstraintsSection(Section):
    grammar = "-", p.attr("name", Restline), nl, attr("items", p.some([nl, Constraint]))


class ExternalSection(Section):
    grammar = "-", p.attr("name", Restline), nl, attr("items", p.some([nl, Metabolite]))


class ObjectiveName(p.Symbol):
    regex = re.compile(".*?(?=[ \t]+{0}[ \t]+{1})".format(re_number_str, re_number_str))

class Objective(p.List):
    grammar = p.attr("reaction", ObjectiveName), p.some(ws), \
              p.attr("coef1", Number), p.some(ws), p.attr("coef2", Number)

    def __repr__(self):
        return "{0} ('{1}', '{2}')".format(self.reaction, self.coef1, self.coef2)

class ObjectiveSection(Section):
    grammar = "-", p.attr("name", Restline), nl, attr("items", p.some([nl, Objective]))



class Bioopt(p.List):
    grammar = attr("sections", p.maybe_some([nl, ReactionsSection, ConstraintsSection, ExternalSection, ObjectiveSection]))

    def __repr__(self):
        return "\n\n".join(repr(s) for s in self.__dict__.get("sections", []))
# Test

f = p.parse("""-REACTIONS 1
R1: A + B -> C
R2: F -> G
R3: A -> B


-REACTIONS 2
R4: A -> B
R5: D -> F
R6: D -> F

-CONSTRAINTS
U214_ [1, 1]
ZYMSTxtI [0, 0]
XANxtI [0, 0]
URIxtI [0, 0]
UREAxtI [0, 0]

-EXTERNAL METABOLITES
acetaldehyde[b]
acetate[b]
adenosine[b]
adenine[b]
2-oxoglutarate[b]

-OBJ
VGRO 1 1
EEEEEEEE 1 1
""", Bioopt, whitespace=None)
print(f)




