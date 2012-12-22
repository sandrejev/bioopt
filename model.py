#!/usr/bin/env python
import sys,re
import thread

class Bounds(object):
    def __init__(self, lb=float("-inf"), ub=float("inf")):
        self.__assert_valid(lb, ub)
        self.__lb = lb
        self.__ub = ub

    @property
    def lb(self):
        return self.__lb

    @lb.setter
    def lb(self, lb):
        self.__assert_valid(lb, self.ub)
        self.__lb = float(lb)

    @property
    def ub(self):
        return self.__ub

    @ub.setter
    def ub(self, value):
        self.__assert_valid(self.lb, value)
        self.__ub = float(value)

    @property
    def direction(self):
        return Direction.forward() if self.lb >= 0 else Direction.reversible()

    def __assert_valid(self, lb, ub):
        if not isinstance(ub, (int, float)):
            raise TypeError("Upper bound is not a number")

        if not isinstance(ub, (int, float)):
            raise TypeError("Upper bound is not a number")

        if not isinstance(lb, (int, float)):
            raise TypeError("Lower bound is not a number")

        if lb > ub:
            raise ValueError("Lower bound is greater than upper bound")

    def __eq__(self, other):
        return self.lb == other.lb and self.ub == other.ub

    def __repr__(self):
        return "[{0}, {1}]".format(self.lb, self.ub)


class Metabolite(object):
    def __init__(self, name, boundary=False):
        self.__assert_name(name)
        self.__assert_boundary(boundary)
        self.__name = name
        self.__boundary = boundary

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, name):
        self.__assert_name(name)
        self.__name = name

    @property
    def boundary(self):
        return self.__boundary

    @boundary.setter
    def boundary(self, boundary):
        self.__assert_boundary(boundary)
        self.__boundary = boundary

    def __eq__(self, other):
        return self.name == other.name and self.boundary == other.boundary

    def __assert_name(self, name):
        if not isinstance(name, str):
            raise TypeError("Metabolite name is not a string")
        if not len(name):
            raise ValueError("Metabolite name is empty string")

    def __assert_boundary(self, boundary):
        if not isinstance(boundary, bool):
            raise TypeError("Metabolite boundary condition is not a boolean")

    def __rmul__(self, other):
        if isinstance(other, (float, int)):
            return ReactionMember(metabolite=self, coefficient=other)

        raise TypeError("Can only multiply by numeric coefficient")

    def __repr__(self):
        return self.name

class ReactionMember(object):
    def __init__(self, metabolite, coefficient):
        self.__assert_valid_metabolite(metabolite)
        self.__assert_valid_coefficient(coefficient)

        self.__metabolite = metabolite
        self.__coefficient = float(coefficient)

    @property
    def metabolite(self):
        return self.__metabolite

    @metabolite.setter
    def metabolite(self, metabolite):
        self.__assert_valid_metabolite(metabolite)
        self.__metabolite = float(metabolite)

    @property
    def coefficient(self):
        return self.__coefficient

    @coefficient.setter
    def coefficient(self, coefficient):
        self.__assert_valid_coefficient(coefficient)
        self.__coefficient = float(coefficient)

    def __add__(self, other):
        if isinstance(other, ReactionMember):
            return ReactionMemberList([self, other])
        elif isinstance(other, ReactionMemberList):
            rml = ReactionMemberList(other)
            rml.insert(0, self)
            return rml
        else:
            raise TypeError("Can only join ReactionMember objects")

    def __repr__(self):
        return "{0:.5g} {1}".format(self.coefficient, self.metabolite)

    def __eq__(self, other):
        return self.metabolite == other.metabolite and self.coefficient == other.coefficient

    def __assert_valid_metabolite(self, metabolite):
        if not isinstance(metabolite, Metabolite):
            raise TypeError("Reaction member is not of type <Metabolite>")

    def __assert_valid_coefficient(self, coefficient):
        if not isinstance(coefficient, (int, float)):
            raise TypeError("Reaction member coefficient is not a number")
        if coefficient <= 0:
            raise ValueError("Reaction member coefficient is not strictly positive")

class Direction(object):
    __lockObj = thread.allocate_lock()
    __forward = None
    __reversible = None

    def __init__(self, type):
        if not type in ["f", "r"]:
            raise ValueError("Invalid reaction direction type")

        self.__type = type

    @staticmethod
    def forward():
        Direction.__lockObj.acquire()
        try:
            if Direction.__forward is None:
                Direction.__forward = Direction("f")
        finally:
            Direction.__lockObj.release()

        return Direction.__forward

    @staticmethod
    def reversible():
        Direction.__lockObj.acquire()
        try:
            if Direction.__reversible is None:
                Direction.__reversible = Direction("r")
        finally:
            Direction.__lockObj.release()

        return Direction.__reversible

    def __repr__(self):
        if self.__type == "f":
            return "->"
        if self.__type == "r":
            return "<->"

class ReactionMemberList(list):
    def __radd__(self, other):
        if isinstance(other, ReactionMember):
            rml = ReactionMemberList()
            rml.append(other)
            rml.extend(self)
            return rml

        return super(ReactionMemberList, self).__radd__(other)

    def __add__(self, other):
        if isinstance(other, ReactionMember):
            rml = ReactionMemberList()
            rml.extend(self)
            rml.append(other)
            return rml

        return super(ReactionMemberList, self).__add__(other)

    def __iadd__(self, other):
        if isinstance(other, ReactionMember):
            self.append(other)
            return self
        elif isinstance(other, ReactionMemberList):
            self.extend(other)
            return self

        return super(ReactionMemberList, self).__iadd__(other)

    def __repr__(self):
        return " + ".join(m.__repr__() for m in self)

class Reaction(object):
    def __init__(self, name, reactants=ReactionMemberList(), products=ReactionMemberList(), direction=None, bounds=Bounds()):
        if direction is None:
            direction = bounds.direction

        self.__assert_name(name)
        self.__assert_members(reactants)
        self.__assert_members(products)
        self.__assert_direction(direction)
        self.__assert_bounds(bounds)

        self.__name = name
        self.__direction = direction
        self.__bounds = bounds

        if isinstance(reactants, ReactionMember):
            self.__reactants = ReactionMemberList([reactants])
        else:
            self.__reactants = reactants

        if isinstance(products, ReactionMember):
            self.__products = ReactionMemberList([products])
        else:
            self.__products = products

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, name):
        self.__assert_name(name)
        self.__name = name

    @property
    def reactants(self):
        return self.__reactants

    @reactants.setter
    def reactants(self, reactants):
        self.__assert_members(reactants)
        if isinstance(reactants, ReactionMember):
            self.__reactants = ReactionMemberList([reactants])
        else:
            self.__reactants = reactants

    @property
    def products(self):
        return self.__products

    @products.setter
    def products(self, products):
        self.__assert_members(self, products)
        if isinstance(products, ReactionMember):
            self.__products = ReactionMemberList([products])
        else:
            self.__products = products

    @property
    def direction(self):
        return self.__direction

    @direction.setter
    def direction(self, direction):
        self.__assert_direction(self, direction)
        self.direction = direction

    @property
    def bounds(self):
        return self.__bounds

    @bounds.setter
    def bounds(self, bounds):
        self.__assert_bounds(self, bounds)
        self.__bounds = bounds

    def find_effective_bounds(self):
        lb = 0 if self.__direction == Direction.forward() and self.__bounds.lb < 0 else self.__bounds.lb
        ub = self.__bounds.ub

        return Bounds(lb, ub)

    def __assert_name(self, name):
        if not isinstance(name, str):
            raise TypeError("Reaction name is not a string")

    def __assert_members(self, reactants):
        if not isinstance(reactants, (ReactionMemberList, ReactionMember)):
            raise TypeError("Reaction reactants is not of type ReactionMemberList or ReactionMember")

    def __assert_direction(self, direction):
        if not isinstance(direction, Direction):
            raise TypeError("Reaction direction is not of type ReactionDirection")

    def __assert_bounds(self, bounds):
        if not isinstance(bounds, Bounds):
            raise TypeError("Reaction bounds is not of type bounds")

    def __repr__(self):
        return "{name}: {lhs} {dir} {rhs}".format(name=self.name, lhs=self.reactants, dir=self.direction, rhs=self.products)

    def __eq__(self, other):
        return self.name == other.name and \
               self.reactants == other.reactants and \
               self.products == other.products and \
               self.bounds == other.bounds and \
               self.direction == other.direction

class Objective(object):
    def __init__(self, metabolites):

class Model(object):
    def __init__(self):
        self.__reactions = list()

    @property
    def reactions(self):
        return self.__reactions

    def find_metabolites(self):
        return set(rm.metabolite for r in self.reactions for rm in r.reactants + r.products)

    def find_boundary_metabolites(self):
        return set(m for m in self.find_metabolites() if m.boundary)

    def __repr__(self):
        ret = "-REACTIONS\n{0}\n\n".format("\n".join(r.__repr__() for r in self.reactions))
        ret += "-CONSTRAINTS\n{0}\n\n".format("\n".join("{0}\t{1}".format(r.name, r.find_effective_bounds()) for r in self.reactions))
        ret += "-EXTERNAL METABOLITES\n{0}".format("\n".join(m.name for m in self.find_boundary_metabolites()))

        return ret