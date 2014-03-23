#!/usr/bin/env python
import sys,re
import thread
import itertools
import warnings
import copy
import math

def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

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

    @staticmethod
    def inf():
        return float("Inf")

    def __assert_valid(self, lb, ub):
        if not isinstance(ub, (int, float)):
            raise TypeError("Upper bound is not a number: {0}".format(type(ub)))

        if not isinstance(lb, (int, float)):
            raise TypeError("Lower bound is not a number: {0}".format(type(lb)))

        if lb > ub:
            raise ValueError("Lower bound is greater than upper bound ({0} > {1})".format(lb, ub))

    def __eq__(self, other):
        return type(self) == type(other) and \
               self.lb == other.lb and \
               self.ub == other.ub

    def __ne__(self, other):
        return not self.__eq__(other)

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
        return type(self) == type(other) and \
               self.name == other.name and \
               self.boundary == other.boundary

    def __ne__(self, other):
        return not self.__eq__(other)

    def __assert_name(self, name):
        if not isinstance(name, str):
            raise TypeError("Metabolite name is not a string: {0}".format(type(name)))
        if not len(name):
            raise ValueError("Metabolite name is empty string")
        if re.search("\s+", name):
            warnings.warn("Metabolite '{0}' contains spaces".format(name), UserWarning)
        if _is_number(name):
            warnings.warn("Metabolite name is a number: '{0}'".format(name), UserWarning)

    def __assert_boundary(self, boundary):
        if not isinstance(boundary, bool):
            raise TypeError("Metabolite boundary condition is not a boolean: {0}".format(type(boundary)))

    def __rmul__(self, other):
        if isinstance(other, (float, int)):
            return ReactionMember(metabolite=self, coefficient=other)

        raise TypeError("Can only multiply by numeric coefficient")

    def __repr__(self):
        b = "*" if self.boundary else ""
        return "{0}{1}".format(self.name, b)

class ReactionMember(object):
    def __init__(self, metabolite, coefficient):
        self.__assert_metabolite(metabolite)
        self.__assert_coefficient(coefficient)

        self.__metabolite = metabolite
        self.__coefficient = float(coefficient)

    @property
    def metabolite(self):
        return self.__metabolite

    @metabolite.setter
    def metabolite(self, metabolite):
        self.__assert_metabolite(metabolite)
        self.__metabolite = metabolite

    @property
    def coefficient(self):
        return self.__coefficient

    @coefficient.setter
    def coefficient(self, coefficient):
        self.__assert_coefficient(coefficient)
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
        return type(self) == type(other) and \
               self.metabolite == other.metabolite and \
               self.coefficient == other.coefficient

    def __ne__(self, other):
        return not self.__eq__(other)

    def __assert_metabolite(self, metabolite):
        if not isinstance(metabolite, Metabolite):
            raise TypeError("Reaction member is not of type <Metabolite>: {0}".format(type(metabolite)))

    def __assert_coefficient(self, coefficient):
        if not isinstance(coefficient, (int, float)):
            raise TypeError("Reaction member coefficient is not a number: {0}".format(type(coefficient)))
        if coefficient <= 0:
            raise ValueError("Reaction member coefficient is not strictly positive: {0}".format(coefficient))

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

    def __eq__(self, other):
        return type(self) == type(other) and self.__type == other.__type

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        if self.__type == "f":
            return "->"
        if self.__type == "r":
            return "<->"

class ReactionMemberList(list):
    def find_member(self, name):
        for mb in self:
            if mb.metabolite.name == name:
                return mb

        return None

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

# TODO: Add class ReactionList
class Reaction(object):
    def __init__(self, name, reactants=ReactionMemberList(), products=ReactionMemberList(), direction=None, bounds=None):
        if bounds is None and direction is None:
            direction = Direction.reversible()
        if direction is None and bounds is not None:
            direction = bounds.direction
        if bounds is None and direction is not None:
            bounds = Bounds(-float('inf'), float('inf')) if direction == Direction.reversible() else Bounds(0, float('inf'))

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

    def find_participants(self):
        return ReactionMemberList(itertools.chain(self.__reactants, self.__products))

    @property
    def products(self):
        return self.__products

    @products.setter
    def products(self, products):
        self.__assert_members(products)
        if isinstance(products, ReactionMember):
            self.__products = ReactionMemberList([products])
        else:
            self.__products = products

    @property
    def direction(self):
        return self.__direction

    @direction.setter
    def direction(self, direction):
        self.__assert_direction(direction)
        self.__direction = direction

    @property
    def bounds(self):
        return self.__bounds

    @bounds.setter
    def bounds(self, bounds):
        self.__assert_bounds(bounds)
        self.__bounds = bounds

    def find_effective_bounds(self):
        lb = 0 if self.__direction == Direction.forward() and self.__bounds.lb < 0 else self.__bounds.lb
        ub = self.__bounds.ub

        return Bounds(lb, ub)

    def reverse(self):
        if self.direction != Direction.reversible():
            raise RuntimeError("Reaction direction is not reversible. Only reversible reactions can be reversed")
        if self.bounds.lb > 0 and self.bounds.ub > 0:
            raise RuntimeError("Reaction effective direction is strictly forward and cannot be reversed")

        tmp = self.__products
        self.__products = self.__reactants
        self.__reactants = tmp
        self.__bounds = Bounds(-self.bounds.ub, -self.bounds.lb)

    def __assert_name(self, name):
        if not isinstance(name, str):
            raise TypeError("Reaction name is not a string: {0}".format(type(name)))
        if not len(name):
            raise ValueError("Reaction name is empty string")
        if re.search("\s+", name):
            warnings.warn("Reaction '{0}' contains spaces".format(name), UserWarning)
        if _is_number(name):
            warnings.warn("Reaction name is a number: '{0}'".format(name), UserWarning)

    def __assert_members(self, reactants):
        if not isinstance(reactants, (ReactionMemberList, ReactionMember)):
            raise TypeError("Reaction reactants is not of type ReactionMemberList or ReactionMember: {0}".format(type(reactants)))

    def __assert_direction(self, direction):
        if not isinstance(direction, Direction):
            raise TypeError("Reaction direction is not of type ReactionDirection: {0}".format(type(direction)))

    def __assert_bounds(self, bounds):
        if not isinstance(bounds, Bounds):
            raise TypeError("Reaction bounds is not of type bounds: {0}".format(type(bounds)))

    def __repr__(self):
        return "{name}{bnds}: {lhs} {dir} {rhs}".format(name=self.name, lhs=self.reactants, dir=self.direction, rhs=self.products, bnds=self.bounds)

    def __eq__(self, other):
        return type(self) == type(other) and \
               self.name == other.name and \
               self.reactants == other.reactants and \
               self.products == other.products and \
               self.bounds == other.bounds and \
               self.direction == other.direction

    def __ne__(self, other):
        return not self.__eq__(other)


class Operation(object):
    __lockObj = thread.allocate_lock()
    __addition = []
    __subtraction = []
    __negation = []
    __multiplication = []
    __division = []

    __unary_priority = ["-"]
    __binary_priority = ["*", "/", "+", "-"]

    def __init__(self, operation, is_unary):
        if not isinstance(is_unary, bool):
            raise TypeError("Parameter is_unary is not of type bool: {0}".format(type(is_unary)))
        if not is_unary and operation not in Operation.__binary_priority:
            raise ValueError("Invalid binary operation: {0}".format(operation))
        if is_unary and operation not in Operation.__unary_priority:
            raise ValueError("Invalid unary operation: {0}".format(operation))

        self.__is_unary = is_unary
        self.__operation = operation

        if is_unary:
            priority = Operation.__unary_priority.index(operation)
        else:
            priority = len(Operation.__unary_priority)
            priority += Operation.__binary_priority.index(operation)

        self.__priority = priority

    @staticmethod
    def __create_singleton(type, operation, instance):
        Operation.__lockObj.acquire()
        try:
            if not instance:
                instance.append(Operation(operation, type))
        finally:
            Operation.__lockObj.release()

        return instance[0]

    @property
    def is_unary(self):
        return self.__is_unary

    @property
    def symbol(self):
        return self.__operation

    @property
    def priority(self):
        return self.__priority

    @staticmethod
    def addition():
        return Operation.__create_singleton(False, "+", Operation.__addition)

    @staticmethod
    def subtraction():
        return Operation.__create_singleton(False, "-", Operation.__subtraction)

    @staticmethod
    def multiplication():
        return Operation.__create_singleton(False, "*", Operation.__multiplication)

    @staticmethod
    def division():
        return Operation.__create_singleton(False, "/", Operation.__division)

    @staticmethod
    def negation():
        return Operation.__create_singleton(True, "-", Operation.__negation)

    def __repr__(self):
        return self.symbol


class MathExpression(object):
    def __init__(self, operation, operands):
        self.__assert_valid(operation, operands)

        self.__operands = operands
        self.__operation = operation

    @property
    def operation(self):
        return self.__operation

    @operation.setter
    def operation(self, operation):
        self.__assert_valid(operation, self.operands)
        self.__operation = operation

    @property
    def operands(self):
        return self.__operands

    @operands.setter
    def operands(self, operands):
        self.__assert_valid(self.operation, operands)
        self.__operands = operands

    def __eq__(self, other):
        return type(self) == type(other) and \
               self.operands == other.operands and \
               self.operation == other.operation

    def __ne__(self, other):
        return not self.__eq__(other)

    def find_variables(self, remove_duplicates=True):
        vars = []
        for o in self.operands:
            if isinstance(o, type(self)):
                vars.extend(o.find_variables(False))
            elif not o is None:
                vars.append(o)

        if remove_duplicates:
            seen = set()
            return [x for x in vars if x not in seen and not seen.add(x)]

        return vars

    def __format_var(self, var):
        var = "({0})".format(var) if isinstance(var, MathExpression) and var.operation.priority > self.operation.priority else var
        var = var.name if isinstance(var, (Reaction, Metabolite)) else var

        return var

    def __indent(self, text, indent="  "):
        return "\n".join([indent + l for l in text.split("\n")])

    def __repr__(self, tree=False):
        if not tree:
            return " {0} ".format(self.operation).join(str(self.__format_var(o)) for o in self.operands)
        else:
            operands = []
            for o in self.operands:
                operands.append(o.__repr__(tree=True) if isinstance(o, MathExpression) else "{0}({1})".format(type(o).__name__, self.__format_var(o)))

            return "{0}(\n{1}\n)".format(
                type(self).__name__,
                self.__indent("\n{0}\n".format(self.operation).join(operands))
            )

    def __assert_valid(self, operation, operands):
        if not isinstance(operation, (Operation, type(None))):
            raise TypeError("Operation is not an instance of class <Operation>: {0}".format(type(operation)))
        if not isinstance(operands, list):
            raise TypeError("Operands are not a list: {0}".format(type(operands)))

        # TODO: test these exceptions
        if operation is None:
            if len(operands) != 1:
                raise ValueError("Math expression not representing any operation (<None>) have number of operands different from one")
        else:
            if operation.is_unary and len(operands) != 1:
                raise ValueError("Unary operation have number of operands different from one")
            elif not operation.is_unary and len(operands) < 2:
                raise ValueError("Binary operation have less than 2 operands")


class Model(object):
    def __init__(self):
        self.__reactions = list()
        self.__objective = None
        self.__design_objective = None

    @property
    def reactions(self):
        return self.__reactions

    @reactions.setter
    def reactions(self, reactions):
        # TODO: assert
        self.__reactions = reactions

    @property
    def objective(self):
        return self.__objective

    @objective.setter
    def objective(self, objective):
        self.__assert_objective(objective)
        self.__objective = objective

    @property
    def design_objective(self):
        return self.__design_objective

    @design_objective.setter
    def design_objective(self, design_objective):
        self.__assert_objective(design_objective)
        self.__design_objective = design_objective

    # TODO: test and test find_metabolites
    # TODO: Test with no reaction section. Reaction present in objective section should be accessible
    def find_reactions(self, names=None):
        if not self.reactions:
            return None

        import collections
        if names is None:
            return [r for r in self.reactions]
        elif isinstance(names, str):
            for r in self.reactions:
                if r.name == names:
                    return r
            return None
        elif isinstance(names, collections.Iterable):
            names = set(names)
            return [r for r in self.reactions if r.name in names]
        else:
            raise TypeError("Names argument should be iterable, string or <None>")

    def unify_metabolite_references(self):
        metabolites = dict((m.name, m) for m in self.find_metabolites())
        for reaction in self.reactions:
            for member in reaction.find_participants():
                member.metabolite = metabolites[member.metabolite.name]

    def __unify_objective_references(self, expression, reactions):
        if isinstance(expression, MathExpression):
            for i, o in enumerate(expression.operands):
                if isinstance(o, MathExpression):
                    self.__unify_objective_references(o, reactions)
                elif isinstance(o, Reaction):
                    if o.name in reactions:
                        expression.operands[i] = reactions[o.name]

    # TODO: unit tests
    def unify_reaction_references(self):
        # TODO: What if more than one reaction with same name (Use first)
        reactions = dict((r.name, r) for r in self.reactions)
        self.__unify_objective_references(self.objective, reactions)
        self.__unify_objective_references(self.design_objective, reactions)

    def unify_references(self):
        self.unify_metabolite_references()
        self.unify_reaction_references()

    def __fix_math_reactions(self, expression, reactions):
        if isinstance(expression, MathExpression):
            for i, o in enumerate(expression.operands):
                if isinstance(o, MathExpression):
                    self.__fix_math_reactions(o, reactions)
                elif isinstance(o, Reaction):
                    expression.operands[i] = reactions[o.name]

    # TODO: Test with multiple instances of the same metabolite (result should contain two instances)
    # TODO: Test with no reaction section. Metabolite present in ext. metabolites should be accessible
    def find_metabolites(self, names=None):
        metabolites_set = set()
        metabolites = []
        for r in self.reactions:
            for rm in r.find_participants():
                if rm.metabolite not in metabolites_set:
                    metabolites_set.add(rm.metabolite)
                    metabolites.append(rm.metabolite)

        import collections
        if names is None:
            return metabolites
        elif isinstance(names, str):
            for m in metabolites:
                if m.name == names:
                    return m
            return None
        elif isinstance(names, collections.Iterable):
            names = set(names)
            return [m for m in metabolites if m.name in names]
        else:
            raise TypeError("Names argument should be iterable, string or <None>")


    def find_boundary_metabolites(self):
        return [m for m in self.find_metabolites() if m.boundary]

    def get_max_bound(self):
        mb = 0
        for r in self.reactions:
            lb = math.fabs(r.bounds.lb)
            if lb > mb:
                mb = lb

            ub = math.fabs(r.bounds.ub)
            if ub > mb:
                mb = ub

        return mb

    @staticmethod
    def commune(models, model_prefix="ML{0:04d}_", env_prefix="ENV_"):
        model = Model()
        for i, mod in enumerate(models):
            mod_new = copy.deepcopy(mod)

            env_reactions = []
            for m in mod.find_boundary_metabolites():
                m_in = copy.deepcopy(m)
                m_in.name = model_prefix.format(i) + m.name
                m_in.boundary = False

                m_out = copy.deepcopy(m)
                m_out.name = env_prefix + m.name
                m_out.boundary = True

                r_out = Reaction(model_prefix.format(i) + 'OUT_' + m.name, ReactionMemberList([ReactionMember(m_in, 1)]), ReactionMemberList([ReactionMember(m_out, 1)]), Direction.forward(), Bounds(0, Bounds.inf()))
                r_in = Reaction(model_prefix.format(i) + 'IN_' + m.name, ReactionMemberList([ReactionMember(m_out, 1)]), ReactionMemberList([ReactionMember(m_in, 1)]), Direction.forward(), Bounds(0, Bounds.inf()))
                env_reactions.extend([r_out, r_in])

            for r in mod_new.reactions:
                r.name = model_prefix.format(i) + r.name

            for m in mod_new.find_metabolites():
                m.name = model_prefix.format(i) + m.name
                m.boundary = False


            model.reactions.extend(mod_new.reactions)
            model.reactions.extend(env_reactions)

        model.unify_references()

        for m_out in model.find_boundary_metabolites():
            m_out.boundary = False

            m_ext = copy.deepcopy(m_out)
            m_ext.name = "{0}xtX".format(m_ext.name)
            m_ext.boundary = True

            r_out = Reaction(m_out.name+"xtO", ReactionMemberList([ReactionMember(m_out, 1)]), ReactionMemberList([ReactionMember(m_ext, 1)]), Direction.forward())
            model.reactions.append(r_out)

            r_in = Reaction(m_out.name+"xtI", ReactionMemberList([ReactionMember(m_ext, 1)]), ReactionMemberList([ReactionMember(m_out, 1)]), Direction.forward())
            model.reactions.append(r_in)

        return model

    def __assert_objective(self, objective):
        if not (objective is None or isinstance(objective, MathExpression)):
            raise TypeError("Objective is not None or <MathExpression>: {0}".format(type(objective)))

    def save(self, path=None, inf=1000):
        ret = "-REACTIONS\n"
        for r in self.reactions:
            reactants = " + ".join("{0:.5g} {1}".format(m.coefficient, m.metabolite.name) for m in r.reactants)
            products = " + ".join("{0:.5g} {1}".format(m.coefficient, m.metabolite.name) for m in r.products)
            dir = "->" if r.direction is Direction.forward() else "<->"
            ret += "{name}\t:\t{lhs} {dir} {rhs}".format(name=r.name, lhs=reactants, dir=dir, rhs=products) + "\n"
        ret += "\n"

        ret += "-CONSTRAINTS\n"
        for r in self.reactions:
            lb = -inf if r.bounds.lb == -Bounds.inf() else r.bounds.lb
            ub = inf if r.bounds.ub == Bounds.inf() else r.bounds.ub
            ret += "{0}\t[{1}, {2}]".format(r.name, lb, ub) + "\n"
        ret += "\n"

        ret += "-EXTERNAL METABOLITES\n"
        for m in self.find_boundary_metabolites():
            ret += m.name + "\n"
        ret += "\n"

        if self.objective:
            ret += "-OBJECTIVE\n"
            ret += str(self.objective)
            ret += "\n\n"

        if self.design_objective:
            ret += "-DESIGN OBJECTIVE\n"
            ret += str(self.design_objective)
            ret += "\n\n"

        if path:
            f = open(path, 'w')
            f.write(ret)
            return f.close()
        else:
            return ret

    def __repr__(self):
        ret = "-REACTIONS\n{0}\n\n".format("\n".join(r.__repr__() for r in self.reactions))
        ret += "-CONSTRAINTS\n{0}\n\n".format("\n".join("{0}\t{1}".format(r.name, r.bounds) for r in self.reactions))
        ret += "-EXTERNAL METABOLITES\n{0}\n\n".format("\n".join(m.__repr__() for m in self.find_boundary_metabolites()))
        ret += "-OBJECTIVE\n{0}\n\n".format(self.objective)
        ret += "-DESIGN OBJECTIVE\n{0}\n\n".format(self.design_objective)

        return ret

    def __eq__(self, other):
        return type(self) == type(other) and \
               self.reactions == other.reactions and \
               self.objective == other.objective and \
               self.design_objective == other.design_objective

    def __ne__(self, other):
        return not self.__eq__(other)