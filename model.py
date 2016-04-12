#!/usr/bin/env python
import re
import thread
import itertools
import warnings
import math

def _is_number(s):
    if s in ['0', '1', '2', '1000']:
        return True

    try:
        float(s)
        return True
    except ValueError:
        return False

def _starts_with_number(s):
    return s[0] in ['-', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0']

class Bounds(object):
    """
    :class:`Bounds` holds description of reactions constraints

    :param lb: Minimal amount of of flux that can go through a reaction (Lower bound). Negative numbers denote reverse direction.
    :param ub: Maximal amount of of flux that can go through a reaction (Upper bound). Negative numbers denote reverse direction.
    :return: :class:`Bounds`
    """

    def __init__(self, lb=float("-inf"), ub=float("inf")):
        self.__assert_valid(lb, ub)
        self.__lb = lb
        self.__ub = ub

    def copy(self):
        """
        Create a deep copy of current object

        :rtype: :class:`Bounds`
        """
        return Bounds(self.__lb, self.__ub)

    @property
    def lb_is_finite(self):
        """
        Returns inf False if lower bound is -infinity or +infinity
        """
        return self.__lb != self.inf() and self.__lb != -self.inf()

    @property
    def lb(self):
        """
        Minimal amount of of flux that can go through a reaction (Lower bound). Negative numbers denote reverse direction.
        """
        return self.__lb

    @lb.setter
    def lb(self, lb):
        self.__assert_valid(lb, self.ub)
        self.__lb = float(lb)

    @property
    def ub_is_finite(self):
        """
        Returns inf False if upper bound is -infinity or +infinity
        """
        return self.__ub != self.inf() and self.__ub != -self.inf()

    @property
    def ub(self):
        """
        Maximal amount of of flux that can go through a reaction (Upper bound). Negative numbers denote reverse direction.
        """
        return self.__ub

    @ub.setter
    def ub(self, value):
        self.__assert_valid(self.lb, value)
        self.__ub = float(value)

    @property
    def direction(self):
        """
        Suggested reaction direction. If lower bound is negative the reaction is suggested to be reversible. Otherwise
        irreversibility is implied. This is only a suggestion to :class:`Reaction` class and this suggestion can be
        broken by enabling reaction reversibility through :attr:`Reaction.direction`

        :rtype: :class:`Direction`
        """
        return Direction.forward() if self.lb >= 0 else Direction.reversible()

    @staticmethod
    def inf():
        """
        Returns infinity

        :return: :class:`float`
        """
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
    """
    :class:`Metabolite` holds information about metabolite. Currently only supported information is metabolite name
    and whether metabolite satisfies boundary condition (imported/exported)

    :param name: Metabolite name. Metabolite name should be a non-empty string.
    :param boundary: Boundary condition (imported/exported - True; otherwise - false)

    :return: :class:`Metabolite`
    """
    def __init__(self, name, boundary=False):
        self.__assert_name(name)
        self.__assert_boundary(boundary)
        self.__name = name
        self.__boundary = boundary
        self.__order_boundary = 0

    def copy(self):
        """
        Create a deep copy of current object

        :rtype: :class:`Metabolite`
        """
        m = Metabolite(self.__name, self.__boundary)
        m.order_boundary = self.__order_boundary

        return m

    @property
    def name(self):
        """
        Metabolite name. Metabolite name should be a non-empty string.
        """
        return self.__name

    @name.setter
    def name(self, name):
        self.__assert_name(name)
        self.__name = name

    @property
    def boundary(self):
        """
        Boundary condition (imported/exported - True; otherwise - false)
        """
        return self.__boundary

    @boundary.setter
    def boundary(self, boundary):
        self.__assert_boundary(boundary)
        self.__boundary = boundary

    @property
    def order_boundary(self):
        """
        Priority of displaying this metabolite. Metabolites with higher priority (lower numbers) are displayed earlier
        in "-External Metabolites" section
        """
        return self.__order_boundary

    @order_boundary.setter
    def order_boundary(self, order_boundary):
        if not _is_number(order_boundary):
            raise TypeError("Display priority should be a number: {0}".format(type(order_boundary)))

        self.__order_boundary = order_boundary

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
        if name == "":
            warnings.warn("Metabolite '{0}' contains spaces".format(name), UserWarning)
        # TODO: Do we need this?
        if _starts_with_number(name) and _is_number(name):
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
    """
    :class:`Bounds` is a wrapper for :class:`Metabolite` object when used in reaction reactants or products. It
    contains reference to the metabolite itself and to it's coefficient in the reaction.

    :param metabolite: Reference to :class:`Metabolite` object
    :param coefficient: Multiplier associated with metabolite
    :return: :class:`ReactionMember`
    """

    def __init__(self, metabolite, coefficient=1):
        self.__assert_metabolite(metabolite)
        self.__assert_coefficient(coefficient)

        self.__metabolite = metabolite
        self.__coefficient = float(coefficient)

    def copy(self):
        """
        Create a deep copy of current object

        :rtype: :class:`ReactionMember`
        """
        rm = ReactionMember(self.__metabolite.copy(), self.__coefficient)

        return rm

    @property
    def metabolite(self):
        """
        Reference to metabolite

        :rtype: :class:`Metabolite`
        """
        return self.__metabolite

    @metabolite.setter
    def metabolite(self, metabolite):
        self.__assert_metabolite(metabolite)
        self.__metabolite = metabolite

    @property
    def coefficient(self):
        """
        Multiplier associated with metabolite
        """
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

class Direction(object):
    """
    Describe reaction directionality. Don't use this class directly! Use factory constructors
    :meth:`Direction.forward` and  :meth:`Direction.reversible`

    :param type: f - Irreversible (**f** orward); r - Reversible (**r** eversible)
    :return: :class:`Direction`
    """

    __lockObj = thread.allocate_lock()
    __forward = None
    __reversible = None

    def __init__(self, type):
        if not type in ["f", "r"]:
            raise ValueError("Invalid reaction direction type. Allowed values are: {0}".format(', '.join(type)))

        self.__type = type

    def copy(self):
        """
        Create a deep copy of current object

        :rtype: :class:`Direction`
        """
        return Direction(self.__type)

    @staticmethod
    def forward():
        """
        Returns irreversible directionality descriptor singleton of type :class:`Direction`

        :return: :class:`Direction`
        """
        Direction.__lockObj.acquire()
        try:
            if Direction.__forward is None:
                Direction.__forward = Direction("f")
        finally:
            Direction.__lockObj.release()

        return Direction.__forward

    @staticmethod
    def reversible():
        """
        Returns reversible directionality descriptor singleton of type :class:`Direction`

        :return: :class:`Direction`
        """
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
    """
    :class:`ReactionMemberList` is a list of :class:`ReactionMember` instances. :class:`ReactionMemberList` inherits
    from :class:`list` all the usual functions to manage a list
    """
    def copy(self):
        """
        Create a deep copy of current object

        :rtype: :class:`ReactionMemberList`
        """
        rml = ReactionMemberList()
        rml.extend(rm.copy() for rm in self)

        return rml

    def find_member(self, name):
        """
        Find metabolite by name in the list

        :rtype: :class:`ReactionMember`
        """
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


class Reaction(object):
    """
    :class:`Reaction` class holds information about reaction including reaction name, members, directionality and constraints

    :param name: Reaction name. Reaction name should be non-empty string
    :param reactants: Reaction left-hand-side. Object of class :class:`ReactionMemberList`.
    :param products: Reaction right-hand-side. Object of class :class:`ReactionMemberList`.
    :param direction: Reaction direction. Object of class :class:`Direction`.
    :param bounds: Reaction constraints. Object of class :class:`Bounds`.
    :rtype: :class:`Reaction`
    """
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

    def copy(self):
        """
        Create a deep copy of current object

        :rtype: :class:`Reaction`
        """
        reactants = self.__reactants.copy()
        products = self.__products.copy()
        bounds = self.__bounds.copy()
        r = Reaction(self.name, reactants, products, direction=self.__direction.copy(), bounds=bounds)
        return r

    @property
    def name(self):
        """
        Reaction name
        """
        return self.__name

    @name.setter
    def name(self, name):
        self.__assert_name(name)
        self.__name = name


    @property
    def reactants(self):
        """
        Reactants

        :rtype: :class:`ReactionMemberList`
        """
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
        """
        Products

        :rtype: :class:`ReactionMemberList`
        """
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
        """
        Reaction direction

        :rtype: :class:`Direction`
        """
        return self.__direction

    @direction.setter
    def direction(self, direction):
        self.__assert_direction(direction)
        self.__direction = direction

    @property
    def bounds(self):
        """
        Reaction constraints

        :rtype: :class:`Bounds`
        """
        return self.__bounds

    @bounds.setter
    def bounds(self, bounds):
        self.__assert_bounds(bounds)
        self.__bounds = bounds

    def bounds_reset(self):
        """
        Reset bounds to predefined default. For reversible reaction defaults bounds are **[-inf, +inf]**. For forward
        reactions it is **[0, +inf]**
        """
        if self.direction == Direction.forward():
            self.bounds = Bounds(0, Bounds.inf())
        else:
            self.bounds = Bounds(-Bounds.inf(), Bounds.inf())


    def find_effective_bounds(self):
        """
        Find effective bounds. For example if reaction is described as irreversible but constraints are [-10, 10] the
        effective bounds would be [0, 10]

        :rtype: :class:`Bounds`
        """
        lb = 0 if self.__direction == Direction.forward() and self.__bounds.lb < 0 else self.__bounds.lb
        ub = self.__bounds.ub

        return Bounds(lb, ub)

    def reverse(self):
        """
        Reverse reactions. This functions change reactants and products places and inverses constraints so that reaction
        would go other way
        """

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
        if name == "":
            warnings.warn("Reaction '{0}' contains spaces".format(name), UserWarning)
        # TODO: Do we need this?
        if _starts_with_number(name) and _is_number(name):
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
    """
    Object describing operation type. **Don't use this class directly! Instead use factory constructors:**

    * :meth:`Operation.addition`
    * :meth:`Operation.subtraction`
    * :meth:`Operation.multiplication`
    * :meth:`Operation.division`
    * :meth:`Operation.negation` (Unary operation)

    :param operation: String describing the operation (binary operations: ``+-/*``, unary operations: ``-``)
    :param is_unary: Describes what type of operation is created. Unary operations support only one argument, while binary support two.
    :return: :class:`Operation`
    """

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
        """
        Is operation unary. True - unary; otherwise False
        :return: :class:`bool`
        """
        return self.__is_unary

    @property
    def symbol(self):
        """
        Short description of operation

        * + Addition
        * - Subtraction
        * / Division
        * * Multiplication
        * - Negation

        :rtype: :class:`Operation`
        """
        return self.__operation

    @property
    def priority(self):
        """
        Operation priority. Operations with higher priority (lower numbers) are executed before lower priority operations

        :rtype: :class:`Operation`
        """
        return self.__priority

    @staticmethod
    def addition():
        """
        Returns addition operation singleton

        :rtype: :class:`Operation`
        """
        return Operation.__create_singleton(False, "+", Operation.__addition)

    @staticmethod
    def subtraction():
        """
        Returns subtraction operation singleton

        :rtype: :class:`Operation`
        """
        return Operation.__create_singleton(False, "-", Operation.__subtraction)

    @staticmethod
    def multiplication():
        """
        Returns multiplication operation singleton

        :rtype: :class:`Operation`
        """
        return Operation.__create_singleton(False, "*", Operation.__multiplication)

    @staticmethod
    def division():
        """
        Returns division operation singleton

        :rtype: :class:`Operation`
        """
        return Operation.__create_singleton(False, "/", Operation.__division)

    @staticmethod
    def negation():
        """
        Returns negation operation (unary) singleton

        :rtype: :class:`Operation`
        """
        return Operation.__create_singleton(True, "-", Operation.__negation)

    def __repr__(self):
        return self.symbol


class MathExpression(object):
    def __init__(self, operation, operands):
        """
        Class describing mathematical expression

        :param operation: Reference to :class:`Operation`
        :param operands: A list of operands (two for binary, one for unary). An operand can be a constant (number) or
        another :class:`MathExpression`
        :return: :class:`MathExpression`
        """
        self.__assert_valid(operation, operands)

        self.__operands = operands
        self.__operation = operation

    @property
    def operation(self):
        """
        Reference to :class:`Operation`

        :rtype: :class:`Operation`
        """
        return self.__operation

    @operation.setter
    def operation(self, operation):
        self.__assert_valid(operation, self.operands)
        self.__operation = operation

    @property
    def operands(self):
        """
        A list of operands (two for binary, one for unary). An operand can be a constant (number) or
        another :class:`MathExpression`

        :rtype: list of :class:`Operation`
        """
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

    @staticmethod
    def format_var(var):
        var = "({0})".format(var) if isinstance(var, MathExpression) else var
        var = var.name if isinstance(var, (Reaction, Metabolite)) else var

        return var

    def __indent(self, text, indent="  "):
        return "\n".join([indent + l for l in text.split("\n")])

    def __repr__(self, tree=False):
        if not tree:
            return " {0} ".format(self.operation).join(str(self.format_var(o)) for o in self.operands)
        else:
            operands = []
            for o in self.operands:
                operands.append(o.__repr__(tree=True) if isinstance(o, MathExpression) else "{0}({1})".format(type(o).__name__, self.format_var(o)))

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
    """
    BioOpt model is a main class in package. It contains list of reactions in the model and other additional information.
    """

    def __init__(self):
        self.__reactions = list()
        self.__objective = None
        self.__design_objective = None

    @property
    def reactions(self):
        """
        List of reactions in the model

        :rtype: list of :class:`Reaction`
        """
        return self.__reactions

    @reactions.setter
    def reactions(self, reactions):
        # TODO: assert
        self.__reactions = reactions

    @property
    def objective(self):
        """
        Optimization target (i.e. Biomass)

        :rtype: class:`MathExpression`
        """
        return self.__objective

    @objective.setter
    def objective(self, objective):
        self.__assert_objective(objective)
        self.__objective = objective

    @staticmethod
    def __extract_expression(expression):
        r = next(r for r in expression.operands if isinstance(r, Reaction))
        c = next(r for r in expression.operands if isinstance(r, (int, float)))

        return r, c

    @staticmethod
    def __extract_objective_dict(objective):
        coefficients = {}
        if objective.operation == Operation.addition():
            for exp in objective.operands:
                r, c = Model.__extract_expression(exp)
                coefficients[r.name] = c
        elif objective.operation == Operation.multiplication():
            r, c = Model.__extract_expression(objective)
            coefficients[r.name] = c

        return coefficients

    @property
    def objective_dict(self):
        return Model.__extract_objective_dict(self.objective)

    @property
    def design_objective_dict(self):
        return Model.__extract_objective_dict(self.objective)

    def get_objective_coefficient(self, reaction):
        # TODO: Implement generic operation handling


        self.objective
        model.objective.operands[2]
        pass

    @property
    def design_objective(self):
        """
        Design optimization target (i.e. Ethanol)

        :rtype: list of :class:`MathExpression`
        """
        return self.__design_objective

    @design_objective.setter
    def design_objective(self, design_objective):
        self.__assert_objective(design_objective)
        self.__design_objective = design_objective

    def find_reaction(self, names=None, regex=False):
        """
        Searches model for reactions with specified names (patterns). If multiple reactions are found this function
        returns only the first one.

        :param names: name or list of names of reactions
        :param regex: If True names argument is assumed to be a regular expression

        :rtype: :class:`Reaction`
        """
        r = self.find_reactions(names, regex=regex)
        if len(r) > 1:
            warnings.warn("Found {0} reactions corresponding to name '{1}'. Returning first!".format(len(r), names), RuntimeWarning)

        return r[0] if len(r) else None


    def find_reactions(self, names=None, regex=False):
        """
        Searches model for reactions with specified names (patterns). This function is capable of returning multiple
        reactions.

        :param names: name or list of names of reactions
        :param regex: If True names argument is assumed to be a regular expression

        :rtype: list of :class:`Reaction`
        """
        if not self.reactions:
            return []

        import collections
        if names is None:
            return [r for r in self.reactions]
        elif isinstance(names, str):
            if regex:
                names = re.compile(names)
                return [r for r in self.reactions if names.search(r.name)]
            else:
                for r in self.reactions:
                    if r.name == names:
                        return [r]
            return []
        elif isinstance(names, collections.Iterable):
            names = set(names)
            if regex:
                names = [re.compile(n) for n in names]
                return [r for r in self.reactions if any(n.search(r.name) for n in names)]
            else:
                return [r for r in self.reactions if r.name in names]
        else:
            raise TypeError("Names argument should be iterable, string or <None>")

    def unify_metabolite_references(self):
        metabolites = dict((m.name, m) for m in self.find_metabolites())
        for reaction in self.reactions:
            for member in reaction.reactants:
                member.metabolite = metabolites[member.metabolite.name]
            for member in reaction.products:
                member.metabolite = metabolites[member.metabolite.name]

    def __unify_objective_references(self, expression, reactions):
        if isinstance(expression, MathExpression):
            for i, o in enumerate(expression.operands):
                if isinstance(o, MathExpression):
                    self.__unify_objective_references(o, reactions)
                elif isinstance(o, Reaction):
                    if o.name in reactions:
                        expression.operands[i] = reactions[o.name]

    def unify_reaction_references(self):
        # TODO: What if more than one reaction with same name (Use first)
        reactions = dict((r.name, r) for r in self.reactions)
        self.__unify_objective_references(self.objective, reactions)
        self.__unify_objective_references(self.design_objective, reactions)

    def unify_references(self):
        """
        Find metabolite with identical names which are stored as different instances and unify instances.
        """
        self.unify_metabolite_references()
        self.unify_reaction_references()

    def __fix_math_reactions(self, expression, reactions):
        if isinstance(expression, MathExpression):
            for i, o in enumerate(expression.operands):
                if isinstance(o, MathExpression):
                    self.__fix_math_reactions(o, reactions)
                elif isinstance(o, Reaction):
                    expression.operands[i] = reactions[o.name]


    def find_metabolite(self, names=None, regex=False):
        """
        Searches model for metabolites with specified names (patterns). If multiple metabolites are found this function
        returns only the first one.

        :param names: name or list of names of metabolites
        :param regex: If True names argument is assumed to be a regular expression

        :rtype: :class:`Metabolite`
        """
        m = self.find_metabolites(names, regex=regex)
        if len(m) > 1:
            warnings.warn("Found {0} metabolites corresponding to name '{1}'. Returning first!".format(len(m), names), RuntimeWarning)

        return m[0] if len(m) else None

    # TODO: Test with multiple instances of the same metabolite (result should contain two instances)
    # TODO: Test with no reaction section. Metabolite present in ext. metabolites should be accessible
    def find_metabolites(self, names=None, regex=False):
        """
        Searches model for metabolites with specified names (patterns). This function is capable of returning multiple
        metabolites.

        :param names: name or list of names of metabolites
        :param regex: If True names argument is assumed to be a regular expression

        :rtype: list of :class:`Metabolite`
        """
        metabolites_set = set()
        metabolites = []
        for r in self.reactions:
            for rm in r.reactants:
                if rm.metabolite.name == "atp_c":
                    pass
                metabolites_set.add(rm.metabolite)

            for rm in r.products:
                if rm.metabolite.name == "atp_c":
                    pass
                metabolites_set.add(rm.metabolite)
            #if rm.metabolite not in metabolites_set:
            #metabolites.append(rm.metabolite)

        metabolites = list(metabolites_set)

        import collections
        if names is None:
            return metabolites
        elif isinstance(names, str):
            if regex:
                names = re.compile(names)
                return [m for m in metabolites if names.search(m.name)]
            else:
                for m in metabolites:
                    if m.name == names:
                        return [m]
            return []
        elif isinstance(names, collections.Iterable):
            names = set(names)
            if regex:
                names = [re.compile(n) for n in names]
                return [m for m in metabolites if any(n.search(m.name) for n in names)]
            else:
                return [m for m in metabolites if m.name in names]
        else:
            raise TypeError("Names argument should be iterable, string or <None>")


    def find_boundary_metabolites(self):
        """
        Searches for imported/exported metabolites

        :rtype: list of :class:`Metabolite`
        """
        return [m for m in self.find_metabolites() if m.boundary]


    def find_boundary_reactions(self):
        """
        Searches for reactions exporting or importing metabolites

        :rtype: list of :class:`Reaction`
        """
        return [r for r in self.reactions if any(m.metabolite.boundary for m in r.reactants) or any(m.metabolite.boundary for m in r.products)]

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
    def commune(models, model_prefix="ML{0:04d}_", env_prefix="ENV_", block=[]):
        """
        Merge two or more models into community model. Community model allows organisms represented by models to share
        metabolites. Briefly, the algorithm first appends reaction and metabolite names in original models with ML****_
        prefix. Then for every metabolite exported in original models boundary condition is set to :class:`False` and new
        export reaction is created. This export reactions allows metabolites to travel from joined models into shared
        environment and back. This setup allows organisms to exchange metabolites through common environment.

        Originally this merging framework was described in `"OptCom: A Multi-Level Optimization Framework for the Metabolic Modeling and Analysis of Microbial Communities" <http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002363>`_
        by Ali R. Zomorrodi and Costas D. Maranas.

        :param models: List of :class:`Model` to join
        :param model_prefix: Model prefix, Model prefix is added to all reaction names to avoid name collision in joined model.
        :param env_prefix: Prefix of metabolites in shared environment.
        :param block: List of names (in original models) of metabolites which should be not allowed to be exchanged between
                organisms. An obvious example of such metabolite is biomass.
        :rtype: :class:`Model`
        """
        block = [block] if not isinstance(block, list) else block
        block = [re.compile(b) for b in block]

        boundary_metabolites = {}

        fwd = Direction.forward()

        model = Model()
        for i, mod in enumerate(models):
            env_reactions = []
            for m in mod.find_boundary_metabolites():
                m_in = Metabolite(model_prefix.format(i) + m.name)
                m_out = Metabolite(env_prefix + m.name)
                boundary_metabolites[m_out.name] = m_out

                r_out = Reaction(model_prefix.format(i) + 'OUT_' + m.name, ReactionMemberList([ReactionMember(m_in, 1)]),
                    ReactionMemberList([ReactionMember(m_out, 1)]), Direction.forward(), Bounds(0, Bounds.inf()))
                r_in = Reaction(model_prefix.format(i) + 'IN_' + m.name, ReactionMemberList([ReactionMember(m_out, 1)]),
                    ReactionMemberList([ReactionMember(m_in, 1)]), Direction.forward(), Bounds(0, Bounds.inf()))
                env_reactions.extend([r_out, r_in])

            reactions = []
            metabolites = set()
            for r in mod.reactions:
                r = r.copy()
                r.name = model_prefix.format(i) + r.name
                for m in r.reactants:
                    metabolites.add(m.metabolite)
                for m in r.products:
                    metabolites.add(m.metabolite)

                reactions.append(r)

            for m in metabolites:
                m.name = model_prefix.format(i) + m.name
                m.boundary = False

            for r in reactions:
                if not any(b.search(r.name) for b in block):
                    model.reactions.append(r)

            for r in env_reactions:
                if not any(b.search(r.name) for b in block):
                    model.reactions.append(r)

        model.unify_references()

        for m_name, m_out in boundary_metabolites.items():
            m_ext = Metabolite("{0}xtX".format(m_name), True)

            r_out = Reaction(m_out.name+"xtO", ReactionMemberList([ReactionMember(m_out)]), ReactionMemberList([ReactionMember(m_ext)]), fwd)
            r_in = Reaction(m_out.name+"xtI", ReactionMemberList([ReactionMember(m_ext)]), ReactionMemberList([ReactionMember(m_out)]), fwd)

            if not any(b.search(r_out.name) for b in block):
                model.reactions.append(r_out)

            if not any(b.search(r_in.name) for b in block):
                model.reactions.append(r_in)

        return model

    def __assert_objective(self, objective):
        if not (objective is None or isinstance(objective, MathExpression)):
            raise TypeError("Objective is not None or <MathExpression>: {0}".format(type(objective)))

    def save(self, path=None, inf=1000):
        """
        Save model on disc in bioopt format

        :param path: The name or full pathname of the file where the BioOpt model is to be written.
        :param inf: Number which would be used for constraints with infinite bounds
        """
        ret = "-REACTIONS\n"
        for r in self.reactions:
            reactants = " + ".join("{0}{1}".format("" if abs(m.coefficient) == 1 else "{0:.5g} ".format(m.coefficient), m.metabolite.name) for m in r.reactants)
            products = " + ".join("{0}{1}".format("" if abs(m.coefficient) == 1 else "{0:.5g} ".format(m.coefficient), m.metabolite.name) for m in r.products)
            dir = "->" if r.direction == Direction.forward() else "<->"
            ret += "{name}\t:\t{lhs} {dir} {rhs}".format(name=r.name, lhs=reactants, dir=dir, rhs=products) + "\n"
        ret += "\n"

        ret += "-CONSTRAINTS\n"
        for r in self.reactions:
            lb = -inf if r.bounds.lb == -Bounds.inf() else r.bounds.lb
            ub = inf if r.bounds.ub == Bounds.inf() else r.bounds.ub
            if not (r.bounds.direction == Direction.forward() and r.bounds == Bounds(0)) and \
                    not (r.bounds.direction == Direction.reversible() and r.bounds == Bounds()):
                ret += "{0}\t[{1:.5g}, {2:.5g}]".format(r.name, lb, ub) + "\n"

        ret += "\n"

        ret += "-EXTERNAL METABOLITES\n"

        b_metabolites = self.find_boundary_metabolites()
        b_metabolites = sorted(b_metabolites, key=lambda x: x.order_boundary)
        for m in b_metabolites:
            ret += m.name + "\n"
        ret += "\n"

        if self.objective:
            ret += "-OBJECTIVE\n"
            ret += " ".join(str(MathExpression.format_var(o)) for o in self.objective.operands)
            ret += "\n\n"

        if self.design_objective:
            ret += "-DESIGN OBJECTIVE\n"
            ret += " ".join(str(MathExpression.format_var(o)) for o in self.design_objective.operands)
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