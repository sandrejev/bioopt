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

    def __assert_valid(self, lb, ub):
        if not isinstance(ub, (int, float)):
            raise TypeError("Upper bound is not a number")

        if not isinstance(ub, (int, float)):
            raise TypeError("Upper bound is not a number")

        if not isinstance(lb, (int, float)):
            raise TypeError("Lower bound is not a number")

        if lb > ub:
            raise ValueError("Lower bound is greater than upper bound")

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
    def __init__(self, name, reactants, products, direction, bounds=Bounds()):
        self.__assert_name(name)
        self.__assert_reactants(reactants)
        self.__assert_products(products)
        self.__assert_direction(direction)
        self.__assert_bounds(bounds)

        self.__name = name
        self.__reactants = reactants
        self.__products = products
        self.__direction = direction
        self.__bounds = bounds

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
        self.__assert_reactants(reactants)
        self.__reactants = reactants

    @property
    def products(self):
        return self.__products

    @products.setter
    def products(self, products):
        self.__assert_products(self, products)
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

    def __assert_reactants(self, reactants):
        if not isinstance(reactants, ReactionMemberList):
            raise TypeError("Reaction reactants is not of type ReactionMemberSet")

    def __assert_products(self, products):
        if not isinstance(products, ReactionMemberList):
            raise TypeError("Reaction products is not of type ReactionMemberSet")

    def __assert_direction(self, direction):
        if not isinstance(direction, Direction):
            raise TypeError("Reaction direction is not of type ReactionDirection")

    def __assert_bounds(self, bounds):
        if not isinstance(bounds, Bounds):
            raise TypeError("Reaction bounds is not of type bounds")

    def __repr__(self):
        return "{name}: {lhs} {dir} {rhs}".format(name=self.name, lhs=self.reactants, dir=self.direction, rhs=self.products)

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

#p = BiooptParser()
#p.parse_file("d:/Users/sandrejev/Desktop/iAG612.bioopt")

#        skip_flag     = 0
#        reactions_flag = 0
#        constraints_flag = 0
#        external_metabolites_flag = 0
#        obj_flag = 0
#        designobj_flag = 0
#
#        fails = 0
#        minor_fails = 0
#
#        reactions = {}
#        constraints = {}
#        external_metabolites = {}
#        obj = {}
#        design_obj = {}
#
#
#        for nr, line in enumerate(f):
#            line = line.strip()
#            c = nr + 1
#
#            #checking flags
#            if not line:
#                continue
#            if re.match('^#',line):
#                continue
#            if re.match('^%',line) and skip_flag == 0:
#                skip_flag = 1
#            elif re.match('^%',line) and skip_flag == 1:
#                skip_flag = 0
#                continue
#            if skip_flag == 1:
#                continue
#
#            if re.match('^-REACTIONS', line):
#                reactions_flag = 1
#                continue
#            if re.match('^-CONSTRAINTS', line) and reactions_flag == 1:
#                reactions_flag = 0
#                constraints_flag = 1
#                continue
#
#            elif re.match('^-CONSTRAINTS', line) and reactions_flag == 0:
#                print parse.__name__, "at line {}, check if -REACTIONS flag was defined before. {}".format(c, error(1))
#                fails += 1
#                continue
#
#            if re.match('^-EXTERNAL METABOLITES', line) and constraints_flag == 1:
#                constraints_flag = 0
#                external_metabolites_flag = 1
#                continue
#            elif re.match('^-EXTERNAL METABOLITES', line) and constraints_flag == 0:
#                print parse.__name__, "at line {}, check if -CONSTRAINTS flag was defined before. {}".format(c, error(1))
#                fails += 1
#                continue
#
#            if re.match('^-OBJ', line) and external_metabolites_flag == 1:
#                external_metabolites_flag = 0
#                obj_flag = 1
#                continue
#            elif re.match('^-OBJ', line) and external_metabolites_flag == 0:
#                print parse.__name__, "at line {}, check if -EXTERNAL METABOLITES flag was defined before. {}".format(c, error(1))
#                fails += 1
#                continue
#
#            if re.match('^-DESIGNOBJ', line) and obj_flag == 1:
#                obj_flag = 0
#                designobj_flag = 1
#                continue
#            elif re.match('^-OBJ', line) and external_metabolites_flag == 0:
#                print parse.__name__, "at line {}, check if -OBJ flag was defined before. {}".format(c, error(1))
#                fails += 1
#                continue
#
#            #colecting data
#            if reactions_flag == 1:
#                #print c, line
#                name, substrates, products, direction = make_reaction(line)
#                if name in reactions:
#                    print "Such reaction {} exists at {}".format(name,c)
#                    fails += 1
#                else:
#                    reactions[name] = [substrates,products,direction]
#                continue
#
#            if constraints_flag == 1:
#                re_number = r"(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?"
#                m = re.match(r'(\w+)\s?\[\s?'+'('+re_number+')'+'\s?,\s?'+'('+re_number+')'+'\s?\]', line)
#                if m:
#                    #re_number if matched returns 3 elements tuple
#                    #m will contain total 7 elements in 0th element will be name +3 for lb and 3 for ub
#                    name = m.groups()[0]
#                    lb   = m.groups()[1]
#                    ub   = m.groups()[4]
#                    if name in reactions:
#                        constraints[name] = [lb,ub]
#                    else:
#                        print "No such reaction {} for constraint {} at line {}".format(name,line,c)
#                else:
#                    print "Wrong format at {}, expected constraints format: reaction[number, number]".format(c)
#                    fails +=1
#                continue
#
#            if external_metabolites_flag == 1:
#                if line in external_metabolites:
#                    print "Duplicate found of external metabolite at {}".format(c)
#                    minor_fails += 1
#                else:
#                    external_metabolites[line] = 1
#                continue
#
#            if obj_flag == 1:
#                re_number = r"(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?"
#                m = re.match(r'(\w+)\s+'+'('+re_number+')'+'\s+'+'('+re_number+')', line)
#                if m:
#                    name = m.groups()[0]
#                    c1   = m.groups()[1]
#                    c2  = m.groups()[4]
#                    obj[name] = [c1,c2]
#                else:
#                    print "Wrong format at {}, expected objective format: reaction number number".format(c)
#                    fails +=1
#                continue
#            if designobj_flag == 1:
#                re_number = r"(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?"
#                m = re.match(r'(\w+)\s+'+'('+re_number+')'+'\s+'+'('+re_number+')', line)
#                if m:
#                    name = m.groups()[0]
#                    c1   = m.groups()[1]
#                    c2  = m.groups()[4]
#                    design_obj[name] = [c1,c2]
#                else:
#                    print "Wrong format at {}, expected design objective format: reaction number number".format(c)
#                    fails +=1
#                continue
#
#        if force == 1:
#
#            if fails:
#                print "There were critical {} fails found and {}".format(fails, minor_fails)
#                quit (-1)
#            else:
#
#                return reactions, constraints, external_metabolites, obj, design_obj
#        else:
#            return reactions, constraints, external_metabolites, obj, design_obj



path="./Scere_master_Kiran.txt"

def error(code):
    error_map = {1: "Check tags order: -REACTIONS -CONSTRAINTS -EXTERNAL METABOLITES -OBJ -DESIGNOBJ",
                2 : "Check format, seems that -REACTIONS is not found of after -CONSTRAINTS",
                3 : "Error3",
                4 : "Direction can be only f, b , r" }
    if code in error_map.keys():
        return error_map[code]
    else:
        return "wrong code"




def make_reaction (line, sep = ":"):
    line = line.strip()
    fails = 0

    parts = line.split(sep)

    if len(parts) != 2:
        print make_reaction.__name__ , "error: delimiter",sep, line
        fails += 1
        exit(-1)

    reaction_name = parts[0].strip()
    d = re.search("(\s+<\->|<\-|\->\s+)", line)
    direction = d.groups()[0]
    parts = parts[1].split(direction)
    #removing white space after splitting into parts
    direction = direction.strip()

    if len(parts) != 2:
        print make_reaction.__name__ , "error: direction",direction, line
        fails += 1
        exit(-1)

    left_part = parts[0]
    right_part = parts[1]

    substrates_part = re.split(r"\s+\+\s+",left_part)
    products_part   = re.split(r"\s+\+\s+",right_part)

    re_number = r"(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?"

    number = re.compile("[\(]?\s?"+"("+re_number+")"+"\s?[\)]?\s+")
    substrates = {}
    for metabolite in substrates_part:
        metabolite = metabolite.strip();
        coefficient = 1
        n = re.search(number, metabolite)
        if n:
            coefficient = float(n.groups()[0])
        metabolite = re.sub(number, "", metabolite)
        substrates[metabolite] = coefficient

    products = {}
    for metabolite in products_part:
        metabolite = metabolite.strip();
        coefficient = 1
        n = re.search(number, metabolite)
        if n:
            coefficient = float(n.groups()[0])
        metabolite = re.sub(number, "", metabolite)
        products[metabolite] = coefficient

    if direction == "->":
        direction = "f"
    elif direction == "<->":
        direction = "r"

    if direction == "<-":
        return reaction_name, products, substrates, "f"
    else:
        return reaction_name, substrates, products, direction













"""



a = parse_file(path)

print a


name, substrates, products, direction = make_reaction(r"reaction: 0.25 2'-bgdfg + 3445 z <-> 0.34354E10 2'-bgdfdfg ")

r = Reaction(name, Metabolites(substrates), Metabolites(products), direction)
"""